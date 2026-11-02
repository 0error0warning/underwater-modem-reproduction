function records=run_link_v2(stage,resume)
%RUN_LINK_V2 Paired time-slot observations, with real acquisition and Doppler.
% Separate development/test files support policy freezing before test access.
% Slot-based delivered rate and counterfactual occupied-time efficiency must
% remain distinct. A separate sequential rollout validates actual scheduling.
if nargin<1, stage='development'; end
if nargin<2, resume=true; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','v2'); if ~exist(out,'dir'), mkdir(out); end
switch stage
    case 'development', recordIds=1:10;
    case 'test'
        frozenPath=fullfile(out,'adaptive','frozen_policy.json');
        assert(exist(frozenPath,'file')==2,'Freeze adaptive policy before final test');
        frozen=jsondecode(fileread(frozenPath));
        assert(strcmp(frozen.status,'FROZEN_BEFORE_FINAL_TEST'),'Policy is not frozen');
        recordIds=31:40;
    case 'smoke', recordIds=1;
    otherwise, error('Unknown stage');
end
metaFile=fullfile(out,['metadata_' stage '.json']);
save_json(metaFile,struct('status','RUNNING','stage',stage));
cfgBase=uac2.default_config();
fs=cfgBase.fs; nBits=128; slotSeconds=6.4; slotSamples=round(slotSeconds*fs); nSlots=5;
guards=[0,4,8,16,32]; snrLevels=[0,10,20]; channels={'NOF1','NCS1'};
% The predeclared policy is fitted only at 20 dB. Other levels belong to
% frozen-policy robustness evaluation, not extra tuning observations.
if strcmp(stage,'development'), snrLevels=20; end
baseSeed=2026110202;
rows=struct([]); rowCount=0; clockRun=tic;
recordsPath=fullfile(out,['records_' stage '.csv']);
if resume && exist(recordsPath,'file')==2
    previous=readtable(recordsPath);
    assert(all(ismember(previous.record,recordIds)) && all(ismember(previous.snr_db,snrLevels)), ...
        'Existing checkpoint does not match this experiment stage');
    rows=table2struct(previous); rowCount=numel(rows);
end
for ci=1:2
    channel=channels{ci};
    for recordNumber=recordIds
        if ~isempty(rows)
            prior=strcmp({rows.channel},channel) & [rows.record]==recordNumber;
            expected=numel(guards)*nSlots*numel(snrLevels)*3;
            if sum(prior)==expected
                fprintf('Resume: complete %s record=%d (%d rows)\n',channel,recordNumber,expected);
                continue;
            elseif any(prior)
                rows=rows(~prior); rowCount=numel(rows);
            end
        end
        channelFile=fullfile(root,'vendor','Watermark','input','channels',channel,'mat',sprintf('%s_%03d.mat',channel,recordNumber));
        data=load(channelFile);
        maximumSeconds=(size(data.h,1)-1)/data.fs_t;
        assert(nSlots*slotSeconds<maximumSeconds,'Sounding too short');
        meanPower=mean(abs(data.h).^2,1); [~,peak]=max(meanPower);
        oracleDelay=(peak-1)/double(data.fs_tau);
        % Independent no-guard BFSK calibration. Receive power is fixed for
        % every guard and packet in this recording, including all pilots.
        calCfg=uac.default_config();
        nCal=floor((maximumSeconds-.14)/calCfg.symbol_s);
        rng(baseSeed+ci*100000+recordNumber,'twister');
        calBits=randi([0,1],nCal,1);
        calTx=[zeros(round(.002*fs),1);uac.modulate(calBits,calCfg);zeros(round(.132*fs),1)];
        calRx=uac.replayfilter_fast(calTx,fs,data.h,data.fs_t,double(data.fs_tau),double(data.fc));
        calRx=uac2.restore_bulk(calRx,data.V0);
        calibrationPower=sum(calRx.^2)/(nCal*calCfg.symbol_s*fs/(1-data.V0/1500));
        for guard=guards
            cfg=cfgBase; cfg.symbol_s=.008; cfg.guard_s=guard/1000;
            train=zeros(nSlots*slotSamples,1);
            bitMatrix=zeros(nBits,nSlots); layouts=cell(1,nSlots); jitters=zeros(1,nSlots); durations=zeros(1,nSlots);
            for slot=1:nSlots
                rng(baseSeed+ci*10000000+recordNumber*10000+slot,'twister');
                bitMatrix(:,slot)=randi([0,1],nBits,1);
                % Same unknown packet start and payload in all guard modes.
                jitter=3*round((.12+.22*rand)*fs/3);
                txCfg=cfg; txCfg.phase_rad=2*pi*rand();
                [packet,layout]=uac2.make_packet(bitMatrix(:,slot),txCfg);
                first=(slot-1)*slotSamples+jitter;
                assert(jitter+numel(packet)+round(.14*fs)<slotSamples,'Frame exceeds fixed measurement slot');
                train(first+(1:numel(packet)))=packet;
                layouts{slot}=layout; jitters(slot)=jitter;
                durations(slot)=numel(packet)/fs;
            end
            cleanTrain=uac.replayfilter_fast(train,fs,data.h,data.fs_t,double(data.fs_tau),double(data.fc));
            for slot=1:nSlots
                clean=cleanTrain((slot-1)*slotSamples+(1:slotSamples));
                [clean,alpha]=uac2.restore_bulk(clean,data.V0);
                layout=layouts{slot}; bits=bitMatrix(:,slot);
                trueStart=1+alpha*(jitters(slot)+oracleDelay*fs);
                oraclePayload=1+alpha*(jitters(slot)+oracleDelay*fs+layout.payload_start_sample-1);
                for snrDb=snrLevels
                    rng(baseSeed+ci*10000000+recordNumber*10000+slot*100+snrDb+50,'twister');
                    variance=calibrationPower/10^(snrDb/10);
                    rx=clean+sqrt(variance)*randn(size(clean));
                    result=uac2.receive_packet(rx,cfg,nBits,struct('compensate',true));
                    oracleBits=uac2.demodulate_scaled(rx,cfg,nBits,oraclePayload,alpha);
                    modes={'sync_compensated','sync_uncompensated','oracle_timing_scale'};
                    for mi=1:numel(modes)
                        rowCount=rowCount+1;
                        detected=result.detected;
                        if mi==1
                            estimate=result.bits_hat_compensated;
                        elseif mi==2
                            estimate=result.bits_hat_uncompensated;
                        else
                            estimate=oracleBits; detected=true;
                        end
                        if detected && numel(estimate)==nBits && all(isfinite(estimate))
                            errors=sum(estimate(:)~=bits(:)); erasures=0;
                        else
                            detected=false; errors=0; erasures=nBits;
                        end
                        rows(rowCount).channel=channel;
                        rows(rowCount).record=recordNumber;
                        if recordNumber<=5, split='train'; elseif recordNumber<=10, split='validation'; else, split='test'; end
                        rows(rowCount).split=split;
                        rows(rowCount).slot=slot;
                        rows(rowCount).slot_s=slotSeconds;
                        rows(rowCount).on_ms=8;
                        rows(rowCount).guard_ms=guard;
                        rows(rowCount).snr_db=snrDb;
                        rows(rowCount).noise_mode='fixed_noise_active_reference';
                        rows(rowCount).receiver_mode=modes{mi};
                        rows(rowCount).detected=double(detected);
                        rows(rowCount).payload_bits=nBits;
                        rows(rowCount).bit_errors=errors;
                        rows(rowCount).bit_erasures=erasures;
                        rows(rowCount).packet_error=double(~detected || errors>0);
                        % 128 ms predeclared tail is common across channels.
                        % Random observation padding/slot idling is separate.
                        rows(rowCount).occupied_s=durations(slot)+.128;
                        rows(rowCount).preamble_s=(layout.payload_start_sample-1)/fs;
                        rows(rowCount).delay_spread_ms=field_or_nan(result,'delay_spread_ms');
                        rows(rowCount).sync_score=field_or_nan(result,'sync_score');
                        rows(rowCount).velocity_hat=field_or_nan(result,'velocity_hat');
                        rows(rowCount).alpha_hat=field_or_nan(result,'alpha_hat');
                        rows(rowCount).alpha_true=alpha;
                        rows(rowCount).velocity_true=data.V0;
                        rows(rowCount).packet_start_error_ms_vs_mean_pdp=(field_or_nan(result,'preamble_start_sample')-trueStart)/fs*1000;
                        rows(rowCount).calibration_power=calibrationPower;
                        rows(rowCount).noise_variance=variance;
                    end
                end
            end
            fprintf('%s record=%d guard=%d ms done (%.1fs)\n',channel,recordNumber,guard,toc(clockRun));
        end
        records=struct2table(rows);
        writetable(records,recordsPath);
    end
end
metadata.status='COMPLETE'; metadata.stage=stage; metadata.record_ids=recordIds;
metadata.channels=channels; metadata.fs=fs; metadata.payload_bits=nBits; metadata.guards_ms=guards;
metadata.snr_db=snrLevels; metadata.seed=baseSeed; metadata.slot_s=slotSeconds;
metadata.slots_per_record=nSlots; metadata.config=cfgBase;
metadata.phase='Independent unknown uniform carrier phase per packet, shared across guard counterfactuals';
metadata.channel='Unmodified measured TVIR plus official sfetch bulk Doppler reinstatement';
metadata.receiver='HFM acquisition and pilot scale estimation; receiver has no h, true V0, packet start or payload bits';
metadata.oracle='Known transmit start, full-record mean-PDP peak delay and true V0; diagnostic, not a guaranteed upper bound';
metadata.noise='One independent 8 ms/no-guard receive-active-power reference per recording; fixed across guard';
metadata.timing='Identical slot start/jitter/data/noise across guard choices. Empty slot duration is counted separately from occupied-time efficiency.';
metadata.erasures='Missed acquisition: no invented bit decisions; all 128 bits marked erased and packet failed. Conditional BER excludes erasures.';
metadata.elapsed_s=toc(clockRun); metadata.matlab_version=version;
records=struct2table(rows);
save_json(metaFile,metadata);
fprintf('COMPLETE %s: %d rows in %.1fs\n',stage,height(records),toc(clockRun));
end

function value=field_or_nan(s,name)
if isfield(s,name) && ~isempty(s.(name)), value=s.(name); else, value=NaN; end
end

function save_json(path,data)
fid=fopen(path,'w'); assert(fid>=0,'Cannot write metadata');
cleanup=onCleanup(@()fclose(fid)); fwrite(fid,jsonencode(data),'char');
end
