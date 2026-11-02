function run_nof1_diagnostic_v3(stage)
%RUN_NOF1_DIAGNOSTIC_V3 Paired DEVELOPMENT-only baseline mechanism diagnosis.
% No frozen production function or adaptive rule is modified by this run.
if nargin < 1, stage = 'full'; end
root = fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out = fullfile(root,'results','v3','nof1_diagnostic');
recordIds = 1:3;
if strcmp(stage,'smoke'), out=fullfile(out,'smoke'); recordIds=1;
elseif strcmp(stage,'confirmation'), out=fullfile(out,'confirmation'); recordIds=4:5;
elseif ~strcmp(stage,'full'), error('Use full, smoke or confirmation'); end
if ~exist(out,'dir'), mkdir(out); end
base = uac2.default_config(); fs=base.fs; nBits=64; nSlots=5;
slotS=6.4; slotN=round(slotS*fs); seed=2026110205;
% Fixed beforehand: four duration/guard controls; two additional tone pairs
% only at 4ms/no guard; each with reset and continuous oscillator phase.
design = [4 0 13875 14125;4 4 13875 14125;8 0 13875 14125;8 8 13875 14125; ...
          4 0 13750 14000;4 0 14000 14250];
phases = {'reset','continuous','payload_reset_zero'}; noiseNames={'clean','paired_sample20'};
meta=struct('status','RUNNING','stage',stage,'record_ids',recordIds,'channel','NOF1', ...
    'design_columns',{{'on_ms','guard_ms','f0_hz','f1_hz'}},'design',design, ...
    'phase_modes',{phases},'noise_modes',{noiseNames},'seed',seed, ...
    'scope','Development 001-003 only; mechanistic paired sensitivity, not untouched final evaluation', ...
    'phase_definition','continuous oscillator phase accumulates current tone over active+guard stride; reset is unchanged V2', ...
    'noise_definition','20dB variance from baseline 13875/14125 RESET frame at same on/guard/record/slot; unchanged across phase/frequency variants', ...
    'timing_sweep','4/8ms no guard, paired_sample20 only: mean full-record PDP peak +/-4 ms in 0.25 ms steps, true bulk scale; diagnostic only, never deployed/selected by payload BER');
if strcmp(stage,'confirmation')
    meta.executed_design_indices=[1,3]; meta.executed_phase_modes={'reset'};
    meta.timing_sweep='none';
    meta.scope='Previously used DEVELOPMENT records004-005, fixed first-versus-main-peak comparison after001-003 diagnosis; no tuning and not an untouched final test';
end
writeJson(fullfile(out,'metadata.json'),meta);
reference=readtable(fullfile(root,'results','v2','paper_alignment','per_packet.csv'));
rows={}; sweeps={}; events={}; clockRun=tic;
for record=recordIds
    channelPath=fullfile(root,'vendor','Watermark','input','channels','NOF1','mat',sprintf('NOF1_%03d.mat',record));
    data=load(channelPath); meanPdp=mean(abs(data.h).^2,1); [~,peak]=max(meanPdp);
    meanDelay=(peak-1)/double(data.fs_tau);
    referenceVariance=nan(4,nSlots);
    for di=1:size(design,1)
        if strcmp(stage,'confirmation') && ~ismember(di,[1,3]), continue; end
        cfg=base; cfg.symbol_s=design(di,1)/1000; cfg.guard_s=design(di,2)/1000;
        cfg.frequencies=design(di,3:4);
        for phaseIndex=1:numel(phases)
            if strcmp(stage,'confirmation') && phaseIndex~=1, continue; end
            if phaseIndex==3 && di~=1, continue; end
            tx=zeros(nSlots*slotN,1); payloads=cell(nSlots,1); jitters=zeros(nSlots,1);
            carrierPhases=zeros(nSlots,1);
            for slot=1:nSlots
                rng(seed+10000000+record*10000+slot,'twister');
                bits=randi([0,1],nBits,1); payloads{slot}=bits;
                jitter=3*round((.12+.22*rand())*fs/3); jitters(slot)=jitter;
                txCfg=cfg; txCfg.phase_rad=2*pi*rand(); carrierPhases(slot)=txCfg.phase_rad;
                [packet,layout]=uac2.make_packet(bits,txCfg);
                if phaseIndex==2, packet(layout.payload_start_sample:end)=continuousPayload(bits,txCfg); end
                if phaseIndex==3
                    payloadCfg=txCfg; payloadCfg.phase_rad=0;
                    packet(layout.payload_start_sample:end)=uac.modulate(bits,payloadCfg);
                end
                first=(slot-1)*slotN+jitter;
                tx(first+(1:numel(packet)))=packet;
            end
            cleanTrain=uac.replayfilter_fast(tx,fs,data.h,data.fs_t,double(data.fs_tau),double(data.fc));
            for slot=1:nSlots
                bits=payloads{slot};
                [clean,alphaTrue]=uac2.restore_bulk(cleanTrain((slot-1)*slotN+(1:slotN)),data.V0);
                energy=sum(clean.^2); powerN=ceil((layout.packet_duration_s+.128)*alphaTrue*fs);
                refIndex=di; if di>4, refIndex=1; end
                if phaseIndex==1 && di<=4, referenceVariance(refIndex,slot)=energy/powerN/100; end
                variance=referenceVariance(refIndex,slot); assert(isfinite(variance));
                rng(seed+10000000+record*10000+slot*100+20+50,'twister'); unitNoise=randn(size(clean));
                for ni=1:2
                    appliedVariance=variance*(ni==2); rx=clean+sqrt(appliedVariance)*unitNoise;
                    result=uac2.receive_packet(rx,cfg,nBits,struct('compensate',true));
                    assert(result.detected,'Unexpected capture failure: diagnose explicitly before extending.');
                    first=result.preamble_start_sample;
                    strong=1+(result.diagnostic.delay.strongest_arrival_sample_at_sync_rate-1)*fs/cfg.sync_fs;
                    pdpStart=1+alphaTrue*(jitters(slot)+meanDelay*fs);
                    methods={'first_est','strong_est','first_true_scale','strong_true_scale','pdp_true','pdp_est_scale','first_est_fft_nearest'};
                    starts=[first,strong,first,strong,pdpStart,pdpStart,first];
                    scales=[result.alpha_hat,result.alpha_hat,alphaTrue,alphaTrue,alphaTrue,result.alpha_hat,result.alpha_hat];
                    for mi=1:numel(methods)
                        detector='exact'; if mi==7, detector='fft_nearest'; end
                        hat=decodeDiagnostic(rx,cfg,nBits,starts(mi),scales(mi),layout,detector);
                        if mi==1, assert(isequal(hat,result.bits_hat_compensated(:)),'Diagnostic decoder differs from frozen receiver'); end
                        r=baseRow(record,slot,di,phaseIndex,ni,design,phases,noiseNames);
                        r.method=methods{mi}; r.payload_bits=nBits; r.bit_errors=sum(hat~=bits);
                        r.packet_error=double(r.bit_errors>0); r.preamble_start_sample=starts(mi);
                        r.alpha_used=scales(mi); r.offset_from_pdp_ms=(starts(mi)-pdpStart)/fs*1000;
                        r.detector=detector; rows{end+1}=r; %#ok<AGROW>
                        if record<=3 && phaseIndex==1 && di<=4 && ni==2 && mi==1
                            old=reference(strcmp(reference.channel,'NOF1') & reference.record==record & reference.slot==slot & ...
                                reference.on_ms==cfg.symbol_s*1000 & reference.guard_ms==cfg.guard_s*1000 & ...
                                reference.snr_label_db==20 & strcmp(reference.noise_mode,'sample_snr_packet'),:);
                            assert(height(old)==1 && old.bit_errors==r.bit_errors,'Original V2 bit-count regression failed');
                            assert(abs(old.noise_variance-variance)<1e-9*max(1,variance),'Original noise variance regression failed');
                            assert(abs(old.carrier_phase_rad-carrierPhases(slot))<1e-12,'Original initial phase regression failed');
                        end
                    end
                    e=baseRow(record,slot,di,phaseIndex,ni,design,phases,noiseNames);
                    e.carrier_phase_rad=carrierPhases(slot); e.injected_start_s=jitters(slot)/fs;
                    e.payload_initial_phase_rad=carrierPhases(slot); if phaseIndex==3, e.payload_initial_phase_rad=0; end
                    e.alpha_true=alphaTrue; e.alpha_hat=result.alpha_hat; e.velocity_true=data.V0;
                    e.velocity_hat=result.velocity_hat; e.first_minus_pdp_ms=(first-pdpStart)/fs*1000;
                    e.strong_minus_pdp_ms=(strong-pdpStart)/fs*1000; e.strong_minus_first_ms=(strong-first)/fs*1000;
                    e.mean_pdp_delay_ms=meanDelay*1000; e.delay_spread_ms=result.delay_spread_ms;
                    e.received_energy=energy; e.noise_variance=appliedVariance;
                    e.sample_snr_actual_db=10*log10(energy/powerN/max(appliedVariance,realmin));
                    e.detected=double(result.detected); e.original_decoder_match=1; events{end+1}=e; %#ok<AGROW>
                    if ~strcmp(stage,'confirmation') && ismember(di,[1,3]) && phaseIndex==1 && ni==2
                        for shiftMs=-4:.25:4
                            hat=decodeDiagnostic(rx,cfg,nBits,pdpStart+shiftMs/1000*fs,alphaTrue,layout,'exact');
                            r=baseRow(record,slot,di,phaseIndex,ni,design,phases,noiseNames);
                            r.shift_ms=shiftMs; r.payload_bits=nBits; r.bit_errors=sum(hat~=bits);
                            sweeps{end+1}=r; %#ok<AGROW>
                        end
                    end
                end
            end
            writetable(struct2table(vertcat(rows{:})),fullfile(out,'decoder_results.csv'));
            writetable(struct2table(vertcat(events{:})),fullfile(out,'events.csv'));
            if ~isempty(sweeps), writetable(struct2table(vertcat(sweeps{:})),fullfile(out,'timing_sweep.csv')); end
            fprintf('NOF1 %03d design=%d phase=%s done; %.1fs\n',record,di,phases{phaseIndex},toc(clockRun));
        end
    end
end
meta.status='COMPLETE'; meta.elapsed_s=toc(clockRun); meta.event_rows=numel(events);
meta.decoder_rows=numel(rows); meta.timing_sweep_rows=numel(sweeps); meta.matlab_version=version;
writeJson(fullfile(out,'metadata.json'),meta);
fprintf('NOF1 DIAGNOSTIC COMPLETE: %d events, %.1fs\n',numel(events),toc(clockRun));
end

function tx=continuousPayload(bits,cfg)
[nActive,nGuard,t]=uac.waveform_parameters(cfg); stride=nActive+nGuard;
freq=cfg.frequencies(double(bits(:)')+1);
% Accumulate with exact stride; phase during guard uses the current tone.
phase=cfg.phase_rad+2*pi*[0,cumsum(freq(1:end-1)*stride/cfg.fs)];
symbols=zeros(stride,numel(bits));
symbols(1:nActive,:)=cfg.amplitude*cos(2*pi*t*freq+phase);
tx=symbols(:);
end

function bits=decodeDiagnostic(rx,cfg,nBits,start,alpha,layout,detector)
[nActive,nGuard,time]=uac.waveform_parameters(cfg);
offsets=bsxfun(@plus,(0:nActive-1)',(0:nBits-1)*(nActive+nGuard));
positions=start+alpha*(layout.payload_start_s*cfg.fs+offsets);
assert(all(positions(:)>=1 & positions(:)<=numel(rx)));
samples=uac2.sinc_sample(rx,positions,cfg.sinc_half_width);
if strcmp(detector,'exact')
    angles=2*pi*time*cfg.frequencies(:)';
    energy=(cos(angles)'*samples).^2+(sin(angles)'*samples).^2;
else
    bins=round(cfg.frequencies*nActive/cfg.fs); spectrum=fft(samples,nActive,1);
    energy=abs(spectrum(bins+1,:)).^2;
end
bits=double(energy(2,:)'>energy(1,:)');
end

function r=baseRow(record,slot,di,pi,ni,design,phases,noiseNames)
r=struct('channel','NOF1','record',record,'slot',slot,'design',di, ...
    'on_ms',design(di,1),'guard_ms',design(di,2),'f0_hz',design(di,3),'f1_hz',design(di,4), ...
    'phase_mode',phases{pi},'noise_mode',noiseNames{ni});
end

function writeJson(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(value),'char');
end
