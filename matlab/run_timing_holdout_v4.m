function run_timing_holdout_v4(stage)
%RUN_TIMING_HOLDOUT_V4 Predeclared paired timing test; no adaptive selection.
% Requires a JVM for portable SHA-256. 'development' reads only record001;
% 'final' verifies the frozen manifest before loading any record011-020.
if nargin < 1, stage = 'development'; end
assert(any(strcmp(stage, {'development','final'})), 'Use development or final.');
assert(usejava('jvm'), 'Run MATLAB with its JVM (do not pass -nojvm).');
root = fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out = fullfile(root,'results','v4','timing_holdout');
channels = {'NOF1','NCS1'}; records = 11:20;
freezePath = fullfile(out,'freeze.json');
freezeHash = '';
if strcmp(stage,'development')
    records = 1; out = fullfile(out,'development');
else
    assert(exist(freezePath,'file')==2, 'Missing preregistration freeze.json.');
    checkFreeze(root,freezePath); freezeHash=hashFile(freezePath);
    assert(exist(fullfile(out,'metadata.json'),'file')~=2, ...
        'uac4:NoOverwrite','Final metadata exists; preserve it before an explicitly documented rerun.');
end
if ~exist(out,'dir'), mkdir(out); end
cfgBase=uac2.default_config(); fs=cfgBase.fs; nBits=64; nSlots=5;
slotS=6.4; slotN=round(slotS*fs); seed=2026110205;
design=[4 0;4 4;8 0;8 8]; methods={'first','strong'};
meta=struct('status','RUNNING','stage',stage,'channels',{channels}, ...
    'record_ids',records,'reserved_unread_record_ids',21:30, ...
    'design_columns',{{'on_ms','guard_ms'}},'design',design,'payload_bits',nBits, ...
    'slots_per_record',nSlots,'slot_duration_s',slotS,'horizon_s',nSlots*slotS, ...
    'seed',seed,'methods',{methods},'noise_mode','sample_snr_packet','snr_label_db',20, ...
    'phase_mode','reset','primary_condition','NOF1; on_ms=4; guard_ms=0', ...
    'freeze_sha256',freezeHash,'rx_hash_encoding','SHA256 of MATLAB double column bytes; little-endian IEEE754', ...
    'scope','Previously unused records of existing channel sequences; not independent environments', ...
    'started_utc',utcNow(),'matlab_version',version);
writeJson(fullfile(out,'metadata.json'),meta);
rows={}; sources={}; clockRun=tic;
if strcmp(stage,'development')
    old=readtable(fullfile(root,'results','v3','nof1_diagnostic','decoder_results.csv'));
    checkSynthetic(cfgBase);
end
for channelIndex=1:numel(channels)
    channel=channels{channelIndex};
    for record=records
        relative=sprintf('vendor/Watermark/input/channels/%s/mat/%s_%03d.mat',channel,channel,record);
        channelPath=fullfile(root,strrep(relative,'/',filesep));
        sourceHash=hashFile(channelPath);
        sources{end+1}=struct('path',relative,'sha256',sourceHash); %#ok<AGROW>
        data=load(channelPath);
        for designIndex=1:size(design,1)
            cfg=cfgBase; cfg.symbol_s=design(designIndex,1)/1000; cfg.guard_s=design(designIndex,2)/1000;
            tx=zeros(nSlots*slotN,1); payloads=cell(nSlots,1);
            jitters=zeros(nSlots,1); carrierPhases=zeros(nSlots,1);
            for slot=1:nSlots
                rng(seed+channelIndex*10000000+record*10000+slot,'twister');
                bits=randi([0,1],nBits,1); payloads{slot}=bits;
                jitter=3*round((.12+.22*rand())*fs/3); jitters(slot)=jitter;
                txCfg=cfg; txCfg.phase_rad=2*pi*rand(); carrierPhases(slot)=txCfg.phase_rad;
                [packet,layout]=uac2.make_packet(bits,txCfg);
                first=(slot-1)*slotN+jitter;
                tx(first+(1:numel(packet)))=packet;
            end
            cleanTrain=uac.replayfilter_fast(tx,fs,data.h,data.fs_t,double(data.fs_tau),double(data.fc));
            for slot=1:nSlots
                bits=payloads{slot};
                [clean,alphaTrue]=uac2.restore_bulk(cleanTrain((slot-1)*slotN+(1:slotN)),data.V0);
                energy=sum(clean.^2); powerN=ceil((layout.packet_duration_s+.128)*alphaTrue*fs);
                variance=energy/powerN/100;
                rng(seed+channelIndex*10000000+record*10000+slot*100+20+50,'twister');
                rx=clean+sqrt(variance)*randn(size(clean));
                rxHash=hashBytes(typecast(double(rx(:)),'uint8'));
                for methodIndex=1:2
                    method=methods{methodIndex}; result=uac4.receive_timing(rx,cfg,nBits,method);
                    if result.detected && numel(result.bits_hat)==nBits
                        errors=sum(result.bits_hat(:)~=bits); erasures=0; decoded=nBits;
                    else
                        errors=0; erasures=nBits; decoded=0;
                    end
                    r=struct('channel',channel,'record',record,'slot',slot,'design',designIndex, ...
                        'on_ms',design(designIndex,1),'guard_ms',design(designIndex,2), ...
                        'method',method,'payload_bits',nBits,'decoded_bits',decoded, ...
                        'bit_errors',errors,'bit_erasures',erasures,'packet_error',double(errors+erasures>0), ...
                        'detected',double(result.detected),'reason',result.diagnostic.reason, ...
                        'preamble_start_sample',result.preamble_start_sample,'alpha_hat',result.alpha_hat, ...
                        'sync_score',result.sync_score,'delay_spread_ms',result.delay_spread_ms, ...
                        'carrier_phase_rad',carrierPhases(slot),'injected_start_s',jitters(slot)/fs, ...
                        'alpha_true_simulation_only',alphaTrue,'noise_variance',variance, ...
                        'received_samples',numel(rx),'rx_sha256',rxHash,'source_sha256',sourceHash);
                    rows{end+1}=r; %#ok<AGROW>
                    if strcmp(stage,'development') && strcmp(channel,'NOF1')
                        oldMethod='first_est'; if methodIndex==2, oldMethod='strong_est'; end
                        reference=old(old.record==record & old.slot==slot & old.design==designIndex & ...
                            strcmp(old.phase_mode,'reset') & strcmp(old.noise_mode,'paired_sample20') & ...
                            strcmp(old.method,oldMethod),:);
                        assert(height(reference)==1 && erasures==0 && errors==reference.bit_errors, ...
                            'Development counts differ from frozen V3 diagnostics.');
                        assert(abs(result.preamble_start_sample-reference.preamble_start_sample)<1e-9);
                        assert(abs(result.alpha_hat-reference.alpha_used)<1e-12);
                    end
                end
            end
            writetable(struct2table(vertcat(rows{:})),fullfile(out,'per_packet.csv'));
            fprintf('%s %03d on/guard=%g/%g ms complete; %.1f s\n',channel,record,design(designIndex,:),toc(clockRun));
        end
    end
end
if strcmp(stage,'final'), checkFreeze(root,freezePath); end
meta.status='COMPLETE'; meta.completed_utc=utcNow(); meta.elapsed_s=toc(clockRun);
meta.waveform_count=numel(rows)/2; meta.method_rows=numel(rows); meta.sources=vertcat(sources{:});
meta.frozen_sources_unchanged=strcmp(stage,'final');
meta.development_v3_checks=0;
if strcmp(stage,'development'), meta.development_v3_checks=40; end
writeJson(fullfile(out,'metadata.json'),meta);
fprintf('TIMING %s COMPLETE: %d method rows, %.1f s\n',stage,numel(rows),toc(clockRun));
end

function checkSynthetic(cfg)
rng(4096,'twister'); bits=randi([0,1],64,1);
cfg.phase_rad=1.23; [packet,~]=uac2.make_packet(bits,cfg);
rx=[zeros(round(.2*cfg.fs),1);packet;zeros(round(.2*cfg.fs),1)];
legacy=uac2.receive_packet(rx,cfg,numel(bits));
first=uac4.receive_timing(rx,cfg,numel(bits),'first');
strong=uac4.receive_timing(rx,cfg,numel(bits),'strong');
assert(isequal(legacy.bits_hat,first.bits_hat) && legacy.detected==first.detected);
assert(first.detected && strong.detected && isequal(strong.bits_hat,bits));
[decoded,available]=uac4.decode_fractional(rx,cfg,numel(bits),legacy.preamble_start_sample,legacy.alpha_hat);
assert(available && isequal(decoded,legacy.bits_hat));
for mode={'first','strong'}
    z=uac4.receive_timing(zeros(size(rx)),cfg,numel(bits),mode{1});
    assert(~z.detected && isempty(z.bits_hat),'Silence must not produce a payload.');
    short=uac4.receive_timing(zeros(32,1),cfg,numel(bits),mode{1});
    assert(~short.detected && isempty(short.bits_hat),'Truncation must not produce a payload.');
end
[~,available]=uac4.decode_fractional(rx,cfg,numel(bits),numel(rx),1);
assert(~available,'Out-of-bounds payload must not be silently decoded.');
end

function checkFreeze(root,path)
freeze=jsondecode(fileread(path));
assert(isfield(freeze,'files') && ~isempty(freeze.files),'Invalid freeze manifest.');
for index=1:numel(freeze.files)
    entry=freeze.files(index);
    assert(strcmp(hashFile(fullfile(root,strrep(entry.path,'/',filesep))),entry.sha256), ...
        'uac4:FreezeMismatch','Frozen source changed: %s',entry.path);
end
end

function hash=hashFile(path)
fid=fopen(path,'rb'); assert(fid>=0,'Cannot open input: %s',path);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
hash=hashBytes(fread(fid,Inf,'*uint8'));
end

function hash=hashBytes(bytes)
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(typecast(uint8(bytes(:)),'int8')); raw=typecast(digest.digest(),'uint8');
hash=lower(reshape(dec2hex(raw,2)',1,[]));
end

function value=utcNow()
value=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
end

function writeJson(path,value)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fwrite(fid,jsonencode(value),'char');
end
