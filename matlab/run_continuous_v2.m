function records=run_continuous_v2()
%RUN_CONTINUOUS_V2 Frozen policies in causal, continuously scheduled TVIR runs.
% The next packet really moves in channel time when the guard changes.
% Feedback is an ideal one-packet-delayed model with explicit time cost.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','v2');
policyPath=fullfile(out,'adaptive','frozen_policy.json');
assert(exist(policyPath,'file')==2,'Freeze policy before opening holdout data');
policy=jsondecode(fileread(policyPath));
assert(strcmp(policy.status,'FROZEN_BEFORE_FINAL_TEST'),'Invalid policy state');
kernelProof=jsondecode(fileread(fullfile(root,'results','replay_at_validation.json')));
assert(strcmp(kernelProof.status,'PASS'),'Offset replay must be verified');
save_json(fullfile(out,'metadata_continuous.json'),struct('status','RUNNING'));
names={'best_fixed','adaptive_stump'};
models={policy.best_fixed,policy.adaptive_stump};
if isempty(models{2}), names=names(1); models=models(1); end
cfgBase=uac2.default_config(); fs=cfgBase.fs;
channels={'NOF1','NCS1'}; recordIds=31:40;
nBits=128; horizon=32; tail=.128; control=policy.base_control_s;
snrDb=20; seed=2026110203; rows=struct([]); rowCount=0; timer=tic;
for ci=1:2
    channel=channels{ci};
    for number=recordIds
        data=load(fullfile(root,'vendor','Watermark','input','channels',channel,'mat',sprintf('%s_%03d.mat',channel,number)));
        % Same calibration construction as paired runs, without policy data.
        maximumSeconds=(size(data.h,1)-1)/data.fs_t;
        calCfg=uac.default_config(); nCal=floor((maximumSeconds-.14)/.008);
        rng(2026110202+ci*100000+number,'twister');
        calTx=[zeros(round(.002*fs),1);uac.modulate(randi([0,1],nCal,1),calCfg);zeros(round(.132*fs),1)];
        calRx=uac.replayfilter_fast(calTx,fs,data.h,data.fs_t,double(data.fs_tau),double(data.fc));
        calRx=uac2.restore_bulk(calRx,data.V0);
        calibrationPower=sum(calRx.^2)/(nCal*.008*fs/(1-data.V0/1500));
        rng(seed+ci*100000+number,'twister');
        % A common noise realization on the absolute observation timeline.
        noise=sqrt(calibrationPower/10^(snrDb/10))*randn(round((horizon+4)*fs),1);
        for pi=1:numel(models)
            model=models{pi}; previous=[]; cursor=.05; packetId=0;
            while true
                guard=choose_guard(model,previous);
                cfg=cfgBase; cfg.symbol_s=.008; cfg.guard_s=guard/1000;
                layout=uac2.packet_layout(cfg,nBits);
                duration=layout.packet_duration_s;
                start=ceil(cursor*double(data.fs_tau))/double(data.fs_tau);
                % The sender uses its known frame duration, not true alpha.
                % Tail and feedback budget cover receive/filter clearance.
                next=start+duration+tail+control;
                if next>horizon, break; end
                packetId=packetId+1;
                rng(seed+ci*10000000+number*10000+packetId,'twister');
                bits=randi([0,1],nBits,1); cfg.phase_rad=2*pi*rand();
                padding=round((.12+.22*rand())*fs);
                tx=uac2.make_packet(bits,cfg);
                clean=uac2.replay_at([tx;zeros(round(tail*fs),1)],fs,data.h,data.fs_t,double(data.fs_tau),double(data.fc),start);
                [clean,alpha]=uac2.restore_bulk(clean,data.V0);
                signal=[zeros(padding,1);clean;zeros(round(.03*fs),1)];
                observationStart=round((start-padding/fs+1)*fs)+1;
                assert(observationStart>=1 && observationStart+numel(signal)-1<=numel(noise));
                rx=signal+noise(observationStart+(0:numel(signal)-1));
                receiverCfg=cfgBase; receiverCfg.symbol_s=.008; receiverCfg.guard_s=guard/1000;
                result=uac2.receive_packet(rx,receiverCfg,nBits,struct('compensate',true));
                detected=result.detected && numel(result.bits_hat)==nBits && all(isfinite(result.bits_hat));
                if detected, errors=sum(result.bits_hat(:)~=bits); erasures=0;
                else, errors=0; erasures=nBits; end
                rowCount=rowCount+1;
                rows(rowCount).channel=channel; rows(rowCount).record=number;
                rows(rowCount).policy=names{pi}; rows(rowCount).packet=packetId;
                rows(rowCount).start_s=start; rows(rowCount).next_start_s=next;
                rows(rowCount).horizon_s=horizon; rows(rowCount).guard_ms=guard;
                rows(rowCount).snr_db=snrDb; rows(rowCount).detected=double(detected);
                rows(rowCount).payload_bits=nBits; rows(rowCount).bit_errors=errors;
                rows(rowCount).bit_erasures=erasures; rows(rowCount).packet_error=double(~detected || errors>0);
                rows(rowCount).occupied_s=duration+tail; rows(rowCount).control_s=control;
                rows(rowCount).delay_spread_ms=result.delay_spread_ms;
                rows(rowCount).sync_score=result.sync_score;
                rows(rowCount).velocity_hat=result.velocity_hat;
                rows(rowCount).alpha_hat=result.alpha_hat; rows(rowCount).alpha_true=alpha;
                rows(rowCount).noise_variance=calibrationPower/10^(snrDb/10);
                previous=result;
                previous.detected=detected;
                cursor=next;
            end
            fprintf('continuous %s record=%d policy=%s %d packets (%.1fs)\n',channel,number,names{pi},packetId,toc(timer));
        end
        records=struct2table(rows); writetable(records,fullfile(out,'records_continuous.csv'));
    end
end
metadata.status='COMPLETE'; metadata.record_ids=recordIds; metadata.horizon_s=horizon;
metadata.policy_file=policyPath; metadata.policy_frozen_content=policy;
metadata.channel='Verified absolute-time TVIR replay; true bulk V0 restored; no channel truth enters receiver or policy';
metadata.schedule='Each next transmitted frame starts after selected frame duration + 128 ms common tail + explicit feedback time; channel time advances on that trajectory';
metadata.control_s=control; metadata.snr_db=snrDb; metadata.seed=seed;
metadata.time_denominator='Full 32 s horizon including startup, ACK/feedback, grid rounding and unused end time';
metadata.feedback='Ideal successful control/ACK within stated budget, one-packet delay; missed data acquisition falls back. Return-link acoustics and propagation are not simulated.';
metadata.phase='Unknown random per-packet phase; common bit sequences and continuous-time noise across policies';
metadata.elapsed_s=toc(timer); metadata.matlab_version=version;
save_json(fullfile(out,'metadata_continuous.json'),metadata);
fprintf('COMPLETE continuous %d rows in %.1fs\n',height(records),toc(timer));
end

function guard=choose_guard(model,previous)
if strcmp(model.kind,'fixed'), guard=model.guard_ms; return; end
guard=model.fallback_guard_ms;
if isempty(previous) || ~previous.detected, return; end
if strcmp(model.feature,'abs_velocity_hat'), value=abs(previous.velocity_hat);
elseif isfield(previous,model.feature), value=previous.(model.feature);
else, return; end
if ~isfinite(value), return; end
if value<=model.threshold, guard=model.left_guard_ms; else, guard=model.right_guard_ms; end
end

function save_json(path,data)
fid=fopen(path,'w'); assert(fid>=0); cleanup=onCleanup(@()fclose(fid));
fwrite(fid,jsonencode(data),'char');
end
