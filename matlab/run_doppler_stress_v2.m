function report=run_doppler_stress_v2()
%RUN_DOPPLER_STRESS_V2 Predeclared long-packet time-scale stress, not TVIR data.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'matlab'));
out=fullfile(root,'results','receiver_v2');
cfg=uac2.default_config(); cfg.symbol_s=.008; cfg.guard_s=0;
velocities=[-1.5,0,1.5]; nBits=512; rng(2026110204,'twister'); rows=struct([]);
for k=1:numel(velocities)
    bits=randi([0,1],nBits,1); txCfg=cfg; txCfg.phase_rad=2*pi*rand();
    tx=uac2.make_packet(bits,txCfg);
    [clean,alpha]=uac2.restore_bulk(tx,velocities(k));
    delay=round((.15+.3*rand())*cfg.fs);
    rx=[zeros(delay,1);clean;zeros(round(.15*cfg.fs),1)];
    % Sample-domain active power SNR=20dB, not Eb/N0.
    rx=rx+sqrt(.5/100)*randn(size(rx));
    r=uac2.receive_packet(rx,cfg,nBits,struct('compensate',true));
    rows(k).velocity_true=velocities(k); rows(k).alpha_true=alpha;
    rows(k).velocity_hat=r.velocity_hat; rows(k).detected=r.detected;
    rows(k).payload_bits=nBits;
    if r.detected
        rows(k).errors_corrected=sum(r.bits_hat_compensated~=bits);
        rows(k).errors_uncorrected=sum(r.bits_hat_uncompensated~=bits);
    else
        rows(k).errors_corrected=NaN; rows(k).errors_uncorrected=NaN;
    end
end
writetable(struct2table(rows),fullfile(out,'long_packet_doppler.csv'));
passed=all([rows.detected]) && all([rows.errors_corrected]==0);
report.status='FAIL'; if passed, report.status='PASS'; end
report.scope='512-bit synthetic single-path packets with unknown phase/start; 20 dB active-sample SNR; no real-channel performance claim';
report.results=rows; report.seed=2026110204;
fid=fopen(fullfile(out,'long_packet_doppler.json'),'w'); fwrite(fid,jsonencode(report),'char'); fclose(fid);
disp(struct2table(rows));
assert(all([rows.detected]),'Acquisition failed on prescribed long-packet test');
assert(all([rows.errors_corrected]==0),'Nonzero errors after correction on prescribed high-SNR single-path test');
end
