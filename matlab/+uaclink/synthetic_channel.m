function rx = synthetic_channel(tx, cfg, seed, velocityMps, sampleSnrDb)
%SYNTHETIC_CHANNEL Test harness only: delay, time scale, AWGN; NOT Watermark.
% Truth parameters remain here, outside the receiving node API.
rng(seed,'twister');
alpha=1/(1-velocityMps/cfg.sound_speed_mps);
[p,q]=rat(alpha,1e-7);
wave=resample(tx,p,q);
offset=round((0.13+0.09*rand())*cfg.fs);
rx=[zeros(offset,1);wave;zeros(round(0.12*cfg.fs),1)];
if isfinite(sampleSnrDb)
    variance=mean(wave.^2)/10^(sampleSnrDb/10);
    rx=rx+sqrt(variance)*randn(size(rx));
end
end
