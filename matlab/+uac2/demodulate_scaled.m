function bits = demodulate_scaled(rx,cfg,nBits,startSample,alpha)
%DEMODULATE_SCALED Known-start, known-scale diagnostic receiver.
% startSample is ONE-BASED in received samples. Public estimated receiver
% is separate; this helper is used for an explicitly labelled oracle control.
nOn=round(cfg.symbol_s*cfg.fs*alpha);
stride=(cfg.symbol_s+cfg.guard_s)*cfg.fs*alpha;
starts=round(startSample-1+(0:nBits-1)*stride);
indices=(1:nOn)'+starts;
if isempty(indices), bits=zeros(0,1); return; end
assert(min(indices(:))>=1 && max(indices(:))<=numel(rx),'uac2:ShortPacket','Insufficient received samples.');
t=(0:nOn-1)'/cfg.fs;
frequency=cfg.frequencies(:)'/alpha;
blocks=rx(indices);
a=cos(2*pi*t*frequency)'*blocks;
b=sin(2*pi*t*frequency)'*blocks;
energy=a.^2+b.^2;
bits=double(energy(2,:)>energy(1,:))';
end
