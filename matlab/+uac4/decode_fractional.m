function [bits, available] = decode_fractional(rx, cfg, nBits, preambleStart, alpha)
%DECODE_FRACTIONAL Exact V2 detector at a supplied, publicly estimated start.
% This is the unchanged fractional sample/exact tone-energy logic from
% uac2.receive_packet. It does not estimate a start or consult payload truth.
validateattributes(nBits, {'numeric'}, {'scalar','integer','nonnegative'});
validateattributes(preambleStart, {'numeric'}, {'scalar','real','finite'});
validateattributes(alpha, {'numeric'}, {'scalar','real','finite','positive'});
layout = uac2.packet_layout(cfg, nBits);
[nActive, nGuard, time] = uac.waveform_parameters(cfg);
offsets = bsxfun(@plus, (0:nActive-1)', (0:nBits-1)*(nActive+nGuard));
positions = preambleStart+alpha*(layout.payload_start_s*cfg.fs+offsets);
available = all(positions(:) >= 1 & positions(:) <= numel(rx));
if ~available, bits = zeros(0, 1); return; end
samples = uac2.sinc_sample(rx, positions, cfg.sinc_half_width);
angles = 2*pi*time*cfg.frequencies(:)';
energy = (cos(angles)'*samples).^2+(sin(angles)'*samples).^2;
bits = double(energy(2, :)' > energy(1, :)');
end
