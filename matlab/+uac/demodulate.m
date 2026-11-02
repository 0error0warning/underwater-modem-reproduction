function [bitsHat, energy] = demodulate(rx, cfg, nBits, offsetSamples)
%DEMODULATE Noncoherent BFSK by two quadrature matched filters per tone.
%   offsetSamples is the ZERO-BASED start of the first active symbol (default
%   zero). Known timing is an explicit assumption; this is not a synchronizer.
%   Only the active interval is integrated. The final guard may be omitted.
%   The energy outputs are squared quadrature projection magnitudes / nActive.
%   Orthogonal equal-energy tones are required for the standard AWGN theory.
if nargin < 4, offsetSamples = 0; end
[nActive, nGuard, time] = uac.waveform_parameters(cfg);
validateattributes(rx, {'numeric'}, {'real', 'finite', 'vector'});
validateattributes(nBits, {'numeric'}, {'real', 'finite', 'scalar', 'integer', 'nonnegative'});
validateattributes(offsetSamples, {'numeric'}, {'real', 'finite', 'scalar', 'integer', 'nonnegative'});
if nBits == 0
    bitsHat = zeros(0, 1);
    energy = zeros(0, 2);
    return;
end
stride = nActive + nGuard;
lastRequired = offsetSamples + (nBits-1)*stride + nActive;
if numel(rx) < lastRequired
    error('uac:ShortReceiveSignal', 'Need at least %d receive samples; got %d.', lastRequired, numel(rx));
end
indices = bsxfun(@plus, (1:nActive)', offsetSamples + (0:nBits-1)*stride);
rx = rx(:);
blocks = rx(indices);
angles = 2*pi*time*cfg.frequencies(:)';
inPhase = cos(angles)' * blocks;
quadrature = sin(angles)' * blocks;
energy = ((inPhase.^2 + quadrature.^2) / nActive)';
bitsHat = double(energy(:, 2) > energy(:, 1));
end
