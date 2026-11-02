function tx = modulate(bits, cfg)
%MODULATE Real-valued BFSK; bit 0/1 selects cfg.frequencies(1/2).
%   The phase resets to cfg.phase_rad (default 0) at each active symbol.
%   Every symbol, including the last, is followed by guard_s seconds of zero.
%   Output is a column vector. Use chunks for large Monte Carlo experiments.
[nActive, nGuard, time] = uac.waveform_parameters(cfg);
validateattributes(bits, {'numeric', 'logical'}, {'real', 'finite', 'vector'});
if any(bits(:) ~= 0 & bits(:) ~= 1)
    error('uac:InvalidBits', 'bits must contain only zeros and ones.');
end
amplitude = 1;
phase = 0;
if isfield(cfg, 'amplitude'), amplitude = cfg.amplitude; end
if isfield(cfg, 'phase_rad'), phase = cfg.phase_rad; end
validateattributes(amplitude, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
validateattributes(phase, {'numeric'}, {'real', 'finite', 'scalar'});
tones = amplitude * cos(2*pi*time*cfg.frequencies(:)' + phase);
symbols = zeros(nActive + nGuard, numel(bits));
symbols(1:nActive, :) = tones(:, double(bits(:)') + 1);
tx = symbols(:);
end
