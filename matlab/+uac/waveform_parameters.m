function [nActive, nGuard, time] = waveform_parameters(cfg)
%WAVEFORM_PARAMETERS Validate waveform settings without optional toolboxes.
required = {'fs', 'frequencies', 'symbol_s', 'guard_s'};
for k = 1:numel(required)
    if ~isfield(cfg, required{k})
        error('uac:InvalidConfig', 'Missing cfg.%s.', required{k});
    end
end
validateattributes(cfg.fs, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
validateattributes(cfg.symbol_s, {'numeric'}, {'real', 'finite', 'scalar', 'positive'});
validateattributes(cfg.guard_s, {'numeric'}, {'real', 'finite', 'scalar', 'nonnegative'});
validateattributes(cfg.frequencies, {'numeric'}, {'real', 'finite', 'vector', 'numel', 2, 'positive'});
if any(cfg.frequencies >= cfg.fs / 2) || cfg.frequencies(1) == cfg.frequencies(2)
    error('uac:InvalidConfig', 'Use two distinct tones below the Nyquist frequency.');
end
samples = cfg.fs * [cfg.symbol_s, cfg.guard_s];
if any(abs(samples - round(samples)) > 1e-8) || samples(1) < 1
    error('uac:InvalidConfig', 'Active and guard durations must contain integer sample counts.');
end
nActive = round(samples(1));
nGuard = round(samples(2));
time = (0:nActive-1)' / cfg.fs;
end
