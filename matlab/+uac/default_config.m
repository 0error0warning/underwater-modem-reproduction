function cfg = default_config()
%DEFAULT_CONFIG Explicit assumptions for the independent BFSK implementation.
%   The paper specifies 250 Hz spacing and a 14 kHz channel centre, but does
%   not specify the exact two BFSK tones. The symmetric tones below are our
%   documented choice, not an exact frequency mapping reported by the paper.
%   symbol_s means ACTIVE transmission duration; guard_s is additional silence.
cfg.fs = 48000;
cfg.frequencies = [13875, 14125];
cfg.symbol_s = 0.008;
cfg.guard_s = 0;
cfg.amplitude = 1;
cfg.phase_rad = 0;
end
