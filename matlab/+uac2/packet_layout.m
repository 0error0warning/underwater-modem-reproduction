function layout = packet_layout(cfg, nBits)
%PACKET_LAYOUT Public packet timing; the receiver may derive this from cfg.
%   Sample indices are MATLAB one-based. No received arrival time is known.
validateattributes(nBits, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
[nActive, nGuard] = uac.waveform_parameters(cfg);
layout.up_start_s = 0;
layout.down_start_s = cfg.hfm_pulse_s + cfg.hfm_gap_s;
layout.pilot0_start_s = layout.down_start_s + cfg.hfm_pulse_s + cfg.pre_pilot_guard_s;
layout.pilot1_start_s = layout.pilot0_start_s + cfg.pilot_s + cfg.pilot_gap_s;
layout.payload_start_s = layout.pilot1_start_s + cfg.pilot_s + cfg.post_pilot_guard_s;
times = [layout.down_start_s, layout.pilot0_start_s, layout.pilot1_start_s, ...
    layout.payload_start_s, cfg.hfm_pulse_s, cfg.pilot_s] * cfg.fs;
assert(all(abs(times-round(times)) < 1e-8), 'Packet timing must use integer transmit samples.');
layout.payload_start_sample = round(layout.payload_start_s*cfg.fs) + 1;
layout.active_samples = nActive;
layout.guard_samples = nGuard;
layout.n_bits = nBits;
layout.packet_samples = layout.payload_start_sample-1 + nBits*(nActive+nGuard);
layout.packet_duration_s = layout.packet_samples/cfg.fs;
end
