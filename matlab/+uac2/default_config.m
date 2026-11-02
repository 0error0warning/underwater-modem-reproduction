function cfg = default_config()
%DEFAULT_CONFIG Receiver-extension assumptions, separate from paper baselines.
cfg = uac.default_config();
cfg.hfm_low_hz = 10000;
cfg.hfm_high_hz = 18000;
cfg.hfm_pulse_s = 0.064;
cfg.hfm_gap_s = 0.064;
cfg.hfm_taper_fraction = 0.25;
cfg.pre_pilot_guard_s = 0.064;
cfg.pilot_s = 0.064;
cfg.pilot_gap_s = 0.016;
cfg.post_pilot_guard_s = 0.128;
cfg.pilot_trim_s = 0.016;
cfg.sync_fs = 12000;
cfg.sound_speed_mps = 1500;
cfg.velocity_range_mps = [-3, 3];
cfg.velocity_grid_step_mps = 0.05;
cfg.false_alarm_probability = 1e-3;
cfg.delay_search_radius_s = 0.064;
cfg.sinc_half_width = 8;
% The 2024 paper does not specify these HFM parameters. Its cited 2023
% synchronization paper motivates the 64 ms pulse/gap example. Frequencies,
% taper, both known-tone pilots, and the pilot guards are OUR design choices.
cfg.receiver_design = 'HFM envelope acquisition; independent known-tone Doppler pilots; fractional-time payload correction';
end
