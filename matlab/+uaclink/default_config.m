function cfg = default_config()
%DEFAULT_CONFIG Public PHY and bounded point-to-point packet protocol settings.
cfg = uac2.default_config();
cfg.link = struct('max_payload_bytes', 64, 'max_attempts', 3, ...
    'ack_wait_s', 2.0, 'timing_mode', 'strong');
end
