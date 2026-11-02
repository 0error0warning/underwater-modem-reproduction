function [tx, layout] = make_packet(bits, cfg)
%MAKE_PACKET HFM pair, known-tone Doppler pilots, then unknown BFSK payload.
%   layout is transmitter metadata. NEVER pass this instance or the bits to
%   receive_packet; its public relative timings are derived from cfg there.
bits = bits(:);
assert(all(bits == 0 | bits == 1), 'Payload must be binary.');
layout = uac2.packet_layout(cfg, numel(bits));
tx = zeros(layout.packet_samples, 1);
phase = 0;
if isfield(cfg, 'phase_rad'), phase = cfg.phase_rad; end
amplitude = 1;
if isfield(cfg, 'amplitude'), amplitude = cfg.amplitude; end
t = (0:round(cfg.hfm_pulse_s*cfg.fs)-1)'/cfg.fs;
up = uac2.hfm_analytic(t, cfg.hfm_low_hz, cfg.hfm_high_hz, cfg.hfm_pulse_s, cfg.hfm_taper_fraction);
down = uac2.hfm_analytic(t, cfg.hfm_high_hz, cfg.hfm_low_hz, cfg.hfm_pulse_s, cfg.hfm_taper_fraction);
tx(1:numel(up)) = amplitude*real(up*exp(1i*phase));
start = round(layout.down_start_s*cfg.fs);
tx(start+(1:numel(down))) = amplitude*real(down*exp(1i*phase));
t = (0:round(cfg.pilot_s*cfg.fs)-1)'/cfg.fs;
starts = [layout.pilot0_start_s, layout.pilot1_start_s];
for k = 1:2
    start = round(starts(k)*cfg.fs);
    tx(start+(1:numel(t))) = amplitude*cos(2*pi*cfg.frequencies(k)*t+phase);
end
payload = uac.modulate(bits, cfg);
tx(layout.payload_start_sample:end) = payload;
end
