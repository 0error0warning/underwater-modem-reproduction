function result = receive_timing(rx, cfg, nBits, mode)
%RECEIVE_TIMING Frozen V2 acquisition with a predeclared timing rule.
% Inputs contain no transmitted bits, channel response, velocity, or start.
% 'first' is the original receiver; 'strong' replaces only the integration
% start by the already-computed, scale-refined HFM joint-statistic main peak.
if nargin < 4, mode = 'strong'; end
assert(ischar(mode) && any(strcmp(mode, {'first','strong'})), ...
    'uac4:TimingMode', 'Timing mode must be first or strong.');
result = uac2.receive_packet(rx, cfg, nBits, struct('compensate', true));
result.diagnostic.timing_mode = mode;
result.diagnostic.original_first_start_sample = result.preamble_start_sample;
if strcmp(mode, 'first'), return; end
% A payload may be out of bounds even after acquisition. Check the HFM
% diagnostic, not result.detected, so each rule gets its own bounds check.
if ~isfield(result.diagnostic, 'delay') || ~isfinite(result.alpha_hat)
    return;
end
strong = result.diagnostic.delay.strongest_arrival_sample_at_sync_rate;
start = 1+(strong-1)*cfg.fs/cfg.sync_fs;
layout = uac2.packet_layout(cfg, nBits);
result.preamble_start_sample = start;
result.payload_start_sample = start+result.alpha_hat*layout.payload_start_s*cfg.fs;
result.payload_start_sample_uncompensated = start+layout.payload_start_s*cfg.fs;
[corrected, available] = uac4.decode_fractional(rx, cfg, nBits, start, result.alpha_hat);
[uncorrected, uncorrectedAvailable] = uac4.decode_fractional(rx, cfg, nBits, start, 1);
result.bits_hat_compensated = corrected;
result.bits_hat_uncompensated = uncorrected;
result.bits_hat = corrected;
result.detected = available;
result.diagnostic.corrected_payload_available = available;
result.diagnostic.uncorrected_payload_available = uncorrectedAvailable;
result.diagnostic.reason = 'acquired';
if ~available, result.diagnostic.reason = 'payload_outside_record'; end
end
