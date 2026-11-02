function result = receive(waveform, cfg)
%RECEIVE Blind waveform-to-packet receiver: only waveform and public cfg.
% Acquire/decode the 64-bit header first, then use its bounded LENGTH field.
% Transmitted bits, length, offset, velocity, channel and layout never enter.
result = struct('valid', false, 'reason', 'not_acquired', 'type', NaN, ...
    'sequence', NaN, 'payload', zeros(1,0,'uint8'), 'payload_length', NaN, ...
    'crc_received', NaN, 'crc_computed', NaN, 'header_bits', 64, ...
    'decoded_frame_bits', 0, 'alpha_hat', NaN, 'preamble_start_sample', NaN);
if isempty(waveform), result.reason='empty_waveform'; return; end
acquired=uac4.receive_timing(waveform,cfg,64,cfg.link.timing_mode);
result.alpha_hat=acquired.alpha_hat;
result.preamble_start_sample=acquired.preamble_start_sample;
if ~acquired.detected, result.reason=acquired.diagnostic.reason; return; end
header=uaclink.parse_header(uaclink.bits_to_bytes(acquired.bits_hat),cfg);
result.reason=header.reason;
if ~header.valid, return; end
result.payload_length=header.payload_length;
nBits=8*header.frame_bytes;
[bits,available]=uac4.decode_fractional(waveform,cfg,nBits, ...
    acquired.preamble_start_sample,acquired.alpha_hat);
if ~available, result.reason='truncated_frame'; return; end
decoded=uaclink.decode_frame(uaclink.bits_to_bytes(bits),cfg);
names=fieldnames(decoded);
for k=1:numel(names), result.(names{k})=decoded.(names{k}); end
result.decoded_frame_bits=nBits;
end
