function result = decode_frame(bytes, cfg)
%DECODE_FRAME Bounded framing and CRC validation; invalid frames deliver no data.
result = struct('valid', false, 'reason', 'short_header', 'type', NaN, ...
    'sequence', NaN, 'payload', zeros(1,0,'uint8'), 'payload_length', NaN, ...
    'crc_received', NaN, 'crc_computed', NaN);
header = uaclink.parse_header(bytes,cfg);
result.reason = header.reason;
if ~header.valid, return; end
result.type=header.type; result.sequence=header.sequence;
result.payload_length=header.payload_length;
if numel(bytes) ~= header.frame_bytes, result.reason='frame_length_mismatch'; return; end
bytes=uint8(bytes(:)');
result.crc_received=double(bytes(end-1))*256+double(bytes(end));
result.crc_computed=double(uaclink.crc16(bytes(1:end-2)));
if result.crc_received ~= result.crc_computed, result.reason='crc_mismatch'; return; end
result.valid=true; result.reason='crc_valid'; result.payload=bytes(9:end-2);
end
