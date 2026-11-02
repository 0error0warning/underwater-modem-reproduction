function header = parse_header(bytes, cfg)
%PARSE_HEADER Parse exactly the public fixed 8-byte wire header; no truth.
header = struct('valid', false, 'reason', 'short_header', 'type', NaN, ...
    'sequence', NaN, 'payload_length', NaN, 'frame_bytes', NaN);
if numel(bytes) < 8, return; end
bytes = double(bytes(:)');
if ~isequal(bytes(1:2), [211,145]), header.reason='bad_magic'; return; end
if bytes(3) ~= 1, header.reason='unsupported_version'; return; end
header.type = bytes(4); header.sequence = bytes(5)*256+bytes(6);
header.payload_length = bytes(7)*256+bytes(8);
if ~ismember(header.type,[1,2]), header.reason='unknown_type'; return; end
if header.payload_length > cfg.link.max_payload_bytes
    header.reason='length_exceeds_limit'; return;
end
if header.type == 2 && header.payload_length ~= 0
    header.reason='ack_with_payload'; return;
end
header.frame_bytes = 8+header.payload_length+2;
header.valid = true; header.reason = 'header_valid';
end
