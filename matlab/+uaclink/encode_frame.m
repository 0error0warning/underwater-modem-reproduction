function [bytes, bits] = encode_frame(payload, sequence, frameType, cfg)
%ENCODE_FRAME Serialize DATA(1)/ACK(2); CRC covers header and payload.
validateattributes(sequence, {'numeric'}, {'scalar','integer','>=',0,'<=',65535});
validateattributes(frameType, {'numeric'}, {'scalar','integer','>=',1,'<=',2});
if ~isempty(payload), validateattributes(payload, {'numeric'}, {'vector','integer','>=',0,'<=',255}); end
payload = uint8(payload(:)');
assert(numel(payload) <= cfg.link.max_payload_bytes, ...
    'uaclink:PayloadTooLong', 'Payload exceeds the configured byte bound.');
assert(frameType ~= 2 || isempty(payload), 'uaclink:InvalidAck', 'ACK must have zero payload.');
lengthBytes = [floor(numel(payload)/256),mod(numel(payload),256)];
header = uint8([211,145,1,frameType,floor(sequence/256),mod(sequence,256),lengthBytes]);
covered = [header,payload]; crc = double(uaclink.crc16(covered));
bytes = [covered,uint8([floor(crc/256),mod(crc,256)])];
bits = uaclink.bytes_to_bits(bytes);
end
