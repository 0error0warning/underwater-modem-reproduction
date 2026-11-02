function value = crc16(bytes)
%CRC16 CRC-16/CCITT-FALSE: poly1021 initFFFF refin=false refout=false xorout0.
if ~isempty(bytes), validateattributes(bytes, {'numeric'}, {'vector', 'integer', '>=', 0, '<=', 255}); end
value = uint16(65535);
for index = 1:numel(bytes)
    value = bitxor(value, bitshift(uint16(bytes(index)), 8));
    for bit = 1:8
        high = bitget(value, 16);
        value = bitshift(value, 1);
        if high, value = bitxor(value, uint16(hex2dec('1021'))); end
    end
end
end
