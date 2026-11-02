function bits = bytes_to_bits(bytes)
%BYTES_TO_BITS Network order: most significant bit first within each byte.
if ~isempty(bytes), validateattributes(bytes, {'numeric'}, {'vector', 'integer', '>=', 0, '<=', 255}); end
bytes = uint8(bytes(:));
matrix = zeros(8, numel(bytes));
for k = 1:8, matrix(k,:) = double(bitget(bytes, 9-k)); end
bits = matrix(:);
end
