function bytes = bits_to_bytes(bits)
%BITS_TO_BYTES Inverse of bytes_to_bits, rejecting incomplete bytes.
bits = bits(:);
assert(all(bits == 0 | bits == 1) && mod(numel(bits),8) == 0, ...
    'uaclink:InvalidBits', 'Expected whole bytes of binary bits.');
bytes = uint8((2.^(7:-1:0))*reshape(double(bits),8,[]));
end
