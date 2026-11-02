function values = sinc_sample(x, positions, halfWidth)
%SINC_SAMPLE Windowed-sinc samples at fractional one-based input positions.
%   Lanczos interpolation avoids the large carrier attenuation and timing
%   bias of linear interpolation at 14 kHz with a 48 kHz sample rate.
shape = size(positions);
q = positions(:);
x = x(:);
centre = floor(q);
values = zeros(size(q), 'like', x);
weightSum = zeros(size(q));
for offset = -halfWidth+1:halfWidth
    index = centre+offset;
    distance = q-index;
    weight = localSinc(distance).*localSinc(distance/halfWidth);
    good = index >= 1 & index <= numel(x);
    values(good) = values(good) + weight(good).*x(index(good));
    weightSum(good) = weightSum(good) + weight(good);
end
values = values ./ max(weightSum, realmin);
values = reshape(values, shape);
end

function y = localSinc(x)
y = ones(size(x));
nonzero = abs(x) > 1e-12;
y(nonzero) = sin(pi*x(nonzero))./(pi*x(nonzero));
end
