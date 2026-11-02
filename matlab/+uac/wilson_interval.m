function [lower, upper] = wilson_interval(nErrors, nTrials, confidence)
%WILSON_INTERVAL Binomial Wilson score interval; no Statistics Toolbox.
%   Scalar or same-sized array inputs are accepted. For zero errors the
%   upper endpoint is nonzero; zero observed errors is not zero true BER.
if nargin < 3, confidence = 0.95; end
validateattributes(nErrors, {'numeric'}, {'real', 'finite', 'integer', 'nonnegative'});
validateattributes(nTrials, {'numeric'}, {'real', 'finite', 'integer', 'positive'});
validateattributes(confidence, {'numeric'}, {'real', 'finite', 'scalar', '>', 0, '<', 1});
if ~(isequal(size(nErrors), size(nTrials)) || isscalar(nErrors) || isscalar(nTrials))
    error('uac:InvalidCounts', 'Counts must have compatible scalar or equal sizes.');
end
if any(nErrors(:) > nTrials(:))
    error('uac:InvalidCounts', 'Error counts cannot exceed trial counts.');
end
z = sqrt(2) * erfinv(confidence);
p = nErrors ./ nTrials;
denominator = 1 + z^2 ./ nTrials;
centre = (p + z^2 ./ (2*nTrials)) ./ denominator;
halfWidth = z ./ denominator .* sqrt(p.*(1-p)./nTrials + z^2./(4*nTrials.^2));
lower = max(0, centre - halfWidth);
upper = min(1, centre + halfWidth);
lower(p == 0) = 0;
upper(p == 1) = 1;
end
