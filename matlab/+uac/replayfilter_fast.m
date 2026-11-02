function y = replayfilter_fast(x, fs_x, h, fs_t, fs_tau, fc)
%REPLAYFILTER_FAST FFT acceleration of Watermark 1.0 direct replay.
%   Same input/output contract and baseband/resampling operations as the
%   official replayfilter.m. h is time-by-delay, with increasing causal lag.
%   Adjacent channel snapshots are linearly interpolated at every sample.
%   This is an implementation optimization, not a different channel model.
%
%   Within each K-sample block, convolution with each endpoint impulse
%   response is evaluated by overlap-save FFT. Linearity allows their outputs
%   to be blended with the original sample-dependent interpolation weights.
%   Complex coefficients are NEVER conjugated. In particular, use .' not '.
validateattributes(x, {'numeric'}, {'vector', 'column', 'nonempty', 'finite'});
validateattributes(h, {'numeric'}, {'2d', 'nonempty', 'finite'});
validateattributes(fs_x, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(fs_t, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(fs_tau, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(fc, {'numeric'}, {'scalar', 'real', 'finite'});
K = size(h, 2);
assert(abs(fs_t*K-fs_tau) <= 1e-10*max(1, abs(fs_tau)), ...
    'uac:ChannelSamplingMismatch', 'This replay convention requires fs_t*size(h,2) == fs_tau.');

% Keep these transformations identical to the official implementation.
t = (0:length(x)-1)' / fs_x;
x = x .* exp(-2*pi*fc*1i*t);
[N, D] = rat(fs_tau/fs_x);
x = resample(x, N, D);

Ly = length(x);
nBlocks = ceil(Ly/K);
assert(size(h, 1) >= nBlocks+1, 'uac:ShortChannelRecord', ...
    'The input requires %d channel snapshots, but h has only %d.', nBlocks+1, size(h, 1));
y = 0*x;
x = [zeros(K-1, 1); x];

% The FFT spans the full history-plus-new-data block. Circular aliasing is
% confined to its first K-1 outputs, which overlap-save discards. This also
% works for non-power-of-two K and a partial final block.
nFFT = 2^nextpow2(2*K-1);
channelFFT = fft(h(1:nBlocks+1, :).', nFFT, 1);
for block = 0:nBlocks-1
    first = block*K;
    count = min(K, Ly-first);
    segment = x(first + (1:K-1+count));
    signalFFT = fft(segment, nFFT, 1);
    left = ifft(signalFFT .* channelFFT(:, block+1), [], 1);
    right = ifft(signalFFT .* channelFFT(:, block+2), [], 1);
    keep = K:K+count-1;
    % Use the same expression as the original for identical floating-point
    % interpolation weights, including non-power-of-two tap counts.
    sample = (first:first+count-1)';
    fraction = sample/K - floor(sample/K);
    y(first+(1:count)) = (1-fraction).*left(keep) + fraction.*right(keep);
end

[N, D] = rat(fs_x/fs_tau);
y = resample(y, N, D);
t = (0:length(y)-1)' / fs_x;
y = real(y .* exp(2*pi*fc*1i*t));
end
