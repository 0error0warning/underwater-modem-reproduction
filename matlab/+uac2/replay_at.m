function y = replay_at(x, fs_x, h, fs_t, fs_tau, fc, start_s)
%REPLAY_AT Watermark direct replay at an absolute measured-channel time.
%   Equivalent to replaying [zeros(start_s*fs_x,1); x] and retaining the
%   output from start_s onward, without calculating the long zero prefix.
%   This function assumes zero input before this local packet. The caller
%   must include channel-tail padding if preceding packets could overlap.
%
%   Supported clock relation: fs_x/fs_tau is a positive integer (3 for the
%   supplied NOF1/NCS1 records), and start_s is on the 1/fs_tau time grid.
%   Channel interpolation and carrier phases use absolute sample indices.
%   Default MATLAB resample filters are preserved, including their small
%   pre-ringing before the local packet, by retaining 64 baseband zero cells.
validateattributes(x, {'numeric'}, {'column','nonempty','finite'});
validateattributes(h, {'numeric'}, {'2d','nonempty','finite'});
validateattributes(fs_x, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(fs_t, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(fs_tau, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(fc, {'numeric'}, {'scalar','real','finite'});
validateattributes(start_s, {'numeric'}, {'scalar','real','finite','nonnegative'});
fs_x = double(fs_x); fs_t = double(fs_t); fs_tau = double(fs_tau); fc = double(fc);
ratio = round(fs_x/fs_tau);
assert(ratio >= 1 && abs(fs_x/fs_tau-ratio) < 1e-12, ...
    'uac2:NonIntegerSampleRatio','fs_x/fs_tau must be a positive integer.');
startBase = round(start_s*fs_tau);
assert(abs(start_s*fs_tau-startBase) <= 1e-7, ...
    'uac2:OffGridStart','start_s must lie on the channel delay-sampling grid.');
K = size(h,2);
assert(abs(fs_t*K-fs_tau) <= 1e-10*max(1,abs(fs_tau)), ...
    'uac2:ChannelSamplingMismatch','Requires fs_t*size(h,2) == fs_tau.');

% A centered downsampling FIR can produce nonzero values before the packet;
% those values can contribute to the retained output after channel filtering.
% Preserve that prefix rather than silently treating local resample edges as
% the beginning of the full transmission. 64 cells exceed the support of
% MATLAB's default integer-ratio resample filters on both conversion steps.
extraBase = min(startBase,64);
windowStartBase = startBase-extraBase;
windowStartPass = windowStartBase*ratio;
cropPass = extraBase*ratio;
inputLength = numel(x);
x = [zeros(cropPass,1); x];
absoluteTime = (windowStartPass+(0:numel(x)-1)')/fs_x;
x = x.*exp(-2*pi*fc*1i*absoluteTime);
[N,D] = rat(fs_tau/fs_x);
x = resample(x,N,D);
Ly = numel(x);
firstBlock = floor(windowStartBase/K);
lastBlock = floor((windowStartBase+Ly-1)/K);
assert(size(h,1) >= lastBlock+2, 'uac2:ShortChannelRecord', ...
    'Absolute replay needs %d channel snapshots; available %d.',lastBlock+2,size(h,1));

yBase = 0*x;
x = [zeros(K-1,1); x];
nFFT = 2^nextpow2(2*K-1);
channelFFT = fft(h(firstBlock+1:lastBlock+2,:).',nFFT,1);
localFirst = 0;
while localFirst < Ly
    globalFirst = windowStartBase+localFirst;
    block = floor(globalFirst/K);
    count = min(K-mod(globalFirst,K),Ly-localFirst);
    signalFFT = fft(x(localFirst+(1:K-1+count)),nFFT,1);
    channelColumn = block-firstBlock+1;
    left = ifft(signalFFT.*channelFFT(:,channelColumn),[],1);
    right = ifft(signalFFT.*channelFFT(:,channelColumn+1),[],1);
    keep = K:K+count-1;
    globalSample = (globalFirst:globalFirst+count-1)';
    fraction = globalSample/K-floor(globalSample/K);
    yBase(localFirst+(1:count)) = (1-fraction).*left(keep)+fraction.*right(keep);
    localFirst = localFirst+count;
end
[N,D] = rat(fs_x/fs_tau);
y = resample(yBase,N,D);
absoluteTime = (windowStartPass+(0:numel(y)-1)')/fs_x;
y = real(y.*exp(2*pi*fc*1i*absoluteTime));
outputLength = ratio*ceil(inputLength/ratio);
y = y(cropPass+(1:outputLength));
end
