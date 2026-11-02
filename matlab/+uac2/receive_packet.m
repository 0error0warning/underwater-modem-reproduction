function result = receive_packet(rx, cfg, nBits, opts)
%RECEIVE_PACKET Acquire and decode without transmitted bits, h, V0, or offset.
%   Only rx, public waveform settings, and payload length enter this receiver.
%   alpha_hat = received duration / transmitted duration = 1/(1-v/c).
%   Positive velocity stretches the waveform and lowers its received tones.
%   Both corrected/uncorrected bit vectors use the SAME acquired preamble.
%   opts.compensate selects bits_hat; it never changes acquisition/estimation.
if nargin < 4, opts = struct(); end
if ~isfield(opts, 'compensate'), opts.compensate = true; end
validateattributes(rx, {'numeric'}, {'real', 'finite', 'vector', 'nonempty'});
rx = rx(:);
layout = uac2.packet_layout(cfg, nBits);
result = emptyResult();
result.diagnostic = struct('reason', 'not_acquired', 'used_payload_truth', false, ...
    'used_channel_truth', false, 'pilot_is_new_design', true);
fs = cfg.fs;
fsSync = cfg.sync_fs;
fc = mean(cfg.frequencies);
analytic = hilbert(rx, 2^nextpow2(numel(rx)));
analytic = analytic(1:numel(rx));
t = (0:numel(rx)-1)'/fs;
[P, Q] = rat(fsSync/fs);
baseband = resample(analytic.*exp(-2*pi*1i*fc*t), P, Q);

up = template(cfg, 1, true, fc);
down = template(cfg, 1, false, fc);
if numel(baseband) < numel(up)+round(layout.down_start_s*fsSync)+1
    result.diagnostic.reason = 'short_preamble_record';
    return;
end
[~, upRho] = correlate(baseband, up);
[~, downRho] = correlate(baseband, down);
% HFM time/Doppler ambiguity can shift the two peaks differently. Permit a
% small DOWN timing neighborhood during initial acquisition, then use pilots.
tol = ceil(0.002*fsSync);
downMaximum = movmax(downRho, [tol, tol]);
spacing = round(layout.down_start_s*fsSync);
pairedDown = zeros(size(upRho));
available = min(numel(upRho), numel(downMaximum)-spacing);
pairedDown(1:available) = downMaximum(spacing+(1:available));
statistic = -numel(up)*(log1p(-upRho)+log1p(-pairedDown));
statistic(available+1:end) = 0;
[coarseStatistic, coarse] = max(statistic);
nTrials = max(1, available*(2*tol+1));
threshold = gammaThreshold(cfg.false_alarm_probability/nTrials);
result.sync_score = (upRho(coarse)*pairedDown(coarse))^0.25;
result.diagnostic.coarse_statistic = coarseStatistic;
result.diagnostic.acquisition_threshold = threshold;
result.diagnostic.false_alarm_target = cfg.false_alarm_probability;
if coarseStatistic < threshold
    result.diagnostic.reason = 'preamble_below_detection_threshold';
    return;
end

% Frequency search uses ONLY central windows of the two PUBLIC single-tone
% pilots. Unknown channel/initial phases disappear in squared correlations.
[velocity, pilotDiagnostic] = pilotVelocity(baseband, coarse-1, cfg, layout, fc);
result.diagnostic.pilot = pilotDiagnostic;
if ~pilotDiagnostic.available
    result.diagnostic.reason = 'short_pilot_record';
    return;
end
alpha = 1/(1-velocity/cfg.sound_speed_mps);
result.alpha_hat = alpha;
result.velocity_hat = velocity;

% With the estimated time scale, sharpen both known HFM correlations. Search
% around the coarse strong arrival, then estimate earlier echoes from HFM only.
up = template(cfg, alpha, true, fc);
down = template(cfg, alpha, false, fc);
[upCorr, upRho] = correlate(baseband, up);
[downCorr, downRho] = correlate(baseband, down);
positions = (1:numel(upRho))';
downPositions = positions + alpha*layout.down_start_s*fsSync;
downAtUp = interp1((1:numel(downRho))', downRho, downPositions, 'pchip', 0);
downAtUp = min(max(downAtUp, 0), 1-1e-12);
joint = -numel(up)*(log1p(-upRho)+log1p(-downAtUp));
local = max(1, coarse-round(0.005*fsSync)):min(numel(joint), coarse+round(0.005*fsSync));
[~, peakLocal] = max(joint(local));
strong = local(peakLocal);
powerDown = interp1((1:numel(downCorr))', abs(downCorr).^2/sum(abs(down).^2), downPositions, 'pchip', 0);
profile = abs(upCorr).^2/sum(abs(up).^2) + max(powerDown, 0);
[arrival, spreadMs, spreadDiagnostic] = estimateSpread(profile, strong, fsSync, cfg);
result.delay_spread_ms = spreadMs;
result.diagnostic.delay = spreadDiagnostic;
result.diagnostic.refined_statistic = joint(strong);
result.diagnostic.reason = 'acquired';
result.preamble_start_sample = 1+(arrival-1)*fs/fsSync;
result.payload_start_sample = result.preamble_start_sample+alpha*layout.payload_start_s*fs;
result.payload_start_sample_uncompensated = result.preamble_start_sample+layout.payload_start_s*fs;

[corrected, correctedAvailable] = decode(rx, cfg, nBits, result.preamble_start_sample, alpha, layout);
[uncorrected, uncorrectedAvailable] = decode(rx, cfg, nBits, result.preamble_start_sample, 1, layout);
result.bits_hat_compensated = corrected;
result.bits_hat_uncompensated = uncorrected;
if opts.compensate
    result.bits_hat = corrected;
    result.detected = correctedAvailable;
else
    result.bits_hat = uncorrected;
    result.detected = uncorrectedAvailable;
end
result.diagnostic.corrected_payload_available = correctedAvailable;
result.diagnostic.uncorrected_payload_available = uncorrectedAvailable;
if ~result.detected, result.diagnostic.reason = 'payload_outside_record'; end
end

function result = emptyResult()
result = struct('detected', false, 'bits_hat', zeros(0, 1), ...
    'bits_hat_compensated', zeros(0, 1), 'bits_hat_uncompensated', zeros(0, 1), ...
    'preamble_start_sample', NaN, 'payload_start_sample', NaN, ...
    'payload_start_sample_uncompensated', NaN, 'alpha_hat', NaN, ...
    'velocity_hat', NaN, 'sync_score', 0, 'delay_spread_ms', NaN, 'diagnostic', struct());
end

function h = template(cfg, alpha, isUp, fc)
t = (0:ceil(cfg.hfm_pulse_s*alpha*cfg.sync_fs)-1)'/cfg.sync_fs;
first = cfg.hfm_low_hz;
last = cfg.hfm_high_hz;
if ~isUp, first = cfg.hfm_high_hz; last = cfg.hfm_low_hz; end
h = uac2.hfm_analytic(t/alpha, first, last, cfg.hfm_pulse_s, cfg.hfm_taper_fraction);
h = h.*exp(-2*pi*1i*fc*t);
end

function [correlation, rhoSquared] = correlate(x, h)
N = numel(x);
L = numel(h);
nFFT = 2^nextpow2(N+L-1);
full = ifft(fft(x, nFFT).*fft(flipud(conj(h)), nFFT));
correlation = full(L:N);
cumulative = [0; cumsum(abs(x).^2)];
energy = cumulative(L+1:end)-cumulative(1:end-L);
rhoSquared = abs(correlation).^2 ./ max(energy*sum(abs(h).^2), realmin);
rhoSquared = min(max(real(rhoSquared), 0), 1-1e-12);
rhoSquared(energy <= max(max(energy)*1e-12, realmin)) = 0;
end

function threshold = gammaThreshold(tailProbability)
% Gamma(shape=2,scale=1) survival is exp(-s)*(1+s); Bonferroni over searches.
target = log(tailProbability);
lo = 0;
hi = max(100, -2*target);
for k = 1:50
    mid = (lo+hi)/2;
    if -mid+log1p(mid) > target, lo = mid; else, hi = mid; end
end
threshold = hi;
end

function [velocity, diagnostic] = pilotVelocity(x, startZero, cfg, layout, fc)
velocity = NaN;
diagnostic = struct('available', false, 'score', NaN, 'at_search_boundary', false);
grid = (cfg.velocity_range_mps(1):cfg.velocity_grid_step_mps:cfg.velocity_range_mps(2))';
scores = zeros(size(grid));
starts = [layout.pilot0_start_s, layout.pilot1_start_s];
for tone = 1:2
    first = round(startZero+(starts(tone)+cfg.pilot_trim_s)*cfg.sync_fs)+1;
    last = round(startZero+(starts(tone)+cfg.pilot_s-cfg.pilot_trim_s)*cfg.sync_fs);
    if first < 1 || last > numel(x) || last <= first, return; end
    segment = x(first:last);
    time = (0:numel(segment)-1)'/cfg.sync_fs;
    frequencies = cfg.frequencies(tone)*(1-grid'/cfg.sound_speed_mps)-fc;
    bank = exp(-2*pi*1i*time*frequencies);
    projection = segment.'*bank;
    scores = scores+(abs(projection).^2/(numel(segment)*max(sum(abs(segment).^2), realmin)))';
end
[score, best] = max(scores);
velocity = grid(best);
if best > 1 && best < numel(grid)
    denominator = scores(best-1)-2*scores(best)+scores(best+1);
    if denominator < 0
        offset = 0.5*(scores(best-1)-scores(best+1))/denominator;
        velocity = velocity+max(-1, min(1, offset))*cfg.velocity_grid_step_mps;
    end
end
diagnostic.available = true;
diagnostic.score = score/2;
diagnostic.at_search_boundary = best == 1 || best == numel(grid);
diagnostic.search_range_mps = cfg.velocity_range_mps;
diagnostic.trim_s = cfg.pilot_trim_s;
end

function [arrival, spreadMs, diagnostic] = estimateSpread(profile, strong, fs, cfg)
radius = round(cfg.delay_search_radius_s*fs);
range = (max(1, strong-radius):min(numel(profile), strong+radius))';
p = profile(range);
noiseFloor = median(p);
% With white complex Gaussian background, the sum of two independently
% matched HFM powers has a Gamma(2,scale) null distribution. Its unit-scale
% median is 1.67834699. A Bonferroni bound avoids labeling a rare isolated
% noise maximum as a long delayed echo merely because the scan is wide.
noiseScale = noiseFloor/1.67834699001666;
noiseCutoff = noiseScale*gammaThreshold(cfg.false_alarm_probability/numel(p));
cutoff = max(0.04*max(p), noiseCutoff);
peaks = find(p(2:end-1) >= p(1:end-2) & p(2:end-1) > p(3:end) & p(2:end-1) >= cutoff)+1;
if isempty(peaks)
    [~, peaks] = max(p);
end
% Resolve arrivals only beyond four inverse swept-bandwidths (0.5 ms at
% 8 kHz). Keep the strongest peak in that neighborhood so the known HFM
% autocorrelation point-spread/sidelobes do not masquerade as channel echoes.
resolutionS = 4/(cfg.hfm_high_hz-cfg.hfm_low_hz);
separation = max(1, ceil(resolutionS*fs));
[~, order] = sort(p(peaks), 'descend');
selected = zeros(0, 1);
for k = 1:numel(order)
    candidate = peaks(order(k));
    if isempty(selected) || all(abs(candidate-selected) >= separation)
        selected(end+1, 1) = candidate; %#ok<AGROW>
    end
end
peaks = sort(selected);
% Only matched-HFM peak powers contribute; payload samples are never used.
weights = max(p(peaks)-noiseFloor, 0);
cumulative = cumsum(weights)/max(sum(weights), realmin);
first90 = peaks(find(cumulative >= 0.05, 1));
last90 = peaks(find(cumulative >= 0.95, 1));
first = range(peaks(1));
last = range(peaks(end));
arrival = first;
if first > 1 && first < numel(profile)
    local = log(max(profile(first-1:first+1), realmin));
    denominator = local(1)-2*local(2)+local(3);
    if denominator < 0
        arrival = first+max(-0.5, min(0.5, 0.5*(local(1)-local(3))/denominator));
    end
end
spreadMs = (last-first)/fs*1000;
diagnostic = struct('source', 'known HFM UP/DOWN correlation peaks only', ...
    'search_radius_ms', cfg.delay_search_radius_s*1000, 'noise_floor', noiseFloor, ...
    'peak_power_threshold', cutoff, 'noise_cfar_threshold', noiseCutoff, 'significant_peaks', numel(peaks), ...
    'arrival_resolution_ms', resolutionS*1000, ...
    'spread_90_percent_ms', (last90-first90)/fs*1000, ...
    'strongest_arrival_sample_at_sync_rate', strong, ...
    'first_arrival_sample_at_sync_rate', arrival, ...
    'near_search_boundary', peaks(1) <= 2 || peaks(end) >= numel(p)-1);
end

function [bits, available] = decode(rx, cfg, nBits, preambleStart, alpha, layout)
[nActive, nGuard, time] = uac.waveform_parameters(cfg);
offsets = bsxfun(@plus, (0:nActive-1)', (0:nBits-1)*(nActive+nGuard));
positions = preambleStart+alpha*(layout.payload_start_s*cfg.fs+offsets);
available = all(positions(:) >= 1 & positions(:) <= numel(rx));
if ~available, bits = zeros(0, 1); return; end
samples = uac2.sinc_sample(rx, positions, cfg.sinc_half_width);
angles = 2*pi*time*cfg.frequencies(:)';
energy = (cos(angles)'*samples).^2+(sin(angles)'*samples).^2;
bits = double(energy(2, :)' > energy(1, :)');
end
