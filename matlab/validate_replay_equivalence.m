function report = validate_replay_equivalence()
%VALIDATE_REPLAY_EQUIVALENCE Compare FFT replay to untouched Watermark code.
%   Tests complex synthetic channels, partial/exact block boundaries, and
%   one-second BFSK signals through the first official NOF1/NCS1 records.
%   Both maximum absolute error and relative L2 error must be below 1e-10.
matlabDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(matlabDir);
watermarkDir = fullfile(projectDir, 'vendor', 'Watermark');
addpath(matlabDir, fullfile(watermarkDir, 'matlab'));
resultDir = fullfile(projectDir, 'results');
if ~exist(resultDir, 'dir'), mkdir(resultDir); end
resultPath = fullfile(resultDir, 'replay_equivalence.json');
seed = 20261103;
rng(seed, 'twister');
report = struct('status', 'RUNNING', 'seed', seed, 'threshold', 1e-10, ...
    'matlab_version', version, 'official_function', which('replayfilter'), ...
    'fast_function', which('uac.replayfilter_fast'));
writeJson(resultPath, report);

cfg = uac.default_config();
tx = uac.modulate(randi([0, 1], 125, 1), cfg);
assert(numel(tx) == cfg.fs, 'Real-channel test signal must span exactly one second.');
cases = cell(4, 1);
K = 37;
lengths = [7001, 3*K*7];
names = {'synthetic_complex_partial_block', 'synthetic_complex_exact_blocks'};
for k = 1:2
    c = struct();
    c.name = names{k};
    c.fs_x = 48000;
    c.fs_tau = 16000;
    c.fs_t = c.fs_tau/K;
    c.fc = 14000;
    c.tx = randn(lengths(k), 1);
    nResponses = ceil(ceil(lengths(k)/3)/K) + 1;
    c.h = (randn(nResponses, K) + 1i*randn(nResponses, K)) / sqrt(2*K);
    c.source = 'fixed-seed synthetic complex impulse responses';
    cases{k} = c;
end
channels = {'NOF1', 'NCS1'};
for k = 1:2
    name = channels{k};
    path = fullfile(watermarkDir, 'input', 'channels', name, 'mat', [name '_001.mat']);
    record = load(path, 'h', 'fs_t', 'fs_tau', 'fc');
    c = record;
    c.name = [name '_001_first_second'];
    c.fs_x = cfg.fs;
    c.tx = tx;
    c.source = path;
    cases{k+2} = c;
end

resultCells = cell(numel(cases), 1);
for k = 1:numel(cases)
    c = cases{k};
    args = {c.tx, c.fs_x, c.h, c.fs_t, c.fs_tau, c.fc};
    % Warm both functions before timing, and measure three runs each.
    official = replayfilter(args{:});
    fast = uac.replayfilter_fast(args{:});
    assert(isequal(size(official), size(fast)), 'Replay output shapes differ.');
    delta = fast-official;
    maxError = max(abs(delta));
    relativeL2 = norm(delta) / max(norm(official), realmin);
    officialSeconds = zeros(3, 1);
    fastSeconds = zeros(3, 1);
    for trial = 1:3
        started = tic;
        replayfilter(args{:});
        officialSeconds(trial) = toc(started);
        started = tic;
        uac.replayfilter_fast(args{:});
        fastSeconds(trial) = toc(started);
    end
    entry = struct();
    entry.name = c.name;
    entry.source = c.source;
    entry.input_samples = numel(c.tx);
    entry.output_samples = numel(official);
    entry.channel_taps = size(c.h, 2);
    entry.channel_snapshots = size(c.h, 1);
    entry.channel_is_complex = ~isreal(c.h);
    entry.fs_t_times_taps_minus_fs_tau = c.fs_t*size(c.h, 2)-c.fs_tau;
    entry.max_absolute_error = maxError;
    entry.relative_l2_error = relativeL2;
    entry.official_median_seconds = median(officialSeconds);
    entry.fast_median_seconds = median(fastSeconds);
    entry.speedup = entry.official_median_seconds/entry.fast_median_seconds;
    entry.passed = maxError < report.threshold && relativeL2 < report.threshold;
    resultCells{k} = entry;
    fprintf('%s: max abs %.3g, relative L2 %.3g, %.3fx speedup (%.4fs -> %.4fs)\n', ...
        entry.name, maxError, relativeL2, entry.speedup, entry.official_median_seconds, entry.fast_median_seconds);
end
results = vertcat(resultCells{:});
report.status = 'PASS';
if ~all([results.passed]), report.status = 'FAIL'; end
report.results = results;
report.created_utc = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
report.scope = 'Same direct-replay operator on tested inputs; float-roundoff agreement, not bitwise identity. Includes original resample and passband transformations.';
writeJson(resultPath, report);
fprintf('Replay equivalence: %s. Results: %s\n', report.status, resultPath);
assert(all([results.passed]), 'Accelerated replay failed the declared equivalence threshold.');
end

function writeJson(path, value)
file = fopen(path, 'w');
if file < 0, error('uac:WriteFailure', 'Cannot write %s.', path); end
cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
fwrite(file, jsonencode(value), 'char');
end
