function summary = run_awgn_validation()
%RUN_AWGN_VALIDATION Validate noncoherent orthogonal BFSK using base MATLAB.
%   No Signal Processing, Communications, or Statistics Toolbox is required.
%   Run from any directory after adding this file's directory to the path.
%   This validates our PHY implementation; it does not reproduce Watermark
%   channel results or establish agreement with the paper's SNR convention.
matlabDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(matlabDir);
addpath(matlabDir);
outputDir = fullfile(projectDir, 'results', 'awgn');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
cfg = uac.default_config();
seed = 20261102;
rng(seed, 'twister');
[nActive, ~, time] = uac.waveform_parameters(cfg);

% Independent quadrature branches and equal-energy tones are assumptions of
% Pb = 0.5*exp(-Eb/N0/2). Check the sampled geometry, not just frequency gaps.
angles = 2*pi*time*cfg.frequencies(:)';
basis = [cos(angles), sin(angles)];
gram = 2/nActive * (basis' * basis);
gramError = max(abs(reshape(gram-eye(4), [], 1)));
assert(gramError < 1e-10, 'Quadrature templates are not mutually orthonormal.');

% Check both guard settings, unknown carrier phase, and a nonzero time offset.
checkBits = randi([0 1], 10000, 1);
guards = [0, 0.008];
phases = [0, pi/7, pi/2];
guardColumn = zeros(6, 1);
phaseColumn = zeros(6, 1);
checkErrors = zeros(6, 1);
checkCount = repmat(numel(checkBits), 6, 1);
row = 0;
for guard = guards
    for phase = phases
        row = row + 1;
        checkCfg = cfg;
        checkCfg.guard_s = guard;
        checkCfg.phase_rad = phase;
        clean = uac.modulate(checkBits, checkCfg);
        delay = 37;
        decoded = uac.demodulate([zeros(delay, 1); clean], checkCfg, numel(checkBits), delay);
        guardColumn(row) = guard;
        phaseColumn(row) = phase;
        checkErrors(row) = sum(decoded ~= checkBits);
    end
end
[checkLo, checkHi] = uac.wilson_interval(checkErrors, checkCount);
zeroUpper = -expm1(log(0.05)./checkCount);
noiseless = table(guardColumn, phaseColumn, checkCount, checkErrors, checkLo, checkHi, zeroUpper, ...
    'VariableNames', {'guard_s', 'phase_rad', 'n_bits', 'n_errors', 'wilson95_low', 'wilson95_high', 'zero_error_exact_one_sided95_upper'});
writetable(noiseless, fullfile(outputDir, 'noiseless_checks.csv'));
assert(all(checkErrors == 0), 'Noiseless BFSK decoding failed.');

% Real passband sampled noise, with two-sided PSD N0/2:
%   Eb = sum(x.^2)/fs, sigma^2 = N0*fs/2
%      = sum(x.^2)/(2*(Eb/N0)).
% There is one bit per BFSK symbol; silence contributes no bit energy.
% Eb/N0 is NOT per-sample SNR: the latter is 2*(Eb/N0)/nActive here.
oneBit = uac.modulate(0, cfg);
energyDiscrete = sum(oneBit.^2);
otherEnergy = sum(uac.modulate(1, cfg).^2);
assert(abs(otherEnergy-energyDiscrete) < 1e-10*energyDiscrete, 'Unequal bit energies.');
ebn0dB = [0; 4; 8; 12];
nPerPoint = 200000;
chunkSize = 2000;
nBits = repmat(nPerPoint, numel(ebn0dB), 1);
nErrors = zeros(size(ebn0dB));
noiseVariance = zeros(size(ebn0dB));
elapsedSeconds = zeros(size(ebn0dB));
for k = 1:numel(ebn0dB)
    started = tic;
    ratio = 10^(ebn0dB(k)/10);
    noiseVariance(k) = energyDiscrete / (2*ratio);
    for first = 1:chunkSize:nPerPoint
        count = min(chunkSize, nPerPoint-first+1);
        bits = randi([0, 1], count, 1);
        tx = uac.modulate(bits, cfg);
        rx = tx + sqrt(noiseVariance(k))*randn(size(tx));
        decoded = uac.demodulate(rx, cfg, count);
        nErrors(k) = nErrors(k) + sum(decoded ~= bits);
    end
    elapsedSeconds(k) = toc(started);
    fprintf('Eb/N0 %2g dB: %d/%d errors, BER %.6g, %.2f s\n', ...
        ebn0dB(k), nErrors(k), nBits(k), nErrors(k)/nBits(k), elapsedSeconds(k));
end
ber = nErrors ./ nBits;
theory = 0.5*exp(-10.^(ebn0dB/10)/2);
[lower, upper] = uac.wilson_interval(nErrors, nBits);
withinWilson = theory >= lower & theory <= upper;
zScore = (ber-theory) ./ sqrt(theory.*(1-theory)./nBits);
% A 95% interval can miss the true value by chance. For the automated smoke
% check use a predeclared 3.5-sigma bound at each of the four points, while
% reporting the requested 95% intervals and their coverage separately.
pointPassed = abs(zScore) <= 3.5;
zeroErrorUpper = NaN(size(nBits));
zeroErrorUpper(nErrors == 0) = -expm1(log(0.05)./nBits(nErrors == 0));
awgn = table(ebn0dB, nBits, nErrors, ber, theory, lower, upper, withinWilson, zScore, ...
    pointPassed, noiseVariance, zeroErrorUpper, elapsedSeconds, ...
    'VariableNames', {'ebn0_db', 'n_bits', 'n_errors', 'ber', 'theory_ber', 'wilson95_low', ...
    'wilson95_high', 'theory_inside_wilson95', 'theory_z_score', 'within_3p5_sigma', ...
    'real_noise_variance', 'zero_error_exact_one_sided95_upper', 'elapsed_s'});
writetable(awgn, fullfile(outputDir, 'awgn_results.csv'));

f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 680]);
fineDb = linspace(0, 13, 200);
semilogy(fineDb, 0.5*exp(-10.^(fineDb/10)/2), 'k-', 'LineWidth', 1.5);
hold on;
errorbar(ebn0dB, ber, ber-lower, upper-ber, 'o', 'LineStyle', 'none', ...
    'LineWidth', 1.3, 'MarkerSize', 7, 'MarkerFaceColor', [0.15, 0.4, 0.8], 'Color', [0.15, 0.4, 0.8]);
set(gca, 'YScale', 'log', 'FontSize', 12);
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit error rate');
title({'Noncoherent orthogonal BFSK in real AWGN', '200,000 bits/point; 95% Wilson intervals'});
legend('0.5 exp(-E_b/(2N_0))', 'Simulation and 95% interval', 'Location', 'southwest');
xlim([-0.5, 13]);
ylim([1e-5, 0.5]);
print(f, fullfile(outputDir, 'awgn_validation.png'), '-dpng', '-r160');
close(f);

summary = struct();
summary.status = 'PASS';
if ~all(pointPassed), summary.status = 'FAIL'; end
summary.created_utc = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
summary.matlab_version = version;
summary.seed = seed;
summary.config = cfg;
summary.bits_per_point = nPerPoint;
summary.chunk_bits = chunkSize;
summary.template_orthogonality_max_abs_error = gramError;
summary.noiseless_cases_passed = all(checkErrors == 0);
summary.all_theory_inside_wilson95 = all(withinWilson);
summary.all_theory_within_3p5_sigma = all(pointPassed);
summary.max_absolute_theory_z_score = max(abs(zScore));
summary.bit_energy_discrete = energyDiscrete;
summary.bit_energy_seconds = energyDiscrete/cfg.fs;
summary.noise_rule = 'real sigma^2 = sum(active_bit.^2)/(2*10^(EbN0_dB/10)); two-sided PSD N0/2';
summary.frequency_assumption = '13875/14125 Hz: our symmetric choice around 14 kHz with 250 Hz spacing; paper does not specify exact BFSK tones';
summary.timing_assumption = 'symbol_s is active duration; guard_s adds silence; known timing; no channel estimation or synchronization';
summary.acceptance_rule = 'zero noiseless errors and absolute theoretical-binomial z score <= 3.5 at each of four points; 95% Wilson intervals reported separately';
summary.scope = 'AWGN implementation validation only, not Watermark or a reproduction of the paper tables';
writeJson(fullfile(outputDir, 'validation_status.json'), summary);
fprintf('AWGN validation: %s. Results: %s\n', summary.status, outputDir);
assert(all(pointPassed), 'AWGN Monte Carlo results exceed the predeclared theory tolerance.');
end

function writeJson(path, value)
file = fopen(path, 'w');
if file < 0, error('uac:WriteFailure', 'Cannot write %s.', path); end
cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
fwrite(file, jsonencode(value), 'char');
end
