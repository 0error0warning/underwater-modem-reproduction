function report = validate_receiver_v2()
%VALIDATE_RECEIVER_V2 Blind receiver tests; truth is confined to this harness.
%   No transmitted bits, physical offset, velocity, or channel taps enter
%   uac2.receive_packet. Noise-only tests measure finite-sample false alarms.
matlabDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(matlabDir);
addpath(matlabDir);
outputDir = fullfile(projectDir, 'results', 'receiver_v2');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
resultPath = fullfile(outputDir, 'validation.json');
seed = 20261104;
rng(seed, 'twister');
cfg = uac2.default_config();
report = struct('status', 'RUNNING', 'seed', seed, 'matlab_version', version);
writeJson(resultPath, report);
rows = {};
definitions = {};
for k = 1:4
    definitions{end+1} = definition(sprintf('clean_random_phase_%02d', k), 'clean', 0, Inf, 0, true); %#ok<AGROW>
end
velocities = [0.0625, -0.158, 0.16, -0.16, 0.5, -0.5, 1.5, -1.5];
for k = 1:numel(velocities)
    definitions{end+1} = definition(sprintf('doppler_%+.4f_mps', velocities(k)), 'doppler', velocities(k), Inf, 0, true); %#ok<AGROW>
end
for snr = [0, 10, 20]
    for k = 1:3
        definitions{end+1} = definition(sprintf('awgn_%02ddb_%02d', snr, k), 'awgn', 0.5, snr, 0, snr >= 10); %#ok<AGROW>
    end
end
for v = [0, 0.5, -0.5]
    definitions{end+1} = definition(sprintf('multipath_%+.1f_mps', v), 'multipath', v, 20, 0.016, true); %#ok<AGROW>
end
for k = 1:numel(definitions)
    d = definitions{k};
    rxCfg = cfg;
    rxCfg.guard_s = d.guard_s;
    txCfg = rxCfg;
    txCfg.phase_rad = 2*pi*rand();
    bits = randi([0, 1], 128, 1);
    [tx, layout] = uac2.make_packet(bits, txCfg);
    if strcmp(d.category, 'multipath')
        delays = round([0, 0.004, 0.011]*cfg.fs);
        taps = complex(zeros(delays(end)+1, 1));
        taps(delays+1) = [1, 0.65*exp(0.8i), 0.45*exp(-1.1i)];
        extendedLength = numel(tx)+numel(taps)-1;
        analytic = hilbert(tx, 2^nextpow2(numel(tx)));
        analytic = analytic(1:numel(tx));
        nFFT = 2^nextpow2(extendedLength);
        wave = ifft(fft(analytic, nFFT).*fft(taps, nFFT));
        tx = real(wave(1:extendedLength));
    end
    alpha = 1/(1-d.velocity_mps/cfg.sound_speed_mps);
    [P, Q] = rat(alpha, 1e-7);
    actualAlpha = P/Q;
    actualVelocity = cfg.sound_speed_mps*(1-1/actualAlpha);
    % Independent MATLAB polyphase resampling simulates physical time scale.
    stretched = resample(tx, P, Q);
    offset = round((0.15+0.45*rand())*cfg.fs);
    rx = [zeros(offset, 1); stretched; zeros(round(0.12*cfg.fs), 1)];
    if isfinite(d.ebn0_db)
        nActive = round(cfg.symbol_s*cfg.fs);
        noiseVariance = (nActive/2)/(2*10^(d.ebn0_db/10));
        rx = rx+sqrt(noiseVariance)*randn(size(rx));
    end
    started = tic;
    result = uac2.receive_packet(rx, rxCfg, numel(bits), struct('compensate', true));
    runtime = toc(started);
    row = emptyRow(d.name, d.category);
    row.required = d.required;
    row.true_velocity_mps = actualVelocity;
    row.ebn0_db = d.ebn0_db;
    row.true_offset_s = offset/cfg.fs;
    row.guard_s = d.guard_s;
    row.n_bits = numel(bits);
    row.detected = result.detected;
    row.velocity_hat_mps = result.velocity_hat;
    row.velocity_error_mps = result.velocity_hat-actualVelocity;
    row.preamble_error_ms = (result.preamble_start_sample-(offset+1))/cfg.fs*1000;
    row.payload_error_ms = (result.payload_start_sample-(offset+1+actualAlpha*(layout.payload_start_sample-1)))/cfg.fs*1000;
    row.delay_spread_ms = result.delay_spread_ms;
    row.sync_score = result.sync_score;
    row.runtime_s = runtime;
    row.reason = result.diagnostic.reason;
    if result.detected
        row.errors_corrected = sum(result.bits_hat_compensated ~= bits);
        row.ber_corrected = row.errors_corrected/numel(bits);
        if numel(result.bits_hat_uncompensated) == numel(bits)
            row.errors_uncorrected = sum(result.bits_hat_uncompensated ~= bits);
            row.ber_uncorrected = row.errors_uncorrected/numel(bits);
        end
    end
    if strcmp(d.category, 'clean') || strcmp(d.category, 'doppler')
        row.passed = row.detected && row.errors_corrected == 0 && abs(row.velocity_error_mps) < 0.05 && abs(row.preamble_error_ms) < 0.5;
    elseif strcmp(d.category, 'multipath')
        row.passed = row.detected && row.ber_corrected <= 0.02 && abs(row.preamble_error_ms) < 0.75 && row.delay_spread_ms >= 8 && row.delay_spread_ms <= 14;
    else
        row.passed = row.detected && row.ber_corrected <= 0.10 && abs(row.velocity_error_mps) < 0.6 && row.delay_spread_ms < 2;
    end
    rows{end+1} = row; %#ok<AGROW>
    fprintf('%s: detected=%d BER=%g vhat=%+.4f timing=%+.3f ms delay=%g ms %.3fs pass=%d\n', ...
        row.name, row.detected, row.ber_corrected, row.velocity_hat_mps, row.preamble_error_ms, row.delay_spread_ms, runtime, row.passed);
end

% Pure noise is never labelled a valid packet. The threshold is fixed in cfg,
% derived from normalized matched-filter null statistics, not payload labels.
nNoiseTests = 50;
for k = 1:nNoiseTests
    duration = 1.8+0.6*rand();
    rx = 10^(-1+2*rand())*randn(round(duration*cfg.fs), 1);
    started = tic;
    result = uac2.receive_packet(rx, cfg, 128, struct('compensate', true));
    row = emptyRow(sprintf('noise_only_%02d', k), 'noise_only');
    row.required = true;
    row.n_bits = 0;
    row.detected = result.detected;
    row.sync_score = result.sync_score;
    row.runtime_s = toc(started);
    row.reason = result.diagnostic.reason;
    row.passed = ~result.detected;
    rows{end+1} = row; %#ok<AGROW>
end
records = vertcat(rows{:});
tableResults = struct2table(records);
writetable(tableResults, fullfile(outputDir, 'validation_cases.csv'));
required = [records.required];
passed = [records.passed];
report.status = 'PASS';
if ~all(passed(required)), report.status = 'FAIL'; end
report.config = cfg;
report.total_cases = numel(records);
report.required_cases = sum(required);
report.required_passed = sum(passed(required));
report.noise_only_cases = nNoiseTests;
noiseRows = strcmp({records.category}, 'noise_only');
report.noise_only_false_alarms = sum([records(noiseRows).detected]);
report.noise_only_zero_count_one_sided95_upper = -expm1(log(0.05)/nNoiseTests);
stress = strcmp({records.category}, 'awgn') & [records.ebn0_db] == 0;
report.ebn0_zero_stress_cases = sum(stress);
report.ebn0_zero_stress_detected = sum([records(stress).detected]);
report.median_runtime_s = median([records.runtime_s]);
report.max_runtime_s = max([records.runtime_s]);
report.detected_median_runtime_s = median([records([records.detected]).runtime_s]);
publicLayout = uac2.packet_layout(cfg, 128);
report.payload_start_s = publicLayout.payload_start_s;
report.truth_isolation = 'Receiver receives only rx, public cfg, nBits, and compensation flag. Truth appears only in this validation harness.';
report.low_snr_scope = '0 dB cases are stress characterization, not required zero-error checks. Detection failures and BER are retained in the CSV.';
report.false_alarm_scope = 'Zero observed false alarms in 50 tests is not proof of zero false-alarm probability; exact one-sided bound is reported.';
report.created_utc = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));
writeJson(resultPath, report);
fprintf('Receiver v2 validation: %s (%d/%d required); noise false alarms %d/%d; median %.3fs\n', ...
    report.status, report.required_passed, report.required_cases, report.noise_only_false_alarms, nNoiseTests, report.median_runtime_s);
assert(all(passed(required)), 'Blind receiver failed a required validation case.');
end

function d = definition(name, category, velocity, ebn0, guard, required)
d = struct('name', name, 'category', category, 'velocity_mps', velocity, ...
    'ebn0_db', ebn0, 'guard_s', guard, 'required', required);
end

function row = emptyRow(name, category)
row = struct('name', name, 'category', category, 'required', true, 'passed', false, ...
    'true_velocity_mps', NaN, 'ebn0_db', NaN, 'true_offset_s', NaN, 'guard_s', NaN, ...
    'n_bits', 0, 'detected', false, 'velocity_hat_mps', NaN, 'velocity_error_mps', NaN, ...
    'preamble_error_ms', NaN, 'payload_error_ms', NaN, 'delay_spread_ms', NaN, ...
    'sync_score', NaN, 'errors_corrected', NaN, 'ber_corrected', NaN, ...
    'errors_uncorrected', NaN, 'ber_uncorrected', NaN, 'runtime_s', NaN, 'reason', '');
end

function writeJson(path, value)
file = fopen(path, 'w');
if file < 0, error('uac:WriteFailure', 'Cannot write %s.', path); end
cleanup = onCleanup(@() fclose(file)); %#ok<NASGU>
fwrite(file, jsonencode(value), 'char');
end
