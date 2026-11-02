function records = run_paper_alignment_v2(stage)
%RUN_PAPER_ALIGNMENT_V2 Disclosed two-noise comparison with paper Tables 2/3.
% All records are development records 001:003. No held-out files are opened.
% The receiver observes an entire 6.4 s slot and receives no channel truth.
if nargin < 1, stage = 'full'; end
root = fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root, 'matlab'));
out = fullfile(root, 'results', 'v2', 'paper_alignment');
if strcmp(stage, 'smoke'), out = fullfile(out, 'smoke'); end
if ~exist(out, 'dir'), mkdir(out); end
validation = jsondecode(fileread(fullfile(root, 'results', 'receiver_v2', 'validation.json')));
assert(strcmp(validation.status, 'PASS'), 'Receiver validation must pass before this run.');
cfgBase = uac2.default_config(); fs = cfgBase.fs;
channels = {'NOF1', 'NCS1'}; recordIds = 1:3; onMs = [4, 8, 16, 32];
nSlots = 5; nBits = 64; slotSeconds = 6.4; slotSamples = round(slotSeconds*fs);
snrLevels = [0, 10, 20]; noiseModes = {'sample_snr_packet', 'ebn0_received'};
seed = 2026110205; tailSeconds = .128;
if strcmp(stage, 'smoke')
    channels = channels(1); recordIds = 1; onMs = 8; nSlots = 1;
elseif ~strcmp(stage, 'full')
    error('stage must be full or smoke');
end
metadata = struct('status', 'RUNNING', 'stage', stage, 'seed', seed);
rows = {}; prior = table(); priorPath = fullfile(out, 'per_packet.csv');
if exist(priorPath, 'file')
    prior = readtable(priorPath, 'TextType', 'char');
    assert(all(ismember(prior.record, recordIds)) && all(prior.payload_bits == nBits), 'Unexpected checkpoint scope.');
    assert(all(ismember(prior.on_ms, onMs)), 'Unexpected checkpoint symbol durations.');
    rows = num2cell(table2struct(prior));
    fprintf('Resuming from %d completed packet observations.\n', height(prior));
end
metadata.resumed_from_rows = numel(rows);
writeJson(fullfile(out, 'metadata.json'), metadata);
started = tic;
for ci = 1:numel(channels)
    channel = channels{ci};
    for recordNumber = recordIds
        expectedPerConfiguration = nSlots*numel(snrLevels)*numel(noiseModes);
        if ~isempty(prior)
            completedRecord = strcmp(prior.channel, channel) & prior.record == recordNumber;
            if sum(completedRecord) == numel(onMs)*2*expectedPerConfiguration, continue; end
        end
        file = fullfile(root, 'vendor', 'Watermark', 'input', 'channels', channel, 'mat', sprintf('%s_%03d.mat', channel, recordNumber));
        data = load(file);
        assert(nSlots*slotSeconds < (size(data.h, 1)-1)/data.fs_t, 'Sounding too short.');
        for durationMs = onMs
            for guardMs = [0, durationMs]
                if ~isempty(prior)
                    completed = strcmp(prior.channel, channel) & prior.record == recordNumber & prior.on_ms == durationMs & prior.guard_ms == guardMs;
                    assert(sum(completed) == 0 || sum(completed) == expectedPerConfiguration, 'Partial configuration checkpoint; inspect before resuming.');
                    if sum(completed) == expectedPerConfiguration, continue; end
                end
                cfg = cfgBase; cfg.symbol_s = durationMs/1000; cfg.guard_s = guardMs/1000;
                tx = zeros(nSlots*slotSamples, 1); payloads = cell(nSlots, 1);
                jitters = zeros(nSlots, 1); frameDurations = zeros(nSlots, 1); phases = zeros(nSlots, 1);
                for slot = 1:nSlots
                    rng(seed+ci*10000000+recordNumber*10000+slot, 'twister');
                    bits = randi([0, 1], nBits, 1); payloads{slot} = bits;
                    jitter = 3*round((.12+.22*rand())*fs/3); jitters(slot) = jitter;
                    txCfg = cfg; txCfg.phase_rad = 2*pi*rand(); phases(slot) = txCfg.phase_rad;
                    packet = uac2.make_packet(bits, txCfg);
                    frameDurations(slot) = numel(packet)/fs+tailSeconds;
                    assert(jitter+numel(packet)+round(tailSeconds*fs) < slotSamples, 'Frame exceeds slot.');
                    first = (slot-1)*slotSamples+jitter;
                    tx(first+(1:numel(packet))) = packet;
                end
                cleanTrain = uac.replayfilter_fast(tx, fs, data.h, data.fs_t, double(data.fs_tau), double(data.fc));
                for slot = 1:nSlots
                    clean = cleanTrain((slot-1)*slotSamples+(1:slotSamples));
                    [clean, alpha] = uac2.restore_bulk(clean, data.V0);
                    % All slot energy includes any measured-channel tail. The
                    % sample-power denominator is ONLY the declared packet +
                    % 128 ms tail, stretched by alpha, not random/listen padding.
                    energySamples = sum(clean.^2);
                    frameSamples = ceil(frameDurations(slot)*alpha*fs);
                    bits = payloads{slot};
                    for snrDb = snrLevels
                        rng(seed+ci*10000000+recordNumber*10000+slot*100+snrDb+50, 'twister');
                        unitNoise = randn(size(clean));
                        for ni = 1:numel(noiseModes)
                            if ni == 1, variance = energySamples/frameSamples/10^(snrDb/10);
                            else, variance = energySamples/(2*nBits*10^(snrDb/10)); end
                            result = uac2.receive_packet(clean+sqrt(variance)*unitNoise, cfg, nBits, struct('compensate', true));
                            detected = result.detected && numel(result.bits_hat_compensated) == nBits && all(isfinite(result.bits_hat_compensated));
                            errors = 0; erasures = nBits; uncompensatedErrors = NaN;
                            if detected
                                errors = sum(result.bits_hat_compensated(:) ~= bits); erasures = 0;
                                if numel(result.bits_hat_uncompensated) == nBits
                                    uncompensatedErrors = sum(result.bits_hat_uncompensated(:) ~= bits);
                                end
                            end
                            row = struct('channel', channel, 'record', recordNumber, 'slot', slot, ...
                                'on_ms', durationMs, 'guard_ms', guardMs, 'snr_label_db', snrDb, ...
                                'noise_mode', noiseModes{ni}, 'receiver_mode', 'sync_compensated', ...
                                'payload_bits', nBits, 'detected', double(detected), 'bit_errors', errors, ...
                                'bit_erasures', erasures, 'packet_error', double(~detected || errors > 0), ...
                                'uncompensated_bit_errors', uncompensatedErrors, 'noise_variance', variance, ...
                                'received_energy_samples', energySamples, 'power_reference_samples', frameSamples, ...
                                'frame_with_tail_s', frameDurations(slot), 'slot_s', slotSeconds, ...
                                'alpha_true', alpha, 'alpha_hat', result.alpha_hat, ...
                                'velocity_true_mps', double(data.V0), 'velocity_hat_mps', result.velocity_hat, ...
                                'sync_score', result.sync_score, 'delay_spread_ms', result.delay_spread_ms, ...
                                'injected_start_s', jitters(slot)/fs, 'carrier_phase_rad', phases(slot), ...
                                'receiver_reason', result.diagnostic.reason);
                            rows{end+1} = row; %#ok<AGROW>
                        end
                    end
                end
                records = struct2table(vertcat(rows{:}));
                writetable(records, fullfile(out, 'per_packet.csv'));
                fprintf('%s %03d Ton=%d Tg=%d: %d rows, %.1fs\n', channel, recordNumber, durationMs, guardMs, height(records), toc(started));
            end
        end
    end
end
records = struct2table(vertcat(rows{:}));
runFftDiagnostic(out, seed, strcmp(stage, 'smoke'));
metadata.status = 'COMPLETE'; metadata.channels = channels; metadata.record_ids = recordIds;
metadata.config = cfgBase; metadata.payload_bits = nBits; metadata.active_ms = onMs;
metadata.guard_rule = '0 or active duration'; metadata.snr_labels_db = snrLevels; metadata.noise_modes = noiseModes;
metadata.slots_per_record = nSlots; metadata.slot_s = slotSeconds; metadata.tail_s = tailSeconds;
metadata.receiver_validation_status = validation.status; metadata.receiver_required_passed = validation.required_passed;
metadata.frame_energy = 'All received clean slot energy includes HFM, new known-tone pilots, payload and multipath tail.';
metadata.sample_snr_packet = 'variance=received_energy_samples/ceil(alpha*(packet_duration+0.128)*fs)/10^(label/10); random listening padding excluded from denominator';
metadata.ebn0_received = 'variance=received_energy_samples/(2*64*10^(label/10)); preamble and pilot energy included per information bit';
metadata.noise_pairing = 'Same standard-normal vector for both noise definitions and all active/guard choices at matched channel/record/slot/label';
metadata.phase = 'One uniform random phase per packet, symbol-reset modulation; phase/payload/jitter matched across active/guard choices';
metadata.capture_loss = 'Receiver reports no complete decoded payload; one packet per slot. This is not a separate verified timing-association metric.';
metadata.truth_isolation = 'Truth only in waveform generation, channel/noise scaling and scoring; receive_packet sees rx, cfg and 64';
metadata.scope = 'Development records 001-003 only; descriptive alignment sensitivity, no exact-replication or independent-test claim';
metadata.elapsed_this_process_s = toc(started); metadata.matlab_version = version;
metadata.resume = 'Deterministic per-configuration CSV checkpoints; completed conditions are preserved across process restarts';
writeJson(fullfile(out, 'metadata.json'), metadata);
fprintf('PAPER ALIGNMENT COMPLETE %s: %d packets, %.1fs\n', stage, height(records), toc(started));
end

function runFftDiagnostic(out, seed, isSmoke)
% Exact correlations and N-point rectangular FFT see IDENTICAL real samples.
fs = 48000; placements = [13875, 14125; 13750, 14000; 14000, 14250];
nBits = 20000; if isSmoke, nBits = 2000; end
chunk = 1000; rows = {};
for onMs = [4, 8]
    N = round(onMs/1000*fs); time = (0:N-1)'/fs;
    for placement = 1:3
        frequencies = placements(placement, :); bins = round(frequencies*N/fs);
        bank = exp(-2*pi*1i*time*frequencies);
        for ebn0 = [0, 10, 20]
            rng(seed+1000+onMs*100+ebn0, 'twister'); errors = [0, 0]; ties = 0;
            for first = 1:chunk:nBits
                count = min(chunk, nBits-first+1); bits = randi([0, 1], 1, count); phases = 2*pi*rand(1, count);
                tone = frequencies(bits+1);
                clean = cos(bsxfun(@plus, 2*pi*time*tone, phases));
                variance = N/(4*10^(ebn0/10));
                received = clean+sqrt(variance)*randn(size(clean));
                exact = abs(bank.'*received).^2;
                spectrum = fft(received, N, 1); nearest = abs(spectrum(bins+1, :)).^2;
                errors(1) = errors(1)+sum((exact(2, :) > exact(1, :)) ~= bits);
                errors(2) = errors(2)+sum((nearest(2, :) > nearest(1, :)) ~= bits);
                ties = ties+sum(abs(nearest(2, :)-nearest(1, :)) < 1e-10);
            end
            row = struct('on_ms', onMs, 'f0_hz', frequencies(1), 'f1_hz', frequencies(2), ...
                'fft_samples', N, 'nearest_bin0', bins(1), 'nearest_bin1', bins(2), ...
                'ebn0_db', ebn0, 'bits', nBits, 'exact_errors', errors(1), 'fft_errors', errors(2), ...
                'exact_ber', errors(1)/nBits, 'fft_ber', errors(2)/nBits, 'numerical_ties', ties);
            rows{end+1} = row; %#ok<AGROW>
        end
    end
end
writetable(struct2table(vertcat(rows{:})), fullfile(out, 'fft_awgn.csv'));
writeJson(fullfile(out, 'fft_metadata.json'), struct('status', 'COMPLETE', 'bits_per_condition', nBits, ...
    'frequencies_hz', placements, 'window', 'rectangular', 'zero_padding', false, ...
    'frequency_rule', 'positive real-passband DFT bins; nearest integer via MATLAB round, half ties toward +infinity', ...
    'decision_rule', 'one bin per tone, maximum squared magnitude; decision ties choose zero', ...
    'waveform_pairing', 'Same received waveform for exact-tone and FFT detectors; same bits, phases and standard noise across frequency placements', ...
    'scope', 'AWGN detector sensitivity only; no Watermark retuning and not the unspecified paper narrow-band integration detector'));
end

function writeJson(path, value)
fid = fopen(path, 'w'); assert(fid >= 0, 'Cannot write metadata.');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(value), 'char');
end
