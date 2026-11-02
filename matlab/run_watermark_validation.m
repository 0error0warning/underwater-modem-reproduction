function summary = run_watermark_validation(nSoundings, runTag)
%RUN_WATERMARK_VALIDATION Declared, finite BFSK reproduction on measured TVIRs.
% Known packet boundaries and mean-PDP peak delay; bulk Doppler removed.
% No synchronization claim. Three noise conventions are kept separate.
if nargin < 1, nSoundings = 10; end
if nargin < 2, runTag = sprintf('first_%02d_soundings', nSoundings); end
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'matlab'));
addpath(fullfile(projectRoot, 'vendor', 'Watermark', 'matlab'));
out = fullfile(projectRoot, 'results', 'watermark', runTag);
if ~exist(out, 'dir'), mkdir(out); end
fid=fopen(fullfile(out,'run_metadata.json'),'w');
fwrite(fid,jsonencode(struct('status','RUNNING','run_tag',runTag)),'char'); fclose(fid);
channels = {'NOF1', 'NCS1'};
onTimes = [8, 4, 16, 32]; % 8 ms/no-guard first: common noise calibration.
noiseModes = {'ebn0_received', 'sample_snr_packet', 'fixed_noise_active_reference'};
levels = [0, 10, 20];
nBits = 128;
fs = 48000;
baseSeed = 20261102;
useFast = false;
proofPath = fullfile(projectRoot, 'results', 'replay_equivalence.json');
if exist(proofPath, 'file') && exist(fullfile(projectRoot,'matlab','+uac','replayfilter_fast.m'), 'file')
    proof = jsondecode(fileread(proofPath));
    useFast = isfield(proof, 'status') && strcmp(proof.status, 'PASS');
end
fprintf('Verified fast replay kernel enabled: %d\n', useFast);
packetRows = struct([]);
fileRows = struct([]);
rowIndex = 0;
fileIndex = 0;
started = tic;
for ci = 1:numel(channels)
    channel = channels{ci};
    for fileNumber = 1:nSoundings
        fileName = sprintf('%s_%03d.mat', channel, fileNumber);
        channelPath = fullfile(projectRoot, 'vendor', 'Watermark', 'input', ...
            'channels', channel, 'mat', fileName);
        assert(exist(channelPath, 'file') == 2, 'Missing channel file: %s', channelPath);
        data = load(channelPath);
        K = size(data.h, 2);
        powerDelay = mean(abs(data.h).^2, 1);
        [~, peakIndex] = max(powerDelay);
        delaySamples = round((peakIndex-1)/double(data.fs_tau)*fs);
        cumulative = cumsum(powerDelay)/sum(powerDelay);
        p05 = find(cumulative >= .05, 1);
        p95 = find(cumulative >= .95, 1);
        delay90ms = (p95-p05)/double(data.fs_tau)*1000;
        maxSamples = floor((size(data.h,1)-1)/data.fs_t*fs);
        lead = ceil(.002*fs);
        tail = ceil(K/double(data.fs_tau)*fs) + lead;
        calibrationPower = NaN;
        for oi = 1:numel(onTimes)
            for guardRatio = [0, 1]
                cfg = uac.default_config();
                cfg.symbol_s = onTimes(oi)/1000;
                cfg.guard_s = cfg.symbol_s*guardRatio;
                nActive = round(cfg.symbol_s*fs);
                nGuard = round(cfg.guard_s*fs);
                coreLength = nBits*(nActive+nGuard);
                packetLength = lead+coreLength+tail;
                nPackets = floor(maxSamples/packetLength);
                assert(nPackets > 0, 'Packet cannot fit in sounding.');
                % Identical per-packet bits across configurations; changing
                % duration changes placement. Store positions explicitly.
                bits = zeros(nBits, nPackets);
                train = zeros(nPackets*packetLength, 1);
                for pi = 1:nPackets
                    rng(baseSeed + ci*100000 + fileNumber*1000 + pi, 'twister');
                    bits(:,pi) = randi([0,1], nBits, 1);
                    first = (pi-1)*packetLength + lead;
                    train(first+(1:coreLength)) = uac.modulate(bits(:,pi), cfg);
                end
                if useFast
                    received = uac.replayfilter_fast(train, fs, data.h, ...
                        data.fs_t, double(data.fs_tau), double(data.fc));
                else
                    received = replayfilter(train, fs, data.h, ...
                        data.fs_t, double(data.fs_tau), double(data.fc));
                end
                if oi == 1 && guardRatio == 0
                    % Same fixed input amplitude; one receive-power reference
                    % for every duration/guard in this channel recording.
                    calibrationPower = sum(received.^2)/(nPackets*nBits*nActive);
                end
                errors = zeros(numel(noiseModes), numel(levels));
                noiseFreeErrors = 0;
                for pi = 1:nPackets
                    first = (pi-1)*packetLength;
                    clean = received(first+(1:packetLength));
                    offset = lead + delaySamples;
                    cleanHat = uac.demodulate(clean, cfg, nBits, offset);
                    cleanErrors = sum(cleanHat ~= bits(:,pi));
                    noiseFreeErrors = noiseFreeErrors + cleanErrors;
                    for mi = 1:numel(noiseModes)
                        for li = 1:numel(levels)
                            gamma = 10^(levels(li)/10);
                            switch noiseModes{mi}
                                case 'ebn0_received'
                                    variance = sum(clean.^2)/(2*nBits*gamma);
                                case 'sample_snr_packet'
                                    % Includes all guard, lead and channel tail.
                                    variance = mean(clean.^2)/gamma;
                                case 'fixed_noise_active_reference'
                                    variance = calibrationPower/gamma;
                            end
                            % Shared draws across conventions and durations
                            % where lengths match. No search against payload.
                            rng(baseSeed + ci*10000000 + fileNumber*100000 + pi*100 + li, 'twister');
                            noisy = clean + sqrt(variance)*randn(size(clean));
                            decoded = uac.demodulate(noisy, cfg, nBits, offset);
                            nErrors = sum(decoded ~= bits(:,pi));
                            errors(mi,li) = errors(mi,li)+nErrors;
                            rowIndex = rowIndex+1;
                            packetRows(rowIndex).channel = channel;
                            packetRows(rowIndex).sounding = fileNumber;
                            packetRows(rowIndex).on_ms = onTimes(oi);
                            packetRows(rowIndex).guard_ms = onTimes(oi)*guardRatio;
                            packetRows(rowIndex).noise_mode = noiseModes{mi};
                            packetRows(rowIndex).level_db = levels(li);
                            packetRows(rowIndex).packet = pi;
                            packetRows(rowIndex).start_s = first/fs;
                            packetRows(rowIndex).bits = nBits;
                            packetRows(rowIndex).errors = nErrors;
                            packetRows(rowIndex).packet_error = double(nErrors>0);
                            packetRows(rowIndex).noise_variance = variance;
                            packetRows(rowIndex).noiseless_errors = cleanErrors;
                        end
                    end
                end
                for mi = 1:numel(noiseModes)
                    for li = 1:numel(levels)
                        countForConfiguration = nPackets*numel(noiseModes)*numel(levels);
                        selected = packetRows(rowIndex-countForConfiguration+1:rowIndex);
                        selected = selected(strcmp({selected.noise_mode}, noiseModes{mi}) & [selected.level_db] == levels(li));
                        fileIndex = fileIndex+1;
                        fileRows(fileIndex).channel = channel;
                        fileRows(fileIndex).sounding = fileNumber;
                        fileRows(fileIndex).on_ms = onTimes(oi);
                        fileRows(fileIndex).guard_ms = onTimes(oi)*guardRatio;
                        fileRows(fileIndex).noise_mode = noiseModes{mi};
                        fileRows(fileIndex).level_db = levels(li);
                        fileRows(fileIndex).bits = nPackets*nBits;
                        fileRows(fileIndex).errors = errors(mi,li);
                        fileRows(fileIndex).ber = errors(mi,li)/(nPackets*nBits);
                        fileRows(fileIndex).packets = nPackets;
                        fileRows(fileIndex).packet_errors = sum([selected.packet_error]);
                        fileRows(fileIndex).per = sum([selected.packet_error])/nPackets;
                        fileRows(fileIndex).noiseless_errors = noiseFreeErrors;
                        fileRows(fileIndex).nominal_bit_rate = 1/(cfg.symbol_s+cfg.guard_s);
                        fileRows(fileIndex).payload_rate_with_padding = nBits/(packetLength/fs);
                        fileRows(fileIndex).error_free_packet_goodput = nBits*sum([selected.packet_error]==0)/(nPackets*packetLength/fs);
                        fileRows(fileIndex).delay_peak_ms = delaySamples/fs*1000;
                        fileRows(fileIndex).delay_90_ms = delay90ms;
                        fileRows(fileIndex).bulk_velocity_mps = data.V0;
                        fileRows(fileIndex).calibration_power = calibrationPower;
                    end
                end
                fprintf('%s on=%d ms guard=%d ms: %d bits; clean BER=%.5g; elapsed=%.1fs\n', ...
                    fileName,onTimes(oi),onTimes(oi)*guardRatio,nPackets*nBits,noiseFreeErrors/(nPackets*nBits),toc(started));
            end
        end
        % Recoverable progress after each complete measured recording.
        writetable(struct2table(fileRows), fullfile(out,'per_sounding.csv'));
    end
end
writetable(struct2table(packetRows), fullfile(out,'per_packet.csv'));
summaryRows = struct([]);
si=0;
for ci=1:numel(channels)
    for onTime=sort(onTimes)
        for guardRatio=[0,1]
            for mi=1:numel(noiseModes)
                for li=1:numel(levels)
                    mask = strcmp({fileRows.channel},channels{ci}) & [fileRows.on_ms]==onTime & ...
                        [fileRows.guard_ms]==onTime*guardRatio & strcmp({fileRows.noise_mode},noiseModes{mi}) & [fileRows.level_db]==levels(li);
                    rows=fileRows(mask);
                    totalBits=sum([rows.bits]); totalErrors=sum([rows.errors]);
                    si=si+1;
                    summaryRows(si).channel=channels{ci};
                    summaryRows(si).on_ms=onTime;
                    summaryRows(si).guard_ms=onTime*guardRatio;
                    summaryRows(si).noise_mode=noiseModes{mi};
                    summaryRows(si).level_db=levels(li);
                    summaryRows(si).soundings=numel(rows);
                    summaryRows(si).bits=totalBits;
                    summaryRows(si).errors=totalErrors;
                    summaryRows(si).ber=totalErrors/totalBits;
                    summaryRows(si).per=sum([rows.packet_errors])/sum([rows.packets]);
                    summaryRows(si).min_sounding_ber=min([rows.ber]);
                    summaryRows(si).max_sounding_ber=max([rows.ber]);
                    summaryRows(si).noiseless_ber=sum([rows.noiseless_errors])/totalBits;
                    summaryRows(si).nominal_bit_rate=rows(1).nominal_bit_rate;
                    summaryRows(si).payload_rate_with_padding=rows(1).payload_rate_with_padding;
                    summaryRows(si).error_free_packet_goodput=sum([rows.error_free_packet_goodput].*[rows.packets])/sum([rows.packets]);
                end
            end
        end
    end
end
summary=struct2table(summaryRows);
writetable(summary,fullfile(out,'summary.csv'));
metadata.status='COMPLETE_LIMITED_SCOPE';
metadata.n_soundings_per_channel=nSoundings;
metadata.channels=channels;
metadata.matlab_version=version;
metadata.seed=baseSeed;
metadata.tones_hz=[13875,14125];
metadata.fs=fs;
metadata.packet_bits=nBits;
metadata.phase='Carrier phase reset to zero for each active symbol; implementation assumption';
metadata.timing='Known transmit boundary plus peak of mean measured power-delay profile; no synchronizer';
metadata.doppler='Bulk V0 remains removed, as in supplied TVIRs; time-varying residual channel retained. Not the default sfetch Doppler reinstatement.';
metadata.detector='Exact-tone quadrature matched energies, not the unspecified paper FFT-band detector';
metadata.noise_modes=noiseModes;
metadata.fixed_noise_reference='Received energy divided by active samples in 8 ms/no-guard reference, per sounding; fixed across guard and duration';
metadata.packet_placement='Each configuration fills the record with its own packet duration; endpoints and packet times differ. Not paired packet observations.';
metadata.goodput='Truth-known error-free 128-bit packets divided by transmit+guard+padding time; idealized discard model without CRC, preamble, ACK or retransmission';
metadata.uncertainty='Per-sounding min/max are descriptive, not confidence intervals. Adjacent soundings and bits may be dependent.';
metadata.used_verified_fast_kernel=useFast;
metadata.elapsed_seconds=toc(started);
fid=fopen(fullfile(out,'run_metadata.json'),'w');
fwrite(fid,jsonencode(metadata),'char'); fclose(fid);
fprintf('COMPLETE: %s (%.1fs)\n',out,toc(started));
end
