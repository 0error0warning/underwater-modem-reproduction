function result = validate_replay_at()
%VALIDATE_REPLAY_AT Compare local absolute-time replay to a long zero prefix.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab'));
addpath(fullfile(root,'vendor','Watermark','matlab'));
proof = jsondecode(fileread(fullfile(root,'results','replay_equivalence.json')));
assert(strcmp(proof.status,'PASS'),'Validated reference kernel proof is required.');
rng(20261104,'twister');
fs = 48000;
threshold = 1e-10;
starts = [0,5.0125];
rows = struct([]); index = 0;
for kind = 1:3
    if kind == 1
        fsTau = 16000; K = 37; fsT = fsTau/K; fc = 14000;
        h = (randn(2250,K)+1i*randn(2250,K))/sqrt(2*K);
        % Non-multiple-of-three input exercises resampling length rounding;
        % non-power-of-two K and off-block starts exercise both partial ends.
        x = randn(7001,1);
        source = 'Synthetic complex TVIR; no leading padding; partial final block';
        name = 'synthetic_complex_partial';
    else
        names = {'NOF1','NCS1'};
        channel = names{kind-1};
        path = fullfile(root,'vendor','Watermark','input','channels',channel,'mat',[channel '_001.mat']);
        data = load(path);
        h = data.h; fsTau = double(data.fs_tau); fsT = data.fs_t; fc = double(data.fc);
        x = [zeros(96,1); randn(fs,1); zeros(96,1)];
        source = path;
        name = [channel '_001_one_second'];
    end
    for start = starts
        prefix = round(start*fs);
        fullInput = [zeros(prefix,1); x];
        tic; fullOutput = uac.replayfilter_fast(fullInput,fs,h,fsT,fsTau,fc); fullTime=toc;
        reference = fullOutput(prefix+1:end);
        tic; actual = uac2.replay_at(x,fs,h,fsT,fsTau,fc,start); localTime=toc;
        assert(numel(reference)==numel(actual),'Output length mismatch.');
        relative = norm(actual-reference)/max(norm(reference),realmin);
        absolute = max(abs(actual-reference));
        fullTimes = zeros(1,3); localTimes = zeros(1,3);
        for repeat=1:3
            tic; uac.replayfilter_fast(fullInput,fs,h,fsT,fsTau,fc); fullTimes(repeat)=toc;
            tic; uac2.replay_at(x,fs,h,fsT,fsTau,fc,start); localTimes(repeat)=toc;
        end
        index=index+1;
        rows(index).name=name;
        rows(index).source=source;
        rows(index).start_s=start;
        rows(index).input_samples=numel(x);
        rows(index).output_samples=numel(actual);
        rows(index).prefix_samples=prefix;
        rows(index).relative_l2_error=relative;
        rows(index).max_absolute_error=absolute;
        rows(index).reference_first_seconds=fullTime;
        rows(index).local_first_seconds=localTime;
        rows(index).reference_median_seconds=median(fullTimes);
        rows(index).local_median_seconds=median(localTimes);
        rows(index).speedup=median(fullTimes)/median(localTimes);
        rows(index).passed=relative<threshold;
        fprintf('%s start %.4f: relative %.3g, max %.3g, speedup %.2fx\n', ...
            name,start,relative,absolute,rows(index).speedup);
    end
end
% An independent short direct-kernel comparison also checks the composition
% of equivalence proofs without paying for a five-second direct zero prefix.
rng(20261105,'twister');
h = (randn(40,37)+1i*randn(40,37))/sqrt(74);
x = randn(1201,1);
start = 0.010125; prefix=round(start*fs);
direct = replayfilter([zeros(prefix,1);x],fs,h,16000/37,16000,14000);
actual = uac2.replay_at(x,fs,h,16000/37,16000,14000,start);
direct = direct(prefix+1:end);
directRelative=norm(actual-direct)/max(norm(direct),realmin);
offGridRejected=false;
try
    uac2.replay_at(x,fs,h,16000/37,16000,14000,0.010126);
catch exception
    offGridRejected=strcmp(exception.identifier,'uac2:OffGridStart');
end
passed=all([rows.passed]) && directRelative<threshold && offGridRejected;
if passed, result.status='PASS'; else, result.status='FAIL'; end
result.threshold=threshold;
result.matlab_version=version;
result.results=rows;
result.direct_kernel_relative_l2_error=directRelative;
result.off_grid_start_rejected=offGridRejected;
result.scope='Zero input before each local packet; absolute channel interpolation and carrier phase; integer fs_x/fs_tau; no bulk Doppler reinstatement in this operator.';
result.reference='Zero prefix plus independently validated uac.replayfilter_fast, with an additional original direct-kernel check.';
result.timestamp_utc=char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd''T''HH:mm:ss''Z'''));
fid=fopen(fullfile(root,'results','replay_at_validation.json'),'w');
fwrite(fid,jsonencode(result),'char'); fclose(fid);
assert(passed,'uac2:ReplayAtValidationFailed','Absolute-time replay failed validation.');
fprintf('REPLAY_AT %s; direct-kernel relative error %.3g\n',result.status,directRelative);
end
