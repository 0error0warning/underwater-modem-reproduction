function [y, alpha] = restore_bulk(x, velocity)
%RESTORE_BULK Reinstate Watermark's removed bulk Doppler, as in sfetch.
% Positive velocity means increasing range and a longer received waveform.
alpha = 1/(1-velocity/1500);
[p,q] = rat(alpha);
persistent keys filters
if isempty(keys), keys={}; filters={}; end
key=sprintf('%d/%d',p,q);
index=find(strcmp(keys,key),1);
if isempty(index)
    [y,b]=resample(x,p,q);
    keys{end+1}=key; filters{end+1}=b;
else
    y=resample(x,p,q,filters{index});
end
end
