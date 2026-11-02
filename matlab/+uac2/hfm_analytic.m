function signal = hfm_analytic(time, firstHz, lastHz, duration, taperFraction)
%HFM_ANALYTIC Analytic hyperbolic chirp evaluated at arbitrary transmit times.
%   Instantaneous frequency varies from firstHz to lastHz. A symmetric Tukey
%   envelope (implemented explicitly) limits correlation sidelobes.
time = time(:);
signal = complex(zeros(size(time)));
inside = time >= 0 & time < duration;
t = time(inside);
if firstHz == lastHz
    phase = 2*pi*firstHz*t;
else
    scale = firstHz*lastHz*duration/(firstHz-lastHz);
    phase = 2*pi*scale*log1p((firstHz-lastHz)*t/(lastHz*duration));
end
envelope = ones(size(t));
edge = taperFraction*duration/2;
if edge > 0
    left = t < edge;
    right = t > duration-edge;
    envelope(left) = 0.5*(1-cos(pi*t(left)/edge));
    envelope(right) = 0.5*(1-cos(pi*(duration-t(right))/edge));
end
signal(inside) = envelope .* exp(1i*phase);
end
