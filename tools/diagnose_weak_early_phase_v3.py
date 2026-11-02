"""Synthetic weak-first/strong-late phase-by-timing interaction, no retuning."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.signal import hilbert

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/v3/nof1_diagnostic/toy'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(2026110210)
    bits = rng.integers(0, 2, 1024)
    fs = 48000
    frequencies = np.array([13875., 14125.])
    phase = .7
    rows = []
    for on_ms in [4, 8]:
        n_on = round(on_ms/1000*fs)
        t = np.arange(n_on)/fs
        bank = np.exp(-2j*np.pi*t[:, None]*frequencies)
        local = 2*np.pi*frequencies[bits, None]*t+phase
        starts = np.r_[0., np.cumsum(frequencies[bits[:-1]]*n_on/fs)]
        signals = {'reset': np.cos(local).ravel(),
                   'continuous': np.cos(local+2*np.pi*(starts % 1)[:, None]).ravel()}
        pad = 2*n_on
        used = np.arange(8, len(bits)-8)
        truth = bits[used]
        for phase_mode, x in signals.items():
            analytic = hilbert(np.pad(x, (pad, pad)))
            for delay_ms in [.75, 1., 1.25]:
                delay = round(delay_ms/1000*fs)
                for timing, offset in [('first', 0), ('strong', delay)]:
                    start = pad+offset
                    direct_samples = analytic[start:start+len(bits)*n_on].real.reshape(-1, n_on)
                    late_samples = analytic[start-delay:start-delay+len(bits)*n_on].reshape(-1, n_on)
                    direct = (direct_samples @ bank)[used]
                    late_real = (late_samples.real @ bank)[used]
                    late_imag = (late_samples.imag @ bank)[used]
                    for early_amplitude in [.2, .4, .6, .8]:
                        for late_phase_index in range(8):
                            late_phase = late_phase_index*np.pi/4
                            late = np.cos(late_phase)*late_real-np.sin(late_phase)*late_imag
                            early = early_amplitude*direct
                            combined = early+late
                            energies = np.abs(combined)**2
                            errors = (energies[:, 1] > energies[:, 0]) != truth
                            correct_early = early[np.arange(len(truth)), truth]
                            correct_late = late[np.arange(len(truth)), truth]
                            isolated = np.abs(correct_early)**2+np.abs(correct_late)**2
                            cross = 2*np.real(correct_early*np.conj(correct_late))
                            rows.append({'on_ms': on_ms, 'phase_mode': phase_mode, 'timing': timing,
                                         'late_delay_ms': delay_ms, 'early_amplitude': early_amplitude,
                                         'late_amplitude': 1., 'late_phase_pi': late_phase/np.pi,
                                         'bits': len(truth), 'errors': int(errors.sum()), 'ber': errors.mean(),
                                         'mean_correct_tone_cross_over_isolated_energy': cross.sum()/isolated.sum()})
    data = pd.DataFrame(rows)
    data.to_csv(OUT/'weak_early_phase_timing.csv', index=False)
    keys = ['on_ms', 'late_delay_ms', 'early_amplitude', 'late_phase_pi']
    wide = data.pivot(index=keys, columns=['phase_mode', 'timing'], values='errors')
    wide.columns = ['_'.join(x) for x in wide.columns]
    wide = wide.reset_index()
    wide['reset_first_minus_strong'] = wide.reset_first-wide.reset_strong
    wide['first_reset_minus_continuous'] = wide.reset_first-wide.continuous_first
    wide.to_csv(OUT/'weak_early_phase_timing_paired.csv', index=False)
    eight = wide[wide.on_ms == 8]
    assert (eight.reset_first == eight.continuous_first).all()
    assert (eight.reset_strong == eight.continuous_strong).all()
    four = wide[wide.on_ms == 4]
    example = four[(four.first_reset_minus_continuous > 0) & (four.reset_first_minus_strong > 0)].sort_values(
        ['first_reset_minus_continuous', 'reset_first_minus_strong'], ascending=False)
    result = {'status': 'PASS', 'seed': 2026110210, 'bits_per_condition': int(data.bits.iloc[0]),
              'rows': len(data), 'four_ms_path_configurations': len(four),
              'four_ms_reset_first_worse_than_strong': int((four.reset_first_minus_strong > 0).sum()),
              'four_ms_continuous_first_better_than_reset': int((four.first_reset_minus_continuous > 0).sum()),
              'eight_ms_reset_continuous_decisions_equal': True,
              'example_selection': 'Largest synthetic error-count difference within the declared 3 delays x 4 gains x 8 phases; illustrative only, not fitted to NOF1.',
              'examples': example.head(3).to_dict(orient='records'),
              'scope': 'Two synthetic complex-gain rays, no noise, exact known timing choices and true unit scale. This tests an interaction mechanism; it is not a fit to the measured NOF1 impulse response.'}
    (OUT/'weak_early_phase_timing_validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))
    if len(example):
        row = example.iloc[0]
        selected = data[(data.on_ms == row.on_ms) & (data.late_delay_ms == row.late_delay_ms) &
                        (data.early_amplitude == row.early_amplitude) & (data.late_phase_pi == row.late_phase_pi)]
        print(selected.to_string(index=False))


if __name__ == '__main__':
    main()
