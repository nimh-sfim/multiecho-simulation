import os.path as op
import sys
import numpy as np
from scipy import signal

sys.path.append(op.join(op.dirname(op.abspath(__file__)), "..", "simlib"))

import generate_timeseries as gt


def test_gen_randfreq_timeseries():
    n_reps = 3
    n_timepoints = 1000
    n_freq = 2
    seed = 42
    fs = 1
    timeseries = gt.gen_randfreq_timeseries(n_reps, n_timepoints, n_freq, seed=seed)
    assert timeseries.shape == (n_reps, n_timepoints)

    freqs, pspectrum = signal.periodogram(timeseries, fs=fs, axis=1)

    peak_vals = np.where(pspectrum > 100)
    assert np.allclose(peak_vals[0], np.array([0, 0, 0, 1, 2, 2, 2]))
    assert np.allclose(peak_vals[1], np.array([121, 167, 168, 32, 32, 37, 38]))


def test_gen_randn_timeseries():
    n_reps = 2
    n_timepoints = 3
    seed = 42
    timeseries = gt.gen_randn_timeseries(n_reps, n_timepoints, seed=seed)
    assert timeseries.shape == (n_reps, n_timepoints)

    expected_result = np.array(
        [[0.49671415, -0.1382643, 0.64768854], [1.52302986, -0.23415337, -0.23413696]]
    )
    assert np.allclose(timeseries, expected_result, atol=1e-4)


def test_gen_bandpass_randn_timeseries():
    n_reps = 2
    n_timepoints = 1000
    seed = 42
    passband = [0.1, 0.2]
    fs = 1.0
    timeseries = gt.gen_bandpass_randn_timeseries(
        n_reps, n_timepoints, passband=passband, fs=fs, seed=seed
    )
    assert timeseries.shape == (n_reps, n_timepoints)

    freqs, pspectrum = signal.periodogram(timeseries, fs=fs, axis=1)

    zero_bands = np.where((freqs < 0.1) | (freqs > 0.2))
    assert np.mean(pspectrum[:, zero_bands]) < 0.07
    pass_bands = np.where((freqs >= 0.1) & (freqs <= 0.2))
    assert np.mean(pspectrum[:, pass_bands]) > 1.8
