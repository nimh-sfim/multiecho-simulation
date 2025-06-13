""" Generate time series with different characteristics to use as input to generate_sims at different echoes"""

from typing import List

import numpy as np
import numpy.typing as npt
from scipy import signal


def gen_randfreq_timeseries(n_reps: int, n_timepoints: int, n_freq: int) -> npt.NDArray:
    """
    time series that is a combination of sine waves with uniform random freq & phase and amplitude

    This should be similar to band pass filtering random time series, but with a bit more control over frequencies

    Parameters
    ----------
    n_reps: (int) Number of repetitions (voxels) to create time series for
    n_timepoints: (int) Number of time points for each repetition
    n_freq: (int) number of sin waves to

    Returns
    -------
    A 2D matrix n_reps time series of n_timepoints durations.
    Each time series has is a combination of n_freq sinusoids each with a random frequency (between 1 cycle for the entire time course & 1 cycle every 5 timepoints),
    a random phase, and a random amplitude (0.1-1.0).
    After the sinusoids are added together they are scaled to have a standard deviation of 1.

    Note
    ----
    This function is used as one of several options for generating simulated time series and should be fairly similar to bandpass filtering random noise (another option)
    """

    freq_random = np.zeros((n_reps, n_timepoints))
    time_scale = np.matlib.repmat(
        np.linspace(0, 2 * np.pi, num=n_timepoints), n_reps, 1
    )  # one cycle of a sin wave with freq=1
    min_freq = 1
    max_freq = 5 / n_timepoints  # a cycle repeating every 5 time points
    for freq_idx in range(n_freq):
        random_amp = np.matlib.repmat(
            0.1 + np.random.rand(n_reps) / 0.9, n_timepoints, 1
        ).T
        random_freq = np.matlib.repmat(
            min_freq + np.random.rand(n_reps) / max_freq, n_timepoints, 1
        ).T
        random_phase = np.matlib.repmat(
            2 * np.pi * np.random.rand(n_reps), n_timepoints, 1
        ).T
        freq_random = freq_random + random_amp * np.sin(
            random_freq * time_scale + random_phase
        )
    tmp_mean = np.matlib.repmat(np.mean(freq_random, axis=1), n_timepoints, 1).T
    tmp_std = np.matlib.repmat(np.std(freq_random, axis=1), n_timepoints, 1).T
    return (freq_random - tmp_mean) / tmp_std


def gen_randn_timeseries(n_reps: int, n_vals: int) -> npt.NDArray:
    """
    time series that are Gaussian random noise

    Parameters
    ----------
    n_reps: (int) Number of repetitions (voxels) to create time series for
    n_timepoints: (int) Number of time points for each repetition

    Returns
    -------
    A 2D matrix n_reps time series of n_timepoints durations of Gaussian random noise
    """
    return np.random.randn(n_reps, n_vals)


def gen_bandpass_randn_timeseries(
    n_reps: int,
    n_timepoints: int,
    passband: List[float] = [1 / 100, 1 / 20],
    fs: float = 0.5,
) -> npt.NDArray:
    """
    time series that are Gaussian random noise and then bandpass filtered

    Parameters
    ----------
    n_reps: (int) Number of repetitions (voxels) to create time series for
    n_timepoints: (int) Number of time points for each repetition
    passband: List(float) Two numbers that are the top and bottom frequencies of the pass band
    fs: (float) The sampling frequency of the digital system

    Returns
    -------
    A 2D matrix n_reps time series of n_timepoints durations of Gaussian random noise
    that is then bandpass filtered with the passband
    """

    sos = signal.butter(10, passband, "bandpass", fs=fs, output="sos")
    return signal.sosfilt(sos, gen_randn_timeseries(n_reps, n_timepoints), axis=1)
