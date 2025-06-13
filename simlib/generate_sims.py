"""Tools to generate simulated data based on desired inputs."""

from typing import Dict, List, Tuple, Union

import numpy as np
import numpy.typing as npt


def monoexponential(
    tes: Union[int, float, List[int], List[float], npt.NDArray],
    s0: Union[float, List[float], npt.NDArray],
    r2star: Union[float, List[float], npt.NDArray],
) -> npt.NDArray:
    """Specify a monoexponential model for MRI signal calculations.

    (Based on from https://github.com/ME-ICA/tedana/blob/main/tedana/decay.py, but flipped for R2*)

    Parameters
    ----------
    tes : (E,) :obj:`list` :obj:`npt.NDArray`
        Echo times
    s0 : :obj:`list` :obj:`npt.NDArray`
        S0 parameter
    r2star ::obj:`list` :obj:`npt.NDArray`
        R2* parameter

    Returns
    -------
    S : obj:`float`
        Calculated signal based on S = s0 * np.exp(-tes * r2star)

    Notes
    -----
    Units for tes, s0, and r2star must be the same (i.e. ms OR sec).

    This function is a simple multiplication.
    For the 3 parameters:
    2-3 need to be a single value
    s0 and r2star must be the same size and tes is a single value
    tes is multiple values and s0 and r2star are a single value
    TODO Add tests to give descriptive errors instead of just having this function crash
    """
    return s0 * np.exp(-tes * r2star)


def calc_delta_r2s_s0_given_s_pchange_proportion(
    te: float,
    s0_baseline: float,
    r2s_baseline: float,
    inputted_s: Union[float, npt.NDArray],
    proportion_s0_r2s: float,
) -> Union[Tuple[float, float], Tuple[npt.NDArray, npt.NDArray]]:
    """Finds relationship between delta R2* and S0 for a given signal proportion of the two.

    Parameters
    ----------
    te : :obj:`float`
        Echo time
        This function can be paired with monoexponential.
        The ratio of delta S0 and R2* can be calculated with the te given here.
        Those same S0 and R2* values can generate time series at any TE using monoexponential.
    s0 : :obj:`float`
        Initial signal parameter
    r2star : :obj:`float`
        R2* parameter
    inputted_s : :obj:`float` :obj:`npt.NDArray`
        The percent change of the overall signal above the baseline value (100==100% change)
    proportion_s0_r2s : :obj:`float`
        The proportion of delta S0 / delta R2* used to generate inputted_s.
        1 is  pure delta S0 changes. 0 is pure delta R2* changes.

    Returns
    -------
    delta_s0: obj:`float`
        The value(s) for delta S0 for all values in inputted_s given baseline parameters
        Will be the same size as inputted_s
    delta_r2s: obj:`float`
        The value(s) for delta R2* for all values in inputted_s given baseline parameters
        Will be the same size as inputted_s

    Math
    ----
    This function is based on:
    (s0_baseline + proportion_s0_r2s*delta_s0_scale) * np.exp(-te(r2s_baseline + (1-proportion_s0_r2s)*delta_r2s_scale)) = S
    That is, for a constant S, one can calculate delta_s0 and delta_r2s for any ratio of delta s0 and r2s signal.
    By using this formulation, one can take a simulated fMRI time series, S, and calculate the delta s0 and r2s values
    if the fluctuations are pure S0, pure R2s or anything in between.
    Delta s0 and r2s can be calculated at one echo time and then applied to other echoes to show how the signal evolves across the echoes.

    In this formulation S is inputted_s different from the value when delta s0 and r2s are both 0 (i.e. the baseline).
    That is, if the baseline is 200 and one is calculating for a 10% signal change, S=220
    -te*r2s_baseline -te*(1-proportion_s0_r2s)*delta_r2s_scale = np.log(S) - np.log(s0_baseline + proportion_s0_r2s*delta_s0_scale)
    (1-proportion_s0_r2s)*delta_r2s_scale = (np.log(s0_baseline + proportion_s0_r2s*delta_s0_scale) - te*r2s_baseline- np.log(S)) / TE
    For proportion_s0_r2s==1:
        0 = (np.log(s0_baseline + delta_s0_scale) - te*r2s_baseline - np.log(S)) / TE
        np.log(s0_baseline + delta_s0_scale) = te*r2s_baseline + np.log(S)
        delta_s0_scale = np.exp(te*r2s_baseline + np.log(S)) - s0_baseline
        The above eq is used to calculate delta_s0_scale in the function
    With a known delta_s0_scale, the above eq with proportion_s0_r2s can be used to solve for delta_r2s_scale
        (1-proportion_s0_r2s)*delta_r2s_scale = (np.log(s0_baseline + proportion_s0_r2s*delta_s0_scale) - te*r2s_baseline- np.log(S)) / TE
        delta_r2s_scale = (np.log(s0_baseline+proportion_s0_r2s*delta_s0_scale) - te*r2s_baseline - np.log(S))/((1-proportion_s0_r2s)*te)
        The above eq is used to calculate delta_r2s_scale in the function, but it fails for 1-proportion_s0_r2s==0, but,
        in that case (1-proportion_s0_r2s)*delta_r2s_scale should be 0 so it's just set to 0
    The returned delta_s0 and delta_r2s values are the scaled values multiplied by their proportions
    """
    # calculate expected S (raw signal value)
    # baseline_value is the value for the decay curve with baseline s0 and r2s values
    baseline_value = monoexponential(te, s0_baseline, r2s_baseline)
    # S is the signal that's percent_change larger or smaller than baseline_value
    s = (1 + (inputted_s / 100)) * baseline_value

    # Calculate delta_s0 so that the decay curve results in S if changing signal is purely S0 (i.e. delta_r2s==0)
    delta_s0_scale = np.exp((te * r2s_baseline) + np.log(s)) - s0_baseline

    if proportion_s0_r2s < 1:
        delta_r2s_scale = (
            np.log(s0_baseline + proportion_s0_r2s * delta_s0_scale)
            - te * r2s_baseline
            - np.log(s)
        ) / ((1 - proportion_s0_r2s) * te)
    else:
        delta_r2s_scale = np.zeros((inputted_s.shape))

    # Returning the delta values scaled by proportions
    return proportion_s0_r2s * delta_s0_scale, (1 - proportion_s0_r2s) * delta_r2s_scale


def generate_s_with_noise(
    te_baseline: float,
    s0_baseline: float,
    r2s_baseline: float,
    inputted_s: npt.NDArray,
    proportion_s0_r2s: float,
    te_output: float,
    noise_var_ratio: float,
    noise_seed: float = None,
    pchange: bool = True,
) -> npt.NDArray:
    """Generate time series for a given echo time.

    Parameters
    ----------
    te_baseline : :obj:`float`
        Echo time for which the delta S0 and R2* values calcualted to match proportion_s0_r2s
    s0_baseline : :obj:`float`
        S0 parameter
    r2star_baseline : :obj:`float`
        R2* parameter
    inputted_s (samples x time): :obj:`npt.NDArray`
        The percent change of the overall signal above the baseline value (100==100% change)
    proportion_s0_r2s : :obj:`float`
        The proportion of delta S0 / delta R2* used to generate inputted_s.
        1 is  pure delta S0 changes. 0 is pure delta R2* changes.
    te_output : :obj:`float`
        Echo time to use for generating the outputted time series
    noise_var_ratio : :obj:`float`
        The ratio variance explained by the noise vs signal
        0 is no noise and 0.5 is 50% noise
    noise_seed: :obj:`int`
        Gaussian noise will be randomly generated each time the function is run.
        If this is set to a number, that will be used as a seed for the noise generator
        so that consistent results are generated across repeated runs
        Default=None
    pchange: :obj:`bool`
        If True, convert time series to percent change from mean (1=100%)
        Default = True

    Returns
    -------
    sig_and_noise: :obj:`npt.NDArray`
        An array the same size as inputted_s with the time series and the specified echo
        with the specified delta R2* and S0 proportions
    """

    # number of samples (voxels) in the inputted data
    n_samples, n_timepoints = inputted_s.shape

    delta_s0, delta_r2s = calc_delta_r2s_s0_given_s_pchange_proportion(
        te_baseline, s0_baseline, r2s_baseline, inputted_s, proportion_s0_r2s
    )
    s_timeseries = monoexponential(te_output, s0_baseline + delta_s0, delta_r2s)

    if noise_seed:
        rng = np.random.default_rng(seed=noise_seed)
    else:
        rng = np.random.default_rng()

    noise_vals = rng.standard_normal(s_timeseries.shape)

    # axis 0 is samples or voxels and variance is calculated over axis 1 (time)
    noise_var = np.var(noise_vals, axis=1)
    sig_var = np.var(s_timeseries, axis=1)

    noise_scale = np.sqrt(noise_var_ratio * sig_var / noise_var)

    # To match the desired signal vs noise ratio, this should be the scaling factor for the signal
    sig_var_scaling = np.sqrt(1 - noise_var_ratio)

    sig_and_noise = (sig_var_scaling * s_timeseries) + (
        np.matlib.repmat(noise_scale, n_timepoints, 1).T * noise_vals
    )

    if pchange:
        total_variance = np.var(sig_and_noise, axis=1)
        total_mean = np.mean(sig_and_noise, axis=1)
        return (sig_and_noise - total_mean[:, None]) / total_variance[:, None]

    else:
        return sig_and_noise
