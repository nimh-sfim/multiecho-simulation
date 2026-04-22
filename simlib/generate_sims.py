"""Tools to generate simulated data based on desired inputs."""

from turtle import pd
from typing import Dict, List, Tuple, Union

import numpy as np
import numpy.typing as npt
import nibabel as nb
import pandas as pd


def tile_variables(
    vars_to_tile: List[Union[float, int, npt.NDArray]],
    target_shape: Union[Tuple[int], None] = None,
) -> Tuple[npt.NDArray]:
    """Tile variables to match a target shape.

    Each variable must be a single value, have matching initial dimensions,
    or dimensions of size 1.
    If a variable is (1 x N) then it is squeezed to be size N

    Parameters
    ----------
    vars_to_tile : List[float, int, or npt.NDArray]
        A list of variables to tile.
        Each variable can be a single value or an array.
    target_shape : Tuple[int] None
        The shape to tile the variables to.
        If None, then the target shape will be the variable
        with the maximum number of non-single dimensions.

    Returns
    -------
    vars_to_tile : List[npt.NDArray]
        A list of tiled variables.
        Each variable will have the same shape as target_shape
        with repeated values in the dimensions to match target_shape
    """
    var_shapes = []
    if target_shape is None:
        tmp_target_shape = np.array([])
    for varidx, var in enumerate(vars_to_tile):
        if isinstance(var, (int, float, List)):
            # Make each variable an np.array
            vars_to_tile[varidx] = np.array([var])
        var_shapes.append(vars_to_tile[varidx].shape)
        if len(var_shapes[varidx]) == 2 and var_shapes[varidx][0] == 1:
            # squeeze if var is a 2D array with the first dimension of size 1
            vars_to_tile[varidx] = np.squeeze(vars_to_tile[varidx])
            var_shapes[varidx] = vars_to_tile[varidx].shape
        if target_shape is None:
            # If no predefined target shape, track the largest size input var
            non_single_dims = [dim for dim in var_shapes[varidx] if dim > 1]
            if len(non_single_dims) >= len(tmp_target_shape):
                tmp_target_shape = tuple(non_single_dims)
    if target_shape is None:
        target_shape = tmp_target_shape

    # print(f"Before {var_shapes=}")

    dims_good = True
    for varidx, var in enumerate(vars_to_tile):
        for dim_idx in range(len(target_shape)):
            if (
                len(var_shapes[varidx]) > dim_idx
                and var_shapes[varidx][dim_idx] != 1
                and var_shapes[varidx][dim_idx] != target_shape[dim_idx]
            ):
                # Raise an error if dimensions>1 aren't the same
                dims_good = False
            elif len(var_shapes[varidx]) <= dim_idx or var_shapes[varidx][0] == 1:
                vars_to_tile[varidx] = np.tile(
                    vars_to_tile[varidx],
                    ((target_shape[dim_idx],) + tuple(np.ones(dim_idx, dtype=int))),
                ).transpose(tuple(range(1, dim_idx + 1)) + (0,))
                var_shapes[varidx] = vars_to_tile[varidx].shape

    if not dims_good:
        raise ValueError(f"Dimensions not compatiable. {target_shape=}, {var_shapes=}")
    # print(f"After {var_shapes=}")

    return vars_to_tile


def monoexponential(
    tes: Union[int, float, List[int], List[float], npt.NDArray],
    s0: Union[float, List[float], npt.NDArray],
    r2star: Union[float, List[float], npt.NDArray],
    t2star_input: bool = False,
) -> npt.NDArray:
    """Specify a monoexponential model for MRI signal calculations.

    (Based on from https://github.com/ME-ICA/tedana/blob/main/tedana/decay.py,
    but flipped to inputting R2* instead of T2*)

    Parameters
    ----------
    tes : (E,) :obj:`list` :obj:`npt.NDArray`
        Echo times
    s0 : :obj:`list` :obj:`npt.NDArray`
        S0 parameter
    r2star : :obj:`list` :obj:`npt.NDArray`
        R2* parameter
    t2star_input : bool
        If True, the r2star parameter is actually T2* and will be inverted to R2* before calculating the signal
        Default is False

    Returns
    -------
    S : :obj:`float` :obj:`npt.NDArray`
        Calculate signal based on S = s0 * np.exp(-tes * r2star)
        If tes is a list or array of echo times and s0 and r2star are not single values,
        then the output will be an array with dimensions of s0 and r2star with tes as the final dimension.
        If tes is 2D or higher, then s0 and r2star should be single values that will be applied to all tes values.

    Notes
    -----
    Units for tes, s0, and r2star must be the same (i.e. ms OR sec).
    """
    if t2star_input:
        r2star = 1 / r2star

    tes, s0, r2star = tile_variables([tes, s0, r2star])

    return s0 * np.exp(-tes * r2star)


def first_order_taylor_approx(
    tes: Union[int, float, List[int], List[float], npt.NDArray],
    delta_s0: Union[float, List[float], npt.NDArray],
    mean_s0: Union[float, List[float], npt.NDArray],
    delta_r2star: Union[float, List[float], npt.NDArray],
    mean_r2star: Union[float, List[float], npt.NDArray],
    t2star_input: bool = False,
) -> npt.NDArray:
    if t2star_input:
        delta_r2star = 1 / delta_r2star
        mean_r2star = 1 / mean_r2star

    tes, delta_s0, mean_s0, delta_r2star, mean_r2star = tile_variables(
        [tes, delta_s0, mean_s0, delta_r2star, mean_r2star]
    )

    s_percent_change = (
        delta_s0 / mean_s0 - tes * delta_r2star - tes * delta_r2star * delta_s0 / mean_s0
    )

    s_mean = monoexponential(tes=tes, s0=mean_s0, r2star=mean_r2star)

    return s_mean * (1 + s_percent_change)


def linear_approx(
    tes: Union[int, float, List[int], List[float], npt.NDArray],
    delta_s0: Union[float, List[float], npt.NDArray],
    mean_s0: Union[float, List[float], npt.NDArray],
    delta_r2star: Union[float, List[float], npt.NDArray],
    mean_r2star: Union[float, List[float], npt.NDArray],
    t2star_input: bool = False,
) -> npt.NDArray:
    if t2star_input:
        delta_r2star = 1 / delta_r2star
        mean_r2star = 1 / mean_r2star

    tes, delta_s0, mean_s0, delta_r2star, mean_r2star = tile_variables(
        [tes, delta_s0, mean_s0, delta_r2star, mean_r2star]
    )

    s_percent_change = delta_s0 / mean_s0 - tes * delta_r2star

    s_mean = monoexponential(tes=tes, s0=mean_s0, r2star=mean_r2star)

    return s_mean * (1 + s_percent_change)


def all_decay_models(
    tes: Union[int, float, List[int], List[float], npt.NDArray],
    delta_s0: Union[float, List[float], npt.NDArray],
    mean_s0: Union[float, List[float], npt.NDArray],
    delta_r2star: Union[float, List[float], npt.NDArray],
    mean_r2star: Union[float, List[float], npt.NDArray],
    t2star_input: bool = False,
) -> Dict[str, npt.NDArray]:
    return {
        "monoexponential": monoexponential(
            tes=tes,
            s0=mean_s0 + delta_s0,
            r2star=mean_r2star + delta_r2star,
            t2star_input=t2star_input,
        ),
        "taylor1_approx": first_order_taylor_approx(
            tes=tes,
            delta_s0=delta_s0,
            mean_s0=mean_s0,
            delta_r2star=delta_r2star,
            mean_r2star=mean_r2star,
            t2star_input=t2star_input,
        ),
        "linear_approx": linear_approx(
            tes=tes,
            delta_s0=delta_s0,
            mean_s0=mean_s0,
            delta_r2star=delta_r2star,
            mean_r2star=mean_r2star,
            t2star_input=t2star_input,
        ),
    }


def calc_delta_r2s_s0_given_s_pchange_proportion(
    inputted_s: Union[float, npt.NDArray],
    s0_mean: Union[float, npt.NDArray],
    r2s_mean: Union[float, npt.NDArray],
    te_baseline: float,
    proportion_s0_r2s: float,
    prop_to_scale: str = "signal",
) -> Union[Tuple[float, float], Tuple[npt.NDArray, npt.NDArray]]:
    """Finds relationship between delta R2* and S0 for a given signal proportion of the two.

    Parameters
    ----------
    inputted_s : :obj:`float` :obj:`npt.NDArray`
        A single value, a 1D array, or a 2D array of values
        The percent change of the overall signal above the baseline value
        S_spc = (S-mean(S))/mean(S)
        (1==100% change)
    s0_mean : :obj:`float` :obj:`npt.NDArray`
        mean s0 signal parameter across time
        If this is a single value, it will be applied to all values in inputted_s.
        If this is 1D array, the length should match the 0th dimension of inputted_s.
        If this is a 2D array, the shape should match the shape of inputted_s.
    r2s_mean : :obj:`float` :obj:`npt.NDArray`
        mean R2* parameter across time
        If this is a single value, it will be applied to all values in inputted_s.
        If this is 1D array, the length should match the 0th dimension of inputted_s.
        If this is a 2D array, the shape should match the shape of inputted_s.
    te_baseline : :obj:`float`
        Echo time
        The delta S0 and R2* values with the prespecified proportionwill be calculated for this echo time (TE).
        At other echo times, the proportional relationship between delta S0 and R2* may be different.
    proportion_s0_r2s : :obj:`float` :obj:`npt.NDArray`
        The proportion of delta S0 / delta R2* used to generate inputted_s.
        1 is  pure delta S0 changes. 0 is pure delta R2* changes.
        If this is a single value, it will be applied to all values in inputted_s.
        If this is 1D array, the length should match the 0th dimension of inputted_s.
        If this is a 2D array, the shape should match the shape of inputted_s.
    prop_to_scale : str
        If "signal" then delta_s0 and delta_r2s are scaled so proportion_s0_r2s represents
        the proportion of delta s0 and delta r2s contribution to the overall signal change (inputted_s).
        If "variance" then delta_s0 and delta_r2s are scaled so that the variance of the signal change from delta s0 and delta r2s matches the proportion_s0_r2s.
        If variance is kept constant, the signal will change.
        Default is "signal" but that might change.
    Returns
    -------
    delta_s0: obj:`float` :obj:`npt.NDArray`
        The value(s) for delta S0 for all values in inputted_s given baseline parameters
        Will be the same size as inputted_s
    delta_r2s: obj:`float` :obj:`npt.NDArray`
        The value(s) for delta R2* for all values in inputted_s given baseline parameters
        Will be the same size as inputted_s

    Math
    ----
    Full math explanation is currently in the README.md for this repository.
    The math is based on a first order Taylor expansion of the monoexponential function.
    There is an additional approximation in that an interaction term between delta s0 and delta r2s is removed.
    The removal makes the relationship linear.
    The assumption is that this interaction term is small
    delta_s0 is proportion_s0_r2s * inputted_s * s0_mean
    delta_r2s is -(1-proportion_s0_r2s) * inputted_s / te_baseline
    """
    if isinstance(inputted_s, (int, float)):
        # if inputted_s is a single value, convert to array for consistent processing
        inputted_s = np.array([inputted_s])

    s_dims = inputted_s.shape
    s0_mean, r2s_mean, proportion_s0_r2s = tile_variables(
        [s0_mean, r2s_mean, proportion_s0_r2s], target_shape=s_dims
    )

    if prop_to_scale.lower() == "signal":
        delta_s0 = proportion_s0_r2s * inputted_s * s0_mean
        delta_r2s = -(1 - proportion_s0_r2s) * inputted_s / te_baseline
    elif prop_to_scale.lower() == "variance":
        delta_s0 = np.sqrt(proportion_s0_r2s) * inputted_s * s0_mean
        delta_r2s = -np.sqrt(1 - proportion_s0_r2s) * inputted_s / te_baseline
    return delta_s0, delta_r2s


def add_noise_to_signal(
    signal: npt.NDArray, noise_var_ratio: Union[float, npt.NDArray], noise_seed: float = None
) -> npt.NDArray:
    """Add Gaussian noise to a signal with a specified variance ratio.

    Parameters
    ----------
    signal : :obj:`npt.NDArray`
        The input signal to which noise will be added.
        Can be an ND array and the variance will be calculated over the last axis (i.e. time).
    noise_var_ratio : :obj:`float` :obj:`npt.NDArray`
        The ratio of variance explained by the noise vs signal.
        0 means no noise and 0.5 means 50% noise.
            If this is a single value, it will be applied to all samples in the signal.
            If this is a ND array, the length should match the the first ND dimensions of the signal
            and it will be tiled to match all dimensions of the signal except the last dimension.
    noise_seed : :obj:`int`
        Gaussian noise will be randomly generated each time the function is run.
        If this is set to a number, that will be used as a seed for the noise generator
        so that consistent results are generated across repeated runs.
        Default=None

    Returns
    -------
    signal_with_noise : :obj:`npt.NDArray`
        An array the same size as the input signal with Gaussian noise added according to the specified variance ratio.
        The variance of each time series will be scaled to make the input, but will now be a combination of the inputted
        signal and the generated noise according to the specified noise variance ratio.
    """

    if isinstance(noise_var_ratio, (int, float)):
        noise_var_ratio = np.array([noise_var_ratio])
    if (
        noise_var_ratio.ndim >= 1
        and noise_var_ratio.shape[0] > 1
        and noise_var_ratio.shape != signal.shape
    ):
        raise ValueError(
            "If noise_var_ratio is an array, the dimensions should match signal's dimensions. "
            f"{noise_var_ratio.shape=}, {signal.shape=}"
        )

    if np.sum(noise_var_ratio < 0) or np.sum(noise_var_ratio >= 1):
        raise ValueError("noise_var_ratio values must be >= 0 and < 1.")

    if np.sum(noise_var_ratio == 0):
        # If noise_var_ratio is 0, return the original signal without adding noise
        return signal

    noise_var_ratio = tile_variables([noise_var_ratio], target_shape=signal.shape[:-1])[0]
    if noise_seed:
        rng = np.random.default_rng(seed=noise_seed)
    else:
        rng = np.random.default_rng()

    noise_vals = rng.standard_normal(signal.shape)

    # variance is calculated over the last axis (i.e.time)
    # Total variance should match the variance of the initial signal,
    # with both signal and noise scaled to match the specified noise variance ratio
    sig_var_scaling = np.sqrt(1 - noise_var_ratio)
    noise_var_scaling = (
        np.sqrt(noise_var_ratio) * np.std(signal, axis=-1) / np.std(noise_vals, axis=-1)
    )

    # Tile to include time points as the last dimension
    sig_var_scaling = tile_variables([sig_var_scaling], target_shape=signal.shape)[0]
    noise_var_scaling = tile_variables([noise_var_scaling], target_shape=signal.shape)[0]

    signal_with_noise = (sig_var_scaling * signal) + (noise_var_scaling * noise_vals)

    return signal_with_noise


def generate_echoes_from_ica(
    spatial_comp_file: str,
    mixing_matrix: Union[str, npt.NDArray],
    proportion_s0_r2s: Union[str, npt.NDArray, List[float], List[int]],
    s0_mean: Union[str, float, npt.NDArray],
    t2s_mean: Union[str, float, npt.NDArray],
    te_baseline: Union[float, int],
    tes: Union[int, float, List[int], List[float], npt.NDArray],
    output_prefix: str,
    comp_scaling: float = 0.1,
    prop_to_scale: str = "signal",
    noise_scaling: float = 0.0,
    noise_seed: Union[float, int] = 42,
    verbose: bool = False,
) -> npt.NDArray:
    """Generate echo signals based on ICA components.

    Parameters
    ----------
    spatial_comp_file : :obj:`str`
        A nifti file containing the spatial ICA components.
        Like the desc-ICA_components.nii.gz file output from tedana.
    mixing_matrix : :obj:`str` :obj:`npt.NDArray`
        Mixing matrix for ICA components.
        Can be a file like the desc-ICA_mixing.tsv file output from tedana or
        a (time x components) np.array.
    proportion_s0_r2s : :obj:`str` :obj:`npt.NDArray` :obj:`List[float]` :obj:`List[int]`
        Proportion of S0 and R2* signals.
        1 is pure S0, 0 is pure R2*.
        Can be a vector with a value for each component,
        or a file like the desc-tedana_metrics.tsv file output from tedana
        where the values are read from the "classification" column
        (1 for rejected, 0 for accepted).
    s0_mean : :obj:`str` :obj:`float` :obj:`int` :obj:`npt.NDArray`
        Mean S0 parameter
        Can be a Nifti file like the S0map.nii output from tedana,
        a single value, or an array.
        The file or array should have the same dimensions as the spatial_comp_file.
        A single value will be applied to all components and voxels.
    t2s_mean : :obj:`str` :obj:`float` :obj:`int` :obj:`npt.NDArray`
        Mean T2* parameter (NOT R2*, T2* = 1/R2*) (seconds, not milliseconds)
        Can be a Nifti file like the S0map.nii output from tedana,
        a single value, or an array.
        The file or array should have the same dimensions as the spatial_comp_file.
        A single value will be applied to all components and voxels.
    te_baseline : :obj:`float` :obj:`int`
        Baseline echo time (ms)
        To proportion of S0 and R2* changes for a pre-specfied signal change is
        calculated for this echo time.
        The entire point of this simulation to model how the signal changes at other echo times
        given a certain proportion of S0 and R2* changes at this baseline echo time.
    tes : (E,) :obj:`list` :obj:`npt.NDArray`
        Echo times (ms) to write out simulated data for
    output_prefix : :obj:`str`
        Prefix for output files.
        Can be a path + prefix or just a prefix to be written to the current directory.
    comp_scaling : :obj:`float`
        Scaling factor for components and the mixing matrix.
        This controls the overall size of the signal changes in the simulated data.
        The spatial components and mixing matrix are scaled relative to their global max value,
        but, given the ways they are combined,
        this doesn't control the precise max signal change in the output data.
        Default is 0.1 to ensure that the simulated signal changes have a rough maximum of 10%.
        If the source data has a larger linear drift, a larger scaling factor may be appropriate.
    prop_to_scale : :obj:`str`
        Property to scale
        If "signal" then the delta_s0 and delta_r2s are scaled so that the proportion_s0_r2s
        represents the proportion of S0 and R2* contribution to the overall signal change.
        If "variance" then the delta_s0 and delta_r2s are scaled so that the variance of the
        signal change from delta s0 and delta r2s matches the proportion_s0_r2s.
        The equations are the same expect s0 and r2s are multiplied by the proportion_s0_r2s
        for signal and the square root for variance.
        If "variance" the signal will change based on proporation even at TE=te_baseline.
        Default is "signal" but that might change.
    noise_scaling : :obj:`float` :obj:`int`
        Scaling factor for added Gaussian noise.
        0 means no noise and 0.5 means 50% noise.
        The time series variance is kept constant so that the signal is scaled down
        to keep the overall variance the same when noise is added.
        Note: Each voxel's time series is scaled separately, so the SNR is similar across voxels.
            In real fMRI data, the baseline thermal/Gaussian noise may be constant across voxels
            so the SNR is higher in voxels with higher signal.
            May tweak this in the future to calculate a globally constant noise floor
            for a pre-specific total SNR across voxels.
        Default is 0.
    noise_seed : :obj:`float` :obj:`int`
        Seed for the random number generator used to generate Gaussian noise.
        Default is 42.
    verbose : :obj:`bool`
        If True, return additional variables from the function call
        Default is False

    Returns
    -------
    sim_echo_data : :obj:`npt.NDArray`
        A (voxels, timepoints, echo_times) array with the
        simulated echo data based on the inputted parameters and the monoexponential decay model.

    If verbose is True, also returns a dictionary with the following variables:
    delta_s0 & delta_r2s:
        (voxels, timepoints, components) npt.NDArray with the delta S0 and delta R2* values
        for each component and time point in each voxel.
    voxelwise_delta_s0 & voxelwise_delta_r2s,
        (voxels, timepoints) npt.NDArray with the overall delta S0 and delta R2* values for each voxel and time point
        based on the combination of the component values and the mixing matrix values for each voxel and time point
    dimensions:
        A tuple with the x, y, and z dimension for the full volume
            non_zero_vox,
    non_zero_vox:
        A boolean array of length x*y*z with True for voxels
        with non-zero values in the spatial components.
        Can be used with dimensions to reshape the output data back into a 3D or 4D volume.

    Notes
    -----
    Currently assumes tes, and te_baseline are in milliseconds and t2s_mean is in seconds.
    TODO: Add checks and conversions to ensure that the units are consistent across parameters.
    """
    # comp_scaling = 0.1
    # mixing_matrix = "desc-ICA_mixing.tsv"
    # spatial_comp_file = "tedana_out_for_sim/desc-ICA_components.nii.gz"
    # proportion_s0_r2s = "tedana_out_for_sim/desc-tedana_metrics.tsv"
    # s0_mean = "/S0map.nii"
    # t2s_mean = "T2starmap.nii"
    # te_baseline = 28
    # prop_to_scale = "signal"
    # tes = [10, 20, 30, 40, 50]
    # output_file = "test"

    # Load ICA mixing matrix from a file (tedana-style output)
    # Will also work with a directly inputted mixing matrix as an np.array
    # but validity checks not currently implemented
    if isinstance(mixing_matrix, str):
        mixing_matrix = np.loadtxt(mixing_matrix, delimiter="\t", skiprows=1)

    n_timepoints, n_comps = mixing_matrix.shape

    # Scale the mixing matrix so the max val across all components and time points is comp_scaling
    # (e.g. 0.1) to ensure that the simulated signal changes are on the order of 10% or less,
    # which is more realistic for fMRI data.
    mixing_matrix = comp_scaling * mixing_matrix / np.max(mixing_matrix)

    # Can read the classification row from a tedana metrics file and
    # assign 1 to rejected (pure S0) and 0 (pure R2*) to accepted.
    if isinstance(proportion_s0_r2s, str):
        print(f"Loading metrics from {proportion_s0_r2s}")
        metrics_df = pd.read_csv(proportion_s0_r2s, sep="\t")
        proportion_s0_r2s = np.zeros(len(metrics_df))
        proportion_s0_r2s[metrics_df["classification"].to_numpy() == "rejected"] = 1
    if isinstance(proportion_s0_r2s, list):
        proportion_s0_r2s = np.array(proportion_s0_r2s)

    if len(metrics_df) != n_comps:
        raise ValueError(
            "# of components in metrics file does not match # of components in mixing matrix. "
            f"{len(metrics_df)=}, {n_comps=}"
        )

    # spatial components must be a nii file to get the header info to write out new simulated data
    comps_img = nb.load(spatial_comp_file)
    nii_header = comps_img.header
    comps_data = comps_img.get_fdata()
    x, y, z, c = comps_data.shape
    spatial_comps = comps_data.reshape(x * y * z, c)

    if c != n_comps:
        raise ValueError(
            "Number of components in spatial components does not match number of "
            f"components in mixing matrix. {c=}, {n_comps=}"
        )

    # reduce comps to only those with non-zero values
    non_zero_vox = np.sum(spatial_comps != 0, axis=1) > 0
    n_vox = np.sum(non_zero_vox)
    spatial_comps = spatial_comps[non_zero_vox, :]

    # Scale comps so that the maximum absolute value of the sum across components is 0.1
    # This is to ensure that the simulated signal changes roughtly are on the order of 10% or less,
    # which is more realistic for fMRI data.
    spatial_comps = spatial_comps - np.mean(spatial_comps)
    spatial_comps = spatial_comps / (comp_scaling * np.max(np.abs(np.sum(spatial_comps, axis=1))))

    # Load S0 and T2* maps to get mean values for each parameter across the brain
    # No checks for accuracy yet, but can also give an (x, y, z) array
    if isinstance(s0_mean, str):
        s0_input = nb.load(s0_mean)
        s0_mean = s0_input.get_fdata()
    # if it's a single value, tile to all voxels
    if isinstance(s0_mean, (float, int)):
        s0_mean = np.tile(s0_mean, (x, y, z))
    s0_mean = s0_mean.reshape(x * y * z)
    s0_mean = s0_mean[non_zero_vox]

    if isinstance(t2s_mean, str):
        t2s_input = nb.load(t2s_mean)
        t2s_mean = t2s_input.get_fdata()
    # if it's a single value, tile to all voxels
    if isinstance(t2s_mean, (float, int)):
        t2s_mean = np.tile(t2s_mean, (x, y, z))
    t2s_mean = t2s_mean.reshape(x * y * z)
    r2s_mean = np.zeros_like(t2s_mean)
    # Avoid dividing by zero for voxels outside of the T2* map mask.
    r2s_mean[t2s_mean != 0] = 0.001 / t2s_mean[t2s_mean != 0]
    r2s_mean = r2s_mean[non_zero_vox]

    # tile everything to (voxels, time points, components)
    s0_mean_tiled = np.tile(s0_mean, [n_timepoints, n_comps, 1]).transpose([2, 0, 1])
    r2s_mean_tiled = np.tile(r2s_mean, [n_timepoints, n_comps, 1]).transpose([2, 0, 1])
    mixing_matrix_tiled = np.tile(mixing_matrix, [n_vox, 1, 1])
    proportion_s0_r2s_tiled = np.tile(proportion_s0_r2s, [n_vox, 1, 1])
    spatial_comps_tiled = np.tile(spatial_comps, [n_timepoints, 1, 1]).transpose([1, 0, 2])
    print(
        f"{s0_mean_tiled.shape=}, {r2s_mean_tiled.shape=}, "
        f"{proportion_s0_r2s_tiled.shape=}, {spatial_comps_tiled.shape=}"
    )

    # Calculate the proportional delta s0 and delta r2s values
    delta_s0, delta_r2s = calc_delta_r2s_s0_given_s_pchange_proportion(
        inputted_s=mixing_matrix_tiled,
        s0_mean=s0_mean_tiled,
        r2s_mean=r2s_mean_tiled,
        proportion_s0_r2s=proportion_s0_r2s_tiled,
        te_baseline=te_baseline,
        prop_to_scale=prop_to_scale,
    )

    voxelwise_delta_s0 = np.mean(spatial_comps_tiled * delta_s0, axis=2)
    voxelwise_delta_r2s = np.mean(spatial_comps_tiled * delta_r2s, axis=2)
    print(f"{voxelwise_delta_s0.shape=} {voxelwise_delta_r2s.shape=}")
    print(
        f"{np.mean(voxelwise_delta_s0)=} {np.std(voxelwise_delta_s0)=}, "
        f"{np.max(np.abs(voxelwise_delta_s0 / s0_mean_tiled[:, :, 0]))=}"
    )
    print(
        f"{np.mean(voxelwise_delta_r2s)=} {np.std(voxelwise_delta_r2s)=}, "
        f"{np.max(np.abs(voxelwise_delta_r2s / r2s_mean_tiled[:, :, 0]))=}"
    )

    sim_echo_data = np.full((n_vox, n_timepoints, len(tes)), np.nan)
    for te_idx, te in enumerate(tes):
        print(f"Calculating {te=}ms echo data")
        sim_echo_data[:, :, te_idx] = monoexponential(
            tes=te,
            s0=s0_mean_tiled[:, :, 0] + voxelwise_delta_s0,
            r2star=r2s_mean_tiled[:, :, 0] + voxelwise_delta_r2s,
        )

        sim_echo_data[:, :, te_idx] = add_noise_to_signal(
            sim_echo_data[:, :, te_idx], noise_var_ratio=noise_scaling, noise_seed=noise_seed
        )

        sim_echo_data_full_shape = np.zeros((x * y * z, n_timepoints))
        sim_echo_data_full_shape[non_zero_vox, :] = sim_echo_data[:, :, te_idx]
        sim_echo_data_full_shape = np.reshape(sim_echo_data_full_shape, (x, y, z, n_timepoints))

        sim_echo_data_img = nb.Nifti1Image(
            sim_echo_data_full_shape, affine=nii_header.get_best_affine(), header=nii_header
        )

        nb.save(sim_echo_data_img, f"{output_prefix}_TE{te}")

    if verbose:
        return (
            sim_echo_data,
            {
                "voxelwise_delta_s0": voxelwise_delta_s0,
                "voxelwise_delta_r2s": voxelwise_delta_r2s,
                "delta_s0": delta_s0,
                "delta_r2s": delta_r2s,
                "dimensions": (x, y, z),
                "non_zero_vox": non_zero_vox,
            },
        )
    else:
        return sim_echo_data
