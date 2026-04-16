import os
from glob import glob
import os.path as op

import warnings
# from IPython.display import HTML

# import math
# from math import ceil
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
# import numpy.typing as npt
import numpy as np
import nibabel as nb
from numpy import polyfit

from nilearn.image import check_niimg, new_img_like
from scipy.stats import zscore, linregress
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA, FastICA

from tedana import io, stats, utils
import tensorly as tl
from tensorly.decomposition import parafac
from tensorly.cp_tensor import cp_normalize
import tensorly.contrib.sparse as stl
from tensorly.random import random_cp

# ============================================================
# Functions for tensor tools, etc
# file_to_2D, mask_data, plot_time_and_echoes_brain, plot_echo_profiles,
# compare_rsqared, pca_dim_reduced, get_new_slopes, initialize_cp, parafac_iter, ica_to_parafac_iter 
# ============================================================


def file_to_2d(drive_loc, file_name):
    """
    Reshape input data.

    Parameters
    ----------
    drive_loc : string
    file_name : string

    Returns
    -------
    data_img_2d : 2D data
    shape : shape of data_img_2d
    """
    file_path = op.join(drive_loc, file_name)
    data_img = nb.load(file_path).get_fdata()
    shape = data_img.shape
    data_img_2d = data_img.reshape(shape[0]*shape[1]*shape[2], shape[3])
    return data_img_2d, shape


def mask_data(mask, data):
    """
    Mask data.

    Parameters
    ----------
    mask : ndarray
    data : ndarray

    Returns
    -------
    masked_data : masked data
    """
    masked_data = data[mask, :]
    return masked_data


def plot_time_and_echoes_brain(
        factors, weights, rank, mask, vol_dims, n_slice=20, slice_rot=1):
    """
    Plot spatial maps, time series, and echo profiles for all components.

    Parameters
    ----------
    factors : list
    weights: list
    rank : int
    mask : ndarray
    vol_dims : ndarray
    n_slice : int
    slice_rot : int
    """
    spatial_maps = factors[0]
    echo_profiles = factors[2]
    time_series = factors[1]

    total_weights = np.sum(weights)

    fig = plt.figure(figsize=(10, 2*rank), constrained_layout=False)
    fig.suptitle(f"Spatial Maps, Time Series, Echo Profiles: {rank} Components", y=0.89)
    outer_grid = fig.add_gridspec(rank, 2)

    for component in range(rank):
        left_cell = outer_grid[component, 0]
        right_cell = outer_grid[component, 1]

        inner_grid_left = left_cell.subgridspec(1, 3)
        inner_grid_right = right_cell.subgridspec(2, 1)
        
        ax_spatial_1 = fig.add_subplot(inner_grid_left[0, 0])
        ax_spatial_2 = fig.add_subplot(inner_grid_left[0, 1])
        ax_spatial_3 = fig.add_subplot(inner_grid_left[0, 2])

        ax_time = fig.add_subplot(inner_grid_right[0, 0])
        ax_echo = fig.add_subplot(inner_grid_right[1, 0])

        spatial_data = spatial_maps[:, component]
        unmasked_data = utils.unmask(spatial_data.T, mask > 0)
        brain = np.reshape(unmasked_data,
                           (vol_dims[0], vol_dims[1], vol_dims[2]))
        
        vmin = np.percentile(spatial_data, 5)
        vmax = np.percentile(spatial_data, 95)

        ax_spatial_1.imshow(np.rot90(brain[n_slice, :, :], slice_rot),
                            vmin=vmin, vmax=vmax, aspect="auto")
        ax_spatial_2.imshow(np.rot90(brain[:, n_slice, :], slice_rot),
                            vmin=vmin, vmax=vmax, aspect="auto")
        ax_spatial_3.imshow(np.rot90(brain[:, :, n_slice], slice_rot),
                            vmin=vmin, vmax=vmax, aspect="auto")

        ax_spatial_1.axis('off')
        ax_spatial_1.annotate(component, xy=(0, 1), xycoords='axes fraction',
                              xytext=(+0.5, -0.5),
                              textcoords='offset fontsize',
                              fontsize='medium', verticalalignment='top',
                              fontfamily='serif',
                              bbox=dict(facecolor='0.7',
                                        edgecolor='none', pad=3.0))

        ax_spatial_2.axis('off')
        ax_spatial_3.axis('off')
        
        ax_time.plot(time_series[:, component], color="limegreen")
        ax_time.tick_params(bottom=False, labelbottom=False)
        ax_time.annotate(f"{round((weights[component]/total_weights)*100, 1)}%", xy=(0, 1),
                         xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',
                         fontsize='medium', verticalalignment='top', fontfamily='serif',
                         bbox=dict(facecolor='0.7', edgecolor='none', pad=3.0))

        ax_echo.plot(echo_profiles[:, component], color="purple")
        ax_echo.tick_params(bottom=False, labelbottom=False)

        if component == 0:
            ax_spatial_1.set_title("Spatial Maps")
            ax_time.set_title("Time Series, Echo Profiles")

        ax_time.margins(x=0)
        ax_echo.margins(x=0)
    
    plt.show()
    # plt.savefig(filepath, bbox_inches='tight')


def plot_echo_profiles(all_factors, n_comps, echo_times, iter_num, deg=1):
    """
    Plot plot echo profiles.

    Parameters
    ----------
    all_factors : list
    n_comps : int
    echo_times : ndarray
    iter_num : int

    Returns
    -------
    r_squared_vals : ndarray
    """
    fig_rows = int(n_comps/10)
    r_squared_vals = [[[] for _ in range(n_comps)] for _ in range(iter_num)]

    if deg == 2:
        def func_1(x, a, b):
            return a*(x**2) + b
        
        def func_2(x, a, b):
            return a*x + b
        
    y_vmin = np.max(all_factors[0][2])
    y_vmax = np.min(all_factors[0][2])

    for iter in range(iter_num):
        factors = all_factors[iter]
        figure, axs = plt.subplots(fig_rows, 10, figsize=(25, 8))
        figure.suptitle("Echo Profiles with Best Fit Line and Statistics")
        figure.tight_layout()

        for component in range(n_comps):
            col = component % 10
            ax = axs[int(component/10), col]
            y = factors[2][:, component]

            ax.set_xlim(echo_times[0], echo_times[len(echo_times)-1])
            ax.set_ylim(y_vmin, y_vmax)
            #ax.set_ylim(-0.6, 1.0)
            ax.plot(echo_times, y)

            if deg == 1:
                slope, intercept, r, p, se = linregress(echo_times, y)
                r_squared_vals[iter][component] = r**2

                ax.plot(echo_times, slope*echo_times + intercept)
                intercept = np.sign(slope)*intercept
                slope = np.abs(slope)*10
                se = se*100

                textstr = f"slope/se: {round(slope/se, 2)} \n intercept/se: {round(intercept/se, 2)}"
    
            elif deg == 2:
                popt_1, _, infodict_1, _, _ = curve_fit(func_1, echo_times, y,
                                                        full_output=True)
                popt_2, _, infodict_2, _, _ = curve_fit(func_2, echo_times, y,
                                                        full_output=True)

                sst = np.sum(y**2)
                r_squared_1 = 1 - np.sum(infodict_1["fvec"]**2)/sst
                r_squared_2 = 1 - np.sum(infodict_2["fvec"]**2)/sst

                if np.sum(infodict_1["fvec"]**2) > np.sum(infodict_2["fvec"]**2):
                    ax.plot(echo_times,
                            func_1(echo_times, popt_1[0], popt_1[1]))
                    r_squared_vals[iter][component] = r_squared_1

                else:
                    ax.plot(echo_times, func_2(echo_times, popt_2[0], popt_2[1]))
                    r_squared_vals[iter][component] = r_squared_2
                
                textstr = f"slope/se: {round(slope/se, 2)} \n intercept/se: {round(intercept/se, 2)}"

            else:
                print("Can't plot for degree {deg}")
            
            # these are matplotlib.patch.Patch properties
            props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.5)
            num_label = f"{component}"
            props2 = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
            
            # place a text box in upper left in axes coords
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)
            ax.text(0.02, 0.97, num_label, transform=ax.transAxes, fontsize=10,
                    verticalalignment='bottom', bbox=props2)
        plt.show()

    return r_squared_vals


def compare_rsqared(r_squared_vals, n_comps, iter_num):
    """
    Compare values of R^2 across iterations.

    Parameters
    ----------
    r_squared_vals : ndarray
    n_comps : int
    iter_num : int
    """
    figure, axs = plt.subplots(2, 1)
    x_vals = np.linspace(0, n_comps, n_comps)
    for iter in range(iter_num):
        ax1 = axs[0]
        ax1.scatter(x_vals, r_squared_vals[iter],
                    label=f"R^2 for iteration {iter}")
        ax1.set_title("R^2 Values for each component")
        ax1.legend()

    ax2 = axs[1]
    ax2.plot(np.array(r_squared_vals))
    ax2.set_title("R^2 Values over iterations")

    plt.legend(bbox_to_anchor=(1, 0.5))
    plt.show()


def pca_dim_reduced(comp_ts, voxel_comp_weights, varex, n_comps): 
    """
    Use a subset of PCA components to create v x t data with reduced dimensionality.

    Parameters
    ----------
    comp_ts : ndarray
    voxel_comp_weights : ndarray
    varex : ndarray
    n_comps : int

    Returns
    -------
    data_reduced : ndarray
    """
    # Create a boolean array where the first n_comps values are true
    keep_comps = np.zeros(comp_ts.shape[1], dtype=bool)
    keep_comps[:n_comps] = True
    voxel_kept_comp_weighted = voxel_comp_weights[:, keep_comps] * varex[None, keep_comps]
    data_reduced = np.dot(voxel_kept_comp_weighted, comp_ts[:, keep_comps].T)
    
    return data_reduced


def get_new_slopes(factors, rank, echo_times, deg=1, plot=False):
    """
    Get new slopes to initialize tensor with.

    Parameters
    ----------
    factors : list
    rank : int
    echo_times : ndarray
    deg : int
    plot : bool

    Returns
    -------
    new_echo_slopes.T : ndarray
    """
    if deg == 2:
        def func_1(x, a, b):
            return a*(x**2) + b
        
        def func_2(x, a, b):
            return a*x + b
    
    for component in range(rank):
        y = factors[2][:, component]
        
        if deg == 1:
            slope, intercept, r, p, se = linregress(echo_times, y)

            if component == 0:
                new_echo_slopes = slope*echo_times + intercept
            else:
                new_slope = slope*echo_times + intercept
                new_echo_slopes = np.vstack((new_echo_slopes, new_slope))

        elif deg == 2:
            popt_1, _, infodict_1, _, _ = curve_fit(func_1, echo_times,
                                                    y, full_output=True)
            popt_2, _, infodict_2, _, _ = curve_fit(func_2, echo_times,
                                                    y, full_output=True)

            if np.sum(infodict_1["fvec"]**2) > np.sum(infodict_2["fvec"]**2):
                if component == 0:
                    new_echo_slopes = func_1(echo_times, popt_1[0], popt_1[1])
                else:
                    new_slope = func_1(echo_times, popt_1[0], popt_1[1])
                    new_echo_slopes = np.vstack((new_echo_slopes, new_slope))

            else:
                if component == 0:
                    new_echo_slopes = func_2(echo_times, popt_2[0], popt_2[1])
                else:
                    new_slope = func_1(echo_times, popt_2[0], popt_2[1])
                    new_echo_slopes = np.vstack((new_echo_slopes, new_slope))

        else:
            print("Can't get slopes for degree {deg}")

        if plot is True:
            for row in range(new_echo_slopes.shape[0]):
                plt.plot(new_echo_slopes[row, :])

        plt.show()

    return new_echo_slopes.T


def initialize_cp_with_ica(data, rank, mixing, space, echos=None):
    """
    Initialize a tensor with the desired properties.

    Parameters
    ----------
    data : ndarray
    rank : int
    mixing : ndarray
    space : ndarray
    echos : bool

    Returns
    -------
    init_tensor : cp.tensor
    """
    # choose between random init and fixed modes init
    init_tensor = random_cp(data.shape, rank, full=False, random_state=0)
    print(f"initialization tensor is shape: {init_tensor.shape}")
    # mixing_mc = mixing - np.mean(np.mean(mixing,axis=1, keepdims=True),
    # axis=0, keepdims=True)
    # mixing_scaled = mixing_mc/np.std(mixing_mc)

    init_tensor[1][1][:, :] = mixing
    init_tensor[1][0][:, :] = space

    if isinstance(echos, np.ndarray):
        init_tensor[1][2] = echos

    return cp_normalize(init_tensor)


def pca_preproc(oc_data, echo_data, n_comps, verbose =True):
    """
    Preprocess echo data and optimally combined data and run PCA on optimally combined data.

    Parameters
    ----------
    oc_data : ndarray
    echo_data : int
    n_comps : int
    verbose : bool

    Returns
    -------
    data_reduced : ndarray
    percent_change_echo_data : ndarray
    comp_ts : ndarray
    """
    n_vols = oc_data.shape[1]

    # preprocess optimally combined data
    if verbose:
        print("Preprocessing data...")
    oc_data_z = ((oc_data - oc_data.mean(axis=1, keepdims=True)) / oc_data.std(axis=1, keepdims=True))
    oc_data_z = (oc_data_z - oc_data_z.mean()) / oc_data_z.std()  # var normalize everything

    mean_centered_echo_data = echo_data - np.mean(echo_data, axis=1, keepdims=True)
    percent_change_echo_data = mean_centered_echo_data /np.mean(echo_data, axis=1, keepdims=True)
    #scaled_by_mean_echo_data = echo_data/ np.mean(echo_data, axis=1, keepdims=True)

    if verbose:
        print(f"OC data preprocessed by demeaning and dividing by standard deviation acorss both axes")
        print(f"Echo data preprocessed by:")

    # Run PCA using SVD
    if verbose:
        print("Running PCA on optimally combined data...")
        # print(f"Dimension of preprocessed, optimally combined data for PCA is {oc_data_z.shape}, voxels x time")

    ppca = PCA(copy=False, n_components=n_vols - 1, svd_solver="full")
    ppca.fit(oc_data_z)

    comp_ts = ppca.components_.T  # PCA components 
    varex = ppca.explained_variance_  # this is the variance explained by each principle component
    voxel_comp_weights = np.dot(np.dot(oc_data_z, comp_ts), np.diag(1.0 / varex))  # multiplies normalized data with pca components 
    # varex_norm = ppca.explained_variance_ratio_ # this is the percentage of total variance explained by each principle component

    # dimensionality reduced data
    data_reduced = pca_dim_reduced(comp_ts, voxel_comp_weights, varex, n_comps)  # run pca_dim_reduced on data

    return percent_change_echo_data, data_reduced, comp_ts


def ica_preproc(oc_data, echo_data, n_comps, verbose=True):
    """
    Preprocess echo data and optimally combined data and run ICA on optimally combined data.

    Parameters
    ----------
    oc_data : ndarray
    echo_data : int
    n_comps : int
    verbose : bool

    Returns
    -------
    percent_change_echo_data : ndarray
    data_reduced : ndarray
    mixing : ndarray
    space : ndarray
    """
    n_vols = oc_data.shape[1]

    # preprocess optimally combined data
    if verbose:
        print("Preprocessing data...")
    oc_data_z = ((oc_data - oc_data.mean(axis=1, keepdims=True)) / oc_data.std(axis=1, keepdims=True))
    oc_data_z = (oc_data_z - oc_data_z.mean()) / oc_data_z.std()  # var normalize everything

    mean_centered_echo_data = echo_data - np.mean(echo_data, axis=1, keepdims=True)
    percent_change_echo_data = mean_centered_echo_data / np.mean(echo_data, axis=1, keepdims=True)
    # scaled_by_mean_echo_data = echo_data/ np.mean(echo_data, axis=1, keepdims=True)

    if verbose:
        print(f"OC data preprocessed by demeaning and dividing by standard deviation acorss both axes")
        print(f"Echo data preprocessed by:")

    # Run PCA using SVD
    if verbose:
        print("Running PCA on optimally combined data...")
        # print(f"Dimension of preprocessed, optimally combined data for PCA is {oc_data_z.shape}, voxels x time")

    ppca = PCA(copy=False, n_components=n_vols - 1, svd_solver="full")
    ppca.fit(oc_data_z)

    comp_ts = ppca.components_.T  # PCA components 
    varex = ppca.explained_variance_  # this is the variance explained by each principle component
    voxel_comp_weights = np.dot(np.dot(oc_data_z, comp_ts), np.diag(1.0 / varex))  # multiplies normalized data with pca components 
    # varex_norm = ppca.explained_variance_ratio_ # this is the percentage of total variance explained by each principle component

    # dimensionality reduced data
    data_reduced = pca_dim_reduced(comp_ts, voxel_comp_weights, varex, n_comps) #run pca_dim_reduced on data

    # run ICA
    if verbose:
        print(f"Running ICA with {n_comps} components on optimally combined data ...")
    ica = FastICA(n_comps, algorithm="parallel", whiten="arbitrary-variance", fun="logcosh", max_iter=200, tol=1e-4, random_state=42)
    # ICA_spatial = ica.fit_transform(data_reduced)
    ica.fit_transform(data_reduced)

    # save to set initialization tensor
    mixing = ica.mixing_
    space = np.dot(data_reduced, ica.components_.T)

    return percent_change_echo_data, data_reduced, space, mixing


def parafac_iter(pre_proc_data, rank, mixing, space, echo_times, iter_num, ortho=True, deg=1, verbose=True): 
    """
    Iterate PARAFAC decomposition.

    Parameters
    ----------
    pre_proc_data : ndarray
    rank : int
    mixing : ndarray
    space : ndarray
    echos_times : ndarray
    iter_num : int
    deg : int
    verbose : bool

    Returns
    -------
    all_weights : list
    all_factors : list
    """
    all_weights = [[] for _ in range(iter_num)]
    all_factors = [[[] for _ in range(3)] for _ in range(iter_num)]

    # init = initialize_cp(pre_proc_data, rank, mixing, space, echos=None)
    init = initialize_cp_with_ica(pre_proc_data, rank, mixing, space, echos=None)

    # might need to change normalize_factors back to false
    for iter in range(iter_num):
        if verbose:
            print(f"PARAFAC Iteration {iter+1}")
        weights, factors = parafac(pre_proc_data, rank, init=init,
                                   svd="truncated_svd",
                                   normalize_factors=True,
                                   orthogonalise=ortho,
                                   tol=1e-5,
                                   verbose=True)

        if deg == 1:
            new_echo_slopes = get_new_slopes(factors, rank, echo_times,
                                             deg=1, plot=False)
        elif deg == 2:
            new_echo_slopes = get_new_slopes(factors, rank, echo_times,
                                             deg=2, plot=False)
        else:
            print("Can't iterate PARAFAC with echo slopes of degree {deg}")

        mixing = factors[1]
        space = factors[0]

        init = initialize_cp_with_ica(pre_proc_data, 
                                      rank, mixing, space, 
                                      new_echo_slopes)

        # init = initialize_cp(pre_proc_data, rank, mixing, space,
        #                      new_echo_slopes)

        all_weights[iter] = weights
        all_factors[iter] = factors

    return all_weights, all_factors


def ica_to_parafac_iter(oc_data, echo_data, n_comps, echo_times,
                        iter_num, ortho=True, deg=1, verbose=True):
    """
    Run ICA and iterate PARAFAC decomposition.

    Parameters
    ----------
    pre_proc_data : ndarray
    rank : int
    mixing : ndarray
    space : ndarray
    echos_times : ndarray
    iter_num : int
    deg : int
    verbose : bool

    Returns
    -------
    all_weights : list
    all_factors : list
    """
    n_vols = oc_data.shape[1]

    # preprocess optimally combined data
    if verbose:
        print("Preprocessing data...")
    oc_data_z = ((oc_data - oc_data.mean(axis=1, keepdims=True)) /
                 oc_data.std(axis=1, keepdims=True))
    oc_data_z = (oc_data_z - oc_data_z.mean()) / oc_data_z.std()

    mean_centered_echo_data = echo_data - np.mean(echo_data, axis=1,
                                                  keepdims=True)
    percent_change_echo_data = mean_centered_echo_data / np.mean(echo_data,
                                                                 axis=1,
                                                                 keepdims=True)
    # scaled_by_mean_echo_data =
    # echo_data/ np.mean(echo_data, axis=1, keepdims=True)

    if verbose:
        print("OC data preprocessed by demeaning and dividing by standard "
              "deviation acorss both axes")
        print("Echo data preprocessed by:")

    # Run PCA using SVD
    if verbose:
        print("Running PCA on optimally combined data...")
        # print(f"Dimension of preprocessed, optimally combined data for PCA 
        # is {oc_data_z.shape}, voxels x time")

    ppca = PCA(copy=False, n_components=n_vols - 1, svd_solver="full")
    ppca.fit(oc_data_z)

    comp_ts = ppca.components_.T  # PCA components
    varex = ppca.explained_variance_  # variance explained by pcs
    voxel_comp_weights = np.dot(np.dot(oc_data_z, comp_ts), 
                                np.diag(1.0 / varex))
    # dimensionality reduced data
    data_reduced = pca_dim_reduced(comp_ts, voxel_comp_weights, varex, n_comps)

    # run ICA
    if verbose:
        print(f"Running ICA with {n_comps}"
              "components on optimally combined data ...")
    ica = FastICA(n_comps, algorithm="parallel", whiten="arbitrary-variance",
                  fun="logcosh", max_iter=200, tol=1e-4, random_state=42)
    # ICA_spatial = ica.fit_transform(data_reduced)
    ica.fit_transform(data_reduced)

    # save to set initialization tensor
    mixing = ica.mixing_
    space = np.dot(data_reduced, ica.components_.T)

    if verbose:
        print("Setting up for PARAFAC decomposition...")
        print(f"Spatial map matrix is size: {space.shape}")
        print(f"Time series matrix is size: {mixing.shape}")

    if verbose:
        print("Running PARAFAC decomposition...")

    echo_data_preproc = percent_change_echo_data

    if deg == 1:
        all_weights, all_factors = parafac_iter(echo_data_preproc, n_comps,
                                                mixing, space, echo_times,
                                                iter_num, ortho=ortho, deg=1,
                                                verbose=verbose)
    elif deg == 2:
        all_weights, all_factors = parafac_iter(echo_data_preproc, n_comps,
                                                mixing, space, echo_times,
                                                iter_num, ortho=ortho, deg=2,
                                                verbose=verbose)
    else:
        print("Can't run PARAFAC with echo slopes of degree {deg}")

    return mixing, all_weights, all_factors


def save_with_headers(matrix_to_save, drive_loc, folder, filename,
                      n_components, preprocess=False):
    """
    Save mixing matrix with headers to run tedana.

    Parameters
    ----------
    matrix_to_save : ndarray
    drive_loc : string
    folder : string
    filename : string
    n_components : int
    preprocess : bool
    """
    # zscore data
    if preprocess:
        matrix_to_save = (matrix_to_save - np.mean(matrix_to_save,
                                                   axis=0, keepdims=True)) / np.std(matrix_to_save, axis=0, keepdims=True)
        print(f"{np.mean(matrix_to_save, axis=0)}")
        print(f"{np.std(matrix_to_save, axis=0)}")

    # Add header and save file as testmixing.tsv to current folder
    header_str = ""
    for i in range(n_components):
        header_str += f"ICA_{(i+1):03d}\t"
    print(len(header_str))
    header_str = header_str.strip()
    print(len(header_str))
    print(header_str)
    np.savetxt(f"{drive_loc}/{folder}/{filename}.tsv", matrix_to_save,
               header=header_str, fmt="%.6f",
               delimiter="\t", comments='')
    
    print(f"File saved to {drive_loc}/{folder}!")