import os.path as op
from unittest import result
import pytest
import sys
import numpy as np
from scipy import signal

sys.path.append(op.join(op.dirname(op.abspath(__file__)), "..", "simlib"))

import generate_sims as gen_sims

def test_tile_variables():
    # test with int + array and no target shape
    vars_to_tile = [4, np.array([1,7,8])]
    expected_result = [np.array([4,4,4]), np.array([1,7,8])]
    result = gen_sims.tile_variables(vars_to_tile)
    assert np.array_equal(result, expected_result)

    # test with a float + array
    vars_to_tile = [4.5, np.array([1,7,8])]
    expected_result = [np.array([4.5, 4.5, 4.5]), np.array([1,7,8])]
    result = gen_sims.tile_variables(vars_to_tile)
    assert np.array_equal(result, expected_result)

    # test with int + array and target shape
    vars_to_tile = [
        2,
        np.array([1,2,3]),
        np.array([[2,6,3],
                  [42,5,16],
                  [7,7,2]]),
        4
    ]
    target_shape = (3,3)
    result = gen_sims.tile_variables(vars_to_tile, target_shape)
    assert result[0].shape == target_shape

    # test 3D arrays
    arr1 = np.arange(24).reshape((2, 3, 4))
    arr2 = np.ones((2, 3, 4))
    vars_to_tile = [arr1, arr2]
    target_shape = (2,3,4)
    result = gen_sims.tile_variables(vars_to_tile, target_shape)
    assert result[0].shape == target_shape

    # test without target shape
    arr1 = np.arange(24).reshape((2, 3, 4))
    arr2 = np.ones((2, 3, 4))
    vars_to_tile = [arr1, arr2]
    largest = max(vars_to_tile, key=lambda x: np.asarray(x).size)
    expected_shape = largest.shape
    result = gen_sims.tile_variables(vars_to_tile)
    assert result[0].shape == expected_shape

    # mix of arrays and int
    arr1 = np.arange(24).reshape((2, 3, 4))
    arr2 = np.zeros((2, 3, 4))
    vars_to_tile = [arr1,4,5.5]
    largest = max(vars_to_tile, key=lambda x: np.asarray(x).size)
    expected_shape = largest.shape
    result = gen_sims.tile_variables(vars_to_tile)
    assert result[0].shape == expected_shape

def test_monoexponential():
    tes = [10, 20]
    s0 = 100
    r2star = 0.05
    # basic test
    expected_result = s0 * np.exp(-np.array(tes) * r2star)
    result = gen_sims.monoexponential(tes, s0, r2star)
    assert np.allclose(result, expected_result)

    # test with t2star input
    result = gen_sims.monoexponential(tes, s0, 1 / r2star, t2star_input=True)
    assert np.allclose(result, expected_result)

    # Test with 1D tes and 1D s0 and r2star with tiling
    s0 = [100, 200]
    r2star = [0.04, 0.05]
    expected_result = np.zeros((2))
    for idx in range(2):
        expected_result[idx] = s0[idx] * np.exp(-np.array(tes[idx]) * r2star[idx])
    result = gen_sims.monoexponential(tes, s0, r2star)
    assert np.allclose(result, expected_result)

    # Separate scaling for tes is currently removed
    # Test with 1D tes and 2D s0 and r2star with tiling
    # s0 = np.array([[100, 125, 175, 200], [300, 325, 375, 400]]).T
    # r2star = np.array([[0.04, 0.043, 0.047, 0.05], [0.02, 0.023, 0.027, 0.03]]).T
    # expected_result = np.zeros((4, 2, 3))
    # for idx1 in range(4):
    #     for idx2 in range(2):
    #         expected_result[idx1, idx2, :] = s0[idx1][idx2] * np.exp(
    #             -np.array(tes) * r2star[idx1][idx2]
    #         )
    # result = gen_sims.monoexponential(tes, s0, r2star)
    # assert np.allclose(result, expected_result)

    # Separate scaling for tes is currently removed
    # tes = np.array([[10, 15, 20], [30, 35, 40]])
    # with pytest.raises(
    #     ValueError,
    #     match="If tes has more than 1 dimension, then s0 and r2star should be single values.",
    # ):
    #     gen_sims.monoexponential(tes, s0, r2star)

    # tes = [10, 20, 30]
    # s0 = [100, 200]
    # with pytest.raises(ValueError, match="s0 and r2star must be the same shape."):
    #     gen_sims.monoexponential(tes, s0, r2star)


def test_calc_delta_r2s_s0_given_s_pchange_proportion():

    inputted_s = 0.1
    delta_s0, delta_r2s = gen_sims.calc_delta_r2s_s0_given_s_pchange_proportion(
        inputted_s=inputted_s,
        s0_mean=100,
        r2s_mean=1 / 28,
        proportion_s0_r2s=0.5,
        te_baseline=28,
    )
    assert len(delta_s0) == 1
    assert len(delta_r2s) == 1

    inputted_s = np.random.randn(3, 4, 5) / 2
    s0_mean = np.array([100, 200, 300])
    r2s_mean = np.array([1 / 28])
    proportion_s0_r2s = np.random.rand(3, 4)
    te_baseline = 28
    delta_s0, delta_r2s = gen_sims.calc_delta_r2s_s0_given_s_pchange_proportion(
        inputted_s=inputted_s,
        s0_mean=s0_mean,
        r2s_mean=r2s_mean,
        proportion_s0_r2s=proportion_s0_r2s,
        te_baseline=te_baseline,
    )
    assert delta_s0.shape == inputted_s.shape
    assert delta_r2s.shape == inputted_s.shape
