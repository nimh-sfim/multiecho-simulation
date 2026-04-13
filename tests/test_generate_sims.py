import os.path as op
import pytest
import sys
import numpy as np
from scipy import signal

sys.path.append(op.join(op.dirname(op.abspath(__file__)), "..", "simlib"))

import generate_sims as gs


def test_monoexponential():
    tes = [10, 20, 30]
    s0 = 100
    r2star = 0.05
    # basic test
    expected_result = s0 * np.exp(-np.array(tes) * r2star)
    result = gs.monoexponential(tes, s0, r2star)
    assert np.allclose(result, expected_result)

    # test with t2star input
    result = gs.monoexponential(tes, s0, 1 / r2star, t2star_input=True)
    assert np.allclose(result, expected_result)

    # Test with 1D tes and 1D s0 and r2star with tiling
    s0 = [100, 200]
    r2star = [0.04, 0.05]
    expected_result = np.zeros((2, 3))
    for idx in range(2):
        expected_result[idx, :] = s0[idx] * np.exp(-np.array(tes) * r2star[idx])
    result = gs.monoexponential(tes, s0, r2star)
    assert np.allclose(result, expected_result)

    # Test with 1D tes and 2D s0 and r2star with tiling
    s0 = np.array([[100, 125, 175, 200], [300, 325, 375, 400]]).T
    r2star = np.array([[0.04, 0.043, 0.047, 0.05], [0.02, 0.023, 0.027, 0.03]]).T
    expected_result = np.zeros((4, 2, 3))
    for idx1 in range(4):
        for idx2 in range(2):
            expected_result[idx1, idx2, :] = s0[idx1][idx2] * np.exp(
                -np.array(tes) * r2star[idx1][idx2]
            )
    result = gs.monoexponential(tes, s0, r2star)
    assert np.allclose(result, expected_result)

    tes = np.array([[10, 15, 20], [30, 35, 40]])
    with pytest.raises(
        ValueError,
        match="If tes has more than 1 dimension, then s0 and r2star should be single values.",
    ):
        gs.monoexponential(tes, s0, r2star)

    tes = [10, 20, 30]
    s0 = [100, 200]
    with pytest.raises(ValueError, match="s0 and r2star must be the same shape."):
        gs.monoexponential(tes, s0, r2star)
