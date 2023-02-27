import pytest
from scipy import io
import numpy as np

def check_two_matrices(m1, m2, eps=1e-6, transpose=True):
    if (transpose):
        m2 = m2.T
    
    assert m1.shape == m2.shape
    assert m1.dtype == m2.dtype
    assert np.amax(np.abs(m1-m2)) <= eps
    return True


def test_py_vs_R():
    m1 = io.mmread('checkpoints/checkpoint_1_vp.mtx')
    m2 = io.mmread('checkpoints/checkpoint_1_vr.mtx')
    assert check_two_matrices(m1, m2)