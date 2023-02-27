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




def check_pca_matrices(m1, m2, eps=1e-6):
    
    assert m1.shape == m2.shape
    assert m1.dtype == m2.dtype
    
    for i in range(m1.shape[1]):
        n1 = np.linalg.norm(m1[:,i])
        n2 = np.linalg.norm(m2[:,i])

        d = np.dot(m1[:,i], m2[:,i])/(n1*n2)

        assert(abs(d)>= 1-eps)

        for j in range(i):       
            n2 = np.linalg.norm(m2[:,j])    
            d = np.dot(m1[:,i], m2[:,j])/(n1*n2)
            assert(abs(d)<= eps)

    return True

def test_py_vs_R():
    m1 = io.mmread('checkpoints/checkpoint_1_vp.mtx')
    m2 = io.mmread('checkpoints/checkpoint_1_vr.mtx')
    assert check_two_matrices(m1, m2)

def test_py_vs_R_pca():
    m1 = io.mmread('checkpoints/checkpoint_2_vp.mtx')
    m2 = io.mmread('checkpoints/checkpoint_2_vr.mtx').todense()



    assert check_pca_matrices(m1, m2)    