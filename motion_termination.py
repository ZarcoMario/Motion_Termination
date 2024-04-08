import numpy as np
from scipy.signal import find_peaks


def termination_detection(x, x_fin, vx, vy, vz, v_th=-0.1):
    '''
    TODO: Comment this method and explain it better
    This function evaluates whether a section of the trajectory is to be discarded
        and calculates the last index if this is the case
        NOTE: This function assumes that the location of the targets
            is different along the x dimension in UNITY
    :param x: dimension tha differs
    :param x_fin:
    :param vx: first derivative of x.
    :param vy: first derivative of y
    :param vz: first derivative of z
    :param v_th: speed threshold below which the function assumes there are minima of speed
    :return:
    '''
    idx_i = np.argwhere(np.abs(vx) == np.max(np.abs(vx))).T[0][0]

    speed = np.sqrt(vz ** 2 + vy ** 2)
    peaks, _ = find_peaks(-speed[idx_i:])
    peaks += idx_i

    if np.any(vz[idx_i:] < v_th) and np.any(vy[idx_i:] < v_th) and peaks.size >= 2:

        idx_max = np.argwhere(np.abs(vx) == np.max(np.abs(vx))).T[0][0]

        indexes = np.argwhere(np.diff(np.sign(vx[idx_max:]))).T[0] + 1
        if indexes.size != 0:
            indexes += idx_max
            idx_termination = indexes[0]
            if np.sign(x[idx_termination]) == np.sign(x_fin):
                return idx_termination
            else:
                return vx.size
        else:
            return vx.size  # - 1

    else:
        return vx.size
