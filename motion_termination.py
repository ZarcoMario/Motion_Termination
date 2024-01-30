import numpy as np
from scipy.signal import find_peaks


def termination_detection(vx, vy, vz, v_th=-0.1):
    '''
    
    :param vx:
    :param vy:
    :param vz:
    :param v_th:
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
            return idx_termination
        else:
            # print(indexes)
            return vx.size # - 1

    else:

        idx_termination = vx.size # - 1
        return idx_termination
