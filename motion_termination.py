import numpy as np
from scipy.signal import find_peaks


def termination_detection(vx, vy, vz, t, v_th=-0.2):
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

    if np.all(speed[peaks] < np.abs(v_th)) and peaks.size >= 2:

        indexes = np.argwhere(np.diff(np.sign(vx[idx_i:]))).T[0] + 1

        if indexes.size != 0:
            indexes += idx_i
            idx_termination = indexes[0]

            # import matplotlib.pyplot as plt
            # from matplotlib.gridspec import GridSpec
            #
            # fig = plt.figure(figsize=(16, 8))
            # gs = GridSpec(5, 1)
            #
            # ax = fig.add_subplot(gs[0, 0])
            # ax.grid(True)
            # ax.plot(t[idx_i:], vx[idx_i:], '.')
            #
            # ax = fig.add_subplot(gs[1, 0])
            # ax.grid(True)
            # ax.plot(t[idx_i:], np.sign(vx[idx_i:]), '.')
            #
            # ax = fig.add_subplot(gs[2, 0])
            # ax.grid(True)
            # ax.plot(np.diff(np.sign(vx[idx_i:])), '.')
            #
            # ax = fig.add_subplot(gs[3, 0])
            # ax.grid(True)
            # ax.plot(t, vx, 'b.')
            # ax.plot(t[indexes], vx[indexes], 'r.')
            #
            # ax = fig.add_subplot(gs[4, 0])
            # ax.grid(True)
            # ax.plot(t, vy, '.', label='vy')
            # ax.plot(t, vz, '.', label='vz')
            # ax.plot(t, speed, '.', label='s')
            # ax.plot(t[indexes], speed[indexes], '.')
            # ax.plot(t[peaks], speed[peaks], 'ro')
            # ax.legend()
            # # plt.show()

            return idx_termination
        else:
            # print(indexes)
            return vx.size # - 1

    else:

        idx_termination = vx.size # - 1
        return idx_termination
