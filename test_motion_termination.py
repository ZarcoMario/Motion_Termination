'''
Test 'maximum deviation' using real data
'''
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from motion_termination import termination_detection
from resample import resample_splines
from derivative import calculate_velocity
from filter import butter_lowpass_filter
from movement_onset_detection import onset_detection


# range_ = [i for i in range(1, 44 + 1) if i != 17]

# PARTICIPANTS
for i, p_ in enumerate(range(1, 4 + 1)):

    # VR-S1 Data
    # p_ = 3
    participant_ = r"\P" + str(p_).zfill(2)
    path_ = os.path.dirname(os.getcwd()) + r"\VR-F1" + participant_ + r"\S001"

    path_results = path_ + r"\trial_results.csv"
    results = pd.read_csv(path_results, usecols=['start_time', 'initial_time'])
    start_time = results['start_time'].to_numpy()
    initiation_time = results['initial_time'].to_numpy()
    t_threshold = initiation_time - start_time

    path_fig = os.path.dirname(os.getcwd()) + r"\VR-F1" + r"\Adjusted" + participant_

    if not os.path.exists(path_fig):
        os.mkdir(path_fig)

    print(participant_)

    for trial_number in range(17, 284 + 1, 1):

        path_trial = path_ + r"\trackers" + r"\controllertracker_movement_T" + str(trial_number).zfill(3) + ".csv"

        # Load Raw Data
        raw_data = pd.read_csv(path_trial, usecols=['time', 'pos_x', 'pos_y', 'pos_z'])

        # Adjust to Zero
        t_raw = raw_data['time'].to_numpy() - start_time[trial_number - 1]
        x_raw = raw_data['pos_x'].to_numpy()
        y_raw = raw_data['pos_y'].to_numpy()
        z_raw = raw_data['pos_z'].to_numpy()

        # Resampling
        resampled_data = resample_splines(t_raw, x_raw, y_raw, z_raw)

        t = resampled_data['t'].to_numpy()
        x_res = resampled_data['x'].to_numpy()
        y_res = resampled_data['y'].to_numpy()
        z_res = resampled_data['z'].to_numpy()

        # Filtering.
        # This is a test. A filter might be helpful but this step is not necessary
        cutoff_fq = 10
        x = butter_lowpass_filter(x_res, cutoff_fq, 90, 2)
        y = butter_lowpass_filter(y_res, cutoff_fq, 90, 2)
        z = butter_lowpass_filter(z_res, cutoff_fq, 90, 2)
        # x, y, z = x_res, y_res, z_res

        # Note: step is typically the same for all trials (e.g. if 90 Hz, step=1/90)
        # Although step is similar across trials, step is quickly calculated here
        step = t[1] - t[0]
        vx = calculate_velocity(step, x)
        vy = calculate_velocity(step, y)
        vz = calculate_velocity(step, z)

        # Movement Onset Time Detection
        delta_T = 0.1  # 100 ms.
        Ts = step
        m = int(delta_T / Ts) - 1
        tm = m * Ts

        res = onset_detection(m, x, z, t, vx, vz, t_th=t_threshold[trial_number - 1], vel_th=0.6)
        to = res[0]

        idx_ub = np.argwhere(t > to).T[0][0]
        idx_lb = np.argwhere(t < to).T[0][-1]

        # Find index
        if abs(t[idx_lb] - to) < abs(t[idx_ub] - to):
            idx_i = idx_lb
        else:
            idx_i = idx_ub

        idx_f = termination_detection(vx, vy, vz, t, v_th=-0.2)

        if idx_f != x.size:

            print("T", str(trial_number).zfill(3))

            fig = plt.figure(figsize=(16, 8))
            gs = GridSpec(1, 2)

            ax = fig.add_subplot(gs[0, 0], projection='3d')
            ax.plot(x, z, y, '.', color='grey', alpha=0.5, label='original')
            ax.plot(x[idx_i:idx_f], z[idx_i:idx_f], y[idx_i:idx_f], 'b.', label='processed')
            ax.legend()

            ax = fig.add_subplot(gs[0, 1])
            ax.grid(True)
            ax.plot(x, z, '.', color='grey', alpha=0.5, label='original')
            ax.plot(x[idx_i:idx_f], z[idx_i:idx_f], 'b.', label='processed')
            ax.legend()

            plt.savefig(path_fig + r"\T" + str(trial_number).zfill(3) + ".png")
            # plt.show()
            plt.close()

