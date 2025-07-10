import numpy as np
import matplotlib.pyplot as plt
from gsd import fl
from scipy.ndimage import gaussian_filter


def untangle_1D(theta):
    phase = 0
    theta_new = np.zeros_like(theta)
    theta_new[0] = theta[0]
    for i in range(1, theta.size):
        dtheta = theta[i] - theta[i-1]
        if dtheta < -np.pi:
            phase += 2 * np.pi
        elif dtheta >= np.pi:
            phase -= 2 * np.pi
        theta_new[i] = theta[i] + phase
    return theta_new


def untangle_2D(theta):
    phase_y = 0
    theta_new = np.zeros_like(theta)

    nrows, ncols = theta.shape
    theta_new[0] = untangle_1D(theta[0])
    for j in range(1, nrows):
        dtheta_y = theta[j , 0] - theta[j-1, 0]
        if dtheta_y < -np.pi:
            phase_y += 2 * np.pi
        elif dtheta_y >= np.pi:
            phase_y -= 2 * np.pi
        theta_new[j] = untangle_1D(theta[j]) + phase_y
    return theta_new


def verify_untangled_angles(theta_untangled):
    nrows, ncols = theta_untangled.shape
    x = np.arange(ncols) + 0.5

    for j in range(nrows):
        dtheta = theta_untangled[j, -1] - theta_untangled[j, 0]
        if dtheta < -np.pi or dtheta >= np.pi:
            plt.plot(x, theta_untangled[j])
            plt.show()
            plt.close()


def get_untangled_angles(ux, uy, sigma=0):
    if sigma == 0:
        theta = np.arctan2(uy, ux)
        theta_untangled = untangle_2D(theta)
    return theta_untangled


def get_theta_fields(fname, i_frame=-1, ret_untangled=True):
    with fl.open(name=fname, mode="r") as fin:
        nframes = fin.nframes
        if i_frame < 0:
            i_frame += nframes
        box = fin.read_chunk(frame=0, name="configuration/box")
        Lx, Ly = int(box[0]), int(box[1])
        theta = fin.read_chunk(frame=i_frame, name="particles/charge").reshape(Ly, Lx)
    if ret_untangled:
        theta_untangled = untangle_2D(theta)
        return theta, theta_untangled
    else:
        return theta


if __name__ == "__main__":
    # folder = "/mnt/sda/LatticeXY/KM/snap"
    # fname = f"{folder}/L1024_1024_T0.1_s0.25_h0.1_S3002.gsd"

    folder = "/mnt/sda/LatticeXY/KM/L4096"
    fname = f"{folder}/L4096_4096_T0.1_s0.25_h0.1_S3000.gsd"
    with fl.open(name=fname, mode="r") as fin:
        nframes = fin.nframes
        print("found", nframes, "frames")
        box = fin.read_chunk(frame=0, name="configuration/box")
        Lx, Ly = int(box[0]), int(box[1])
        
        for i_frame in range(nframes-3, nframes):
            theta = fin.read_chunk(frame=i_frame, name="particles/charge").reshape(Ly, Lx)
            theta_untangled = untangle_2D(theta)

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
            ax1.imshow(theta, origin="lower", cmap="hsv")
            ax2.imshow(theta_untangled, origin="lower")

            plt.show()
            plt.close()
