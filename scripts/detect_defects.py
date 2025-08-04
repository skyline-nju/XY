import numpy as np
import matplotlib.pyplot as plt
from gsd import fl
import os
import glob
import sys


def detect_defects(theta):
    def get_diff_angle(i1, j1, i2, j2):
        theta1 = theta[j1 % ny, i1 % nx]
        theta2 = theta[j2 % ny, i2 % nx]
        dtheta = theta2 - theta1
        if dtheta >= np.pi:
            dtheta -= 2 * np.pi
        elif dtheta < -np.pi:
            dtheta += 2 * np.pi
        return dtheta

    ny, nx = theta.shape

    x = []
    y = []
    charge = []

    for j in range(ny):
        for i in range(nx):
            delta_theta = get_diff_angle(i, j, i+1, j) + get_diff_angle(i+1, j, i+1, j+1) + get_diff_angle(i+1, j+1, i, j+1) + get_diff_angle(i, j+1, i, j)
            my_charge = int(np.round(delta_theta / (2 * np.pi)))
            if my_charge != 0:
                charge.append(my_charge)
                x.append(i + 0.5)
                y.append(j + 0.5)
    return np.array(x), np.array(y), np.array(charge)


def detect_defects_fast(theta):
    def get_diff_angle(theta1, theta2):
        dtheta = theta2 - theta1
        mask = dtheta  >= np.pi
        dtheta[mask] -= 2 * np.pi
        mask = dtheta < -np.pi
        dtheta[mask] += 2 * np.pi
        return dtheta

    theta_left = np.roll(theta, (1, 0), axis=(1, 0))
    theta_upper_left = np.roll(theta, (1, 1), axis=(1, 0))
    theta_upper = np.roll(theta, (0, 1), axis=(1, 0))

    delta_theta = get_diff_angle(theta, theta_left) + get_diff_angle(theta_left, theta_upper_left) + \
        get_diff_angle(theta_upper_left, theta_upper) + get_diff_angle(theta_upper, theta)
    charge_field = np.round(delta_theta / (2 * np.pi)).astype(int)

    ny, nx = charge_field.shape
    x_1D = np.arange(nx) + 0.5
    y_1D = np.arange(ny) + 0.5
    xx, yy = np.meshgrid(x_1D, y_1D)

    mask = charge_field != 0
    x = xx[mask]
    y = yy[mask]
    charge = charge_field[mask]
    return x, y, charge


def show_defects(L, sigma=None, T=None, seed=None, fname=None, beg_frame=None, savefig=False, fmt=".jpg"):
    folder = f"/mnt/sda/LatticeXY/KM/L{L:d}"
    if fname is None:
        if sigma is None or T is None or seed is None:
            print(L, sigma, T, seed)
            sys.exit(1)
        else:
            fname = f"{folder}/L{L}_{L}_T{T:g}_s{sigma:g}_h0.1_S{seed:d}.gsd"
    if savefig == True:
        outdir = f"{folder}/{os.path.basename(fname).rstrip(".gsd")}"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    with fl.open(name=fname, mode="r") as fin:
        nframes = fin.nframes
        print("found", nframes, "frames")
        box = fin.read_chunk(frame=0, name="configuration/box")
        Lx, Ly = int(box[0]), int(box[1])
        extent = [-0.5, Lx-0.5, -0.5, Ly-0.5]
        
        if beg_frame is None:
            snaps = glob.glob(f"{outdir}/*{fmt}")
            beg_frame = len(snaps)
        elif beg_frame < 0:
            beg_frame += nframes

        if L <= 1024:
            figsize = (9, 8)
        else:
            figsize = (11, 10)

        for i_frame in range(beg_frame, nframes):
            print("frame %d/%d" %(i_frame, nframes))
            theta = fin.read_chunk(frame=i_frame, name="particles/charge").reshape(Ly, Lx)
            try:
                t = fin.read_chunk(frame=i_frame, name="configuration/step")[0] * 0.1
            except KeyError:
                t = 0
            x, y, charge = detect_defects_fast(theta)

            fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)
            theta[theta < 0] += np.pi * 2
            im = ax.imshow(theta, cmap="hsv", origin="lower", vmin=0, vmax=2 * np.pi, extent=extent)

            # print(x, y, charge)
            for j in range(charge.size):
                if charge[j] == 1:
                    mk = "ko"
                elif charge[j] == -1:
                    mk = "ks"
                ax.plot(x[j], y[j], mk, fillstyle="none", ms=8, mew=2)
            
            if charge.size > 0:
                n_tot = charge.size
                n_plus = np.sum(charge == 1)
                n_mins = np.sum(charge == -1)
            else:
                n_tot = n_plus = n_mins = 0
            ax.set_title(r"Defect number= %d with $n_{+1}=%d, n_{-1}=%g, t=%g$" % (n_tot, n_plus, n_mins, t), fontsize="xx-large")
            plt.colorbar(im)

            if savefig:
                if L >= 4096:
                    dpi = 150
                else:
                    dpi = 100
                plt.savefig(f"{outdir}/{i_frame:06d}{fmt}", dpi=dpi)
            else:
                plt.show()
            plt.close()

    
def update_imgs(L):
    folder = f"/mnt/sda/LatticeXY/KM/L{L:d}"
    pat = f"{folder}/*.gsd"
    files = glob.glob(pat)

    for f in files:
        show_defects(L, fname=f, beg_frame=None, savefig=True, fmt=".jpg")

if __name__ == "__main__":
    update_imgs(L=4096)