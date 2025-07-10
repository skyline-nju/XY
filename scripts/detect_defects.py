import numpy as np
import matplotlib.pyplot as plt
from gsd import fl


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
    return x, y, charge

if __name__ == "__main__":

    folder = "/mnt/sda/LatticeXY/KM/L4096"
    fname = f"{folder}/L4096_4096_T0.1_s0.2_h0.1_S3000.gsd"
    with fl.open(name=fname, mode="r") as fin:
        nframes = fin.nframes
        print("found", nframes, "frames")
        box = fin.read_chunk(frame=0, name="configuration/box")
        Lx, Ly = int(box[0]), int(box[1])
        extent = [-0.5, Lx-0.5, -0.5, Ly-0.5]
                
        for i_frame in range(nframes-2, nframes):
            theta = fin.read_chunk(frame=i_frame, name="particles/charge").reshape(Ly, Lx)

            x, y, charge = detect_defects(theta)

            fig, ax = plt.subplots(1, 1, figsize=(10, 10), constrained_layout=True)
            ax.imshow(theta, cmap="hsv", origin="lower", vmin=-np.pi, vmax=np.pi, extent=extent)

            for j in range(len(charge)):
                if charge[j] == 1:
                    mk = "ko"
                elif charge[j] == -1:
                    mk = "ks"
                ax.plot(x[j], y[j], mk, fillstyle="none", ms=8, mew=2)
            
            if len(charge) > 0:
                charge = np.array(charge)
                n_tot = charge.size
                n_plus = np.sum(charge == 1)
                n_mins = np.sum(charge == -1)
            else:
                n_tot = n_plus = n_mins = 0
            ax.set_title(r"Defect number= %d with $n_{+1}=%d, n_{-1}=%g$" % (n_tot, n_plus, n_mins))

            plt.show()
            plt.close()