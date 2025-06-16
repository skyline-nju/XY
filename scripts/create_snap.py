import sys
import numpy as np
import matplotlib.pyplot as plt
from gsd import hoomd, fl
import os
import sys


def get_frame(fname_in, i_frame=-1, flag_show=False):
    i_frame = -1
    s = hoomd.Frame()
    with hoomd.open(name=fname_in, mode='r') as fin:
        nframes = len(fin)
        if i_frame < 0:
            i_frame += nframes
        s = fin[i_frame]

    if flag_show:
        # omega = s.particles.position[:, 2]
        theta = s.particles.charge
        plt.imshow(theta.reshape(512, 512), origin="lower")
        plt.show()
        plt.close()
    return s


def reset_theta(snap, theta0="rand"):
    """ reset the theta array to theta0.
        Caution: error when theta0 = 0.

    Args:
        snap (hoomd.Frame): Hoomd snapshot
        theta0 (str, optional): _description_. Defaults to "rand", such that
                                theta is uniformly distributed from -pi to pi. 

    Returns:
        hoomd.Frame: new snapshot 
    """
    s2 = hoomd.Frame()
    s2.configuration.box = snap.configuration.box
    s2.configuration.step = 0
    s2.particles.N = snap.particles.N
    s2.particles.position = snap.particles.position

    N = snap.particles.N
    charge = np.zeros(N, dtype=np.float32)
    if theta0 == "rand":
        rng = np.random.default_rng()
        charge = (rng.random(N) - 0.5) * np.pi * 2
    elif isinstance(theta0, float):
        for i in range(charge.size):
            charge[i] = theta0
    else:
        print("Error, theta0 must be 'rand' or a float number")
        sys.exit(1)
    s2.particles.charge = charge
    return s2


if __name__ == "__main__":
    folder = "/mnt/sda/LatticeXY/KM/snap"
    seed = 1001
    basename = f"L512_512_T0.1_s0.2_h0.1_S{seed:d}.gsd"
    fname_in = f"{folder}/{basename}"

    snap = get_frame(fname_in)
    snap2 = reset_theta(snap, theta0=np.pi/2)

    # plt.imshow(snap2.particles.charge.reshape(512, 512))
    # plt.show()
    # plt.close()
    fname_out = f"{folder}/L512_512_T0.1_s0.2_h0.1_S10013000.gsd"
    print(snap2.particles.charge[0])
    with hoomd.open(name=fname_out, mode="w") as fout:
        fout.append(snap2)
