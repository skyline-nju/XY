import numpy as np
import matplotlib.pyplot as plt
from gsd import hoomd
from plot_fields import untangle_2D
from create_snap import get_frame
from add_line import add_line


def cal_Sq(snap:hoomd.Frame, dx:float=1.):
    box = snap.configuration.box
    Lx, Ly = int(box[0]), int(box[1])

    n = int(Lx / dx)
    qx = np.fft.fftfreq(n, d=dx/(2 * np.pi))
    qy = np.fft.fftfreq(n, d=dx/(2 * np.pi))
    q_radius = qx[1:n//2]
    dq = qx[1] - qx[0]
    qx = np.fft.fftshift(qx)
    qy = np.fft.fftshift(qy)
    qx_ij, qy_ij = np.meshgrid(qx, qy)
    q_module = np.sqrt(qx_ij **2 + qy_ij ** 2)

    theta = snap.particles.charge.reshape(Ly, Lx)
    ux = np.cos(theta)
    uy = np.sin(theta)
    theta_untangled = untangle_2D(theta)

    ux_q = np.fft.fft2(ux, norm="ortho")
    uy_q = np.fft.fft2(uy, norm="ortho")
    Sq_u = ux_q * ux_q.conj() + uy_q * uy_q.conj()

    theta_q = np.fft.fft2(theta_untangled, norm="ortho")
    Sq_theta = theta_q * theta_q.conj()

    Sq_u = np.fft.fftshift(Sq_u).real
    Sq_theta = np.fft.fftshift(Sq_theta).real
    Sq_u_r = np.zeros(q_radius.size)
    Sq_theta_r = np.zeros(q_radius.size)

    for j in range(Sq_u_r.size):
        mask = np.logical_and(q_module > q_radius[j]-dq/2, q_module <= q_radius[j]+dq/2)
        Sq_u_r[j] = np.mean(Sq_u[mask])
        Sq_theta_r[j] = np.mean(Sq_theta[mask])
    return q_radius, Sq_u_r, Sq_theta_r


if __name__ == "__main__":
    folder = "/mnt/sda/LatticeXY/KM/snap"
    seed = 1000
    fname = f"{folder}/L512_512_T0.1_s0.2_h0.1_S{seed:d}.gsd"

    snap = get_frame(fname)

    q, Sq_u, Sq_theta = cal_Sq(snap)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True)
    ax.plot(q, Sq_u, label=r"$\langle |\tilde{u}(\mathbf{k})|^2\rangle$")
    ax.plot(q, Sq_theta, label=r"$\langle |\tilde{\theta}(\mathbf{k})|^2\rangle$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$k$", fontsize="x-large")
    ax.set_ylabel(r"$S_k$", fontsize="x-large")
    ax.legend(fontsize="x-large")
    add_line(ax, 0, 1, 1, -4, label=r"$-4$", yl=0.63)
    plt.show()
    plt.close()