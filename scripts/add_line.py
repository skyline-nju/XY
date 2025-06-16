import numpy as np


def add_line(ax,
             x_beg,
             y_beg,
             x_end,
             slope,
             label=None,
             xl=None,
             yl=None,
             fontsize="x-large",
             scale="log",
             c="#7f7f7f",
             lw=None,
             deg=None):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if scale == "lin":
        slope_new = slope * (xmax - xmin) / (ymax - ymin)
    else:
        slope_new = slope * (np.log10(xmax / xmin) / np.log10(ymax / ymin))
    x = np.linspace(x_beg, x_end, 100)
    y = slope_new * (x - x_beg) + y_beg
    ax.plot(x, y, "-.", transform=ax.transAxes, color=c, lw=lw)
    if label is not None:
        if deg is None:
            width = ax.bbox.width
            height = ax.bbox.height
            deg = np.arctan(slope_new * height / width) * 180 / np.pi
        dx = x_end - x_beg
        if xl is None:
            xl = x_beg + dx * 0.3
        if yl is None:
            yl = y_beg + dx * 0.6 * slope_new
        ax.text(
            xl,
            yl,
            label,
            transform=ax.transAxes,
            rotation=deg,
            color=c,
            fontsize=fontsize)