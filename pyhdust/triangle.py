# -*- coding:utf-8 -*-

"""PyHdust auxiliary module: third-part MCMC plotting tools.

:co-author:Dan Foreman-Mackey
:license: GNU GPL v3.0 https://github.com/danmoser/pyhdust/blob/master/LICENSE
"""
from __future__ import print_function, absolute_import, unicode_literals
import numpy as _np
import warnings as _warn
import pyhdust.phc as _phc
import pyhdust.stats as _stt
from six import string_types as _strtypes

try:
    import matplotlib as _mpl
    import matplotlib.pyplot as _plt
    import matplotlib.cm as _cm
    import matplotlib.gridspec as _gridspec
    from matplotlib.ticker import MaxNLocator as _MaxNLoc
    from scipy.ndimage import gaussian_filter as _gf
    from scipy import stats as _stats
except ImportError:
    gaussian_filter = None
    _warn.warn("matplotlib or scipy not installed!!!")

__author__ = "Daniel Moser"
__email__ = "dmfaes@gmail.com"


def kerncorner(xs, cmapn="gray_r", ncl=3, bestvals="perc", labels=[]):
    """xs = (niteractions, ndim)
    ncl = number os contour lines
    bestvals = None, 'peak' or 'perc'
    """
    ndim = _np.shape(xs)[1]
    if len(labels) < ndim:
        labels += [""] * (ndim - len(labels))

    fig = _plt.figure(figsize=(9, 9))
    gs = _gridspec.GridSpec(ndim, ndim)
    gs.update(hspace=0.01)
    fs = dict(_mpl.rcParams.viewitems())["font.size"] - 4

    axs = []
    xslim = []
    for i in range(ndim):
        ax = []
        xslim += [[_np.min(xs[:, i]), _np.max(xs[:, i])]]
        for j in range(i + 1):
            ax += [_plt.subplot(gs[i, j])]
        axs += [ax]

    for i in range(ndim):
        for j in range(i + 1):
            axs[i][j].locator_params(nbins=5)
            xmin, xmax = xslim[j]
            if i == j:
                kernel = _stats.gaussian_kde(xs[:, i])
                x = _np.linspace(xmin, xmax, 101)
                y = kernel(x)
                axs[i][j].plot(x, y / _np.max(y))
                axs[i][j].set_xlim([xmin, xmax])
                ymin, ymax = axs[i][j].get_ylim()
                if isinstance(bestvals, _strtypes):
                    if bestvals.startswith("perc"):
                        pc = _np.percentile(xs[:, i], [15.9, 50.0, 84.1])
                    elif bestvals.startswith("peak"):
                        mode = x[_phc.find_nearest(y, _np.max(y), idx=True)]
                        mad = _stt.mad(xs[:, i])
                        pc = mode + _np.array([-mad, 0, mad])
                    else:
                        _warn.warn("Invalid `bestvals` in kerncorner", stacklevel=2)
                    if "pc" in locals():
                        dx = xmax - xmin
                        for p in pc:
                            if p < xmin + 0.015 * dx:
                                p = xmin + 0.015 * dx
                            if p > xmax - 0.015 * dx:
                                p = xmax - 0.015 * dx
                            axs[i][j].plot([p, p], [ymin, ymax], "k--")
            else:
                y = xs[:, i]
                x = xs[:, j]
                ymin, ymax = xslim[i]
                #
                values = _np.vstack([x, y])
                kernel = _stats.gaussian_kde(values)
                X, Y = _np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                positions = _np.vstack([X.ravel(), Y.ravel()])
                Z = _np.reshape(kernel(positions), X.shape)
                #
                # axs[i][j].plot(x, y, 'k.')
                axs[i][j].imshow(
                    Z.T,
                    cmap=_plt.get_cmap(cmapn),
                    extent=[xmin, xmax, ymin, ymax],
                    origin="lower",
                )
                if ncl > 0:
                    CS = axs[i][j].contour(
                        X, Y, Z, ncl, colors=_phc.gradColor(range(ncl), cmapn="copper")
                    )
                # axs[i][j].clabel(CS, inline=1, fontsize=fs)
            if j == 0:
                axs[i][j].set_ylabel(labels[i])
            if i == ndim - 1:
                axs[i][j].set_xlabel(labels[j])
            if (i == 0 and j == 0) or (j != 0):
                axs[i][j].set_yticklabels([])
            axs[i][j].set_aspect(abs(xmax - xmin) / abs(ymax - ymin))
            if i != ndim - 1:
                # axs[i][j].get_xaxis().set_visible(False)
                # axs[i][j].xaxis.set_visible(False)
                axs[i][j].set_xticklabels([])
            else:
                _plt.setp(axs[i][j].xaxis.get_majorticklabels(), rotation=50)

    _phc.savefig(fig)  # figname='outname')
    return


def corner(
    xs,
    bins=20,
    range=None,
    weights=None,
    color="k",
    smooth=None,
    smooth1d=None,
    labels=None,
    label_kwargs=None,
    show_titles=False,
    title_fmt=".2f",
    title_kwargs=None,
    truths=None,
    truth_color="#4682b4",
    scale_hist=False,
    quantiles=None,
    verbose=False,
    fig=None,
    max_n_ticks=5,
    top_ticks=False,
    hist_kwargs=None,
    qtdiff=None,
    **hist2d_kwargs
):
    """
    Make a *sick* corner plot showing the projections of a data set in a
    multi-dimensional space. kwargs are passed to hist2d() or used for
    `matplotlib` styling.

    Parameters
    ----------
    xs : array_like (nsamples, ndim)
        The samples. This should be a 1- or 2-dimensional array. For a 1-D
        array this results in a simple histogram. For a 2-D array, the zeroth
        axis is the list of samples and the next axis are the dimensions of
        the space.

    weights : array_like (nsamples,)
        The weight of each sample. If `None` (default), samples are given
        equal weight.

    labels : iterable (ndim,) (optional)
        A list of names for the dimensions. If a ``xs`` is a
        ``pandas.DataFrame``, labels will default to column names.

    show_titles : bool (optional)
        Displays a title above each 1-D histogram showing the 0.5 quantile
        with the upper and lower errors supplied by the quantiles argument.

    title_fmt : string (optional)
        The format string for the quantiles given in titles.
        (default: `.2f`)

    title_args : dict (optional)
        Any extra keyword arguments to send to the `add_title` command.

    extents : iterable (ndim,) (optional)
        A list where each element is either a length 2 tuple containing
        lower and upper bounds (extents) or a float in range (0., 1.)
        giving the fraction of samples to include in bounds, e.g.,
        [(0.,10.), (1.,5), 0.999, etc.].
        If a fraction, the bounds are chosen to be equal-tailed.

    truths : iterable (ndim,) (optional)
        A list of reference values to indicate on the plots.  Individual
        values can be omitted by using ``None``.

    truth_color : str (optional)
        A ``matplotlib`` style color for the ``truths`` makers.

    scale_hist : bool (optional)
        Should the 1-D histograms be scaled in such a way that the zero line
        is visible?

    quantiles : iterable (optional)
        A list of fractional quantiles to show on the 1-D histograms as
        vertical dashed lines.

    verbose : bool (optional)
        If true, print the values of the computed quantiles.

    plot_contours : bool (optional)
        Draw contours for dense regions of the plot.

    plot_datapoints : bool (optional)
        Draw the individual data points.

    max_n_ticks: int (optional)
        maximum number of ticks to try to use

    fig : matplotlib.Figure (optional)
        Overplot onto the provided figure object.

    """
    if quantiles is None:
        quantiles = []
    if title_kwargs is None:
        title_kwargs = dict()
    if label_kwargs is None:
        label_kwargs = dict()
    if qtdiff is None:
        qtdiff = _np.zeros(xs.shape[1])

    # Try filling in labels from pandas.DataFrame columns.
    if labels is None:
        try:
            labels = xs.columns
        except AttributeError:
            pass

    # Deal with 1D sample lists.
    xs = _np.atleast_1d(xs)
    if len(xs.shape) == 1:
        xs = _np.atleast_2d(xs)
    else:
        assert len(xs.shape) == 2, "The input sample array must be 1- or 2-D."
        xs = xs.T
    assert xs.shape[0] <= xs.shape[1], (
        "I don't believe that you want more " "dimensions than samples!"
    )

    # Parse the weight array.
    if weights is not None:
        weights = _np.asarray(weights)
        if weights.ndim != 1:
            raise ValueError("Weights must be 1-D")
        if xs.shape[1] != weights.shape[0]:
            raise ValueError("Lengths of weights must match number of samples")

    # Parse the parameter ranges.
    if range is None:
        if "extents" in hist2d_kwargs:
            _warn.warn("Deprecated keyword argument 'extents'. " "Use 'range' instead.")
            range = hist2d_kwargs.pop("extents")
        else:
            range = [[x.min(), x.max()] for x in xs]
            # Check for parameters that never change.
            m = _np.array([e[0] == e[1] for e in range], dtype=bool)
            if _np.any(m):
                raise ValueError(
                    (
                        "It looks like the parameter(s) in "
                        "column(s) {0} have no dynamic range. "
                        "Please provide a `range` argument."
                    ).format(", ".join(map("{0}".format, _np.arange(len(m))[m])))
                )

    else:
        # If any of the extents are percentiles, convert them to ranges.
        for i, _ in enumerate(range):
            try:
                emin, emax = range[i]
            except TypeError:
                q = [0.5 - 0.5 * range[i], 0.5 + 0.5 * range[i]]
                range[i] = quantile(xs[i], q, weights=weights)

    if len(range) != xs.shape[0]:
        raise ValueError("Dimension mismatch between samples and range")

    # Parse the bin specifications.
    try:
        bins = [float(bins) for _ in range]
    except TypeError:
        if len(bins) != len(range):
            raise ValueError("Dimension mismatch between bins and range")

    # Some magic numbers for pretty axis layout.
    K = len(xs)
    factor = 2.0  # size of one side of one panel
    lbdim = 0.5 * factor  # size of left/bottom margin
    trdim = 0.2 * factor  # size of top/right margin
    whspace = 0.05  # w/hspace size
    plotdim = factor * K + factor * (K - 1.0) * whspace
    dim = lbdim + plotdim + trdim

    # Create a new figure if one wasn't provided.
    if fig is None:
        fig, axes = _plt.subplots(K, K, figsize=(dim, dim), tight_layout=False)
    else:
        try:
            axes = _np.array(fig.axes).reshape((K, K))
        except:
            raise ValueError(
                "Provided figure has {0} axes, but data has "
                "dimensions K={1}".format(len(fig.axes), K)
            )
        fig.set_tight_layout(False)

    # Format the figure.
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(
        left=lb, bottom=lb, right=tr, top=tr, wspace=whspace, hspace=whspace
    )

    # Set up the default histogram keywords.
    if hist_kwargs is None:
        hist_kwargs = dict()
    hist_kwargs["color"] = hist_kwargs.get("color", color)
    if smooth1d is None:
        hist_kwargs["histtype"] = hist_kwargs.get("histtype", "step")

    for i, x in enumerate(xs):
        # Deal with masked arrays.
        if hasattr(x, "compressed"):
            x = x.compressed()

        if _np.shape(xs)[0] == 1:
            ax = axes
        else:
            ax = axes[i, i]
        # Plot the histograms.
        if smooth1d is None:
            n, _, _ = ax.hist(
                x,
                bins=int(round(bins[i])),
                weights=weights,
                range=range[i],
                **hist_kwargs
            )
        else:
            if gaussian_filter is None:
                raise ImportError("Please install scipy for smoothing")
            n, b = _np.histogram(
                x, bins=int(round(bins[i])), weights=weights, range=range[i]
            )
            n = _gf(n, smooth1d)
            x0 = _np.array(zip(b[:-1], b[1:])).flatten()
            y0 = _np.array(zip(n, n)).flatten()
            ax.plot(x0, y0, **hist_kwargs)

        if truths is not None and truths[i] is not None:
            ax.axvline(truths[i], color=truth_color)

        # Plot quantiles if wanted.
        if len(quantiles) > 0:
            qts = _np.array(quantiles) + qtdiff[i]
            qts = _np.where(qts <= 0, 0.01, qts)
            qts = _np.where(qts >= 1, 0.99, qts)
            qts = list(qts)
            qvalues = quantile(x, qts, weights=weights)
            for q in qvalues:
                ax.axvline(q, ls="dashed", color=color)

            if verbose:
                print("Quantiles:")
                print([item for item in zip(qts, qvalues)])

        if show_titles:
            # Compute the quantiles for the title. This might redo
            # unneeded computation but who cares.
            qts = _np.array(quantiles) + qtdiff[i]
            qts = _np.where(qts <= 0, 0.01, qts)
            qts = _np.where(qts >= 1, 0.99, qts)
            qts = list(qts)
            q_16, q_50, q_84 = quantile(x, qts, weights=weights)
            q_m, q_p = q_50 - q_16, q_84 - q_50

            # Format the quantile display.
            if q_50 < 1e6:
                fmt = "{{0:{0}}}".format(title_fmt).format
                title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
                title = title.format(fmt(q_50), fmt(q_m), fmt(q_p))
            else:
                expn = int(_np.log10(q_50))
                fmt = "{{0:{0}}}".format(title_fmt).format
                title = r"${{{0}}}_{{-{1}}}^{{+{2}}}e{3}$"
                title = title.format(
                    fmt(q_50 * 10**-expn),
                    fmt(q_m * 10**-expn),
                    fmt(q_p * 10**-expn),
                    expn,
                )

            # Add in the column name if it's given.
            if labels is not None:
                if len(labels) <= i:
                    labels = list(labels) + [""]
                j = labels[i].find(" ")
                if j == -1:
                    j = len(labels[i])
                title = "{0} = {1}".format(labels[i][:j], title)

            # Add the title to the axis.
            ax.set_title(title, **title_kwargs)

        # Set up the axes.
        ax.set_xlim(range[i])
        if scale_hist:
            maxn = _np.max(n)
            ax.set_ylim(-0.1 * maxn, 1.1 * maxn)
        else:
            ax.set_ylim(0, 1.1 * _np.max(n))
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(_MaxNLoc(max_n_ticks, prune="lower"))

        if i < K - 1:
            if top_ticks:
                ax.xaxis.set_ticks_position("top")
                [l.set_rotation(45) for l in ax.get_xticklabels()]
            else:
                ax.set_xticklabels([])
        else:
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if labels is not None:
                ax.set_xlabel(labels[i], **label_kwargs)
                ax.xaxis.set_label_coords(0.5, -0.3)

        for j, y in enumerate(xs):
            if _np.shape(xs)[0] == 1:
                ax = axes
            else:
                ax = axes[i, j]
            if j > i:
                ax.set_visible(False)
                ax.set_frame_on(False)
                continue
            elif j == i:
                continue

            # Deal with masked arrays.
            if hasattr(y, "compressed"):
                y = y.compressed()
            hist2d(
                y,
                x,
                ax=ax,
                range=[range[j], range[i]],
                weights=weights,
                color=color,
                smooth=smooth,
                bins=[bins[j], bins[i]],
                **hist2d_kwargs
            )

            if truths is not None:
                if truths[i] is not None and truths[j] is not None:
                    ax.plot(truths[j], truths[i], "s", color=truth_color)
                if truths[j] is not None:
                    ax.axvline(truths[j], color=truth_color)
                if truths[i] is not None:
                    ax.axhline(truths[i], color=truth_color)

            ax.xaxis.set_major_locator(_MaxNLoc(max_n_ticks, prune="lower"))
            ax.yaxis.set_major_locator(_MaxNLoc(max_n_ticks, prune="lower"))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if labels is not None:
                    ax.set_xlabel(labels[j], **label_kwargs)
                    ax.xaxis.set_label_coords(0.5, -0.3)

            if j > 0:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if labels is not None:
                    ax.set_ylabel(labels[i], **label_kwargs)
                    ax.yaxis.set_label_coords(-0.3, 0.5)

    return fig


def quantile(x, q, weights=None):
    """
    Like numpy.percentile, but:

    * Values of q are quantiles [0., 1.] rather than percentiles [0., 100.]
    * scalar q not supported (q must be iterable)
    * optional weights on x

    """
    if weights is None:
        return _np.percentile(x, [100.0 * qi for qi in q])
    else:
        idx = _np.argsort(x)
        xsorted = x[idx]
        cdf = _np.add.accumulate(weights[idx])
        cdf /= cdf[-1]
        return _np.interp(q, cdf, xsorted).tolist()


def hist2d(
    x,
    y,
    bins=20,
    range=None,
    weights=None,
    levels=None,
    smooth=None,
    ax=None,
    color=None,
    plot_datapoints=True,
    plot_density=True,
    plot_contours=True,
    fill_contours=False,
    contour_kwargs=None,
    contourf_kwargs=None,
    data_kwargs=None,
    **kwargs
):
    """
    Plot a 2-D histogram of samples.

    """
    if ax is None:
        ax = _plt.gca()

    # Set the default range based on the data range if not provided.
    if range is None:
        if "extent" in kwargs:
            _warn.warn("Deprecated keyword argument 'extent'. " "Use 'range' instead.")
            range = kwargs["extent"]
        else:
            range = [[x.min(), x.max()], [y.min(), y.max()]]

    # Set up the default plotting arguments.
    # if color is None:
    # color = "k"
    #
    # Choose the default "sigma" contour levels.
    # if levels is None:
    # levels = 1.0 - _np.exp(-0.5 * _np.arange(0.5, 2.1, 0.5) ** 2)
    #
    # This is the color map for the density plot, over-plotted to indicate the
    # density of the points near the center.
    # density_cmap = LinearSegmentedColormap.from_list(
    # "density_cmap", [color, (1, 1, 1, 0)])
    #
    # This color map is used to hide the points at the high density areas.
    # white_cmap = LinearSegmentedColormap.from_list(
    # "white_cmap", [(1, 1, 1), (1, 1, 1)], N=2)
    #
    # This "color map" is the list of colors for the contour levels if the
    # contours are filled.
    # rgba_color = colorConverter.to_rgba(color)
    # contour_cmap = [rgba_color] + [list(rgba_color) for l in levels]
    # for i, l in enumerate(levels):
    # contour_cmap[i+1][-1] *= float(len(levels) - i) / (len(levels)+1)
    #
    # We'll make the 2D histogram to directly estimate the density.
    try:
        H, X, Y = _np.histogram2d(
            x.flatten(), y.flatten(), bins=bins, range=range, weights=weights
        )
    except ValueError:
        raise ValueError(
            "It looks like at least one of your sample columns "
            "have no dynamic range. You could try using the "
            "'range' argument."
        )
    #
    # if smooth is not None:
    # if gaussian_filter is None:
    # raise ImportError("Please install scipy for smoothing")
    # H = _gf(H, smooth)
    #
    # Compute the density levels.
    # Hflat = H.flatten()
    # inds = _np.argsort(Hflat)[::-1]
    # Hflat = Hflat[inds]
    # sm = _np.cumsum(Hflat)
    # sm /= sm[-1]
    # V = _np.empty(len(levels))
    # for i, v0 in enumerate(levels):
    # try:
    # V[i] = Hflat[sm <= v0][-1]
    # except:
    # V[i] = Hflat[0]
    #
    # Compute the bin centers.
    # X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])
    #
    # Extend the array for the sake of the contours at the plot edges.
    # H2 = H.min() + _np.zeros((H.shape[0] + 4, H.shape[1] + 4))
    # H2[2:-2, 2:-2] = H
    # H2[2:-2, 1] = H[:, 0]
    # H2[2:-2, -2] = H[:, -1]
    # H2[1, 2:-2] = H[0]
    # H2[-2, 2:-2] = H[-1]
    # H2[1, 1] = H[0, 0]
    # H2[1, -2] = H[0, -1]
    # H2[-2, 1] = H[-1, 0]
    # H2[-2, -2] = H[-1, -1]
    # X2 = _np.concatenate([
    # X1[0] + _np.array([-2, -1]) * _np.diff(X1[:2]),
    # X1,
    # X1[-1] + _np.array([1, 2]) * _np.diff(X1[-2:]),
    # ])
    # Y2 = _np.concatenate([
    # Y1[0] + _np.array([-2, -1]) * _np.diff(Y1[:2]),
    # Y1,
    # Y1[-1] + _np.array([1, 2]) * _np.diff(Y1[-2:]),
    # ])
    #
    # if plot_datapoints:
    # if data_kwargs is None:
    # data_kwargs = dict()
    # data_kwargs["color"] = data_kwargs.get("color", color)
    # data_kwargs["ms"] = data_kwargs.get("ms", 2.0)
    # data_kwargs["mec"] = data_kwargs.get("mec", "none")
    # data_kwargs["alpha"] = data_kwargs.get("alpha", 0.1)
    # ax.plot(x, y, "o", zorder=-1, rasterized=True, **data_kwargs)
    #
    # Plot the base fill to hide the densest data points.
    # if plot_contours or plot_density:
    # ax.contourf(X2, Y2, H2.T, [V[-1], H.max()],
    # cmap=white_cmap, antialiased=False)
    #
    # if plot_contours and fill_contours:
    # if contourf_kwargs is None:
    # contourf_kwargs = dict()
    # contourf_kwargs["colors"] = contourf_kwargs.get("colors", contour_
    # cmap)
    # contourf_kwargs["antialiased"] = contourf_kwargs.get("antialiased",
    # False)
    # ax.contourf(X2, Y2, H2.T, _np.concatenate([[H.max()], V, [0]]),
    # **contourf_kwargs)
    #
    # Plot the density map. This can't be plotted at the same time as the
    # contour fills.
    # elif plot_density:
    # ax.pcolor(X, Y, H.max() - H.T, cmap=density_cmap)
    if plot_density:
        ax.pcolor(X, Y, H.max() - H.T, cmap=_cm.get_cmap("gist_heat"))
    #
    # Plot the contour edge colors.
    # if plot_contours:
    # if contour_kwargs is None:
    # contour_kwargs = dict()
    # contour_kwargs["colors"] = contour_kwargs.get("colors", color)
    # ax.contour(X2, Y2, H2.T, V, **contour_kwargs)

    ax.set_xlim(range[0])
    ax.set_ylim(range[1])
    return


# Main
if __name__ == "__main__":
    pass
