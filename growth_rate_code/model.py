import numpy as np
import matplotlib.pyplot as plt
from lmfit import Minimizer, Parameters
import pandas as pd


def read_csv(csv_path):
    return pd.read_csv(csv_path, header=None, float_precision='high').to_numpy(dtype='float64').flatten()


def piecewise_linear(x_values, TERc, TERl, ORIc, ORIl):
    # parameter to fit are TER and ORI coverage (y) and TER and ORI location (x)
    a = (TERc - ORIc) / (TERl - ORIl)
    x1 = min(TERl, ORIl)
    x2 = max(TERl, ORIl)
    y1 = TERc if x1 == TERl else ORIc
    y2 = TERc if x2 == TERl else ORIc

    return np.piecewise(x_values, [x_values <= x1, (x_values > x1) & (x_values < x2), x_values >= x2],
                        [lambda x:-a*x + y1 + a*x1, lambda x:a*x + y1 - a*x1, lambda x:-a*x + y2 + a*x2])


def cost_func_lmfit_circlepen(params, x, y, ori_ter_dist_dev):
    Tc = params['TERc']
    Tl = params['TERl']
    Oc = params['ORIc']
    Ol = params['ORIl']
    model = piecewise_linear(x, Tc, Tl, Oc, Ol)
    circle_consistency = np.abs(model[0] - model[-1]) * 1.0
    model = model[1:-1]
    # enforce ori and ter distance:
    penalization = max(np.abs(np.abs(Tl - Ol) - 0.5) - ori_ter_dist_dev, 0.0) * 100.0
    return np.abs(model - y) + penalization + circle_consistency


def cost_func_lmfit(params, x, y, ori_ter_dist_dev):
    Tc = params['TERc']
    Tl = params['TERl']
    Oc = params['ORIc']
    Ol = params['ORIl']
    model = piecewise_linear(x, Tc, Tl, Oc, Ol)
    penalization = max(np.abs(np.abs(Tl - Ol) - 0.5) - ori_ter_dist_dev, 0.0) * 100.0
    return np.abs(model - y) + penalization


def compute_ptr(csv_path, n_windows=50, save_plot=True, save_name='test', high_res_plot=False,
                subsample_size=10000, use_all_reads=True):
    """
    computes the peak to trough ratio (growth rate) for the sample at the given csv_path
    :param csv_path: path to sample csv file
    :param n_windows: initial number of bins used
    :param save_plot: whether to save the plots
    :param save_name: plots saved under this name
    :param high_res_plot: save plots in high resolution
    :param subsample_size: maximum number of reads to use
    :param use_all_reads: ignore subsample size and use all reads
    :return: growth rate, number of all reads, used bins, error codes, fit error
    """
    csv_data = read_csv(csv_path)
    number_all_reads = len(csv_data)
    if number_all_reads < subsample_size or use_all_reads:
        subsample_size = number_all_reads

    if high_res_plot:
        img_dpi = 70
    else:
        img_dpi = 40

    window_size = 1 / n_windows

    # normally, ori and ter would be 0.5 apart; because of incorrect alignments and artefacts,
    # ori and ter can be apart [0.5-ori_ter_dist_dev, 0.5+ori_ter_dist_dev]
    ori_ter_dist_dev = 1/6  # => [1/3, 2/3]

    growth_rate_th = 20  # maximal valid growth rate
    fit_err_th = 9.0  # maximal valid fit error
    median_coverage_th_factor = 0.1
    read_coverage_th = 500
    bins_after_filter_th = int(0.5 * n_windows)

    # median filter parameters:
    filter_size = 5
    filter_th = 1.5

    error_codes = []
    """
    [] : no error
    -1 : insufficient median coverage
    -2 : not enough reads left after filtering
    -3 : not enough bins left after filtering
    -4 : growth rate way too high, probably corrupt fitting
    -5 : fit error too large
    """

    # take random subset
    data = np.random.choice(csv_data, subsample_size, replace=False)
    # compute bin coverage
    coverage, borders = np.histogram(data, bins=n_windows, range=(0.0, 1.0))
    borders_track = borders.copy()  # keep track of removed bins in filtering

    # x values for plot are midpoints of bin intervals
    x_values = [(b + borders[i+1]) / 2 for i, b in enumerate(borders) if i < len(borders) - 1]

    # classify as invalid if median coverage is insufficient
    expected_coverage = subsample_size / n_windows
    median_coverage = np.median(coverage)
    if median_coverage < expected_coverage * median_coverage_th_factor:
        error_codes.append(-1)

    # median filter:
    # bins in close neighborhood should have approx the same values
    # => remove outliers (compare against median of values in filter range)
    delete_indices = []
    deleted_bins = 0
    length = len(coverage)
    for i in range(length):
        indices = np.arange(i, i+filter_size) % length
        median_value = np.median(coverage[indices])
        for k in range(filter_size):
            idx = (i+k) % length
            diff = median_value*filter_th - median_value
            if coverage[idx] > median_value+diff or coverage[idx] < median_value-diff or \
               coverage[idx] > median_value+expected_coverage/2 or \
               coverage[idx] < median_value-expected_coverage/2:
                delete_indices.append(idx)

    delete_indices = np.unique(delete_indices)
    deleted_bins += len(delete_indices)
    all_deleted_indices = []
    if len(delete_indices) > 0:
        coverage = np.delete(coverage, delete_indices)
        x_values = np.delete(x_values, delete_indices)
        all_deleted_indices.append(borders_track[delete_indices])

    if deleted_bins > 0:
        all_deleted_indices = np.concatenate(all_deleted_indices, axis=0)
        all_deleted_indices.sort()

    if np.sum(coverage) < read_coverage_th:
        # reject if not enough reads left
        error_codes.append(-2)
    
    if len(x_values) < bins_after_filter_th:
        # reject if not enough bins left
        error_codes.append(-3)

    x = x_values
    y = coverage

    # cannot compute growth rate if there are too few bins for fitting algo
    if len(x) < 4 and save_plot:
        # save plots without fitted line segments:
        plt.figure(figsize=(31, 7), dpi=img_dpi)
        plt.rc('font', size=15)

        plt.subplot(1, 3, 1)
        title = "REJECT! Error codes: %s. %d reads (%d bins)" % (error_codes, subsample_size, n_windows)
        plt.title(title)
        plt.hist(data, bins=n_windows, density=False, range=(0.0, 1.0))

        plt.subplot(1, 3, 2)
        title = "%d reads (%d bins filtered)" % (subsample_size, len(all_deleted_indices))
        plt.title(title)

        if len(all_deleted_indices) > 0:
            # remove reads within filtered bins:
            delete_borders = np.zeros((len(all_deleted_indices), 2))
            delete_borders[:, 0] = all_deleted_indices
            delete_borders[:, 1] = delete_borders[:, 0] + window_size
            masks = np.ones_like(data, dtype=bool)
            for i, r in enumerate(data):
                for k in range(delete_borders.shape[0]):
                    if delete_borders[k, 0] <= r < delete_borders[k, 1]:
                        masks[i] = False
                        break
            data = data[masks]

        plt.hist(data, bins=n_windows, density=False, range=(0.0, 1.0))

        all_data_bins = 100
        plt.subplot(1, 3, 3)
        title = "all %d reads (%d bins)" % (number_all_reads, all_data_bins)
        plt.title(title)
        plt.hist(csv_data, bins=all_data_bins, density=False, range=(0.0, 1.0))
        plt.savefig(save_name, format='png')
        plt.close()

    if len(x) < 4:  # cannot fit function with 4 parameters with less than 4 data points
        return 0, number_all_reads, n_windows - deleted_bins, error_codes, -1

    # initialize parameters
    init_TERc = int((subsample_size / n_windows) * 0.8)
    init_ORIc = int((subsample_size / n_windows) * 1.2)
    init_ORIl = 0.3
    init_TERl = 0.8

    pfit = Parameters()
    pfit.add(name='TERc', value=init_TERc, min=0.1, max=subsample_size)
    pfit.add(name='ORIc', value=init_ORIc, min=0.1, max=subsample_size)
    pfit.add(name='TERl', value=init_TERl, min=0.0, max=1.0)
    pfit.add(name='ORIl', value=init_ORIl, min=0.0, max=1.0)

    myfit = Minimizer(cost_func_lmfit, pfit, fcn_args=(x, y, ori_ter_dist_dev))
    out = myfit.minimize(method='leastsq')
    p = [np.float32(out.params['TERc'].value), np.float32(out.params['TERl'].value),
         np.float32(out.params['ORIc'].value), np.float32(out.params['ORIl'].value)]

    # enforce circularity
    x_circle = np.empty((len(x)+2,))
    # need points at 0.0 and 1.0 for circularity
    x_circle[1:-1] = x
    x_circle[0] = 0.0
    x_circle[-1] = 1.0
    pfit = Parameters()
    # initialize parameters with values of first fit
    pfit.add(name='TERc', value=p[0], min=0.1, max=subsample_size)
    pfit.add(name='ORIc', value=p[2], min=0.1, max=subsample_size)
    pfit.add(name='TERl', value=p[1], min=0.0, max=1.0)
    pfit.add(name='ORIl', value=p[3], min=0.0, max=1.0)

    myfit = Minimizer(cost_func_lmfit_circlepen, pfit, fcn_args=(x_circle, y, ori_ter_dist_dev))
    out = myfit.minimize(method='leastsq')
    p = [np.float32(out.params['TERc'].value), np.float32(out.params['TERl'].value),
         np.float32(out.params['ORIc'].value), np.float32(out.params['ORIl'].value)]

    growth_rate = p[2] / p[0] if p[2] > p[0] else p[0] / p[2]
    if growth_rate > growth_rate_th:
        error_codes.append(-4)

    # compute fit error
    preds = piecewise_linear(x, *p)
    if np.mean(y) == 0.0:
        bins_fit_error = -1  # faulty sample
    else:
        # scale by mean coverage to make it comparable to samples with different number of reads
        bins_fit_error = (np.mean((y - preds)**2) / np.mean(y)**2) * 100.0
    
    if bins_fit_error > fit_err_th or bins_fit_error < 0:
        error_codes.append(-5)

    if save_plot:
        # plot 4 subplots next to each other
        plt.figure(figsize=(31, 7), dpi=img_dpi)
        plt.rc('font', size=15)

        # 1) plot coverage of bins used for fitting + fit result
        plt.subplot(1, 4, 1)
        xd = np.linspace(0, 1.0, n_windows)
        plt.xlabel('Position', fontsize=20)
        plt.ylabel('Coverage (# reads)', fontsize=20)
        plt.plot(x, y, "o")
        plt.plot(xd, piecewise_linear(xd, *p))
        title_text = 'pred rate: %.1f,  fit err: %.1f' % (growth_rate, bins_fit_error)
        plt.title(title_text)

        # 2) plot fit result in all data histogram
        plt.subplot(1, 4, 2)
        if len(error_codes) > 0:
            title = "REJECT! Error codes: %s. %d reads (%d bins)" % (error_codes, subsample_size, n_windows)
        else:
            title = "%d reads (%d bins)" % (subsample_size, n_windows)
        plt.title(title)
        plt.hist(data, bins=n_windows, density=False, range=(0.0, 1.0))
        plt.plot(xd, piecewise_linear(xd, *p), color='red', lw=5)

        # 3) plot fit result in filtered data histogram
        plt.subplot(1, 4, 3)
        title = "%d reads (%d bins filtered)" % (subsample_size, len(all_deleted_indices))
        plt.title(title)

        if len(all_deleted_indices) > 0:
            # remove reads within filtered bins:
            delete_borders = np.zeros((len(all_deleted_indices), 2))
            delete_borders[:, 0] = all_deleted_indices
            delete_borders[:, 1] = delete_borders[:, 0] + window_size
            masks = np.ones_like(data, dtype=bool)
            for i, r in enumerate(data):
                for k in range(delete_borders.shape[0]):
                    if delete_borders[k, 0] <= r < delete_borders[k, 1]:
                        masks[i] = False
                        break
            data = data[masks]

        plt.hist(data, bins=n_windows, density=False, range=(0.0, 1.0))
        plt.plot(xd, piecewise_linear(xd, *p), color='red', lw=5)

        # 4) plot fit result in all data histogram with more number of bins
        all_data_bins = 100
        conversion_fac = (n_windows / all_data_bins) * (number_all_reads/subsample_size)
        plt.subplot(1, 4, 4)
        title = "all %d reads (%d bins)" % (number_all_reads, all_data_bins)
        plt.title(title)
        plt.hist(csv_data, bins=all_data_bins, density=False, range=(0.0, 1.0))
        plt.plot(xd, piecewise_linear(xd, *p) * conversion_fac, color='red', lw=5)
        plt.savefig(save_name, format='png')
        plt.close()

    return growth_rate, number_all_reads, n_windows - deleted_bins, error_codes, bins_fit_error
