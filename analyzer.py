import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit


# Define a Gaussian function
def gaussian(x, amp, cen, wid, bck):
    return amp * np.exp(-((x - cen) ** 2) / (2 * wid**2)) + bck


def load_data(data_dir, startswith: str = "Fe", extension: str = ".dat"):
    file_list = os.listdir(data_dir)
    file_list = [f for f in file_list if f.startswith(startswith)]
    file_list.sort(key=lambda x: float(x.split("_")[2].replace(extension, "")))

    temp_list = []
    chi2_list = []

    for i, file in enumerate(file_list):
        if file.startswith(startswith):
            temp = float(file.split("_")[2].replace(extension, ""))
            tof, counts, err = np.loadtxt(os.path.join(data_dir, file), comments="'").T
            # Identify peaks in the counts data
            peaks, _ = find_peaks(counts, height=0)
            if i == 0:
                tof0 = tof
                counts0 = counts
                err0 = err
                chi2 = 0
            else:
                non_nan_indices = ~np.isnan(counts)
                non_nan_indices &= ~np.isnan(counts0)
                non_nan_indices &= ~np.isnan(err)
                non_nan_indices &= ~np.isnan(err0)
                non_nan_indices &= err > 0
                non_nan_indices &= err0 > 0
                chi2 = np.mean(
                    (counts[non_nan_indices] - counts0[non_nan_indices]) ** 2
                    / err[non_nan_indices] ** 2
                )
                temp_list.append(temp)
            chi2_list.append(chi2)

    return file_list, temp_list, chi2_list


def peak_finder(
    file_path: str, label: str = "", peak_prominence: float = 5, peak_width: float = 10
):
    # TODO: Fit peak in Q space to avoid weird widths

    _data = np.loadtxt(file_path, comments="'")

    tof_range = [0, len(_data) - 100]
    tof, counts, err = _data[tof_range[0] : tof_range[1]].T

    # Identify peaks in the counts data
    peaks, _ = find_peaks(counts, prominence=peak_prominence)

    print(peaks)
    plt.figure(figsize=(10, 6))
    plt.plot(tof, counts, label=label, linestyle="--")

    integral = []
    center = []
    width = []

    for peak in tof[peaks]:
        # Fit a Gaussian to each peak
        delta = 0.04 * tof
        mask = (tof > peak - delta) & (tof < peak + delta)
        try:
            popt, pcov = curve_fit(
                gaussian,
                tof[mask],
                counts[mask],
                p0=[counts[mask].max(), peak, peak_width, 0],
            )
            perr = np.sqrt(np.diag(pcov))

            a = popt[0] * popt[2] * np.sqrt(2 * np.pi)
            integral.append(a)
            center.append(popt[1])
            width.append(popt[2])
        except Exception as e:
            print(f"Failed to fit peak at {peak}")
            continue
        plt.plot(tof[mask], gaussian(tof[mask], *popt), color="green", linewidth=2)
        plt.axvline(x=peak, color="r", linestyle="--", linewidth=0.5)

    plt.xlabel("TOF")
    plt.ylabel("Counts")
    plt.legend()
    plt.show()
    return peaks, center


def analyzer_series(file_list, peaks, peak_width: float = 10, extension: str = ".dat"):
    integral = []
    err_int = []
    center = []
    temperature = []

    print("Number of peaks", len(peaks))
    for i, file in enumerate(file_list):
        temp = float(file.split("_")[2].replace(extension, ""))
        tof, counts, err = np.loadtxt(file, comments="'").T
        integral_for_temp = []
        err_for_temp = []

        for peak in tof[peaks]:
            # Fit a Gaussian to each peak
            delta = 0.04 * tof
            mask = (tof > peak - delta) & (tof < peak + delta)
            try:
                popt, pcov = curve_fit(
                    gaussian,
                    tof[mask],
                    counts[mask],
                    p0=[counts[mask].max(), peak, peak_width, 0],
                )
                perr = np.sqrt(np.diag(pcov))
                a = np.fabs(popt[0] * popt[2] * np.sqrt(2 * np.pi))
                err_a = np.sqrt(2 * np.pi) * np.sqrt(
                    popt[0] ** 2 * perr[2] ** 2 + popt[2] ** 2 * perr[0] ** 2
                )
                integral_for_temp.append(a)
                err_for_temp.append(err_a)
                # center.append(popt[1])
            except:
                print(f"Failed to fit peak at {peak}")
                integral_for_temp.append(0)
                err_for_temp.append(0)
                continue
        temperature.append(temp)
        integral.append(integral_for_temp)
        err_int.append(err_for_temp)

    integral = np.asarray(integral).T
    err_int = np.asarray(err_int).T

    n_tot = len(peaks)
    ysize = 2 * n_tot
    fig, axs = plt.subplots(n_tot, 1, dpi=100, figsize=(7, ysize), sharex=True)

    for i in range(len(peaks)):
        plt.subplot(n_tot, 1, i + 1)
        # plt.plot(temperature, integral[i], label='%g' % tof[peaks[i]])
        plt.errorbar(
            temperature, integral[i], yerr=err_int[i], label="%g" % tof[peaks[i]]
        )
        plt.legend()
    plt.xlabel("Temperature")
    plt.ylabel("Peak integral")
    # plt.yscale('log')
    plt.legend(frameon=False)
    plt.show()
    return temperature, integral, err_int
