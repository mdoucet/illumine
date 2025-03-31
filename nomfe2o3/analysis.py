import os
import json
import numpy as np

import matplotlib.pyplot as plt

from scipy.signal import find_peaks, find_peaks_cwt
from scipy.optimize import curve_fit


TEMPLATE = "template.json"
OUTPUT_DIR = "/SNS/NOM/IPTS-34537/shared/sliced"


def reduce_all():
    import process_nom_Fe2O3 as pnf
    import mantid
    import mantid.simpleapi as api

    for run_number in range(208518, 208532):
        ws_sc = api.Load("NOM_%s" % run_number)

        with open(TEMPLATE) as fd:
            template = json.load(fd)

        ws_names = pnf.slice_data(ws_sc, 30)

        pnf.process_slices(ws_names, "nom_%s" % run_number, OUTPUT_DIR, template)


def load_data(data_dir: str, bank: int = 2) -> list:
    ramp = []
    temps = []

    for j, data_file in enumerate(sorted(os.listdir(data_dir))):
        if not data_file.startswith("nom"):
            continue
        with open(f"{data_dir}/{data_file}", "r") as f:
            toks = data_file.split("_")

            tof = []
            counts = []
            errors = []
            data = f.readlines()
            started = False
            for i, line in enumerate(data):
                if (
                    line.startswith("#")
                    or line.startswith("BANK")
                    or line.startswith("Monitor")
                    or line.startswith("Sample")
                ):
                    if line.startswith("BANK %d" % bank):
                        started = True
                    elif line.startswith("BANK"):
                        started = False
                    continue
                else:
                    _data = line.split()
                    if started and len(_data) == 3:
                        tof.append(float(_data[0]))
                        counts.append(float(_data[1]))
                        errors.append(float(_data[2]))

            # temp = temps[j]
            # np.savetxt(f'{output_dir}/Fe2O3_BANK{bank}_{temp}.dat', np.asarray([tof, counts, errors]).T)
            ramp.append(np.asarray([tof, counts, errors]))
            temp = float(toks[2].replace(".h5.gsa", ""))
            temps.append(temp)

    return temps, ramp


# Define a Gaussian function
def gaussian(x, amp, cen, wid, bck):
    return amp * np.exp(-((x - cen) ** 2) / (2 * wid**2)) + bck


def peak_finder(tof, counts, peak_prominence: float = 2):
    # Identify peaks in the counts data
    peaks, _ = find_peaks(counts, prominence=2)

    plt.figure(figsize=(10, 6))
    plt.plot(tof, counts, linestyle="--")

    integral = []
    center = []

    for peak in tof[peaks]:
        # Fit a Gaussian to each peak
        delta = tof / 25
        mask = (tof > peak - delta) & (tof < peak + delta)
        try:
            popt, pcov = curve_fit(
                gaussian, tof[mask], counts[mask], p0=[counts[mask].max(), peak, 10, 0]
            )
            perr = np.sqrt(np.diag(pcov))

            a = popt[0] * popt[2] * np.sqrt(2 * np.pi)
            integral.append(a)
            center.append(popt[1])
        except:
            print(f"Failed to fit peak at {peak}")
            continue
        plt.plot(
            tof[mask],
            gaussian(tof[mask], *popt),
            color="green",
            linestyle="-",
            linewidth=3,
        )
        plt.axvline(x=peak, color="r", linestyle="--", linewidth=1)

    plt.xlabel("TOF")
    plt.ylabel("Counts")
    plt.show()
    return peaks


def series_analyzer(data: list, peaks: list, temp_list: list):
    integral = []
    err_int = []
    center = []
    temperature = []
    chi2 = []

    print("Number of peaks", len(peaks))
    for i, temp_data in enumerate(data):
        tof = temp_data[0]
        counts = temp_data[1]
        err_counts = temp_data[2]
        integral_for_temp = []
        err_for_temp = []
        skip_point = False
        
        #print(np.mean(err_counts/counts))
        if np.mean(err_counts/counts) > .8:
            continue

        for peak in tof[peaks]:
            # Fit a Gaussian to each peak
            delta = tof / 20 - 5
            mask = (tof > peak - delta) & (tof < peak + delta)
            
            try:
                popt, pcov = curve_fit(
                    gaussian,
                    tof[mask],
                    counts[mask],
                    p0=[counts[mask].max(), peak, 10, counts[mask].min()],
                )
                perr = np.sqrt(np.diag(pcov))
                a = np.fabs(popt[0] * popt[2] * np.sqrt(2 * np.pi))
                err_a = np.sqrt(2 * np.pi) * np.sqrt(
                    popt[0] ** 2 * perr[2] ** 2 + popt[2] ** 2 * perr[0] ** 2
                )

                _chi2 = np.mean(np.asarray(gaussian(tof[mask], *popt) - counts[mask])**2)
            except:
                print(f"Failed to fit peak at {peak}")
                a = 0
                err_a = 0
                _chi2 = 0
                continue
    
            # If the error is too large, we will skip this point
            #if err_a > a / 10.0:
            #    skip_point = True

            integral_for_temp.append(a)
            err_for_temp.append(err_a)
            # center.append(popt[1])

        if not skip_point:
            temperature.append(temp_list[i])
            integral.append(integral_for_temp)
            err_int.append(err_for_temp)
            chi2.append(_chi2)

    integral = np.asarray(integral).T
    err_int = np.asarray(err_int).T
    temperature = np.asarray(temperature)

    n_tot = len(peaks)
    ysize = 2 * n_tot
    fig, axs = plt.subplots(n_tot, 1, dpi=100, figsize=(7, ysize), sharex=True)

    for i in range(len(peaks)):
        plt.subplot(n_tot, 1, i + 1)
        # plt.plot(temperature, integral[i], label='%g' % tof[peaks[i]])

        mask = err_int[i] < integral[i]/2.0
        if len(temperature[mask]) > 0:
            plt.errorbar(
                temperature[mask], integral[i][mask], yerr=err_int[i][mask], label="%g" % tof[peaks[i]]
            )
            plt.legend(frameon=False)
    plt.xlabel("Temperature")
    plt.ylabel("Peak integral")
    # plt.yscale('log')
    plt.legend(frameon=False)
    plt.show()

    return temperature, integral, err_int
