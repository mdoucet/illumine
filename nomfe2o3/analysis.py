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


def peak_finder(tof, counts, peak_prominence: float = 2):
    # Identify peaks in the counts data
    peaks, _ = find_peaks(counts, prominence=2)

    plt.figure(figsize=(10, 6))
    plt.plot(tof, counts, linestyle="--")

    # Define a Gaussian function
    def gaussian(x, amp, cen, wid, bck):
        return amp * np.exp(-((x - cen) ** 2) / (2 * wid**2)) + bck

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
