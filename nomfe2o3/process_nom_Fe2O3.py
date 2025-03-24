"""
Process NOMAD data files, slice them, and produce reduced data
"""

import os
import json
import argparse
import time

import mantid
import mantid.simpleapi as api

mantid.kernel.config.setLogLevel(3)


def slice_data(ws, time_interval: int) -> list:
    """
    Split the provided workspace into slices of the given time interval (seconds).
    """
    splitws, infows = api.GenerateEventsFilter(
        InputWorkspace=ws, TimeInterval=time_interval
    )

    api.FilterEvents(
        InputWorkspace=ws,
        SplitterWorkspace=splitws,
        InformationWorkspace=infows,
        OutputWorkspaceBaseName="time_ws",
        GroupWorkspaces=True,
        FilterByPulseTime=True,
        OutputWorkspaceIndexedFrom1=True,
        CorrectionToSample="None",
        SpectrumWithoutDetector="Skip",
        SplitSampleLogs=False,
        OutputTOFCorrectionWorkspace="mock",
    )
    wsgroup = api.mtd["time_ws"]
    return wsgroup


def process_slices(
    ws_list: list,
    prepend: str,
    output_dir: str,
    template: str,
    variable: str = "BL1B:SE:SampleTemp",
) -> None:
    """
    Process the list of workspaces, save them to disk, and run the MTS reduction.

    ws_list: list of workspaces
    prepend: string to prepend to the file name
    output_dir: directory to save the files
    template: MTS configuration template
    variable: the varying PV for the measurements (for example temperature)

    """
    # Create the output directory for intermediate files if it does not exist
    h5_dir = os.path.join(output_dir, "h5")
    if not os.path.exists(h5_dir):
        os.makedirs(h5_dir)

    # Create the output directory for reduced files if it does not exist
    reduced_dir = os.path.join(output_dir, "reduced")
    if not os.path.exists(reduced_dir):
        os.makedirs(reduced_dir)

    for ws in ws_list:
        # Get the mean value of the variable we are interested in
        t = ws.getRun()[variable]
        value = t.getStatistics().mean

        # Save the workspace to a Nexus file so we can run the MTS reduction
        file_name = "%s_%3.0f.h5" % (prepend, value)
        file_path = os.path.join(h5_dir, file_name)
        api.SaveNexus(ws, Filename=file_path)

        # Update the MTS configuration template
        template["Title"] = file_name
        template["Sample"]["Filenames"] = [file_path]
        template["OutputDir"] = reduced_dir

        with open("mts_config.json", "w") as fd:
            json.dump(template, fd)

        # Run the Mantid Total Scattering reduction
        try:
            os.system("mts mts_config.json")
        except Exception as e:
            print(f"Error processing {file_name}: {e}")
            continue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process NOMAD data files.")
    parser.add_argument(
        "--run_number", type=int, required=True, help="Run number to process"
    )
    parser.add_argument(
        "--time_interval",
        type=int,
        default=30,
        help="Time interval for slicing the data",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=os.path.expanduser("/SNS/NOM/IPTS-34537/shared/sliced"),
        help="Directory to save the output files",
    )
    parser.add_argument(
        "--template",
        type=str,
        default="template.json",
        help="Path to the MTS configuration template",
    )

    args = parser.parse_args()

    ws_sc = api.Load("NOM_%s" % args.run_number)

    duration = ws_sc.getRun()["duration"].value
    print(f"Run {args.run_number}: Duration = {duration} seconds")

    with open(args.template) as fd:
        template = json.load(fd)

    t0 = time.time()
    ws_names = slice_data(ws_sc, args.time_interval)
    number_of_slices = len(ws_names)
    t1 = time.time()
    process_slices(ws_names, "nom_%s" % args.run_number, args.output_dir, template)
    t2 = time.time()

    print(f"\n\n\nTime to slice data: {t1 - t0:.2f} s")
    print(f"Time to process {number_of_slices} slices: {t2 - t1:.2f} s")
    print(f"Total time: {t2 - t0:.2f} s")
