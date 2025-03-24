import json
import process_nom_Fe2O3 as pnf

import mantid
import mantid.simpleapi as api

TEMPLATE = "template.json"
OUTPUT_DIR = "/SNS/NOM/IPTS-34537/shared/sliced"


def reduce_all():
    for run_number in range(208521, 208532):
        ws_sc = api.Load("NOM_%s" % run_number)

        with open(TEMPLATE) as fd:
            template = json.load(fd)

        ws_names = pnf.slice_data(ws_sc, 30)

        pnf.process_slices(ws_names, "nom_%s" % run_number, OUTPUT_DIR, template)


reduce_all()
