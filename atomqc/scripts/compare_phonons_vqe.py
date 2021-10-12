from jarvis.io.qiskit.inputs import HermitianSolver
from dask import delayed, compute
from jarvis.db.figshare import data
from jarvis.db.figshare import (
    get_wann_phonon,
    get_hk_tb,
    get_wann_electron,
)
from jarvis.db.jsonutils import dumpjson
import os

fls = data("raw_files")


def get_phonon(jid=""):
    try:
        w, atoms = get_wann_phonon(jid=jid)
        nwan = w.nwan
        if nwan <= 32:
            hk = get_hk_tb(w=w)
            info = {}
            herm = HermitianSolver(hk)
            np_vals = list(herm.run_numpy()[0])
            max_vqe, _, _ = herm.run_vqe(mode="max_val")
            min_vqe, _, _ = herm.run_vqe(mode="min_val")
            info["jid"] = jid
            info["np_max_vals"] = max(np_vals)
            info["np_min_vals"] = min(np_vals)
            info["vqe_max_val"] = max_vqe
            info["vqe_min_val"] = min_vqe
            info["formula"] = atoms.composition.reduced_formula
            info["nwan"] = nwan
            name = os.getcwd() + "/" + jid + "_phonon.json"
            print(info)
            if not os.path.exists(name):
                dumpjson(data=info, filename=name)
    except Exception as exp:
        print(jid, exp)
        pass


all_jids = []
for i in fls["FD-ELAST"]:
    if isinstance(i, dict):
        jid = i["name"].split(".zip")[0]
        all_jids.append(jid)
#all_jids=all_jids[0:20]
if __name__ == "__main__":
    x=[get_phonon(i) for i in all_jids]
    #values = [delayed(get_phonon)(i) for i in all_jids]
    #resultsDask = compute(*values)#, scheduler="processes")
