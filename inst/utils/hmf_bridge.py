# -*- coding: utf-8 -*-
"""
# hmf_bridge.py  — PHANTOM bridge to the HMF Python package.
# Accepts a JSON task file and writes a text output file.

# Usage (called by MATLAB via system()):
#     python hmf_bridge.py task.json output.txt
# """
# import sys
# import json
# import numpy as np

# def run(task_file, output_file):
#     with open(task_file, 'r') as f:
#         task = json.load(f)

#     from hmf import MassFunction

#     M_arr   = np.array(task['M'])        # [Msun/h]
#     z       = float(task.get('z', 0.0))
#     model   = task.get('model', 'ST')    # hmf model string
#     cosmo_p = task.get('cosmo_params', {})

#     # Default Planck18-like params if not provided
#     H0      = cosmo_p.get('H0',      67.66)
#     Om0     = cosmo_p.get('Om0',     0.3111)
#     Ob0     = cosmo_p.get('Ob0',     0.0490)
#     sigma8  = cosmo_p.get('sigma8',  0.8102)
#     ns      = cosmo_p.get('ns',      0.9665)

#     from hmf import MassFunction

#     cp = task["cosmo_params"]

#     hmf_obj = MassFunction(
#         Mmin=np.log10(M_arr.min()),
#         Mmax=np.log10(M_arr.max()),
#         dlog10m=np.log10(M_arr[1] / M_arr[0]),
#         z=task["z"],
#         hmf_model=task["model"],
#         transfer_model="EH",
#         cosmo_params={
#             "H0": cp["H0"],
#             "Om0": cp["Om0"],
#             "Ob0": cp["Ob0"],
#             "sigma8": cp["sigma8"],
#             "ns": cp["ns"],
#         },
#     )
    
    
#     # Interpolate onto requested M grid
#     from scipy.interpolate import interp1d
#     # hmf gives dndlnM = m * dndm
#     dndlnM_hmf = hmf_obj.m * hmf_obj.dndm   # [h^3 Mpc^-3]
#     interp_fn  = interp1d(np.log10(hmf_obj.m), np.log10(dndlnM_hmf),
#                           kind='linear', bounds_error=False,
#                           fill_value=-np.inf)
#     dndlnM_out = 10.0 ** interp_fn(np.log10(M_arr))

#     np.savetxt(output_file,
#                np.column_stack([M_arr, dndlnM_out]),
#                header='M_Msun_h dndlnM_h3Mpc3', comments='')

# if __name__ == '__main__':
#     run(sys.argv[1], sys.argv[2])
import sys
import json
import numpy as np
from scipy.interpolate import interp1d

def run(task_file, output_file):
    with open(task_file, "r") as f:
        task = json.load(f)

    from hmf import MassFunction

    M_arr = np.asarray(task["M"], dtype=float)
    z = float(task.get("z", 0.0))
    model = task.get("model", "ST")
    cp = task.get("cosmo_params", {})

    H0 = float(cp.get("H0", 67.66))
    Om0 = float(cp.get("Om0", 0.3111))
    Ob0 = float(cp.get("Ob0", 0.0490))
    sigma8 = float(cp.get("sigma8", 0.8102))
    ns = float(cp.get("ns", 0.9665))

    hmf_obj = MassFunction(
        Mmin=np.log10(M_arr.min()),
        Mmax=np.log10(M_arr.max()),
        dlog10m=np.log10(M_arr[1] / M_arr[0]),
        z=z,
        hmf_model=model,
        transfer_model="EH",
        sigma_8=sigma8,
        n=ns,
        cosmo_params={
            "H0": H0,
            "Om0": Om0,
            "Ob0": Ob0,
            "Tcmb0": 2.7255,
        },
    )

    dndlnM_hmf = hmf_obj.m * hmf_obj.dndm
    interp_fn = interp1d(
        np.log10(hmf_obj.m),
        np.log10(dndlnM_hmf),
        kind="linear",
        bounds_error=False,
        fill_value=-np.inf,
    )
    dndlnM_out = 10.0 ** interp_fn(np.log10(M_arr))

    np.savetxt(
        output_file,
        np.column_stack([M_arr, dndlnM_out]),
        header="M_Msun_h dndlnM_h3Mpc3",
        comments="",
    )

if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])