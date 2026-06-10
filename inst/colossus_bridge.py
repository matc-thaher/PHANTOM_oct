# -*- coding: utf-8 -*-
"""
Created on Sun May 24 17:48:51 2026

@author: atcna
"""

import sys
import json
import numpy as np
from colossus.cosmology import cosmology as colossus_cosmo
from colossus.halo import concentration as colossus_conc
from colossus.lss import peaks

def run(task_file, output_file):
    with open(task_file, 'r') as f:
        task = json.load(f)

    # Set cosmology
    cosmo_id = task.get('cosmo_id', 'planck15')
    cosmo = colossus_cosmo.setCosmology(cosmo_id)

    quantity = task['quantity']
    z        = task.get('z', 0.0)

    if quantity == 'power_spectrum':
        from colossus.cosmology import power_spectrum as ps
        k     = np.array(task['k'])
        model = task.get('model', 'eisenstein98')
        Pk    = cosmo.matterPowerSpectrum(k, z, model=model)
        np.savetxt(output_file, np.column_stack([k, Pk]),
                   header='k_h_Mpc P_Mpc_h3', comments='')

    elif quantity == 'variance':
        R     = np.array(task['R'])
        model = task.get('model', 'eisenstein98')
        sig   = cosmo.sigma(R, z, ps_args={'model': model})
        np.savetxt(output_file, np.column_stack([R, sig]),
                   header='R_Mpc_h sigma', comments='')

    elif quantity == 'correlation_function':
        r     = np.array(task['r'])
        model = task.get('model', 'eisenstein98')
        xi    = cosmo.correlationFunction(r, z, ps_args={'model': model})
        np.savetxt(output_file, np.column_stack([r, xi]),
                   header='r_Mpc_h xi', comments='')
    
    elif quantity == 'concentration':
        M     = np.array(task['M'])
        mdef  = task.get('mdef', '200c')
        model = task.get('model', 'ishiyama21')
        c = colossus_conc.concentration(M, mdef, z, model=model,
                                        range_return=False,
                                        range_warning=False)
        np.savetxt(output_file, np.column_stack([M, np.atleast_1d(c)]),
                   header='M_Msun_h concentration', comments='')

    elif quantity == 'peakheight':
        M = np.array(task['M'])        # <-- must assign M here explicitly
        # z = float(task['z'])
        result = peaks.peakHeight(M, z)
        np.savetxt(output_file,
               np.column_stack([M, np.atleast_1d(result)]),
               header='M_Msun_h peakheight', comments='')

    elif quantity == 'neff':
        M   = np.array(task['M'])      # <-- and here
        kappa = 1.00
        R     = peaks.lagrangianR(M)
        k_R   = 2.0 * np.pi / R * kappa
        result = colossus_cosmo.getCurrent().matterPowerSpectrum(
                 k_R, model='eisenstein98_zb', derivative=True)
        np.savetxt(output_file,
               np.column_stack([M, np.atleast_1d(result)]),
               header='M_Msun_h neff', comments='')
    
    elif quantity == 'profile_batch':
        from colossus.halo import profile_nfw, profile_einasto, \
                              profile_hernquist, profile_dk14
        r   = np.array(task['r'])
        M   = float(task['M']);  c = float(task['c'])
        mdef = task.get('mdef', '200c')
        models = {'nfw': profile_nfw.NFWProfile,
              'einasto': profile_einasto.EinastoProfile,
              'hernquist': profile_hernquist.HernquistProfile,
              'dk14': profile_dk14.DK14Profile}
        out = np.column_stack([r] + [
              models[k](M=M, c=c, z=z, mdef=mdef).density(r)
              for k in ['nfw','einasto','hernquist','dk14']])
        np.savetxt(output_file, out, delimiter=',',
               header='r_kpch,nfw,einasto,hernquist,dk14', comments='')

    elif quantity == 'hmf':
        from colossus.lss import mass_function as col_mf
        M     = np.array(task['M'])
        model = task.get('model', 'tinker08')

        # mdef defaults depend on calibration: Tinker08 uses '200m',
        # all FOF-based models (PS, ST, Crocce10, Watson13) use 'fof'.
        # fof_models = {'press74', 'sheth99', 'crocce10', 'watson13',
        #               'bhattacharya11', 'courtin11', 'angulo12', 'reed03', 'reed07'}
        # default_mdef = 'fof' if model in fof_models else '200m'
        mdef_overrides = {
            'press74': 'fof',
            'sheth99': 'fof',
            'angulo12': 'fof',
            'bhattacharya11': 'fof',
            'crocce10': 'fof',
            'watson13': 'fof',
            'reed03': 'fof',
            'reed07': 'fof',
            'courtin11': 'fof',
            'tinker08': '200m',
            'comparat17': 'vir',
            'rodriguezpuebla16': 'vir',
            'yung24': 'vir',
            'yung25': 'vir',
        }
        mdef = task.get('mdef', mdef_overrides.get(model, '200m'))
        # mdef  = task.get('mdef', default_mdef)

        mdef_overrides = {
            'yung24': 'vir',
            'yung25': 'vir',
        }
        if model in mdef_overrides:
            mdef = mdef_overrides[model]

        dndlnM = col_mf.massFunction(M, z, mdef=mdef, model=model,
                                     q_out='dndlnM')
        np.savetxt(output_file,
                   np.column_stack([M, dndlnM]),
                   header='M_Msun_h dndlnM_h3Mpc3', comments='')

    elif quantity == 'delta_c':
        from colossus.lss import peaks
        corrections = task.get('corrections', True)
        dc = peaks.collapseOverdensity(z=z, corrections=corrections)
        np.savetxt(output_file,
               np.array([[z, dc]]),
               header='z delta_c', comments='')

    elif quantity == 'sigma_fof':
        from colossus.lss import peaks
        M     = np.array(task['M'])
        model = task.get('model', 'eisenstein98')
        # Colossus sigma uses Lagrangian radius from M directly
        R     = peaks.lagrangianR(M)
        sig   = cosmo.sigma(R, z, ps_args={'model': model})
        np.savetxt(output_file,
               np.column_stack([M, sig]),
               header='M_Msun_h sigma', comments='')

    else:
        raise ValueError(f"Unknown quantity: {quantity}")

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])