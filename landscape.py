# not using pmut_scan

import os, csv, re
import pyrosetta
from variant import Variant

# ======== pose =============================================================
pdb_filename = 'MEF2C_long.rx.pdb'

param_filenames = []

model = Variant(filename=pdb_filename)


# ======== score =============================================================

if not os.path.exists('landscape'):
    os.mkdir('landscape')

for chain in ('A', 'B'):
    interfaces = (('DNA', 'AB_KLW'), ('homodimer_B', 'AKLW_B'), ('homodimer_A', 'BKLW_A'))


    modelname = f'MEF2C_DNA_hydrate_{chain}'

    with open(modelname+'_scan.csv', 'w') as w:
        fieldnames = ['model',
                    'mutation',
                    'complex_ddG',
                    'complex_native_dG',
                    'complex_mutant_dG',
                     'FA_RMSD',
                     'CA_RMSD'] + [f'{interface_name}_{suffix}' for suffix in ('interface_native_dG',
                                                                                   'interface_mutant_dG',
                                                                                   'interface_ddG')
                                                                    for interface_name, interface_scheme in interfaces]
        out = csv.DictWriter(w, fieldnames=fieldnames)
        out.writeheader()
        scorefxn = pyrosetta.get_fa_scorefxn()
        ## wt
        n = scorefxn(model.pose)
        ref_interface_dG = {}
        for interface_name, interface_scheme in interfaces:
            ref_interface_dG[interface_name] = model.score_interface(model.pose, interface_scheme)['interface_dG']
        ## muts
        seq='MGRKKIQITRIMDERNRQVTFTKRKFGLMKKAYELSVLCDCEIALIIFNSTNKLFQYASTDMDKVLLKYTEYNEPHESRTNSDIVETLRKKG'
        for resi in range(2, 93):
            for to_resn in 'IVLFCMAGTSWYPHNDEQKR': # order take from https://direvo.mutanalyst.com/main/landscape
                mutation = f'{seq[resi-1]}{resi}{to_resn}'
                print(mutation)
                if not model.does_contain(mutation, chain):
                    print('Absent')
                    continue
                variant = model.make_mutant(model.pose,
                                            mutation=mutation,
                                            chain=chain,
                                            distance=6,
                                            cycles=3)
                variant.dump_pdb(f'variants/{modelname}.{mutation}.pdb')
                m = scorefxn(variant)
                data = {'model': modelname,
                          'mutation': mutation,
                          'complex_ddG': m - n,
                          'complex_native_dG': n,
                          'complex_mutant_dG': m,
                           'FA_RMSD': model.FA_RMSD(model.pose, variant, resi=int(mutation[1:-1]), chain=chain, distance=10),
                           'CA_RMSD': model.CA_RMSD(model.pose, variant, resi=int(mutation[1:-1]), chain=chain, distance=10)
                          }
                for interface_name, interface_scheme in interfaces:
                    if model.has_interface(variant, interface_scheme):
                        print(f'{interface_name} ({interface_scheme}) applicable to {modelname}')
                        i = model.score_interface(variant, interface_scheme)['interface_dG']
                    else:
                        print(f'{interface_name} ({interface_scheme}) not applicable to {modelname}')
                        i = float('nan')
                    data[f'{interface_name}_interface_native_dG'] = ref_interface_dG[interface_name]
                    data[f'{interface_name}_interface_mutant_dG'] = i
                    data[f'{interface_name}_interface_ddG']= i - ref_interface_dG[interface_name]
                out.writerow(data)
