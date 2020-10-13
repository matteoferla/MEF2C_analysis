# original name = mutagenesis.ipynb

import pyrosetta
from variant import Variant



# ======== init pyrosetta =============================================================
# https://www.rosettacommons.org/docs/latest/full-options-list
def format_options(**options):
    def format(v):
        if v is None:
            return ''
        elif v in (False, True):
            return ' ' + str(v).lower()
        else:
            return ' ' + str(v)

    args = ['-' + k + format(options[k]) for k in options]
    return ' '.join(args)


pyrosetta.init(extra_options=format_options(no_optH=False,
                                            ex1=None,
                                            ex2=None,
                                            mute='all',
                                            ignore_unrecognized_res=True,
                                            load_PDB_components=False,
                                            ignore_waters=False)
               )

# ======== variants =============================================================
from Bio.SeqUtils import seq3, seq1

case_list = '''p.Arg3Ser
p.Arg15His
p.Lys23Arg
p.Arg24Lys
p.Gly27Ala
p.Leu38Gln
p.Val65Gly
p.Arg3Gly
p.Lys5Met
p.Arg15Cys
p.Arg15Pro
p.Lys30Asn
p.Tyr33Asp
p.Ser36Arg
p.Ser36Arg
p.Arg79Pro
p.Pro193Ser
p.Ala309Thr
p.Arg411Cys'''.replace('p.','').split()

cases = [seq1(case[:3])+case[3:-3]+seq1(case[-3:]) for case in case_list]
MEF2C_gnomAD = ['R89K','S50G','V37L', 'V37M']
MEF2A_gnomAD = ['R10C', 'R10S', 'M12L', 'L28F', 'M29L', 'L35F', 'D40E', 'F48V', 'S50T', 'N52D', 'F55V', 'Q56E',
           'D61N', 'D61H', 'M62T', 'D63G', 'N73S', 'S82L', 'V85A', 'N89K']
MEF2B_gnomAD = ['R3K', 'K4Q', 'Q7R', 'R10C', 'R10P', 'R10H', 'R17W', 'R17Q', 'R24W', 'G27R', 'L28V', 'A32D',
           'Y33C', 'E34D', 'V37M','A44G', 'I46T', 'N49S', 'A51T', 'N52K', 'N52S', 'R53C', 'R53H',
           'T60M', 'M62V', 'M62I', 'R64C', 'R64H', 'S73G', 'S73R', 'E74K', 'E77K', 'S78N',
           'R79C', 'R79H', 'T82S', 'I84M', 'E86K', 'T87M', 'L88M', 'K89R', 'R90W', 'R90N', 'R91K']
MEF2D_gnomAD = ['D13N', 'R15W', 'R15Q', 'E34Q', 'N49S', 'H50Q', 'N52S', 'D61N', 'N73S', 'A82T', 'K80R']

variants = cases + MEF2C_gnomAD + MEF2A_gnomAD + MEF2B_gnomAD + MEF2D_gnomAD

# ======== pose =============================================================
pdb_filename = 'MEF2C_long.rx.pdb'

param_filenames = []

model = Variant(filename=pdb_filename)


# ======== score =============================================================
model.strict_about_starting_residue = False

import os, csv, re

if not os.path.exists('variants'):
    os.mkdir('variants')

for chain in ('A', 'B'):
    interfaces = (('DNA', 'AB_KLW'), ('homodimer_B', 'AKLW_B'), ('homodimer_A', 'BKLW_A'))


    modelname = f'MEF2C_DNA_hydrate_{chain}'

    with open(modelname+'_scores.csv', 'w') as w:
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
        for mutation in variants:
            print(mutation)
            if not model.does_contain(mutation, chain):
                print('Absent')
                continue
            variant = model.make_mutant(model.pose,
                                        mutation=mutation,
                                        chain=chain,
                                        distance=10,
                                        cycles=5)
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

