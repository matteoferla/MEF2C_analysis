# original filename = model_prep.ipynb

import pyrosetta


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


# ======== pose =============================================================
pdb_filename = 'MEF2C.pdb'
ccp4_filename = '6byy.ccp4'
constraint_file = None
param_filenames = []

# load
pose = pyrosetta.rosetta.core.pose.Pose()
if param_filenames:
    params_paths = pyrosetta.rosetta.utility.vector1_string()
    params_paths.extend(param_filenames)
    resiset = pyrosetta.generate_nonstandard_residue_set(pose, params_paths) # pdb6bq1.pdb
pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdb_filename) # pose_from_pdbstring or pose_from_file

# ======== ED setup =============================================================
scorefxnED = pyrosetta.get_fa_scorefxn()
ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(ccp4_filename)
sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
sdsm.apply(pose)
## Set ED constraint
elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
scorefxnED.set_weight(elec_dens_fast, 30)
## Set generic constraints
if constraint_file:
    stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    for contype_name in ("atom_pair_constraint", "angle_constraint", "dihedral_constraint"):
        contype = stm.score_type_from_name(contype_name)
        scorefxnED.set_weight(contype, 5)
    setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    setup.constraint_file(constraint_file)
    setup.apply(pose)

# ======== Swissmodel madness fix =============================================================

movemap = pyrosetta.MoveMap()
movemap.set_bb(False)
movemap.set_chi(True)

scorefxnED.set_weight(elec_dens_fast, 30)
relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 5)
relax.set_movemap(movemap)
# relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
relax.apply(pose)

# ======== Relax proper =============================================================

for w in (30, 20, 10):
    scorefxnED.set_weight(elec_dens_fast, w)
    relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 5)
    relax.apply(pose)

pose.dump_pdb(pdb_filename.replace('.pdb', '.r.pdb'))

# ======== Check =============================================================
# The following is Jupyter notebook only:

# import nglview
#
# view = nglview.show_rosetta(pose)
# view