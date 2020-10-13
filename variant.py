import pyrosetta, re
from typing import Optional, List, Dict


class Variant:
    """
    Copy pasted from PI4KA <- GNB2 <- SnoopCatcher
    """
    strict_about_starting_residue = True

    _name3 = {'A': 'ALA',
              'C': 'CYS',
              'D': 'ASP',
              'E': 'GLU',
              'F': 'PHE',
              'G': 'GLY',
              'H': 'HIS',
              'I': 'ILE',
              'L': 'LEU',
              'K': 'LYS',
              'M': 'MET',
              'N': 'ASN',
              'P': 'PRO',
              'Q': 'GLN',
              'R': 'ARG',
              'S': 'SER',
              'T': 'THR',
              'V': 'VAL',
              'W': 'TRP',
              'Y': 'TYR'}

    def __init__(self, filename: str, params_filenames: Optional[List[str]] = None):
        self.pose = self.load_pose_from_file(filename, params_filenames)

    def load_pose_from_file(self, filename: str, params_filenames: Optional[List[str]] = None) -> pyrosetta.Pose:
        """
        Loads a pose from filename with the params in the params_folder
        :param filename:
        :param params_filenames:
        :return:
        """
        pose = pyrosetta.Pose()
        if params_filenames:
            params_paths = pyrosetta.rosetta.utility.vector1_string()
            params_paths.extend(params_filenames)
            pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)
        return pose

    def relax_around_mover(self,
                           pose: pyrosetta.Pose,
                           resi: int, chain: str,
                           scorefxn=None, cycles=5, distance=5,
                           cartesian=False, own_chain_only=False) -> None:
        """
        Relaxes pose ``distance`` around resi:chain.

        :param resi: PDB residue number.
        :param chain:
        :param pose:
        :param scorefxn:
        :param cycles: of relax (3 quick, 15 thorough)
        :param distance:
        :param cartesian:
        :param own_chain_only:
        :return:
        """
        if pose is None:
            pose = self.pose
        if scorefxn is None:
            scorefxn = pyrosetta.get_fa_scorefxn()
            # self._cst_score(scorefxn)
        movemap = pyrosetta.MoveMap()
        ####
        n = self.get_neighbour_vector(pose=pose, resi=resi, chain=chain, distance=distance,
                                      own_chain_only=own_chain_only)
        print(pyrosetta.rosetta.core.select.residue_selector.ResidueVector(n))
        movemap.set_bb(False)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(False)
        movemap.set_chi(allow_chi=n)
        movemap.set_jump(False)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap(movemap)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.cartesian(cartesian)
        relax.apply(pose)

    def get_neighbour_vector(self, pose: pyrosetta.Pose, resi: int, chain: str, distance: int,
                             include_focus_in_subset: bool = True,
                             own_chain_only: bool = False) -> pyrosetta.rosetta.utility.vector1_bool:
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        resi_sele.set_index(pose.pdb_info().pdb2pose(chain=chain, res=resi))
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(resi_sele, distance=distance,
                                                 include_focus_in_subset=include_focus_in_subset)
        if own_chain_only:
            chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain)
            and_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(neigh_sele, chain_sele)
            return and_sele.apply(pose)
        else:
            return neigh_sele.apply(pose)

    def make_mutant(self,
                    pose: pyrosetta.Pose,
                    mutation: str,
                    chain='A',
                    distance: int = 10,
                    cycles: int = 5
                    ) -> pyrosetta.Pose:
        """
        Make a point mutant (``A23D``).
        :param pose: pose
        :param mutation:
        :param chain:
        :return:
        """
        if pose is None:
            mutant = self.pose.clone()
        else:
            mutant = pose.clone()
        pose2pdb = pose.pdb_info().pdb2pose
        rex = re.match('(\w)(\d+)(\w)', mutation)
        r = pose2pdb(res=int(rex.group(2)), chain=chain)
        rn = pose.residue(r).name1()
        if self.strict_about_starting_residue:
            assert rn == rex.group(1), f'residue {r}(pose)/{rex.group(2)}(pdb) is a {rn}, not a {rex.group()}'
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res=self._name3[rex.group(3)]).apply(mutant)
        self.relax_around_mover(mutant, int(rex.group(2)), chain, distance=distance, cycles=cycles,
                                own_chain_only=False)
        return mutant

    def score_interface(self, pose: pyrosetta.Pose, interface: str) -> Dict[str, float]:
        if pose is None:
            pose = self.pose
        assert self.has_interface(pose, interface), f'There is no {interface}'
        ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
        ia.apply(pose)
        return {'complex_energy': ia.get_complex_energy(),
                'separated_interface_energy': ia.get_separated_interface_energy(),
                'complexed_sasa': ia.get_complexed_sasa(),
                'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
                'interface_dG': ia.get_interface_dG(),
                'interface_delta_sasa': ia.get_interface_delta_sasa()}

    def has_interface(self, pose: pyrosetta.Pose, interface: str) -> bool:
        if pose is None:
            pose = self.pose
        pose2pdb = pose.pdb_info().pose2pdb
        have_chains = {pose2pdb(r).split()[1] for r in range(1, pose.total_residue() + 1)}
        want_chains = set(interface.replace('_', ''))
        return have_chains == want_chains

    def has_residue(self, pose: pyrosetta.Pose, resi: int, chain: str) -> bool:
        if pose is None:
            pose = self.pose
        pdb2pose = pose.pdb_info().pdb2pose
        r = pdb2pose(res=resi, chain=chain)
        return r != 0

    def vector2list(self, vector: pyrosetta.rosetta.utility.vector1_bool) -> pyrosetta.rosetta.std.list_unsigned_long_t:
        rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(vector)
        x = pyrosetta.rosetta.std.list_unsigned_long_t()
        assert len(rv) > 0, 'Vector is empty!'
        for w in rv:
            x.append(w)
        return x

    def CA_RMSD(self, poseA: pyrosetta.Pose, poseB: pyrosetta.Pose, resi: int, chain: str, distance: int) -> float:
        """
        Carbon alpha atoms only RMSD.

        :param poseA:
        :param poseB:
        :param resi:
        :param chain:
        :param distance:
        :return:
        """
        n = self.get_neighbour_vector(pose=poseA, resi=resi, chain=chain, distance=distance)
        residues = self.vector2list(n)
        return pyrosetta.rosetta.core.scoring.CA_rmsd(poseA, poseB, residues)

    def FA_RMSD(self, poseA: pyrosetta.Pose, poseB: pyrosetta.Pose, resi: int, chain: str, distance: int) -> float:
        """
        Full-atom, all atom RMSD.
        :param poseA:
        :param poseB:
        :param resi:
        :param chain:
        :param distance:
        :return:
        """
        n = self.get_neighbour_vector(pose=poseA, resi=resi, chain=chain, distance=distance,
                                      include_focus_in_subset=False)
        residues = self.vector2list(n)
        # pyrosetta.rosetta.core.scoring.automorphic_rmsd(residueA, residueB, False)
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd(poseA, poseB, residues)

    def does_contain(self, mutation: str, chain: str) -> bool:
        r = self.pose.pdb_info().pdb2pose(chain=chain, res=int(mutation[1:-1]))
        if r == 0:
            return False
        if self.pose.residue(r).name1() not in self._name3.keys():
            return False
        else:
            return True
