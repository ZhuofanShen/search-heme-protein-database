import argparse
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamers, \
    IncludeCurrent, OperateOnResidueSubset, RestrictAbsentCanonicalAASRLT, \
    RestrictToRepackingRLT, PreventRepackingRLT
from pyrosetta.rosetta.core.pack import *
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import *
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    ResidueIndexSelector, ResiduePropertySelector, \
    InterGroupInterfaceByVectorSelector, NeighborhoodResidueSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, \
    AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.minimization_packing import \
    PackRotamersMover, MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import \
    MutateResidue
from os.path import basename

parser = argparse.ArgumentParser()
parser.add_argument('pdb', type=str)
parser.add_argument('-ref', type=str)
parser.add_argument('-params', type=str, nargs='*')
parser.add_argument('-chi', type=float)
parser.add_argument('-cst', type=str)
parser.add_argument('-enzdescst', type=str)
parser.add_argument('-sites', type=int, nargs='*')
args = parser.parse_args()

opts = '-ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false'
if args.params:
    opts += ' -extra_res_fa ' + ' '.join(args.params)
if args.cst:
    opts += ' -constraints:cst_fa_file {}'.format(args.cst)
if args.enzdescst:
    opts += ' -enzdes:cstfile {} -run:preserve_header'.format(args.enzdescst)
init(opts)

print("zs251::Building score function")
score_function = create_score_function('ref2015_cst')
score_function.set_weight(ScoreType.fa_intra_rep_nonprotein, 0.545)
score_function.set_weight(ScoreType.fa_intra_atr_nonprotein, 1)

pose = pose_from_pdb(args.pdb)
ref_pose = pose_from_pdb(args.ref)
fold_tree = FoldTree()
fold_tree.add_edge(1, 383, -1)
fold_tree.add_edge(1, 384, 1)
fold_tree.add_edge(384, 385, 'FE', 'N1')
pose.fold_tree(fold_tree)

if args.chi:
    pose.set_chi(2, 385, args.chi)

if args.cst:
    add_fa_constraints_from_cmdline(pose, score_function)
if args.enzdescst:
    enzdes_cst = AddOrRemoveMatchCsts()
    enzdes_cst.set_cst_action(ADD_NEW)
    enzdes_cst.apply(pose)

coord_cst_gen = CoordinateConstraintGenerator()
coord_cst_gen.set_reference_pose(ref_pose)
add_csts = AddConstraints()
add_csts.add_generator(coord_cst_gen)

# RMSD metric
rmsdm = RMSDMetric()
rmsdm.set_comparison_pose(pose)

prefix = basename(args.pdb).replace('.pdb', '').replace('.gz', '')

print("zs251::Setting up and making a mutation")
for site in args.sites:
    mutation_selector = ResidueIndexSelector(site)
    ts_selector = ResidueIndexSelector(385)
    heme_selector = ResidueIndexSelector(384)
    interface_selector = InterGroupInterfaceByVectorSelector()
    focus_selector = OrResidueSelector(mutation_selector, ts_selector)
    interface_selector.group1_selector(focus_selector)
    interface_selector.group2_selector(NotResidueSelector(focus_selector))
    repacking_selector = AndResidueSelector(OrResidueSelector(interface_selector, ts_selector), \
        NotResidueSelector(OrResidueSelector(mutation_selector, heme_selector)))
    static_selector = NotResidueSelector(OrResidueSelector(mutation_selector, repacking_selector))
    # mutation_repacking_selector = NeighborhoodResidueSelector()
    # mutation_repacking_selector.set_focus_selector(mutation_selector)
    # mutation_repacking_selector.set_distance(8.0)
    # mutation_repacking_selector.set_include_focus_in_subset(True)
    # repacking_selector = AndResidueSelector(mutation_repacking_selector, NotResidueSelector(mutation_selector))

    move_map = MoveMap()
    move_map.set_bb(True)
    side_chain_movable_selector = NotResidueSelector(heme_selector)
    side_chain_movable_vector = side_chain_movable_selector.apply(pose)
    move_map.set_chi(side_chain_movable_vector)

    result_str = '"' + str(site) + '":{'
    # for target_res_name1 in ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', 'Y', \
    #     'D', 'E', 'R', 'H', 'K', 'S', 'T', 'C', 'M', 'N', 'Q']:
    for target_res_name1 in ['A', 'G', 'I', 'L', 'V', 'F', 'W', 'Y', 'S', 'T', 'C', 'M']:
        task_factory = TaskFactory()
        task_factory.push_back(IncludeCurrent())

        # Mutated residue
        aa_force = RestrictAbsentCanonicalAASRLT()
        aa_force.aas_to_keep(target_res_name1)
        task_factory.push_back(OperateOnResidueSubset(aa_force, mutation_selector))

        # Repacking residues
        repack = RestrictToRepackingRLT()
        task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))

        # Immobile side chains
        prevent = PreventRepackingRLT()
        task_factory.push_back(OperateOnResidueSubset(prevent, static_selector))

        fast_relax = FastRelax()
        fast_relax.set_scorefxn(score_function)
        fast_relax.set_task_factory(task_factory)
        fast_relax.set_movemap(move_map)

        for decoy in range(1):
            pose_copy = Pose(pose)
            add_csts.apply(pose_copy)
            fast_relax.apply(pose_copy)
            current_energy = score_function(pose_copy)
            if decoy == 0 or current_energy < lowest_energy:
                lowest_energy = current_energy
                point_mutated_pose = pose_copy
        rmsdm.apply(point_mutated_pose)

        energies = point_mutated_pose.energies().total_energies()
        score = energies.get(ScoreType.total_score) - energies.get(ScoreType.coordinate_constraint)
        
        result_str += '"' + target_res_name1 + '":' + str(round(score, 3)) + ','
        
        point_mutated_pose.dump_pdb(prefix + '_' + str(site) + '_' + target_res_name1 + '.pdb')
    result_str = result_str[:-1] + '}\n'
    print("zs251::Printing rotamer and rotamer list to a file")
    with open('_'.join(str(site) for site in args.sites) + '.out', 'a+') as results:
        results.write(result_str)

