from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.toolbox import cleanATOM
from pyrosetta.toolbox import get_secstruct
from pyrosetta.teaching import *
from pyrosetta.toolbox import get_hbonds
from pyrosetta.toolbox import mutate_residue
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.rigid import *
from pyrosetta.rosetta.protocols.docking import *
import sys
init()

def main():
    filename=sys.argv[1]
    pose=pose_from_pdb(filename)
    scorefxn=get_fa_scorefxn()   
    mutate_residue(pose, pose.pdb_info().pdb2pose('A', 51), "P")
    MC(pose, scorefxn, "PRO")

def MC(pose, scorefxn, mutant):
    test=Pose()
    test.assign(pose)
    dumpfile = 'Folding_Output_51_'+str(mutant)+'.pdb'
    txtfile = 'Folding_Output_51_'+str(mutant)+'.txt'
    newfile = open(txtfile, "w")
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
    kT = 1
    mc = MonteCarlo(test, scorefxn, kT)

    residue=int(test.pdb_info().pdb2pose('A', 51))
    residue=test.residue(residue).xyz("CA")
    count=0
    move_list=[]
    for i in range(1, test.total_residue()+1):
        i_residue=test.residue(i).xyz("CA")
        if (residue-i_residue).norm()<10:
            move_list.append(i)
            count +=1

    min_mover = MinMover()
    mm = MoveMap()
    for i in move_list:
        mm.set_bb(i, True)
        mm.set_chi(i, True)
    min_mover.movemap(mm)
    min_mover.score_function(scorefxn)
    min_mover.min_type("dfpmin")
    min_mover.tolerance(0.01)

    task_pack = standard_packer_task(test)
    task_pack.restrict_to_repacking()
    task_pack.temporarily_fix_everything()
    for i in move_list:
        task_pack.temporarily_set_pack_residue(i,True)
    pack_mover=PackRotamersMover(scorefxn, task_pack)    
    pack_mover.apply(test)
    
    for i in range(300000):
        min_mover.apply(test)
        mc.boltzmann(test)
        newfile.write(str(scorefxn(test)))
        newfile.write(' ')
        newfile.write(str(i))
        newfile.write(' ')
        newfile.write(str(CA_rmsd(pose, test)))
        newfile.write('\n')
    
    mc.recover_low(test)
    newfile.write(str(scorefxn(test)))
    newfile.write('\n')
    newfile.write(str(CA_rmsd(pose, test)))
    newfile.close()
    test.dump_pdb(dumpfile)

main()
