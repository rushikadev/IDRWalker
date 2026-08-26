import numpy as np
from IDRWalker import *

chainA = Chain('spe18.fasta', 'spe18noIDR.pdb', chainID = 'A')

print(chainA.gap_list)

box = Box(np.array((100.0, 100.0, 100.0)), 2.0)
box.add_chain(chainA)
box.run()

chainA.write('Spe18_out3.pdb')