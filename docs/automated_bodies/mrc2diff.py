import os, sys
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
#
# chimera --script map.mrc output.mrc thresh1 thresh2
# (where thresh1 > thresh2 )
#
input_mrc  = sys.argv[1]
output_mrc = sys.argv[2]
thresh1    = sys.argv[3]
thresh2    = sys.argv[4]

rc("open " + input_mrc)
rc("vop add #0 ")
rc("vol all step 1")
rc("vop threshold #0 minimum " + thresh1 + " set 0")
rc("vop threshold #1 minimum " + thresh2 + " set 0")
rc("vop subtract #3 #2")
rc("vop gaussian #4 sDev 4")
rc("vol #5 save " + output_mrc)
quit()
