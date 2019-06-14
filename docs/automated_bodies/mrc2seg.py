import os, sys
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages
from Segger import segment_dialog
from chimera import openModels

input_mrc      = sys.argv[1]
output_mrc     = sys.argv[2]
threshold      = sys.argv[3]
numSteps       = sys.argv[4]
stepSize       = sys.argv[5]
minRegionSize  = sys.argv[6]
minContactSize = sys.argv[7]

rc("open " + input_mrc)
rc("vol all step 1")
rc("vol all sdLevel " + threshold)

d = segment_dialog.show_volume_segmentation_dialog()
d.minRegionSize.set(minRegionSize)
d.minContactSize.set(minContactSize)
d.groupMode.set('smooth')
d.numSteps.set(numSteps)
d.stepSize.set(stepSize)
d.Segment()
d.SaveSegmentation()

rc("segment exportmask #1 savePath " + output_mrc)
quit()
