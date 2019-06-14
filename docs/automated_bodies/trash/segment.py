# Set segmentation parameters and segment map.
from Segger import segment_dialog
d = segment_dialog.show_volume_segmentation_dialog()
d.numSteps.set(4)
d.stepSize.set(2)
d.groupMode.set ('smooth')
d.Segment()
# The following line will show a file save dialog.
d.WriteAllRegionsMRCFile()
