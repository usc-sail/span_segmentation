span_segmentation
---
Segmentation tools for real-time magnetic resonance imaging of speech implemented in MATLAB

Bresch, E., & Narayanan, S. (2009). Region segmentation in the frequency domain applied to upper airway real-time magnetic resonance images. *IEEE Trans Med Imaging*, 28(3), 323-338.

Usage instructions:
1. Run `wrap_make_template.m` from the MATLAB Command Window.
  - Navigate between templates using `previous template` and `next template` buttons
  - Switch frames using `previous frame` and `next frame` buttons
  - Select an articulator in the drop-down menu, then modify the contour using `edit` or draw a new contour using `draw`. When done with the contour, click `save` (in the `line segments` menu, not the `templates` menu).
  - Click `delete` in the `templates` menu to delete the current template.
  - Click `replace` in the `templates` menu to update the template with the current contours.
  - Click `add` in the `templates` menu to make a new template with the current contours.
  - When you are all done with your templates, click `convert and exit`.
2. Run `wrap_make_batch.m` to generate the files to be run on the USC HPC Cluster (https://hpcc.usc.edu/). These files appear in the folder `cluster`.
3. Copy `cluster` to the directory on the USC HPC cluster whose path you set in the string variable `hpcFolder` in the script `wrap_make_batch.m`.
4. On the cluster, run the MATLAB script using the command `source /usr/usc/matlab/default/setup.sh; matlab`. This is best done using `screen`. If this is the first time using parallel MATLAB on the cluster, configure Parallel MATLAB and the Distributed Computing Server Toolbox using the instructions in the file `Parallel Matlab on HPC Cluster.docx`.
