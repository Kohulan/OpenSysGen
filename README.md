# OpenSysGen
A tool for systematic conformer generation

This project describes a tool for systematic conformer generation of small molecules without the use of an energy cut-off , which is a user friendly tool, can simply be used by anyone who needs conformers to be generated to their projects. The tool is written in JAVA which incorporates Open Babel 2.4.1 JAVA binding libraries. Rotatable bonds are identified and systematically rotated with a pre-set increment. The pre-set increments for each class of rotatable bonds (from 0 to 10) was derived based on a dataset of 110 approved drugs. The tool was tested on a dataset of 428 drugs and results reveal that in all cases, conformers within 1 Å RMSD of the bioactive conformer were retrieved. 

Licensing

OpenSysGen is distributed under the GNU General Public License (GPL) version 2. This program is free software. A copy of the license file can be found in the folder.
