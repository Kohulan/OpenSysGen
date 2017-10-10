# OpenSysGen : A tool for systematic conformer generation

This project describes a tool for systematic conformer generation of small molecules without the use of an energy cut-off , which is a user friendly tool, can simply be used by anyone who needs conformers to be generated to their projects. The tool is written in JAVA which incorporates Open Babel 2.4.1 JAVA binding libraries. Rotatable bonds are identified and systematically rotated with a pre-set increment. The pre-set increments for each class of rotatable bonds (from 0 to 10) was derived based on a dataset of 110 approved drugs. The tool was tested on a dataset of 428 drugs and results reveal that in all cases, conformers within 1 Å RMSD of the bioactive conformer were retrieved. 

Summary of dependencies:
-------------------------------
OpenSysGen is a tool built using JAVA so to run the tool Java Runtime Environment is thus required.

Also it uses Open Babel JAVA libraries and openbabel.jar is incuded in the OpenSysGen.jar file. To run the tool make sure OpenBabel 2.4.1 is installed on your system and the path is correctly set.

Licensing
---------

OpenSysGen is distributed under the GNU General Public License (GPL) version 2. This program is free software. A copy of the license file can be found in the folder.

Open Babel
> The Open Babel package is a separate open source package.
  Open Babel is licensed under the
  [GNU General Public License, version 2](http://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).
  Only the babel executable is needed (and included) in the chemalot package.
  More information on Open Babel is available at:
  [http://openbabel.org](http://openbabel.org)
