#!/bin/bash

#Source the corect version of root
source /home/camilo/HEPTools/ROOT/root/Root6Python2_7/bin/thisroot.sh

#Path to MG 
export MadGrapgSYS=/home/camilo/HEPTools/MADGRAPH/MG5_aMC_v2_5_5_Root6
export PATH=$PATH:$MadGrapgSYS/bin

#Delphes
export DELPHES=/home/camilo/HEPTools/MADGRAPH/MG5_aMC_v2_5_5_Root6/Delphes
export LD_LIBRARY_PATH=$DELPHES:$LD_LIBRARY_PATH
export EXROOTANALYSYS=/home/camilo/HEPTools/MADGRAPH/MG5_aMC_v2_5_5_Root6/ExRootAnalysis
export LD_LIBRARY_PATH=$EXROOTANALYSYS:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=./:$DELPHES:$DELPHES/external:$EXROOTANALYSYS/ExRootAnalysis:$ROOT_INCLUDE_PATH

#LHAPDF To calculate systematica
export LD_LIBRARY_PATH=/home/camilo/HEPTools/LHAPDF/LHAPDF_6_2_1:$LD_LIBRARY_PATH
export PYTHONPATH=/home/camilo/HEPTools/LHAPDF:$PYTHONPATH
