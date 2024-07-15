#!/bin/sh


listfile=list.dat

cat $listfile | while read line
do
    filename=$line
    #echo root -l -b -q convert_artroot_wf_uvw.C\(\"$filename\"\)
    #root -l -b -q  convert_artroot_wf_uvw.C\(\"$filename\"\)
    root -l -b -q  convert_artroot_DQM_2.C\(\"$filename\"\)

done


