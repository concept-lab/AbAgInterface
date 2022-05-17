#!/bin/bash

for cada in `cat ../../data/pdb.list`
do
    nwat=`grep -c HOH ../../structures/raw/$cada.pdb`
    if [[ $nwat -ne 0 ]]
    then
        echo $cada >> wat_pdb.list
    fi
done
