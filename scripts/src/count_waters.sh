#!/bin/bash

for cada in `cat ../../data/pdb.list`
do
    grep -c HOH ../../structures/raw/$cada.pdb >> wat_count
done
