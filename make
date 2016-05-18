#!/bin/bash
g++ -ggdb -g3 -I /hive/groups/recon/projs/mus_strain_cactus/ragout/variant/include/bamtools/ -L /hive/groups/recon/projs/mus_strain_cactus/ragout/variant/lib/ -o $1 $1.cpp -lz -lbamtools
