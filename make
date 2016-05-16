#!/bin/bash
g++ -I /hive/groups/recon/projs/mus_strain_cactus/ragout/variant/include/bamtools/ -L /hive/groups/recon/projs/mus_strain_cactus/ragout/variant/lib/ -o verify verify.cpp -lz -lbamtools
