### Install :
    git clone --recursive https://github.com/ptdtan/ancestor.git
    cd ancestor/bamtools
    mkdir build && cd build && cmake .. && make
    cd src && make
    export LD_LIBRARY_PATH=$LD_LIBARY_PATH:/path/to/ancestor/bamtools/lib
### Usage:
    ./overlap_graph <bamFile> <ThresHold>


### crocodile:
- ghaGan1 N50: 96,076,944
- croPor2 N50: 84,437,661
- allMis3 N50: 47,204,337
