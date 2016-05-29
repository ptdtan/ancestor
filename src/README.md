### Install :
    git clone --recursive https://github.com/ptdtan/ancestor.git
    cd ancestor/bamtools
    mkdir build && cd build && cmake .. && make
    cd src && make
    export LD_LIBRARY_PATH=$LD_LIBARY_PATH:/path/to/ancestor/bamtools/lib
### Usage:
    ./overlap_graph <bamFile> <ThresHold>
