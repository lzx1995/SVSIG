#!/bin/bash
# module avail 
module load gcc/6.5.0
# g++ --version > g++_version
# gcc --version > gcc_version
# whereis gcc > gcc_version
make clean 
export LD_PRELOAD=../lib/mimalloc/out/release/libmimalloc.so
make -j

./PageRank -numberOfUpdateBatches 5 -nEdges 1000 -streamPath ~/dataset/twitter/twitter-random-op.txt -block_size 4194304 -maxIters 50 -outputFile ../output/correct  ~/dataset/twitter/twitter-half.adj >correct

# ./PageRank -numberOfUpdateBatches 5 -nEdges 1000 -streamPath ~/dataset/twitter/twitter-random-op.txt -block_size 4194304 -maxIters 50 -outputFile ../output/dotzero1001  ~/dataset/twitter/twitter-half.adj >dotzero1001_res
# ./PageRank -numberOfUpdateBatches 1 -nEdges 1000 -streamPath ~/dataset/wikipedia_link_en/wikipedia_link_en-random-op.txt -block_size 1048576 -maxIters 50 -ae 1   ~/dataset/wikipedia_link_en/wikipedia_link_en-half.adj >./activate_edge.txt
# sh ./DZIG_twitter.sh
# sh ./DZIG_twitter_mpi.sh
# sh ./DZIG_friendster.sh

# sh ./DZIG_twitter_iter10.sh
# sh ./DZIG_twitter_mpi_iter10.sh
# sh ./DZIG_friendster_iter10.sh

# sh ./DZIG_twitter_iter20.sh
# sh ./DZIG_twitter_mpi_iter20.sh
# sh ./DZIG_friendster_iter30.sh

# sh ./DZIG_twitter_iter30.sh
# sh ./DZIG_twitter_mpi_iter30.sh
# sh ./DZIG_friendster_iter30.sh