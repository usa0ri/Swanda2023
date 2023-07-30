# !/bin/bash

curdir=`echo $PWD`
datadir=$1

docker run -it \
    --name $2 \
    -v $curdir:/home/result \
    -v $datadir:/home/data \
    -v /home/saori/Git/myRiboSeq_res/ref:/home/ref \
    saori/riboseq:03
docker rm $2
