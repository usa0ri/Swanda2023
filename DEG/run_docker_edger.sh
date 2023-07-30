#!/bin/bash
# $1: container name

host_home=$2
docker_home=/home/rstudio

docker rm $1
docker run -it \
    -p 8787:8787 \
    -v ${host_home}/myRiboSeq/R:${docker_home}/script \
    -v ${host_home}/myRiboSeq_res/result:${docker_home}/result \
    -v ${host_home}/myRiboSeq_res/data:${docker_home}/data \
    -v ${host_home}/myRiboSeq_res/ref:${docker_home}/ref \
    -w ${docker_home} \
    --name $1 \
    -e PASSWORD=1 \
    saori/rstudio:edger
