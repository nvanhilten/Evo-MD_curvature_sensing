#!/bin/bash


if [ $(($OMPI_COMM_WORLD_LOCAL_RANK%2)) -eq 0 ]
then
  echo "even"
  echo $OMPI_COMM_WORLD_LOCAL_RANK
  export CUDA_VISIBLE_DEVICES=0
  echo $CUDA_VISIBLE_DEVICES
else
  echo "odd"
  echo $OMPI_COMM_WORLD_LOCAL_RANK
  export CUDA_VISIBLE_DEVICES=1
  echo $CUDA_VISIBLE_DEVICES
fi

$@