#!/bin/bash

# Allowing for the image to be named
if [ -n "$1" ]; then
    DOCKER_NAME="$1"
else
    DOCKER_NAME="garf_sim:latest"
fi

# Allowing for the image to be named
if [ -n "$2" ]; then
    IMAGE_NAME="$2"
else
    IMAGE_NAME="garf_sim.sif"
fi

apptainer build $IMAGE_NAME docker-daemon:$DOCKER_NAME