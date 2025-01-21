#!/bin/bash

# Allowing for the image to be named
if [ -n "$1" ]; then
    IMAGE_NAME="$1"
else
    IMAGE_NAME="garf_sim:latest"
fi

docker build -t $IMAGE_NAME . 