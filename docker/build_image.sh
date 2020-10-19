#!/bin/bash
TAG=sample-building
IMAGENAME=cspinc/fia:$TAG
docker build --rm -t $IMAGENAME .
