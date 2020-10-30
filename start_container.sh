#!/bin/bash
IMAGENAME=cspinc/fia:sample-building
DATADRIVE=/home/cspadmin/datablob
SECRETS=/home/cspadmin/azure_secrets.env
PORT=8765
#-v $DATADRIVE:/datadrive \
docker run -it --rm --privileged \
    -p $PORT:$PORT \
    -v $('pwd'):/content \
    -w /content \
    --env-file $SECRETS \
    $IMAGENAME \
    /bin/bash
