#!/bin/bash
IMAGENAME=cspinc/fia:sample-building
DATADRIVE=/home/tony/datablob
SECRETS=/home/tony/azure_secrets.env
PORT=8080
#-v $DATADRIVE:/datadrive \
docker run -it --rm --privileged \
    -p $PORT:$PORT \
    -v $('pwd'):/content \
    -w /content \
    --env-file $SECRETS \
    $IMAGENAME \
    /bin/bash
