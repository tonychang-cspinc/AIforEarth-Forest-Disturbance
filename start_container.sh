#!/bin/bash
USER=tony
IMAGENAME=cspinc/fia:sample-building
DATADRIVE=/home/$USER/datablob
SECRETS=/home/$USER/azure_secrets.env
PORT=8888
docker run -it --rm --privileged \
    -p $PORT:$PORT \
    -v $('pwd'):/content \
    -w /content \
    -v /datadrive:/datadrive \
    --env-file $SECRETS \
    $IMAGENAME \
    /bin/bash
