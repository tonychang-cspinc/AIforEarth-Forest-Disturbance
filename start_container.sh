#!/bin/bash
IMAGENAME=cspinc/fia:sample-building
DATADRIVE=/home/tony/datablob
SECRETS=/home/tony/azure_secrets.env
#-v $DATADRIVE:/datadrive \
docker run -it --rm --privileged \
    -p 8888:8888 \
    -v $('pwd'):/content \
    -w /content \
    --env-file $SECRETS \
    $IMAGENAME \
    /bin/bash
