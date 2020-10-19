#!/bin/bash
IMAGENAME=cspinc/fia:sample-building
DATADRIVE=/E
docker run -it --rm \
    -p 8888:8888 \
    -v $('pwd'):/content \
    -v $DATADRIVE:/datadrive \
    -v $(pwd)/authentication/credentials:/root/.config/earthengine/credentials \
    -v $(pwd)/authentication/.boto:/root/.boto \
    -w /content \
    $IMAGENAME \
    /bin/bash
#	 /bin/bash -c "bash authentications.sh;/bin/bash" \
