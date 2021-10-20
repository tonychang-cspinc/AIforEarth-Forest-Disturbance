#!/bin/bash

PORT=8888
IP=0.0.0.0
sh ./blobfuse_container.sh 
jupyter notebook --port $PORT --ip $IP --no-browser --allow-root
