#!/bin/bash
ROOT_DIR="$(pwd)"

for i in SFrameUncerts_$1/*/.
do
    echo "Working on $ROOT_DIR/$i/conf.txt"
    cd "$ROOT_DIR/$i/"
    sframe_batch.py config.xml
done
cd "$ROOT_DIR"
