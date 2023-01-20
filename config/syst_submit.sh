#!/bin/bash
ROOT_DIR="$(pwd)"

for i in SFrameUncerts_$1/*
do
    "Working on $ROOT_DIR/$i/config.xml"
    cd "$ROOT_DIR/$i/"
    sframe_batch.py -s config.xml
done
cd "$ROOT_DIR"
