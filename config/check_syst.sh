ROOT_DIR="$(pwd)"

for i in SFrameUncerts_2017/*/.
do
    echo "Working on $ROOT_DIR/$i/conf.txt"
    cd "$ROOT_DIR/$i/"
    #sframe_batch.py conf.xml
done
cd "$ROOT_DIR"
