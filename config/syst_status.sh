ROOT_DIR="$(pwd)"

for i in SFrameUncerts_2016/*/.
do
    echo "Working on $ROOT_DIR/$i/conf.txt"
    cd "$ROOT_DIR/$i/"
    sframe_batch.py conf.xml
done
cd "$ROOT_DIR"

# for i in SFrameUncerts_2017/*/.
# do
#     echo "Working on $ROOT_DIR/$i/conf.txt"
#     cd "$ROOT_DIR/$i/"
#     sframe_batch.py conf.xml
# done
# cd "$ROOT_DIR"

# for i in SFrameUncerts_2018/*/.
# do
#     echo "Working on $ROOT_DIR/$i/conf.txt"
#     cd "$ROOT_DIR/$i/"
#     sframe_batch.py conf.xml
# done
# cd "$ROOT_DIR"
