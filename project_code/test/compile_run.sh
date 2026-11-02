# variables
PROGRAM="bin/test_simple.exe"

# make dependencies
make -C ../face_addressed_mesh
make -C ../pv_coupling

# make and run program
make $PROGRAM
cat input.txt | $PROGRAM