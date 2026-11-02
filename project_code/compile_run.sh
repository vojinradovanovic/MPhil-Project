make -C face_addressed_mesh
make -C pv_coupling
make -C turbulence
make bin/main.exe
cat input/input.txt | bin/main.exe