make -C ../face_addressed_mesh
make -C ../pv_coupling
make -C ../turbulence
make bin/write_output_data.exe
make bin/compare_data.exe

bin/write_output_data.exe < input.txt
bin/compare_data.exe < input.txt