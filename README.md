# MPhil-Project
RANS k-epsilon turbulence model

Before running a simulation:
- ```cd``` to ```project_code/write_mesh``` and select your desired mesh, initial and boundary conditions based on your simulation case.
- run ```make bin/main.exe```
- run ```bin/main.exe < write_input/my_dir/my_file.txt```

To run the simulation tests:
- ```cd``` to ```project_code/test```. Edit the ```compile_run.sh``` script to select the preferred test form ```test/src```
- run ```bash compile_run.sh```
- afterwards, from the same location, run ```gnu_plot plotting/my_plot_script.gp``` to produce plots in the visualisations file

To run the turbulence simulation:
- ```cd``` to ```project_code```
- run ```bash compile_run.sh```
- afterwards, from the same location, run ```gnu_plot my_plot_script.gp``` to produce plots in the visualisations file.

To generate OpenFOAM results:
- ```cd``` to ```OpenFOAM/my_case```
- run ```blockMesh```, then ```pimpleFoam``` (requires OpenFOAM v2406 installed)

To compare OpenFOAM results with model results:
- ```cd``` to ```dynamicCode```
- run ```bash compile_run.sh```
