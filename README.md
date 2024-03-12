# SonoTraceLab
 
A raytracing-based simulation framework for modelling of active echo-based sensing systems. \
Paper: https://arxiv.org/pdf/2403.06847.pdf \
Youtube Video: https://www.youtube.com/watch?v=_h9p5dmUeeI \

## Publications
We kindly ask to cite our paper if you find this library useful:
```
@inproceedings{sonotracelab2024steckel,
  author={Jansen, Wouter and Steckel, Jan},
  title={SonoTraceLab - A Raytracing-Based Acoustic Modelling System for Simulating Echolocation Behavior of Bats}, 
  year={2024},
  volume={},
  number={},
  pages={},
  doi={}}
```

## Usage

### Dependencies
 - Matlab 2023 or higher
   - Signal Processing Toolbox
   - Deep Learning Toolbox
   - Parallel Computing Toolbox
   - Phased Array System Toolbox
   - Lidar Toolbox
- GPU with CUDA support for parallel processing

### Examples
Two examples are provided with ['exampleSimple.m'](exampleSimple.m) and ['exampleMovingTarget.m'](exampleMovingTarget.m).

### Mex files 
There are some accelerated functions that are implemented in C++/CUDA and compiled to a Mex file. 
Depending on your operating system these may need to be compiled from source. This relates to:
- [patch_normals_double.c](SourceCode/SupportCode/patchnormals_double.c), can be compiled with `mex patchnormals_double.c`
- [mollertrumbore.cu](SourceCode/SupportCode/mollertrumbore.cu), can be compiled with `mexcuda mollertrumbore.cu`

## License
This library is provided as is, will not be actively updated and comes without warranty or support.
Please contact a Cosys-Lab researcher to get more in depth information or if you wish to collaborate.
SonoTraceLab is open source under the MIT license, see the [LICENSE](LICENSE) file.

## Open-Source libraries included in this project
 - Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox by Paul Leopardi for the University of New South Wales [(link)](https://github.com/penguian/eq_sphere_partitions)
 - Curvature Estimation On Triangle Mesh by Itzik Ben Shabat [(link)](https://www.mathworks.com/matlabcentral/fileexchange/47134-curvature-estimationl-on-triangle-mesh)
 - Progress bar by HyunGwang Cho [(link)](https://www.mathworks.com/matlabcentral/fileexchange/121363-progress-bar-cli-gui-parfor?s_tid=srchtitle)
 - Save MAT files more quickly by Tim Holy [(link)](https://www.mathworks.com/matlabcentral/fileexchange/39721-save-mat-files-more-quickly?s_tid=srchtitle)
 - Patch Normals by Dirk-Jan Kroon [(link)](https://www.mathworks.com/matlabcentral/fileexchange/24330-patch-normals)
 - DataHash by Jan [(link)](https://www.mathworks.com/matlabcentral/fileexchange/31272-datahash?s_tid=srchtitle)
