# Computational Design of Custom-Fit PAP Masks

![teaser](/Resources/GitHub/Picture1.png?raw=true "teaser") 

This repository is the official implementation of the paper [Computational Design of Custom-Fit PAP Masks](https://sutd-cgl.github.io/supp/Publication/projects/2024-SMI-MaskDesign/index.html). 
The implementation takes an input human face already aligned with the generic mask, customizes the generic mask cushion and output the customized cushion surface and the mask interface. 

## Usage:
### Input the file path of the human face model:
```bash
cd ./bin/Release
./MaskDesign.exe [path_to_human_face_file]
```
For example, to generate one paper result from Fig.8 :
```bash
./MaskDesign.exe path_to_data_folder/face1.obj
```
To use another human face model, you should align it to the [generice mask](<font color=Red>TODO privide a link</font>) before input. The alignment approximates the human wearing the mask, and a certain amount of interpenetration is desired.

### Adjust design related parameters:
Adjust parameters in the [parameters.txt](<font color=Red>TODO privide a link</font>).
The comments start with # describs the meaning of parameters. The paramters that related to cushion comfort and air leakage measure, and trejectory curve initialization are important ones:
```bash
# cushion optimization parameters 
1e3 # evaluation comfort metric: average pressure 
0.2 # evaluation air leakage metric: force distribution
2e-4 # evaluation air leakage metric: area distribution

# (cushion initialization) trajectory curve adjustment weight:
1e3 # objective function: measuring curve to human face distance
2e-4 # cushion width
5e-6 # curvature; control curve smoothness
0.01 # symmetry
1e8 # alignment
120.0 # angle
```

