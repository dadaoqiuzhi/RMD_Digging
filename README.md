# RMD_Digging: A Toolkit to Remove the Barriers in the Way of Scientific Calculations and Analysis  
# <span style="color: red"> Notice: An upgrade version RMD_Digging_v3.1_dev has been released! This version is more fast and robust, with more practical functions and little bugs. It is also aimed to alleviate the memory consumption in many tasks. Both versions with respective English and Chinese language will be provided. Enjoy it. </span>
# Welcome to join our QQ group: 948210961. 

:rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket:
## Development Target

RMD_Digging is developed by the MATLAB language. It is aimed to provide pre-processing and post-processing tools for the reactive molecular dynamics (ReaxFF) simulations performed on the LAMMPS platform. Its functions involve formatting the reactive force field parameters, statistic anasysis of structures, trajectories and mechanisms and output of the visualization files. Besides, extra modelling by other softwares can be performed by preparing corresponding input files. The following softwares are used to realize these goals, including Materials studio, VMD and Gaussian. Let's make it more versatile, powerful and robust together!

## A Brief Introduction to the RMD_Digging Toolkit

(1)**ReaxFF_Gen toolkit** is used to generate standard file with ReaxFF force field parameters from literatures.<br>
(2)**loglammps toolkit** is mainly used to extract and process data with statistic average method.<br>
(3)**species toolkit** can process the simulation generated species files. It sorts out the formed products with time sequence and can further refine them by molecular weight, elements and so on. Moreover, the evolution of number-average molecular weight and weight-average molecular weight can be calculated.<br>
(4)**bonds_analysis toolkit** can read and handle the bond order (BO) information, further mainly cope with the BO information with the following purpose: a)obtain the BO information of the specific molecule or fragment; b)give the chemical composition of all the molecules and fragments (molecular formula); c)rearrange molecular formula and trace the BO information of the interested species.<br>
(5)**lammpstrj2xyz_arc_pdb toolkit** can output the *.xyz, *.arc and *.pdb file, which can be be used for visualization by virture of Materials studio and VMD software, both in the form of static image and dynamic trajectory. <br>
(6)**chemi_mechanism** toolkit can export the products in specific trajectories and analysis the reaction path/channel. More desired functions like orientation analysis will be added in the furture.<br>
(7)structure_analysis toolkit is aimed to analysis the complex molecular structure, such as benzene ring, phenolic hydroxylic group and other specific chemical bonds or groups. However, its function is incomplete. Much work will be carried out in the furture.

## Please Cite the Following Papers and RMD_Digging ToolKit Accordingly:
(1) **Liu, Q.**, RMD_Digging ToolKit, 2020, https://github.com/dadaoqiuzhi/RMD_Digging, accessed date.<br>
(2) **Liu, Q.**; Liu, S.; Lv, Y.; Hu, P.; Huang, Y.; Kong, M.; Li, G. Atomic-scale insight into the pyrolysis of polycarbonate by ReaxFF-based reactive molecular dynamics simulation. Fuel 2021, 287, 119484, DOI: https://doi.org/10.1016/j.fuel.2020.119484.<br>
(3) **Liu, Q.**; Huang, W.; Liu, B.; Wang, P.-C.; Chen, H.-B. Gamma Radiation Chemistry of Polydimethylsiloxane Foam in Radiation-Thermal Environments: Experiments and Simulations. ACS Appl. Mat. Interfaces 2021, 13 (34), 41287-41302, DOI: https://doi.org/10.1021/acsami.1c10765.<br>
(4) C. Li, **Q. Liu**, W. Gong, Z. Zhou, Z. Yao, X. Meng, Study on the atomic scale of thermal and thermo-oxidative degradation of polylactic acid via reactive molecular dynamics simulation, Thermochim. Acta 709 (2022) 179144, DOI: https://doi.org/10.1016/j.tca.2021.179144<br>
(5) **Liu, Q.**; Huang, W.; Chen, H. Paving the Way to Simulate and Understand the Radiochemical Damage of Porous Polymer Foam. ACS Materials Letters 2023, 2174-2188. DOI: 10.1021/acsmaterialslett.3c00307.<br>

## Disclaimer:
***Using the software and the scripts in this repository is at your own risk***. The software is freeware and it come without any (as in nothing, nada, niente, rien) WARRANTY. Therefore, Neither me nor the Institute of Nuclear Physics and Chemistry is responsible for the loss of data, time, bad results or anything else deriving by the use of my software or this repository. We do not make any warranty, express or implied, or assumes any legal liability or responsibility for the accuracy, completeness, or usefulness of any software, scipts, information, apparatus, product, disclosed process, or represents that its use would not infringe privately owned rights. Reference herein to any specific commercial software, product, process, or service by trade name, trademark, manufacturer, or otherwise, does not necessarily constitute or imply its endorsement, recommendation, or favoring by the People's Republic of China. The views and opinions of authors expressed herein do not necessarily state or reflect those of the People's Republic of China or the Institute of Nuclear Physics and Chemistry, and shall not be used for advertising or product endorsement purposes.   

:rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket::rocket:
