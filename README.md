# RMD_Digging

1 RMD_digging is developed by the MATLAB language. It is aimed to provide pre-processing and post-processing tools for the reactive molecular dynamics (ReaxFF) simulations performed on the LAMMPS  platform. Its functions involve formatting the reactive force field parameters, statistic anasysis of structures, trajectories and mechanisms and output of the visualization files. Besides, extra modelling by other softwares can be performed by preparing corresponding input files. The following softwares are used to realize these goals, including Materials studio, VMD and Gaussian. Let's make it more versatile, powerful and robust together!

2 A brief introduction to the RMD_Digging toolkit:
(1)ReaxFF_Gen toolkit is used to generate standard file with ReaxFF force field parameters from literatures.
(2)loglammps toolkit is mainly used to extract and process data with statistic average method.
(3)species toolkit can process the simulation generated species files. It sort out the formed products with time sequence and can further refine them by molecular weight, elements and so on. Moreover, the evolution of number-average molecular weight and weight-average molecular weight can be calculated.
(4)bonds_analysis toolkit can read and handle the bond order (BO) information, further mainly cope with the BO information with the following purpose: a)obtain the BO information of the specific molecule or fragment; b)give the chemical composition of all the molecules and fragments (molecular formula); c)rearrange molecular formula and trace the BO information of the interested species.
(5)lammpstrj2car_mdf_arc toolkit can output the car and mdf file, which can be be used for visualization by virture of Materials studio. Note: VMD software will be intended as the Fundamentals of Visualization. 
(6)chemi_mechanism toolkit can export the products in specific trajectories and analysis the reaction path/channel. More desired functions will be added in the furture.
(7)structure_analysis toolkit is aimed to analysis the complex molecular structure, such as benzene ring, phenolic hydroxylic group and other specific chemical bonds or groups. However, its function is incomplete. Much work will be carried out in the furture.

3 Please cite the following papers:
(1) Liu, Q.; Liu, S.; Lv, Y.; Hu, P.; Huang, Y.; Kong, M.; Li, G. Atomic-scale insight into the pyrolysis of polycarbonate by ReaxFF-based reactive molecular dynamics simulation. Fuel 2021, 287, 119484, DOI: https://doi.org/10.1016/j.fuel.2020.119484.
(2)  Liu, Q.; Huang, W.; Liu, B.; Wang, P.-C.; Chen, H.-B. Gamma Radiation Chemistry of Polydimethylsiloxane Foam in Radiation-Thermal Environments: Experiments and Simulations. ACS Appl. Mat. Interfaces 2021, DOI: https://doi.org/10.1021/acsami.1c10765.
