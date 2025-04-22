# SSA_Payload_Systems
Contains RF system Sizing and Error propagation. 

# Methodology Step-By-Step
1) Use v3 SAR model in SAR model + ScanSAR folder to define the SAR tradespace - Antenna Dimension array, Power limit, central frequency, resolution, defined by IOD Requirements.
2) Filter the Tradespace through the remaining IOD requirements
3) Compute the Range and Angular Uncertanity for remaining design points. Use worse case scenario - Largest beamwidth dimension. 
4) Select Design Points along the Pareto Front for the SAR tradespace
5) Redo with "ScanSAR" operational SAR mode 
6) Within the Debris_Covariance Folder, Input chosen ScanSAR parameters into IOD_v3 file, and input files
7) Run into line 214. Then Open MC_simulation script to propagte the error ellipsoid to define a revisit time. 

