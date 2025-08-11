## Figure Generation
All relevant figures are generated in plot_new.ipynb.

## Data Files
* **File Format**: <br>
The data are stored in .csv files with the first column corresponding to L or u (corresponding to files with
varying L and u respectively), the second column corresponding to $$\kappa$$, the third column corresponding
to G(A:B:C), and the fourth column corresponding to the (2,2)-reflected entropy.
* **Figure 7**: <br>
kappa_u_0_3_Chern_L_50_APBC.csv <br>
kappa_u_0_3_Chern_L_50_PBC.csv <br>
* **Figure 8**: <br>
8(a). kappa_L_10_70_Chern_u_1.65_APBC.csv, kappa_L_10_70_Chern_u_1.65_PBC.csv <br>
8(b). kappa_L_10_70_Chern_u_1.99_APBC.csv, kappa_L_10_70_Chern_u_1.99_PBC.csv <br>
8(c). kappa_L_10_70_Chern_u_2.01_APBC.csv, kappa_L_10_70_Chern_u_2.01_APBC.csv <br>

## Data Generation
* The data file kappa_u_0_3_Chern_L_50_APBC.csv is generated with Majorana_fermion_calculation_Chern_u_APBC.m <br>
* The data file kappa_u_0_3_Chern_L_50_PBC.csv is generated with Majorana_fermion_calculation_Chern_u_PBC.m <br>
* The data file kappa_L_10_70_Chern_u_1.99_APBC.csv is generated with Majorana_fermion_calculation_Chern_L_APBC.m <br>
* The data file kappa_L_10_70_Chern_u_1.99_PBC.csv is generated with Majorana_fermion_calculation_Chern_L_PBC.m <br>
* The data file kappa_L_10_70_Chern_u_1.65_APBC.csv can be generated with Majorana_fermion_calculation_Chern_L_APBC.m by changing the parameter in the fourth line of the code to u=1.65. <br>
* The data file kappa_L_10_70_Chern_u_1.65_PBC.csv can be generated with Majorana_fermion_calculation_Chern_L_PBC.m by changing the parameter in the fourth line of the code to u=1.65. <br>
* The data file kappa_L_10_70_Chern_u_2.01_APBC.csv can be generated with Majorana_fermion_calculation_Chern_L_APBC.m by changing the parameter in the fourth line of the code to u=2.01. <br>
* The data file kappa_L_10_70_Chern_u_2.01_PBC.csv can be generated with Majorana_fermion_calculation_Chern_L_PBC.m by changing the parameter in the fourth line of the code to u=2.01. <br>
* The codes Chern_insulator_correlation_matrix_2d.m and Chern_insulator_correlation_matrix_2d_APBC.m contain functions necessary to run the above-mentioned codes.
