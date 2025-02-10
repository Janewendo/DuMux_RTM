February 10 Update 2  
<br> 

I simplified the system to a 2-phase model with 4 primary components and 5 secondary components, as detailed in the Du_2p folder. I attempted to modify the code to allow for the input of free ion concentrations and calculate the total concentrations in the mass balance equation. However, this approach has not been successful in running properly.  
<br>  


February 10 Update 1  
<br> 

I simplified the system to a 1-phase, 4 primary components, and 5 secondary components system, see Du_1p folder. I modified the code to be able to input free ion concentration, and calculate the total concentration in the mass balance equation.  
<br>  

The base code are from \dumux\test\porousmediumflow\1pnc\1p2c\isothermal. 
<br> 

The modification are recorded. You may check the 'update **.file' commit on Feb 10 to see changes made.  
<br>  


   

January 23 Update  
<br>  

Following the meeting on January 20, I made changes to the following files:

•	du2024/dumux/porousmediumflow/2picp/volumevariables.hh

•	du2024/dumux/material/chemistry/biogeochemistry/leocarbonicacid.hh

•	du2024/appl/icp/eicp/leocolumnproblem.hh  
<br>  

First, in volumevariables.hh:
   
o	I added lines 481 to 487. These lines are intended to input the total concentration into the mass balance equation.

o	The total concentration is calculated at line 470 using the calculateEquilibriumChemistry function, which is defined in leocarbonicacid.hh. The specific calculation of the total concentration takes place between lines 223 and 243.


However, the simulation is not running successfully, and the log info is saved in du2024/logfile.txt.


Then I made some adjustment in Boundary Condition:
   
o	I suspect the issue could be related to the boundary condition settings for transport in leocolumnproblem.hh. Specifically, in lines 585 to 587, where the volume variable is used as the concentration. Since the volume variable now refers to the total concentration, I attempted to modify the code between lines 589 and 591 to recalculate the concentration for the free ion. 

Unfortunately, this did not resolve the simulation issue.  
<br>  

Question:

•	I’m wondering if other parts of the code need to be adjusted to account for the change to the total concentration. Any guidance on this would be appreciated.
