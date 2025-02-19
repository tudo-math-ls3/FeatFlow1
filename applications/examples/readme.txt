
  FEATFLOW 1.3 - Examples
  -----------------------
  
This document describes how to compute the examples from Chapter 3
of the Featflow manual. Due to the change in the directory structure,
the way of executing the examples has slightly changed.
The following discussions assume that you installed Featflow using
a "make build" command in the main directory of the Featflow repository
as described in /readme.txt; this command has prebuild all
applications.

The following descriptions are restricted to the Unix/Linux command
line for simplicity.


 3.2. The 2D example
---------------------

 3.2.1 - The stationary example
--------------------------------
Our intention is to demonstrate the computation of a stationary
"Flow-around-cylinder" example using the cc2d solver. This solver
is designed for stationary and lower Reynolds-number nonstationary
problems and can be found in the directory

  applications/cc2d

of the featflow repository. Please enter this directory by typing

  cd applications/cc2d
  
in the main Featflow directory.

The basic domain/mesh files for the flow-around-cylinder example
can be found in the following subdirectories:

  applications/cc2d/#adc
    c2d0.tri - The coarse mesh, referred to as "level 1".

  applications/cc2d/#pre
    c2d.prm  - Definition of the domain.

  Local refinement of the coarse mesh
  -----------------------------------
  The basic coarse mesh c2d0.tri was generated using a very old tool
  "Omega2D", which is not maintained anymore. Today, we feature our
  mesh design tool "grid3d" which can be downloaded from the Featflow
  homepage and be used to modify c2d.prm/c2d0.tri. However, this tool
  does also not provide local refinement of a mesh which we want 
  to use in this example.
  
  The computational coarse mesh "c2d2.tri" for this test is generated from the 
  basic mesh "c2d0.tri" in two steps using the external refinement utility
  trigen2d, which can be found in applications/trigen2d. If you are
  in "applications/cc2d", the refinement is applied as follows:
  
  a) Refinement step 1. Enter the commands
  
    cp -f ../examples/trigen2d.dat_0 #data/trigen2d.dat
    ../trigen2d/trigen2d-*
    
  which invokes the trigen2d utility to apply a local refinement
  to the mesh c2d0.tri (with ITYPEL=1, which means that marked elements
  are refined intu 9 subelements). The result is saved to
    
    #adc/c2d1.tri
    
  b) Refinement step 2. In this step, a local refinement around 
  the circle should be applied. Enter the commands
  
    cp -f ../examples/trigen2d.dat_1 #data/trigen2d.dat
    ../trigen2d/trigen2d-*
    
  which again invokes trigen2d to refine around the circle
  (using ITYPEL=3). The result is saved to
  
    #adc/c2d2.tri
    
Next, we have to prepare the Navier-Stokes solver to be applied
to this mesh with the settings we want. We have to specify a data
file "cc2d.dat" which defines the settings of the solver and a
sourcecode file "indat2d.f" which specifies the boundary conditions
and has to be linked into the application. All these files can be
found in applications/examples. Enter the following commands to
use them:

  cp -f ../examples/cc2d.dat #data/cc2d.dat
  cp -f ../examples/indat2d.f_stat ./indat2d.f
  make
  
This overwrites previous files of the cc2d example and invokes the
compiler to compile the boundary conditions into the application.
The Navier-Stokes solver can then be started using

  ./cc2d-*

According to #data/cc2d.dat, a solution is calculated at level "NLMAX=4".
Files for graphical output are saved to the #avs and #gmv directories.
The corresponding solution vector is saved (unformattedly) in the file

  #data/#DX4_stat
  
The chosen discretization scheme for the convective parts is the
Streamline-diffusion scheme, and the stopping criterions are 1E-3
for the maximum relative changes. A further result is the LOG file
#data/cc2d.stat. A detailed explaination of the output can be found
in Section 3.2 of the Featflow manual.


 3.2.2 - The nonstationary example
-----------------------------------
Next, we perform a nonstationary calculation using the nonstationary
solver pp2d which is designed for higher Reynolds number calculations.
In particular, we apply the solver to the flow-around-cylinder example
at Reynolds number Re=100. The solver can be found in the directory

  applications/pp2d
  
of the Featflow repository. Please enter this directory using the
command

  cd applications/pp2d
  
in the main Featflow directory.

Again, we have to prepare the Navier-Stokes solver to be applied
to this mesh with the settings we want. This time, we have to specify 
a data file "pp2d.dat" and, as before, a sourcecode file "indat2d.f" 
for the boundary conditions which has to be linked into the application. 
Again, all these files can be found in applications/examples. 
Enter the following commands to use them:

  cp -f ../examples/pp2d.dat #data/pp2d.dat
  cp -f ../examples/indat2d.f_non ./indat2d.f
  make
  
This overwrites previous files of the pp2d application and invokes the
compiler to compile the boundary conditions into the application.
The differences between indat2d.f_non in comparison to indat2d.f_stat
are as follows:

 * The inflow profile at x=0 is again parabolic, but with a maximum
   velocity Umax=1.5.
 * The parameter UMEAN is changed for the appropriate calculation
   of drag and lift.

In this computation, we want to use the same domain/mesh which we used
in cc2d before, so we copy it into the #adc directory:

  cp -f ../cc2d/#adc/c2d2.tri #adc/
  cp -f ../cc2d/#pre/c2d.prm #pre/
  
We want to use the solution of the stationary solver as initial
solution at time t=0 for the nonstationary solver. Copy the solution
computed with the cc2d example above to the #data directory of
pp2d as follows:

  cp -f ../cc2d/#data/#DX4_stat #data/

Now we are ready to invoke the nonstationary solver using the
command

  ./pp2d-*
  
According to #data/pp2d.dat, a solution is calculated at level "NLMAX=4"
using #data/DX4_stat as initial solution. The underlying discretisation
scheme for the convective part is here the Upwind scheme, and a
Fractional-Step timestepping scheme is applied for the time discretisation.
Timestep length are controlled by a fully adaptive timestepping scheme
until TIMEMX=4.0. 

Files for graphical output are saved to the #avs and #gmv directories
every "1 second". Point value output data (drag/lift coefficients) 
is saved to the directory #points and can be visualised using GNUPLOT.
A further result is the LOG file #data/pp2d.non which is rather similar 
to the previous cc2d.stat. A detailed explaination of the output can be found
in Section 3.2 of the Featflow manual.


 3.3. The 3D example
---------------------

 3.3.1 - The stationary example
--------------------------------
In this example, we want to apply the Navier-Stokes solver cc3d to calculate
a stationary 3D example. This solver was designed for low Reynolds number
calculations, and in particular, stationary simulations. The solver can be
found in the directory

  applications/cc3d
  
of the Featflow repository. Please enter this directory using the
command

  cd applications/cc3d
  
in the main Featflow directory.

In a first step, we want to extend our mesh c2d2.tri from the 2D example
to a 3D counterpart. This is done in two steps:

a) We copy the domain and the mesh from the 2D example with the following 
   commands:

     cp -f ../cc2d/#adc/c2d2.tri #adc/
     cp -rf ../cc2d/#pre .

b) We invke the application "tr2to3" to extend it to a 3D mesh
   using the following commands:

     cp -f ../examples/tr2to3.dat #data/
     ../tr2to3/tr2to3-*
     
   the result is a file "#adc/c3d2.tri" with a 3D mesh.

The next step is to prepare the data file "#data/cc3d.dat" and to
impose the boundary conditions with an appropriate sourcecode file
"indat3d.f". We again take the templates from the examples/
subdirectory and compile the solver.

  cp -f ../examples/cc3d.dat #data/
  cp -f ../examples/indat3d.f_stat ./indat3d.f
  make
  
After the compilation process is finished, the solver can be
invoked with the command

  ./cc3d-*

According to the data file "#data/cc3d.dat" it computes a solution
at level NLMAX=3 and saves the calculated solution to the file

  #data/#DX3_stat
  
similar to the stationary case. A protocol file is saved to
#data/cc3d.stat which is almost identical to the previous cc2d.dat.
Postprocessing output is written to the #avs and #gmv directories.
For the exact definitions of all parameters defined in the cc3d.dat
and indat3d.f, we refer to the Featflow manual.


 3.3.2 - The nonstationary example
-----------------------------------
Similar to the 2D case, we now want to extend the stationary simulation
to a nonstationary one, using the nonstationary solver pp3d.
Similar to pp2d, this is designed for higher Reynolds number calculations.
In particular, we apply the solver to the flow-around-cylinder example
at Reynolds number Re=100. The solver can be found in the directory

  applications/pp3d
  
of the Featflow repository. Please enter this directory using the
command

  cd applications/pp3d
  
in the main Featflow directory.

In a first step, we copy the computational mesh used in cc3d to pp3d.
Furthermore, we copy the solution vector computed with cc3d to pp3d
as well as we want to use it as initial solution at time t=0:

  cp -f ../cc3d/#adc/c3d2.tri #adc
  cp -f ../cc3d/#data/#DX3_stat #data
  
IN a second step, we prepare the data file "pp3d.dat" and the
boundary condition source file "indat3d.f" of pp3d to use the
correct settings for this computation. We take the templates
from the examples/ subdirectory and compile the boundary conditions
into the application.

  cp -f ../examples/pp3d.dat #data/
  cp -f ../examples/indat3d.f_non ./indat3d.f
  make

For description of the parameters in these files and the differences
between indat3d.f_non and indat3d.f_stat, we refer to Section 3.3
of the Featflow manual.

After the compilation process is finished, the solver can be invoked
using the command

  ./pp3d-*
 
It generates graphical output in #avs/ and #gmv/ in "1 second" steps
and point value output data (drag/lift coefficients) in the directory
#points; this data can be visualised using GNUPLOT. A protocol file
is saved to #data/pp3d.non. A detailed explaination of the output can 
be found in Section 3.3 of the Featflow manual.
