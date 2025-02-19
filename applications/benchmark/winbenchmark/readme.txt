Windows-Benchmark
-----------------
This subdirectory contains project files for the benchmark if the
Microsoft Visual Studio 2003 IDE is used together with the Intel
Fortran Compiler. Four directories can be found here:

 winbenchcc2d/ - Contains the cc2d-benchmark project
 winbenchcc3d/ - Contains the cc3d-benchmark project
 winbenchpp2d/ - Contains the pp2d-benchmark project
 winbenchpp3d/ - Contains the pp3d-benchmark project
 
All projects are set up the following way:

1.) Compile the source files in the application/xxxd/src-directory 
    (with xxxd standing for pp2d/3d, cc2d/3d) together with the geometry-
    specific files INDATxD.F/PARQxD.F from the 
    application/benchmark/-subdirectory to an application - either
    cc2d/3d or pp2d/3d.
    
2.) The dependent libraries are set up as sub-projects to the workspace.
    The application is made dependent of the added sub-projects
    (right-click on the application, project-dependencies)
    In the case of the 3D-projects, a more specific dependency is added
    to make sure that the FEAT3D-library has higher priority than the
    FEAT2D-library:

    a) Main project is dependent on:        BLAS, FEAT3D
    b) Sub-project FEAT3D is dependent on:  FEAT2D
    
3.) The project options are modified the following way:
      Check Array and String Bounds = No
      Calling Convention            = C, Reference
      Working Path                  = $(ProjectDir)\..\..
      Additional Include-Path       = $(ProjectDir)\..\..
    The calling convention is changed in advance for future
    enhancements. Although this is not necessary for the current
    FeatFlow version, it's necessary if C-libraries (e.g. UMFPACK-4)
    is later added.

The benchmark in the Visual Studio IDE cannot be executed automatically.
The user has to set up the benchmark configuration in the IDE manually.

For every application (winbenchcc2d/3d, winbenchpp2d/3d) two benchmark-
specific .DAT-files exist in the application/benchmark/#data-
subdirectory: 

  xxxd.dat_upw - A test using upwind stabilization
  xxxd.dat_sdg - A test using streamline-diffusion stabilization
  
To start a specific benchmark, the user has to do the following:

1.) Open and compile the benchmark project file (winbenchxxxd)
2.) Overwrite the .DAT-file #data/xxxd.DAT with the appropriate .DAT-
    file of the test that should be performed (xxxd.dat_upw or
    xxxd.dat_sdg, respectively)
3.) Start the application.

FINAL REMARK:
-------------
The benchmark using the Visual Studio IDE is completely experimental
without any guarantee that it works. It's only provided with this
package for users who "want to give it a try".
Users not familiar with setting up this by hand in the IDE 
should use the automatic benchmark, which is executed by the 
command "make bench" using the Cygwin Unix emulation package instead
of MS Visual Studio. This will directly write out a .log-file with
the results of all benchmarks into the main installation directory.
