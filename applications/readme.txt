________________________________________________________________________
                    FEATFLOW 1.3 - applications/
________________________________________________________________________

1. General structure
2. The Featflow benchmark
3. Applications in the Microsoft Visual Studio IDE
4. Postprocessing data in Windows


                          1. General structure
------------------------------------------------------------------------
This directory contains the main applications for the FEATFLOW
package.

 -  The "cc2d/cc3d" and "pp2d/pp3d" directories contain the main solvers
    for the flow simulation. Each directory contains a couple of .f-
    files that are designed to be changed by the user according
    to the problem to be solved. 
    The main code of these appications is capsuled in the "src"
    subdirectory. The code there is designed to be a "black-box"-solver,
    so normally the user doesn't have to change that. The files
    to be changed by the user are:
    
    indatXd.f  - specification of inflow, outflow,...
    parqXd.f   - eventually specifying special boundary geometries
    YYXD.inc   - include file with definition of workspace;
                 must be changed if reserved memory is not enough!
    
    For changing the code to a specified situation, the corresponding
    solver subdirectory (cc2d/3d, pp2d/3d) should be duplicated
    inside of this applications/-subdirectory and then the .f-files 
    in that directory should be changed to user's needs.
    
    Remark: The file "parpre.f" in the src/-subdirectory is a prototype
     of a parqXd.f-file; if the parqXd.f gets destroyed, it can 
     restored by overwriting it with parpre.f.

    IMPORTANT: Avoid changing the code of ppXd/ccXd in the original 
     directory. The Featflow benchmark depends on the code in these 
     directories (see below), so by changing it you can't be sure
     anymore to obtain benchmark results that are comparable to those
     published on the Featflow homepage (http://www.featflow.de).
   
 -  "trigen2d" contains a ".tri"-file generator for 2D.
    "trigen3d" contains a ".tri"-file generator for 3D.
    "tr2to3" contains a ".tri"-file generator that generates a 3D
             .tri-file by extruding a 2D .tri-file.
    
 -  The "benchmark" subdirectory contains Featflow-benchmark specific
    files. For a detailed description about the benchmark, see below.
    
 -  The "examples" subdirectory contains files for the example
    computations from Chapter 3 of the Featflow manual. See the file
    "examples/readme.txt" for details how to apply the examples
    in the new makefile structure.
    
 -  The "bouss", "bouss_powerlaw", "cc2d_movbc",... directories
    contain experimental code for Boussinesq-model and cc2d with moving
    boundaries. Use these applications with care!


                      2. The Featflow benchmark
------------------------------------------------------------------------
The "benchmark" subdirectory contains Featflow-benchmark specific
files. The benchmark itself is performed the following way:

a) Compile the files in the xxxd/src-directory of each application
   (cc2d/3d and pp2d/3d) together with the benchmark-specific files
   for the geometry inside of the benchmark/-directory and in a couple
   of alternative source files in the altsrc/-directory (see below)
b) Replace the .dat-file inside of the benchmark/#data-directory
   by the appropriate benchmark-specific data file
   (can also be found in benchmark/#data)
c) Execute the application
d) The application generates a .log-file in the #data/-subdirectory
   (xxyd.sdg/xxyd.upw). All thee log files are concatenated together
   to one .log-file in the installation directory of Featflow.
   
So avoid to change the code in the "src/" subdirectories of the
applications since the benchmark is depending on this!

Each application is executed two times with a different
bechmark configuration file. The geometry forming the basis of
the benchmark is the DFG-benchmark geometry ("flow around
cylinder"). This benchmark is executed
 - in 2D as well as in 3D
 - stationary (ccXd) and instationary (ppXd)
 - with streamline diffusion (SDG) as well as with upwind (UPW)
   stabilization of the resulting mathematical system
which gives an overall sum of 8 tests.

To start the benchmark you have the following possibilities
after Featflow is installed:
a) Run a top-level "make bench"
b) In the benchmark/subdirectory run "make bench"

The results of the benchmark are saved in the installation directory
of Featflow. The filename contains the machine architecture
string for easier identification. The filename of the file
containing the results is named as:

  "bench-xxxxxxxxxxxxxxx.log"
  
with "xxxxxxxxxxxxxxx" the identification string of your machine.

If you want to compare the results of your machine to the
"reference" results that are published of the Featflow homepage
(http://www.featflow.de) you have two possibilities:

a) Directly look into the .log file that is generated. The file has
   the following structure:
   
   ...
   Processing benchmark: cc2d with streamline diffusion (cc2d-SDG)
   ------------------------------------------------------------------
            INPUT DATA
   ...
   -------------------------------------------------------------------
      total time :   40.1908912
      appr. time :   40.1848923
      grid  time :   1.71673905
      post  time :   0.0899887085
      lin.  time :   8.88264823
      -> mavec time :   1.4077853
      -> konv. time :   6.76497436
      -> bdry  time :   0.00499820709
      -> LC    time :   0.70489037
      mg    time :   29.4955163
      #substeps  :  1
      #nonlinear :  11
      #mg        :  27
   -------------------------------------------------------------------
   ...

   The first line describes the type of simulation that is executed.
   At the end of each test the time is logged that was necessary
   for the computation. These timings correspond to those published
   on the Featflow homepage.
   
b) Use an extraction script to extract the time-data from the
   log file and print them on screen. This script was formally used
   only internally, but we decided to publish it for easier access
   to the data. The script itself is an AWK-script, so you need
   an installed AWG or GAWK in order to use it.
   
   After the Featflow benchmark is completed, switch to the
   installation directory of Featflow and type:
   
     awk -f bin/extract.awk bench-xxxxxxxxxxxxxxx.log
   
   The extract.awk script will then extract the data from the file
   and print it. 
   
The directory benchmark/results contains a couple of .log files we
computed of our machines that were up-to-date when this Featflow-release
was created. You can use them for comparison.


Remark: Alternative source files
--------------------------------
The Featflow benchmark differs from the standard ccXd/ppXd in some
source files:
 - the geometry/problem specific files (can be found in the benchmark/-
   subdirectory), that adapt the code to a specific problem
 - some "alternative" source files (can be found in the altsrc/-
   subdirectory)
The reason for shipping "alternative" source files for the benchmark
in addition to geometry-specific files is the following: 
The benchmark was created in the early days of Featflow. All these 
alternative source files contain bugs from old versions. 
One example is cc2d: The "old" version that is used for the
benchmark uses an improper defect correction, which makes the solver
having problems when solving for higher accuracy than 10^-6.
The Featflow benchmark solves the system only to an accuracy of 10^-5,
so the problem does not occur here. 
Unfortunately with the "correct" source code the results of the
Featflow benchmark are not directly comparable to those published on
our website: The results are mostly about 10% worse than with the
"old" source. So we decided to ship the old versions of the source
together with the benchmark to allow calculating results that can be
compared to all our previous results.


           3. Applications in the Microsoft Visual Studio IDE
------------------------------------------------------------------------
The standard applications cc2d/3d, pp2d/3d as well as the benchmark
directory also contain projects for the Visual Studio 2003 IDE 
together with the Intel Fortran Compiler for Windows. 
The Visual Studio project can be found in the winxxxx/-subdirectory
of those applications. The applications are assembled from the source 
code in the following way:

 a) All source files are added to the project
 b) The dependent libraries are added as sub-projects to the workspace.
    The application is made dependent of the added sub-projects
    (right-click on the application, project-dependencies)
 c) The project options are modified the following way:
      Check Array and String Bounds = No
      Calling Convention            = C, Reference
      Working Path                  = $(ProjectDir)\..
      Additional Include-Path       = $(ProjectDir)\..
    The calling convention is changed in advance for future
    enhancements. Although this is not necessary for the current
    Featflow version, it's necessary if C-libraries (e.g. UMFPACK-4)
    is later added.

The 3D-versions of cc/pp are a little bit tricky to set up, because
some names of subroutines are the same in the FEAT2D- as well as in
the FEAT3D-library. So the Visual Studio linker has to follow a specific
order how to link the libraries: The FEAT3D-routines must have higher
priority than the FEAT2D-routines. This can be archived by the following
project dependencies:

 a) Main project is dependent on:        BLAS, FEAT3D
 b) Sub-project FEAT3D is dependent on:  FEAT2D

Because the main project is no more directly dependent on FEAT2D,
the linker includes FEAT3D-routines first.

How to setup the benchmark in Visual Studio IDE is described in the
readme.txt-file in the benchmark/winbenchmark/-subdirectory.


                     4. Postprocessing in Windows
------------------------------------------------------------------------
The preprocessing of data (creation of the grid, ...) is normally done
by DeViSoR. As this program is written in Java, it can be used on
every machine. Unfortunately how to postprocess data is not so obvious
in a Windows environment, so we give a short idea here how to do that.

For postprocessing, the program "GMV" can be used. Featflow
applications normally write out a .GMV-file for postprocessing
into the #gmv/ subdirectory of the application. The GMV application
itself to postprocess this data is available on the GMV-homepage:

  http://www-xdiv.lanl.gov/XCM/gmv/GMVHome.html
  
Fortunately there is a Cygwin port of this program available on this
homepage, too. To use it,
 - download and extract it into a directory, rename the file "cygwinogl"
   to "gmv.exe"
 - make sure, you installed the Cygwin Unix emulator package together
   with the X11 and OpenGL extensions
 - start the X11 server of Cygwin by running the "startxwin.bat"
   batch file which can be found in the "usr\X11R6\bin" subdirectory
   of Cygwin (i.e. in the standard installation: 
   "C:\cygwin\usr\X11R6\bin\startxwin.bat")
 - This file launches a graphical bash-shell. Switch to the directory
   you extracted GMV to and start the "gmv.exe" in that directory.
A GMV-window will pop up where you can open GMV-files created by
Featflow.

A final word to MPEG-Movies: For Linux/Unix environments there is a
batch script "gmvmpeg" provided in the utitilies/-subdirectory, you
can use to create an MPEG move from a couple of GMV-files. This
uses an MPEG creation utility which is not present in the Cygwin
emulator package. To create an MPEG from GMV-data in Windows, we found
another way that works quite well, but you need some experience for
it to successfully use it. We can't give a detailed description here,
but we give the rough idea for those users who really want to try -
so don't ask us for support if it does not work :)

At first you need the following additional tools:
 - the "ImageMagick" converter suite; this is an additional package
   of Cygwin, you can install it with the Setup utility
 - the "AVISynth" program/driver package (http://www.avisynth.org/)
 - the "TMPEGEnc" MPEG-encoder program; a freeware version can
   be downloaded from http://www.tmpgenc.net/
   
To convert a couple of MPEG-files into a MPEG movie, you can then try
the following procedure:

a) Rename the GMV-files to allow GMV to read the files consecutively
   and automatically with the "File/Auto Read Simulation" feature.
   Normally the files must match the following way:
   "xxx.gmv.0000", "xxx.gmv.0001", "xxx.gmv.0002",...

b) Open the first GMV-file in GMV. Make the appropriate changes in the
   settings of GMV to create a good-looking view.
   
c) Open the "File/Auto Read - Same Simulation" window. Specify the first
   and last filename-number that should be read. Click on 
   "Auto Snapshots" for activating the automatic generating of snapshots
   of the files that are readed. Specify a "Filename Prefix" for the
   snapshots. Switch the auto-reading to "On" and close the window for
   specifying the filename again.
   
d) Click on "Start" to start auto-reading. GMV will read all your
   GMV-files consecutively and make a snapshot of each data-file into a
   RGB-coded file. All these output files can be found in the directory
   you extracted GMV to.
   
In the following we assume that you specified the filename "output"
in the auto-read feature; in this case GMV named your output files as
"output0000", "output0001", "output0002",...
   
e) Start a Cygwin shell, switch to the directory with the output files
   of GMV and use the "mogrify" command to convert all files into
   .BMP Windows Bitmap files, e.g. in our case:
   
        "mogrify -format bmp output*"

   This will generate "output0000.bmp", "output0001.bmp",
   "output0002.bmp", ...
   
f) Use a text editor to create an .AVS AVISynth-Script that reads
   these BMP-file and converts them on-the-fly to a movie. This script
   file should be placed into the same directory that contains the BMP-
   files. In our case the content of this AVS-Textfile can be e.g.
   something like:
   
       ImageSource("output%04d.bmp",start=0,end=25,fps=15)
       Crop(10,110,-20,-10)

   If you use this, modify the "end"-parameter to the number of frames
   and the parameters of the Crop command according to your needs.
   
g) Open the .avs script in TMPEGEnc like any other .avi or .mpeg file.
   Although the file is a text file that you created on your own,
   TMPGEnc is able to open it! TMPGEnc will not open it directly;
   it will use the AVISynth "driver", which on-the-fly collects the
   BMP's to a movie. Make the appropriate encoding settings in 
   the properties of TMPEGEnc and save the movie into an MPG-file.
