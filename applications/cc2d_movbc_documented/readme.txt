This is the first documented version of CC2D_MOVBC.

The version itself is an inofficial, experimental code. 
There might be arbitrary many errors in the documentation - 
from the spelling point of view as well as in the content.

The code base of CC2D also has slightly changed as the solver can
now be called as a subroutine. The code is cleaned up a little bit
and structured into different subdirectories in the src/ directory.
There is also some experimental test code in the src/notusedgoodies/ 
subdirectory, which is not used by the program but may help in 
developing applications in Fortran 77.

If you find errors, you might like to contact us
(featflow@featflow.de) to inform about that.
You may also write questions to us or me if anything is unclear.
But be aware that (as this version is not really official) there might
be no time for answering it ;)

USE AT YOUR OWN RISK - like everything of Feat(Flow).

P.S.: 1.) To control the memory usage, open the file include/cmem.inc
          and comment in the definition of NNWORK you need!
      2.) To build the executable, simply type "make" in this
          directory (as usual).
      

