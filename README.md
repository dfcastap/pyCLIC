(based on the F77 CLIC)
Written by Diego Castaneda

Based on the Fortran code of Catherine Lovekin (2006ApJ...643..460L), it calculates the observed SED of a rotating stellar model at any inclination. It relies on PHOENIX atmospheric models (more specifically, the *.70 intensity files).
It uses a couple of F90 modules that need to be compiled before it can run: (assuming you have Numpy, run the following)

f2py -c -m fluxtable.f90 fluxtable
gfotran wavelengthinterp.f90 -o wavelengthinterplog.exe

After that, main.py has all the code that needs to be run.
