The goal of this package is to make available the Collider Time Evolution or CTE particle 
tracking simulation code. Originally this code was developed in FORTRAN, making it hard to
read, maintain and extend. This package is the first attempt to improve on that,
where the heavy numerical calculations are kept in FORTRAN routines, distributed in different files with filenames
refering to their functionality. The main program has been shifted out of FORTRAN into PYHON,
increasing readability and extendability (e.g. the code is now extended to have different 
particle types in both beams). 

