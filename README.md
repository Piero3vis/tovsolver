Welcome to the TOV RK4 solver!

To run the code please note:

-You need Python 2.7 or a more recent version.
-From the terminal run the script specifying in the same line:
							- Central density in cgs
							- K (in c=G=Msun=1 units)
							- polytropic exponent gamma 
like in the example below:

%%%%%%%%%%%%%%%%Example%%%%%%%%%%%%%%%%

python TOVsolver.py 5e14 300000 2.75

%%%%%%%%%%%%%%%%Example%%%%%%%%%%%%%%%%%



The solver will provide the pressure, mass, density profile and save it in the same folder.
Furthermore all data are saved in a centraldensityvalue.txt file.

To use the Euler method simply remove the comment marks in the code and comment the RK4 method where is specified .

If there is need to store the final mass/stellar radius to further studies create a .txt file named massradiusgraf.txt
The program automatically will store there the Central density / Total mass / Stellar Radius.


Good luck and happy coding!




