Possible_bound_strange_hypernuclei

calculate the binding energy of three body nuclei in ./main (or ./paralle/main)
calculate the wigner function by the wave function calculated by ./main in ./wigner (or ./paralle/wigner)
in ./blast-wave-threebody, it will give the production of the nuclear;
you can edit nuclear you want in ./include, remember to add the head file in ./RRcoeff, ./main (or ./paralle/main) and ./blast-wave-threebody, then re-make it;
./RRcoeff will produce a RR coefficient matrix; after add a new nuclear, need to calculate it first;
./paralle produce a program to submit on severs
