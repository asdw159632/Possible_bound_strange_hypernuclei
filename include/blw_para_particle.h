
#ifdef E200GeV
//Omega
double Tkin_Omega = 0.1116;//GeV
double rho_0_Omega = 0.9;
const int num_Omega = 100;//need to divide 100*1.8

//nucleon
double Tkin_nucleon = 0.1116;//GeV
double rho_0_nucleon = 0.98;
const int num_nucleon = 17;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//
#elif defined E2_76TeV
//Omega
double Tkin_Omega = 0.122;//GeV
double rho_0_Omega = 1.07;
const int num_Omega = 100;//need to divide 100*1.9785995

//nucleon
double Tkin_nucleon = 0.122;//GeV
double rho_0_nucleon = 1.2;
const int num_nucleon = 35;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//
#endif

