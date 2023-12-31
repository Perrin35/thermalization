OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6808074) q[0];
sx q[0];
rz(2.1587125) q[0];
sx q[0];
rz(8.4150253) q[0];
rz(2.9653964) q[1];
sx q[1];
rz(5.3806452) q[1];
sx q[1];
rz(7.5723958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50531189) q[0];
sx q[0];
rz(-1.8587776) q[0];
sx q[0];
rz(-1.0907537) q[0];
rz(-pi) q[1];
rz(-0.024359811) q[2];
sx q[2];
rz(-1.5640386) q[2];
sx q[2];
rz(0.18505219) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28469052) q[1];
sx q[1];
rz(-0.33387091) q[1];
sx q[1];
rz(2.2536709) q[1];
rz(-pi) q[2];
rz(1.8390158) q[3];
sx q[3];
rz(-1.6383639) q[3];
sx q[3];
rz(-2.4951631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5916799) q[2];
sx q[2];
rz(-2.9788571) q[2];
sx q[2];
rz(3.0060449) q[2];
rz(3.0013951) q[3];
sx q[3];
rz(-1.0547767) q[3];
sx q[3];
rz(-0.31301096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41483375) q[0];
sx q[0];
rz(-0.93678513) q[0];
sx q[0];
rz(2.360789) q[0];
rz(0.26023284) q[1];
sx q[1];
rz(-0.5865016) q[1];
sx q[1];
rz(-1.82812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75269964) q[0];
sx q[0];
rz(-1.7912731) q[0];
sx q[0];
rz(-1.7874784) q[0];
x q[1];
rz(-1.08554) q[2];
sx q[2];
rz(-1.2962356) q[2];
sx q[2];
rz(-2.1833414) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4310415) q[1];
sx q[1];
rz(-1.6574873) q[1];
sx q[1];
rz(-2.980742) q[1];
rz(-pi) q[2];
rz(-1.4649269) q[3];
sx q[3];
rz(-1.6459811) q[3];
sx q[3];
rz(0.74010805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.446622) q[2];
sx q[2];
rz(-1.5335252) q[2];
sx q[2];
rz(-2.1874645) q[2];
rz(-1.702884) q[3];
sx q[3];
rz(-1.8372767) q[3];
sx q[3];
rz(0.71189705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0838098) q[0];
sx q[0];
rz(-0.6018146) q[0];
sx q[0];
rz(2.4457248) q[0];
rz(-0.35481915) q[1];
sx q[1];
rz(-1.4777007) q[1];
sx q[1];
rz(-0.16608873) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21712185) q[0];
sx q[0];
rz(-2.8170601) q[0];
sx q[0];
rz(1.0990418) q[0];
rz(-pi) q[1];
rz(0.66785779) q[2];
sx q[2];
rz(-2.4154426) q[2];
sx q[2];
rz(-2.6454676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3661763) q[1];
sx q[1];
rz(-2.011538) q[1];
sx q[1];
rz(0.52904769) q[1];
rz(-pi) q[2];
x q[2];
rz(2.283964) q[3];
sx q[3];
rz(-1.3751251) q[3];
sx q[3];
rz(-1.5426202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35963905) q[2];
sx q[2];
rz(-2.226318) q[2];
sx q[2];
rz(-2.1598143) q[2];
rz(-1.123547) q[3];
sx q[3];
rz(-2.8765364) q[3];
sx q[3];
rz(-1.8410929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05900255) q[0];
sx q[0];
rz(-2.1868717) q[0];
sx q[0];
rz(0.29378763) q[0];
rz(-2.7000973) q[1];
sx q[1];
rz(-2.0916633) q[1];
sx q[1];
rz(2.3667483) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.144865) q[0];
sx q[0];
rz(-1.5872247) q[0];
sx q[0];
rz(-0.073547151) q[0];
rz(-pi) q[1];
x q[1];
rz(1.61548) q[2];
sx q[2];
rz(-1.2031021) q[2];
sx q[2];
rz(1.9862663) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0171417) q[1];
sx q[1];
rz(-1.972773) q[1];
sx q[1];
rz(0.084852858) q[1];
x q[2];
rz(1.8294228) q[3];
sx q[3];
rz(-1.4050409) q[3];
sx q[3];
rz(-0.37763077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1321156) q[2];
sx q[2];
rz(-2.0602132) q[2];
sx q[2];
rz(0.57007989) q[2];
rz(-1.305497) q[3];
sx q[3];
rz(-0.92731849) q[3];
sx q[3];
rz(-1.6397887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36561361) q[0];
sx q[0];
rz(-1.403084) q[0];
sx q[0];
rz(-0.18138012) q[0];
rz(2.5326305) q[1];
sx q[1];
rz(-0.47157559) q[1];
sx q[1];
rz(0.10770527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90751326) q[0];
sx q[0];
rz(-2.7764752) q[0];
sx q[0];
rz(1.8318729) q[0];
x q[1];
rz(0.55107112) q[2];
sx q[2];
rz(-2.1398009) q[2];
sx q[2];
rz(1.0702733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4957461) q[1];
sx q[1];
rz(-1.1762113) q[1];
sx q[1];
rz(2.2699725) q[1];
x q[2];
rz(2.0884573) q[3];
sx q[3];
rz(-0.62873757) q[3];
sx q[3];
rz(-2.351458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.057377664) q[2];
sx q[2];
rz(-1.6485018) q[2];
sx q[2];
rz(2.8520612) q[2];
rz(0.036751898) q[3];
sx q[3];
rz(-2.3993902) q[3];
sx q[3];
rz(0.010820476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0832131) q[0];
sx q[0];
rz(-0.20861861) q[0];
sx q[0];
rz(1.5104729) q[0];
rz(-2.0026813) q[1];
sx q[1];
rz(-1.4803671) q[1];
sx q[1];
rz(-2.1246134) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4646343) q[0];
sx q[0];
rz(-0.58566739) q[0];
sx q[0];
rz(-0.93924378) q[0];
x q[1];
rz(2.9404042) q[2];
sx q[2];
rz(-1.5474461) q[2];
sx q[2];
rz(-0.77229283) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0000227) q[1];
sx q[1];
rz(-1.1035898) q[1];
sx q[1];
rz(2.1187374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0820886) q[3];
sx q[3];
rz(-1.0923311) q[3];
sx q[3];
rz(2.2389776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5974474) q[2];
sx q[2];
rz(-2.1400673) q[2];
sx q[2];
rz(-0.63465676) q[2];
rz(2.0641816) q[3];
sx q[3];
rz(-2.1118739) q[3];
sx q[3];
rz(-2.6950148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8909661) q[0];
sx q[0];
rz(-0.64193305) q[0];
sx q[0];
rz(-2.2977258) q[0];
rz(-1.5232874) q[1];
sx q[1];
rz(-2.4089676) q[1];
sx q[1];
rz(-1.9218146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6424113) q[0];
sx q[0];
rz(-1.5326323) q[0];
sx q[0];
rz(0.092060815) q[0];
rz(-pi) q[1];
rz(0.90791038) q[2];
sx q[2];
rz(-0.50953509) q[2];
sx q[2];
rz(-2.8223035) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0724832) q[1];
sx q[1];
rz(-1.3246037) q[1];
sx q[1];
rz(-2.2407131) q[1];
rz(-pi) q[2];
rz(-1.1620429) q[3];
sx q[3];
rz(-1.6533274) q[3];
sx q[3];
rz(2.2857034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5119778) q[2];
sx q[2];
rz(-0.93705606) q[2];
sx q[2];
rz(-2.8249595) q[2];
rz(3.1043502) q[3];
sx q[3];
rz(-1.3469632) q[3];
sx q[3];
rz(-2.3855456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1036296) q[0];
sx q[0];
rz(-2.0541971) q[0];
sx q[0];
rz(-2.8102002) q[0];
rz(1.1720852) q[1];
sx q[1];
rz(-2.4786699) q[1];
sx q[1];
rz(-0.40922871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6054571) q[0];
sx q[0];
rz(-0.94263173) q[0];
sx q[0];
rz(-0.34988846) q[0];
rz(-pi) q[1];
rz(0.51034575) q[2];
sx q[2];
rz(-2.5556457) q[2];
sx q[2];
rz(-2.3454391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48620263) q[1];
sx q[1];
rz(-1.7230002) q[1];
sx q[1];
rz(0.57032077) q[1];
x q[2];
rz(0.41142923) q[3];
sx q[3];
rz(-0.55555389) q[3];
sx q[3];
rz(-0.065447741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.99872148) q[2];
sx q[2];
rz(-1.8565535) q[2];
sx q[2];
rz(2.5869353) q[2];
rz(2.5000642) q[3];
sx q[3];
rz(-0.21637622) q[3];
sx q[3];
rz(3.0564814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2239969) q[0];
sx q[0];
rz(-1.6308835) q[0];
sx q[0];
rz(0.0099442033) q[0];
rz(-2.0857281) q[1];
sx q[1];
rz(-1.5751585) q[1];
sx q[1];
rz(-1.508629) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2752339) q[0];
sx q[0];
rz(-2.4320721) q[0];
sx q[0];
rz(-1.8190246) q[0];
rz(-pi) q[1];
rz(-0.091392322) q[2];
sx q[2];
rz(-0.95801991) q[2];
sx q[2];
rz(2.6932655) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.846261) q[1];
sx q[1];
rz(-1.8896855) q[1];
sx q[1];
rz(2.1178513) q[1];
rz(-3.085) q[3];
sx q[3];
rz(-2.2731299) q[3];
sx q[3];
rz(3.1346723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9987954) q[2];
sx q[2];
rz(-1.5441511) q[2];
sx q[2];
rz(0.43668401) q[2];
rz(-1.3302594) q[3];
sx q[3];
rz(-2.4618849) q[3];
sx q[3];
rz(-2.1127624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0325539) q[0];
sx q[0];
rz(-0.80820307) q[0];
sx q[0];
rz(2.57634) q[0];
rz(-2.8887707) q[1];
sx q[1];
rz(-1.0193453) q[1];
sx q[1];
rz(1.3814829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94986445) q[0];
sx q[0];
rz(-1.0690332) q[0];
sx q[0];
rz(-1.6185332) q[0];
rz(-pi) q[1];
rz(-1.6818468) q[2];
sx q[2];
rz(-2.0425218) q[2];
sx q[2];
rz(-1.9327088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6227464) q[1];
sx q[1];
rz(-1.985605) q[1];
sx q[1];
rz(-0.82621375) q[1];
rz(-0.85716565) q[3];
sx q[3];
rz(-1.790159) q[3];
sx q[3];
rz(-2.9227748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4505724) q[2];
sx q[2];
rz(-2.317163) q[2];
sx q[2];
rz(0.59664574) q[2];
rz(0.48554844) q[3];
sx q[3];
rz(-2.2462626) q[3];
sx q[3];
rz(-0.027464494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52994603) q[0];
sx q[0];
rz(-1.2048789) q[0];
sx q[0];
rz(-2.0299029) q[0];
rz(-1.272841) q[1];
sx q[1];
rz(-2.0618989) q[1];
sx q[1];
rz(2.3088574) q[1];
rz(1.3960719) q[2];
sx q[2];
rz(-0.25824418) q[2];
sx q[2];
rz(1.2373409) q[2];
rz(-0.5100526) q[3];
sx q[3];
rz(-1.6418213) q[3];
sx q[3];
rz(2.1591975) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
