OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5481663) q[0];
sx q[0];
rz(3.610008) q[0];
sx q[0];
rz(9.2863291) q[0];
rz(1.8847213) q[1];
sx q[1];
rz(-2.0710019) q[1];
sx q[1];
rz(-0.021477403) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8801418) q[0];
sx q[0];
rz(-2.7316748) q[0];
sx q[0];
rz(-1.0410801) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8201222) q[2];
sx q[2];
rz(-0.73222762) q[2];
sx q[2];
rz(2.2158465) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5925046) q[1];
sx q[1];
rz(-0.23356423) q[1];
sx q[1];
rz(0.84850581) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5455467) q[3];
sx q[3];
rz(-1.2674244) q[3];
sx q[3];
rz(1.9762126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5452177) q[2];
sx q[2];
rz(-2.6863828) q[2];
sx q[2];
rz(1.630416) q[2];
rz(3.0916072) q[3];
sx q[3];
rz(-1.2286011) q[3];
sx q[3];
rz(0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2216457) q[0];
sx q[0];
rz(-1.1476465) q[0];
sx q[0];
rz(3.1044712) q[0];
rz(3.1275753) q[1];
sx q[1];
rz(-0.57828301) q[1];
sx q[1];
rz(0.23606539) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35821298) q[0];
sx q[0];
rz(-1.556408) q[0];
sx q[0];
rz(-3.0667801) q[0];
rz(-pi) q[1];
rz(-2.5949941) q[2];
sx q[2];
rz(-2.3476217) q[2];
sx q[2];
rz(-3.1304718) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3379682) q[1];
sx q[1];
rz(-1.9527862) q[1];
sx q[1];
rz(2.1398786) q[1];
rz(-pi) q[2];
rz(-0.25059741) q[3];
sx q[3];
rz(-1.2211694) q[3];
sx q[3];
rz(-2.6299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0125121) q[2];
sx q[2];
rz(-1.3690288) q[2];
sx q[2];
rz(0.033626076) q[2];
rz(0.29020894) q[3];
sx q[3];
rz(-2.2955194) q[3];
sx q[3];
rz(-0.39053759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.520312) q[0];
sx q[0];
rz(-2.1618167) q[0];
sx q[0];
rz(-2.8424971) q[0];
rz(-1.5215123) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(-0.026195899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55518245) q[0];
sx q[0];
rz(-1.5918709) q[0];
sx q[0];
rz(-0.0064959244) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1059987) q[2];
sx q[2];
rz(-0.31162497) q[2];
sx q[2];
rz(-1.0958042) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4995998) q[1];
sx q[1];
rz(-2.2471554) q[1];
sx q[1];
rz(0.73222668) q[1];
x q[2];
rz(0.44680178) q[3];
sx q[3];
rz(-1.3626771) q[3];
sx q[3];
rz(-1.3835386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4574778) q[2];
sx q[2];
rz(-2.699615) q[2];
sx q[2];
rz(0.6074062) q[2];
rz(-2.7720747) q[3];
sx q[3];
rz(-1.4990467) q[3];
sx q[3];
rz(0.0079689715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5615416) q[0];
sx q[0];
rz(-2.9396368) q[0];
sx q[0];
rz(1.1308905) q[0];
rz(0.35218969) q[1];
sx q[1];
rz(-0.49030855) q[1];
sx q[1];
rz(0.21603781) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433885) q[0];
sx q[0];
rz(-1.8134612) q[0];
sx q[0];
rz(-0.46671861) q[0];
rz(0.23996578) q[2];
sx q[2];
rz(-0.45825175) q[2];
sx q[2];
rz(-1.3461408) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6747848) q[1];
sx q[1];
rz(-2.4744316) q[1];
sx q[1];
rz(-0.902456) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7407577) q[3];
sx q[3];
rz(-1.483184) q[3];
sx q[3];
rz(-0.57249285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5242247) q[2];
sx q[2];
rz(-0.91297954) q[2];
sx q[2];
rz(0.53140223) q[2];
rz(1.6898539) q[3];
sx q[3];
rz(-2.6302591) q[3];
sx q[3];
rz(-2.1112198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78469974) q[0];
sx q[0];
rz(-2.0578616) q[0];
sx q[0];
rz(0.30237958) q[0];
rz(1.1131635) q[1];
sx q[1];
rz(-1.0013564) q[1];
sx q[1];
rz(1.0611634) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36398023) q[0];
sx q[0];
rz(-1.6355231) q[0];
sx q[0];
rz(-1.510468) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9169754) q[2];
sx q[2];
rz(-1.6664522) q[2];
sx q[2];
rz(3.0358853) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50143277) q[1];
sx q[1];
rz(-1.6987015) q[1];
sx q[1];
rz(-0.78269614) q[1];
x q[2];
rz(-1.8228028) q[3];
sx q[3];
rz(-1.8600501) q[3];
sx q[3];
rz(0.94925971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82659668) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(-0.19485168) q[2];
rz(2.3991614) q[3];
sx q[3];
rz(-2.4104379) q[3];
sx q[3];
rz(0.2754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18405296) q[0];
sx q[0];
rz(-2.4803211) q[0];
sx q[0];
rz(1.1141962) q[0];
rz(-0.30955744) q[1];
sx q[1];
rz(-0.93300262) q[1];
sx q[1];
rz(-1.385744) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64460574) q[0];
sx q[0];
rz(-0.33503767) q[0];
sx q[0];
rz(-1.0100288) q[0];
x q[1];
rz(1.0631752) q[2];
sx q[2];
rz(-1.6706027) q[2];
sx q[2];
rz(-0.55528477) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0506865) q[1];
sx q[1];
rz(-0.93877568) q[1];
sx q[1];
rz(1.2463831) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4900804) q[3];
sx q[3];
rz(-2.2527472) q[3];
sx q[3];
rz(1.7891974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4625357) q[2];
sx q[2];
rz(-2.6947196) q[2];
sx q[2];
rz(-0.40979579) q[2];
rz(-3.0637686) q[3];
sx q[3];
rz(-1.9255368) q[3];
sx q[3];
rz(-2.3791544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536612) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(-1.4713564) q[0];
rz(0.81368601) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(-2.241316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48834688) q[0];
sx q[0];
rz(-1.4524967) q[0];
sx q[0];
rz(1.6315881) q[0];
rz(1.8531606) q[2];
sx q[2];
rz(-2.8430364) q[2];
sx q[2];
rz(2.3207842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7309396) q[1];
sx q[1];
rz(-1.5789642) q[1];
sx q[1];
rz(1.3047303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6358709) q[3];
sx q[3];
rz(-1.3378007) q[3];
sx q[3];
rz(0.20291337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6206996) q[2];
sx q[2];
rz(-2.7405379) q[2];
sx q[2];
rz(2.9996784) q[2];
rz(2.9698931) q[3];
sx q[3];
rz(-1.2407691) q[3];
sx q[3];
rz(-2.5920674) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5080344) q[0];
sx q[0];
rz(-1.1548076) q[0];
sx q[0];
rz(-1.5800193) q[0];
rz(-0.70478565) q[1];
sx q[1];
rz(-0.90122688) q[1];
sx q[1];
rz(-2.8660692) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47833035) q[0];
sx q[0];
rz(-1.6323292) q[0];
sx q[0];
rz(1.6986153) q[0];
rz(-pi) q[1];
rz(0.24532206) q[2];
sx q[2];
rz(-1.5333912) q[2];
sx q[2];
rz(-2.268301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2299978) q[1];
sx q[1];
rz(-1.8640638) q[1];
sx q[1];
rz(1.0469584) q[1];
rz(-pi) q[2];
rz(2.4750214) q[3];
sx q[3];
rz(-1.2351994) q[3];
sx q[3];
rz(-0.94387142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0813109) q[2];
sx q[2];
rz(-1.6401446) q[2];
sx q[2];
rz(-2.8509129) q[2];
rz(-2.3447013) q[3];
sx q[3];
rz(-0.47178888) q[3];
sx q[3];
rz(2.1577788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31442916) q[0];
sx q[0];
rz(-1.6570579) q[0];
sx q[0];
rz(2.2737801) q[0];
rz(0.21226352) q[1];
sx q[1];
rz(-2.2607195) q[1];
sx q[1];
rz(2.8689522) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72948706) q[0];
sx q[0];
rz(-0.22815591) q[0];
sx q[0];
rz(1.233393) q[0];
rz(2.3219982) q[2];
sx q[2];
rz(-1.4965881) q[2];
sx q[2];
rz(1.5774262) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8729108) q[1];
sx q[1];
rz(-2.2051393) q[1];
sx q[1];
rz(0.047206248) q[1];
rz(-0.61278559) q[3];
sx q[3];
rz(-0.49344188) q[3];
sx q[3];
rz(2.8145344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0923882) q[2];
sx q[2];
rz(-0.35085756) q[2];
sx q[2];
rz(1.997088) q[2];
rz(2.5149964) q[3];
sx q[3];
rz(-2.383039) q[3];
sx q[3];
rz(0.315061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7476244) q[0];
sx q[0];
rz(-0.67248857) q[0];
sx q[0];
rz(-0.59338635) q[0];
rz(-0.74140948) q[1];
sx q[1];
rz(-2.4948222) q[1];
sx q[1];
rz(-1.5945565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5760959) q[0];
sx q[0];
rz(-1.4617697) q[0];
sx q[0];
rz(-1.7897357) q[0];
x q[1];
rz(-2.6213264) q[2];
sx q[2];
rz(-2.6718585) q[2];
sx q[2];
rz(0.007291468) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50362557) q[1];
sx q[1];
rz(-2.4071998) q[1];
sx q[1];
rz(0.10403688) q[1];
rz(1.3688332) q[3];
sx q[3];
rz(-1.5416577) q[3];
sx q[3];
rz(2.6239708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0934304) q[2];
sx q[2];
rz(-2.5967345) q[2];
sx q[2];
rz(2.1989934) q[2];
rz(-1.4283098) q[3];
sx q[3];
rz(-1.8918248) q[3];
sx q[3];
rz(-1.9863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9485332) q[0];
sx q[0];
rz(-1.5567224) q[0];
sx q[0];
rz(-1.3514883) q[0];
rz(1.1763186) q[1];
sx q[1];
rz(-2.2358924) q[1];
sx q[1];
rz(2.4814503) q[1];
rz(-1.0084441) q[2];
sx q[2];
rz(-2.8165419) q[2];
sx q[2];
rz(-1.4680924) q[2];
rz(0.42849937) q[3];
sx q[3];
rz(-0.57572031) q[3];
sx q[3];
rz(-2.4888103) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
