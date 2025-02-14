OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(-0.88859963) q[0];
rz(-1.3130045) q[1];
sx q[1];
rz(-1.5994025) q[1];
sx q[1];
rz(-1.3808274) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2228969) q[0];
sx q[0];
rz(-1.7641506) q[0];
sx q[0];
rz(-1.3816637) q[0];
rz(-pi) q[1];
rz(1.0384473) q[2];
sx q[2];
rz(-0.84538904) q[2];
sx q[2];
rz(-2.7388245) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4145248) q[1];
sx q[1];
rz(-1.5611851) q[1];
sx q[1];
rz(-2.8060786) q[1];
x q[2];
rz(-0.22486658) q[3];
sx q[3];
rz(-1.9266085) q[3];
sx q[3];
rz(-2.0840621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31972739) q[2];
sx q[2];
rz(-1.3250019) q[2];
sx q[2];
rz(2.3925609) q[2];
rz(0.15394112) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(-0.099893959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(1.9248167) q[0];
rz(-1.0617537) q[1];
sx q[1];
rz(-0.80445015) q[1];
sx q[1];
rz(0.42713508) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0741074) q[0];
sx q[0];
rz(-0.68746131) q[0];
sx q[0];
rz(0.56337728) q[0];
rz(-pi) q[1];
rz(-0.24332328) q[2];
sx q[2];
rz(-2.2017225) q[2];
sx q[2];
rz(2.6532946) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6500351) q[1];
sx q[1];
rz(-1.2833529) q[1];
sx q[1];
rz(-2.9950055) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0892761) q[3];
sx q[3];
rz(-1.7744831) q[3];
sx q[3];
rz(1.2869204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.079166807) q[2];
sx q[2];
rz(-1.4126974) q[2];
sx q[2];
rz(-1.3986577) q[2];
rz(0.22377293) q[3];
sx q[3];
rz(-0.90884915) q[3];
sx q[3];
rz(1.5435425) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021521213) q[0];
sx q[0];
rz(-1.7912309) q[0];
sx q[0];
rz(-3.0960826) q[0];
rz(-1.9728164) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(0.57560903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63059124) q[0];
sx q[0];
rz(-1.5596034) q[0];
sx q[0];
rz(0.63909792) q[0];
rz(3.124442) q[2];
sx q[2];
rz(-0.91020012) q[2];
sx q[2];
rz(-1.2397546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8519046) q[1];
sx q[1];
rz(-0.59483268) q[1];
sx q[1];
rz(1.965305) q[1];
rz(-0.024954114) q[3];
sx q[3];
rz(-0.91812274) q[3];
sx q[3];
rz(-0.11634532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.0098972926) q[2];
sx q[2];
rz(-2.3991149) q[2];
sx q[2];
rz(-2.1843145) q[2];
rz(0.57473985) q[3];
sx q[3];
rz(-1.5301907) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15605536) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(-1.2611457) q[0];
rz(0.082322923) q[1];
sx q[1];
rz(-2.0539093) q[1];
sx q[1];
rz(-1.3538768) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9808423) q[0];
sx q[0];
rz(-0.35355648) q[0];
sx q[0];
rz(-2.4048231) q[0];
rz(-pi) q[1];
rz(-2.8026695) q[2];
sx q[2];
rz(-0.51033516) q[2];
sx q[2];
rz(-1.5938544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30198797) q[1];
sx q[1];
rz(-1.0763554) q[1];
sx q[1];
rz(-1.6322664) q[1];
rz(-pi) q[2];
rz(0.4889826) q[3];
sx q[3];
rz(-0.95150286) q[3];
sx q[3];
rz(1.3448546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7202619) q[2];
sx q[2];
rz(-0.91840863) q[2];
sx q[2];
rz(1.8850231) q[2];
rz(2.4300857) q[3];
sx q[3];
rz(-1.810775) q[3];
sx q[3];
rz(-2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.310815) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(0.34314439) q[0];
rz(-3.0798196) q[1];
sx q[1];
rz(-2.1715178) q[1];
sx q[1];
rz(-1.4168581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34300229) q[0];
sx q[0];
rz(-1.8483642) q[0];
sx q[0];
rz(-2.9647602) q[0];
rz(-pi) q[1];
rz(1.896865) q[2];
sx q[2];
rz(-2.2448501) q[2];
sx q[2];
rz(1.774396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8843536) q[1];
sx q[1];
rz(-1.4352918) q[1];
sx q[1];
rz(1.0865023) q[1];
rz(0.23483488) q[3];
sx q[3];
rz(-1.1560625) q[3];
sx q[3];
rz(2.107634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5477649) q[2];
sx q[2];
rz(-1.4676899) q[2];
sx q[2];
rz(-2.055638) q[2];
rz(-2.2222399) q[3];
sx q[3];
rz(-1.3708401) q[3];
sx q[3];
rz(-2.5277188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3980961) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(0.79291517) q[0];
rz(-0.74053699) q[1];
sx q[1];
rz(-2.1426327) q[1];
sx q[1];
rz(0.83121306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39033) q[0];
sx q[0];
rz(-2.7075276) q[0];
sx q[0];
rz(2.084226) q[0];
rz(-pi) q[1];
rz(1.5969876) q[2];
sx q[2];
rz(-0.27715836) q[2];
sx q[2];
rz(1.8089393) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7523642) q[1];
sx q[1];
rz(-0.78654754) q[1];
sx q[1];
rz(-2.0081372) q[1];
x q[2];
rz(-3.1208349) q[3];
sx q[3];
rz(-1.100173) q[3];
sx q[3];
rz(2.1955459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0703766) q[2];
sx q[2];
rz(-1.2717609) q[2];
sx q[2];
rz(1.1191204) q[2];
rz(-1.2707155) q[3];
sx q[3];
rz(-2.0344574) q[3];
sx q[3];
rz(-1.7601815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0041466) q[0];
sx q[0];
rz(-0.16682145) q[0];
sx q[0];
rz(-1.5420089) q[0];
rz(-2.1265325) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-0.17280811) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5240066) q[0];
sx q[0];
rz(-2.3476971) q[0];
sx q[0];
rz(-2.2909597) q[0];
x q[1];
rz(1.2020993) q[2];
sx q[2];
rz(-2.2208344) q[2];
sx q[2];
rz(-1.4024613) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1041278) q[1];
sx q[1];
rz(-1.3224673) q[1];
sx q[1];
rz(-2.009997) q[1];
x q[2];
rz(2.5739848) q[3];
sx q[3];
rz(-0.85383534) q[3];
sx q[3];
rz(0.53138083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.479594) q[2];
sx q[2];
rz(-2.2981503) q[2];
sx q[2];
rz(-0.97770989) q[2];
rz(0.57502037) q[3];
sx q[3];
rz(-1.9366555) q[3];
sx q[3];
rz(1.2303111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8136895) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(-3.1258702) q[0];
rz(1.9533336) q[1];
sx q[1];
rz(-2.5091722) q[1];
sx q[1];
rz(1.9421008) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6576783) q[0];
sx q[0];
rz(-1.8059314) q[0];
sx q[0];
rz(-1.1170618) q[0];
x q[1];
rz(-2.0616777) q[2];
sx q[2];
rz(-1.617031) q[2];
sx q[2];
rz(0.85109988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3856628) q[1];
sx q[1];
rz(-2.0634868) q[1];
sx q[1];
rz(2.1098968) q[1];
rz(-1.1114798) q[3];
sx q[3];
rz(-2.5190341) q[3];
sx q[3];
rz(0.42296577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99379313) q[2];
sx q[2];
rz(-1.2387929) q[2];
sx q[2];
rz(1.3851059) q[2];
rz(0.6238474) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-2.3210242) q[0];
sx q[0];
rz(0.69325915) q[0];
rz(2.8765826) q[1];
sx q[1];
rz(-2.3131504) q[1];
sx q[1];
rz(-1.5230491) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5833359) q[0];
sx q[0];
rz(-0.77363217) q[0];
sx q[0];
rz(-3.1257854) q[0];
rz(2.5184143) q[2];
sx q[2];
rz(-2.0668536) q[2];
sx q[2];
rz(-3.0219699) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48623006) q[1];
sx q[1];
rz(-1.910348) q[1];
sx q[1];
rz(1.5454253) q[1];
rz(-pi) q[2];
rz(-1.4632439) q[3];
sx q[3];
rz(-0.58379025) q[3];
sx q[3];
rz(-1.8369758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.5134209) q[2];
sx q[2];
rz(-1.4576853) q[2];
rz(2.1863106) q[3];
sx q[3];
rz(-0.49946076) q[3];
sx q[3];
rz(-0.83520755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38462287) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(0.48267522) q[0];
rz(-0.92974281) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(-0.39628705) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9549462) q[0];
sx q[0];
rz(-2.5781693) q[0];
sx q[0];
rz(2.4231829) q[0];
rz(-pi) q[1];
rz(1.7223139) q[2];
sx q[2];
rz(-1.7477525) q[2];
sx q[2];
rz(-2.3529476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0490108) q[1];
sx q[1];
rz(-1.7515469) q[1];
sx q[1];
rz(0.84047079) q[1];
x q[2];
rz(2.4195005) q[3];
sx q[3];
rz(-0.51755899) q[3];
sx q[3];
rz(0.90921569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3489909) q[2];
sx q[2];
rz(-1.7020117) q[2];
sx q[2];
rz(2.8144042) q[2];
rz(0.78091019) q[3];
sx q[3];
rz(-1.9802997) q[3];
sx q[3];
rz(1.4769295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1031716) q[0];
sx q[0];
rz(-0.99089834) q[0];
sx q[0];
rz(0.25767576) q[0];
rz(0.36956638) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(-0.69330458) q[2];
sx q[2];
rz(-1.5360166) q[2];
sx q[2];
rz(1.7770384) q[2];
rz(0.9624858) q[3];
sx q[3];
rz(-1.2506093) q[3];
sx q[3];
rz(0.067230732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
