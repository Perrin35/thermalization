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
rz(-1.5934264) q[0];
sx q[0];
rz(-0.46841535) q[0];
sx q[0];
rz(-3.0031437) q[0];
rz(-1.2568714) q[1];
sx q[1];
rz(-1.0705907) q[1];
sx q[1];
rz(-3.1201153) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3119451) q[0];
sx q[0];
rz(-1.2197131) q[0];
sx q[0];
rz(0.21613516) q[0];
rz(-pi) q[1];
rz(-2.3214705) q[2];
sx q[2];
rz(-2.409365) q[2];
sx q[2];
rz(2.2158465) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8565377) q[1];
sx q[1];
rz(-1.3962588) q[1];
sx q[1];
rz(0.15600295) q[1];
rz(-pi) q[2];
rz(-2.6329166) q[3];
sx q[3];
rz(-2.4812316) q[3];
sx q[3];
rz(2.3213399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5452177) q[2];
sx q[2];
rz(-2.6863828) q[2];
sx q[2];
rz(1.630416) q[2];
rz(0.049985416) q[3];
sx q[3];
rz(-1.2286011) q[3];
sx q[3];
rz(-0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91994691) q[0];
sx q[0];
rz(-1.1476465) q[0];
sx q[0];
rz(-3.1044712) q[0];
rz(3.1275753) q[1];
sx q[1];
rz(-2.5633096) q[1];
sx q[1];
rz(2.9055273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0229313) q[0];
sx q[0];
rz(-0.076181024) q[0];
sx q[0];
rz(-2.9514022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0843956) q[2];
sx q[2];
rz(-0.91569967) q[2];
sx q[2];
rz(2.4156605) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8467847) q[1];
sx q[1];
rz(-0.67344147) q[1];
sx q[1];
rz(2.2114141) q[1];
x q[2];
rz(2.8909952) q[3];
sx q[3];
rz(-1.9204233) q[3];
sx q[3];
rz(2.6299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0125121) q[2];
sx q[2];
rz(-1.7725638) q[2];
sx q[2];
rz(-3.1079666) q[2];
rz(2.8513837) q[3];
sx q[3];
rz(-2.2955194) q[3];
sx q[3];
rz(0.39053759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.520312) q[0];
sx q[0];
rz(-0.97977591) q[0];
sx q[0];
rz(2.8424971) q[0];
rz(1.6200804) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(-0.026195899) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55518245) q[0];
sx q[0];
rz(-1.5497218) q[0];
sx q[0];
rz(3.1350967) q[0];
x q[1];
rz(-2.0355939) q[2];
sx q[2];
rz(-2.8299677) q[2];
sx q[2];
rz(2.0457884) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6419928) q[1];
sx q[1];
rz(-2.2471554) q[1];
sx q[1];
rz(0.73222668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3407767) q[3];
sx q[3];
rz(-1.134308) q[3];
sx q[3];
rz(-3.0530086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4574778) q[2];
sx q[2];
rz(-2.699615) q[2];
sx q[2];
rz(-2.5341865) q[2];
rz(0.36951798) q[3];
sx q[3];
rz(-1.642546) q[3];
sx q[3];
rz(-0.0079689715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5615416) q[0];
sx q[0];
rz(-0.20195584) q[0];
sx q[0];
rz(1.1308905) q[0];
rz(2.789403) q[1];
sx q[1];
rz(-0.49030855) q[1];
sx q[1];
rz(2.9255548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2239922) q[0];
sx q[0];
rz(-2.619714) q[0];
sx q[0];
rz(0.50295575) q[0];
rz(-pi) q[1];
rz(2.9016269) q[2];
sx q[2];
rz(-2.6833409) q[2];
sx q[2];
rz(1.7954519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6591726) q[1];
sx q[1];
rz(-1.9643088) q[1];
sx q[1];
rz(-2.1244786) q[1];
rz(-pi) q[2];
rz(-0.0888865) q[3];
sx q[3];
rz(-1.4014932) q[3];
sx q[3];
rz(2.1282738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.617368) q[2];
sx q[2];
rz(-2.2286131) q[2];
sx q[2];
rz(2.6101904) q[2];
rz(1.4517387) q[3];
sx q[3];
rz(-2.6302591) q[3];
sx q[3];
rz(2.1112198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78469974) q[0];
sx q[0];
rz(-1.0837311) q[0];
sx q[0];
rz(0.30237958) q[0];
rz(1.1131635) q[1];
sx q[1];
rz(-2.1402363) q[1];
sx q[1];
rz(2.0804292) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9386834) q[0];
sx q[0];
rz(-1.6309982) q[0];
sx q[0];
rz(0.06484443) q[0];
x q[1];
rz(0.40675478) q[2];
sx q[2];
rz(-2.8977721) q[2];
sx q[2];
rz(-1.8610473) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94316846) q[1];
sx q[1];
rz(-2.3454002) q[1];
sx q[1];
rz(1.7502341) q[1];
x q[2];
rz(1.3187899) q[3];
sx q[3];
rz(-1.2815426) q[3];
sx q[3];
rz(-0.94925971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.314996) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(-0.19485168) q[2];
rz(-2.3991614) q[3];
sx q[3];
rz(-2.4104379) q[3];
sx q[3];
rz(2.8661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18405296) q[0];
sx q[0];
rz(-2.4803211) q[0];
sx q[0];
rz(1.1141962) q[0];
rz(0.30955744) q[1];
sx q[1];
rz(-0.93300262) q[1];
sx q[1];
rz(1.385744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64460574) q[0];
sx q[0];
rz(-2.806555) q[0];
sx q[0];
rz(-1.0100288) q[0];
rz(2.0784174) q[2];
sx q[2];
rz(-1.6706027) q[2];
sx q[2];
rz(-2.5863079) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0506865) q[1];
sx q[1];
rz(-0.93877568) q[1];
sx q[1];
rz(1.8952096) q[1];
rz(2.4900804) q[3];
sx q[3];
rz(-0.88884547) q[3];
sx q[3];
rz(-1.7891974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.679057) q[2];
sx q[2];
rz(-2.6947196) q[2];
sx q[2];
rz(-2.7317969) q[2];
rz(-0.077824079) q[3];
sx q[3];
rz(-1.9255368) q[3];
sx q[3];
rz(2.3791544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78793144) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(-1.4713564) q[0];
rz(2.3279066) q[1];
sx q[1];
rz(-0.7209456) q[1];
sx q[1];
rz(-2.241316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012205781) q[0];
sx q[0];
rz(-3.0086522) q[0];
sx q[0];
rz(-2.6690527) q[0];
rz(-1.288432) q[2];
sx q[2];
rz(-0.29855628) q[2];
sx q[2];
rz(-2.3207842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7309396) q[1];
sx q[1];
rz(-1.5789642) q[1];
sx q[1];
rz(1.8368624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6358709) q[3];
sx q[3];
rz(-1.8037919) q[3];
sx q[3];
rz(-2.9386793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6206996) q[2];
sx q[2];
rz(-2.7405379) q[2];
sx q[2];
rz(-2.9996784) q[2];
rz(-0.17169954) q[3];
sx q[3];
rz(-1.2407691) q[3];
sx q[3];
rz(-2.5920674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63355821) q[0];
sx q[0];
rz(-1.9867851) q[0];
sx q[0];
rz(1.5800193) q[0];
rz(0.70478565) q[1];
sx q[1];
rz(-0.90122688) q[1];
sx q[1];
rz(-0.27552342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64618212) q[0];
sx q[0];
rz(-0.14178628) q[0];
sx q[0];
rz(1.1205733) q[0];
rz(-pi) q[1];
rz(1.5322379) q[2];
sx q[2];
rz(-1.3256494) q[2];
sx q[2];
rz(2.4347255) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2299978) q[1];
sx q[1];
rz(-1.2775289) q[1];
sx q[1];
rz(1.0469584) q[1];
x q[2];
rz(0.66657127) q[3];
sx q[3];
rz(-1.2351994) q[3];
sx q[3];
rz(-2.1977212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0813109) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(2.8509129) q[2];
rz(-0.79689133) q[3];
sx q[3];
rz(-2.6698038) q[3];
sx q[3];
rz(-0.98381388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8271635) q[0];
sx q[0];
rz(-1.4845347) q[0];
sx q[0];
rz(2.2737801) q[0];
rz(0.21226352) q[1];
sx q[1];
rz(-2.2607195) q[1];
sx q[1];
rz(2.8689522) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51202392) q[0];
sx q[0];
rz(-1.4958515) q[0];
sx q[0];
rz(-1.7864947) q[0];
x q[1];
rz(-3.0402203) q[2];
sx q[2];
rz(-0.82216149) q[2];
sx q[2];
rz(-3.0790975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8729108) q[1];
sx q[1];
rz(-0.93645331) q[1];
sx q[1];
rz(3.0943864) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2708067) q[3];
sx q[3];
rz(-1.172903) q[3];
sx q[3];
rz(2.7949445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.049204443) q[2];
sx q[2];
rz(-0.35085756) q[2];
sx q[2];
rz(1.1445047) q[2];
rz(0.62659621) q[3];
sx q[3];
rz(-2.383039) q[3];
sx q[3];
rz(2.8265317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3939683) q[0];
sx q[0];
rz(-0.67248857) q[0];
sx q[0];
rz(-2.5482063) q[0];
rz(-2.4001832) q[1];
sx q[1];
rz(-0.64677042) q[1];
sx q[1];
rz(-1.5945565) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55047148) q[0];
sx q[0];
rz(-2.897399) q[0];
sx q[0];
rz(-2.0376192) q[0];
rz(0.41489652) q[2];
sx q[2];
rz(-1.7977568) q[2];
sx q[2];
rz(2.0503873) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98982544) q[1];
sx q[1];
rz(-1.5011468) q[1];
sx q[1];
rz(0.73169692) q[1];
x q[2];
rz(-1.4265028) q[3];
sx q[3];
rz(-0.20402597) q[3];
sx q[3];
rz(-0.91183364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0934304) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(0.94259924) q[2];
rz(1.4283098) q[3];
sx q[3];
rz(-1.2497679) q[3];
sx q[3];
rz(1.1552756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1930595) q[0];
sx q[0];
rz(-1.5567224) q[0];
sx q[0];
rz(-1.3514883) q[0];
rz(1.1763186) q[1];
sx q[1];
rz(-2.2358924) q[1];
sx q[1];
rz(2.4814503) q[1];
rz(2.1331486) q[2];
sx q[2];
rz(-2.8165419) q[2];
sx q[2];
rz(-1.4680924) q[2];
rz(-1.3073714) q[3];
sx q[3];
rz(-2.0888804) q[3];
sx q[3];
rz(-2.9874939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
