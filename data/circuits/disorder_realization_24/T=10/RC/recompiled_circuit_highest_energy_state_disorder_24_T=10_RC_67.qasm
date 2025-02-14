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
rz(-0.10676323) q[0];
sx q[0];
rz(-1.8101298) q[0];
sx q[0];
rz(-1.3504299) q[0];
rz(-2.8332233) q[1];
sx q[1];
rz(-2.8928533) q[1];
sx q[1];
rz(-2.1623936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5738663) q[0];
sx q[0];
rz(-1.4970333) q[0];
sx q[0];
rz(0.7598138) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91186422) q[2];
sx q[2];
rz(-1.2871337) q[2];
sx q[2];
rz(-0.58667573) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0268271) q[1];
sx q[1];
rz(-1.2176759) q[1];
sx q[1];
rz(1.4985282) q[1];
x q[2];
rz(0.37110801) q[3];
sx q[3];
rz(-1.733193) q[3];
sx q[3];
rz(1.1134256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8637434) q[2];
sx q[2];
rz(-1.9305482) q[2];
sx q[2];
rz(-2.5085874) q[2];
rz(2.4999319) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(0.7707001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4577456) q[0];
sx q[0];
rz(-2.1774543) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(-1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(0.32018426) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6711222) q[0];
sx q[0];
rz(-1.1129541) q[0];
sx q[0];
rz(1.1127959) q[0];
rz(-pi) q[1];
rz(2.9302017) q[2];
sx q[2];
rz(-2.0703982) q[2];
sx q[2];
rz(1.4818986) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6588179) q[1];
sx q[1];
rz(-1.5482727) q[1];
sx q[1];
rz(0.080959678) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8289964) q[3];
sx q[3];
rz(-0.22919433) q[3];
sx q[3];
rz(-2.3865139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0188633) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(-1.2790206) q[2];
rz(-3.0458798) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(-1.8224645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1529481) q[0];
sx q[0];
rz(-2.4929292) q[0];
sx q[0];
rz(-1.0308107) q[0];
rz(2.5910494) q[1];
sx q[1];
rz(-2.2830453) q[1];
sx q[1];
rz(1.2446838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1079946) q[0];
sx q[0];
rz(-1.4116086) q[0];
sx q[0];
rz(1.6248996) q[0];
rz(-pi) q[1];
rz(-0.16764201) q[2];
sx q[2];
rz(-0.9603921) q[2];
sx q[2];
rz(0.33627015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74725973) q[1];
sx q[1];
rz(-1.366341) q[1];
sx q[1];
rz(-0.89337272) q[1];
x q[2];
rz(1.9774417) q[3];
sx q[3];
rz(-0.79504025) q[3];
sx q[3];
rz(-0.94889489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2391669) q[2];
sx q[2];
rz(-0.5639762) q[2];
sx q[2];
rz(1.172056) q[2];
rz(0.82550448) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(0.66506213) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.091752) q[0];
sx q[0];
rz(-1.4565775) q[0];
sx q[0];
rz(-2.7276584) q[0];
rz(0.44231689) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(-2.5962459) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2065679) q[0];
sx q[0];
rz(-0.27982084) q[0];
sx q[0];
rz(-0.18144515) q[0];
rz(-pi) q[1];
rz(0.92939922) q[2];
sx q[2];
rz(-0.9050194) q[2];
sx q[2];
rz(-2.5432472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87121039) q[1];
sx q[1];
rz(-2.5426513) q[1];
sx q[1];
rz(-2.0763055) q[1];
x q[2];
rz(0.14800565) q[3];
sx q[3];
rz(-1.0665575) q[3];
sx q[3];
rz(-1.3324236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.44117323) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(2.9803661) q[2];
rz(-1.553933) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(2.5509295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081838354) q[0];
sx q[0];
rz(-1.2677544) q[0];
sx q[0];
rz(-2.4181714) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.864121) q[1];
sx q[1];
rz(1.0378729) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53053601) q[0];
sx q[0];
rz(-0.7270455) q[0];
sx q[0];
rz(3.0226991) q[0];
rz(2.4167929) q[2];
sx q[2];
rz(-2.0033859) q[2];
sx q[2];
rz(-3.070433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7241076) q[1];
sx q[1];
rz(-1.8383664) q[1];
sx q[1];
rz(-0.45326434) q[1];
rz(-pi) q[2];
rz(2.9014189) q[3];
sx q[3];
rz(-0.55910149) q[3];
sx q[3];
rz(-0.19776519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7448685) q[2];
sx q[2];
rz(-1.1149656) q[2];
sx q[2];
rz(2.569516) q[2];
rz(-0.62503302) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6201651) q[0];
sx q[0];
rz(-3.1075952) q[0];
sx q[0];
rz(0.52325621) q[0];
rz(-1.322586) q[1];
sx q[1];
rz(-1.4679642) q[1];
sx q[1];
rz(0.62017131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5145743) q[0];
sx q[0];
rz(-0.73643273) q[0];
sx q[0];
rz(-0.66177841) q[0];
rz(-pi) q[1];
rz(-2.6729726) q[2];
sx q[2];
rz(-1.6362564) q[2];
sx q[2];
rz(-2.5602788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4166231) q[1];
sx q[1];
rz(-0.94601224) q[1];
sx q[1];
rz(-3.0333972) q[1];
x q[2];
rz(1.3955388) q[3];
sx q[3];
rz(-1.7850998) q[3];
sx q[3];
rz(-2.8068723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.72914499) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(1.283851) q[2];
rz(-0.9681975) q[3];
sx q[3];
rz(-1.3788393) q[3];
sx q[3];
rz(-1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17484434) q[0];
sx q[0];
rz(-1.1356249) q[0];
sx q[0];
rz(-1.211776) q[0];
rz(0.9696331) q[1];
sx q[1];
rz(-1.3215093) q[1];
sx q[1];
rz(-1.2744354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14718854) q[0];
sx q[0];
rz(-2.4321803) q[0];
sx q[0];
rz(1.8038007) q[0];
rz(-pi) q[1];
rz(1.4757122) q[2];
sx q[2];
rz(-2.3543752) q[2];
sx q[2];
rz(-1.9958391) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3962773) q[1];
sx q[1];
rz(-2.4419129) q[1];
sx q[1];
rz(1.0902283) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46044965) q[3];
sx q[3];
rz(-1.8619974) q[3];
sx q[3];
rz(-1.2429383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6756639) q[2];
sx q[2];
rz(-1.6589386) q[2];
sx q[2];
rz(2.8181804) q[2];
rz(-0.0023500738) q[3];
sx q[3];
rz(-1.6833143) q[3];
sx q[3];
rz(0.6215483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49331409) q[0];
sx q[0];
rz(-1.6596153) q[0];
sx q[0];
rz(-1.7907273) q[0];
rz(1.3510652) q[1];
sx q[1];
rz(-1.2850782) q[1];
sx q[1];
rz(1.8538808) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0072216) q[0];
sx q[0];
rz(-2.8343081) q[0];
sx q[0];
rz(-1.1728806) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9016198) q[2];
sx q[2];
rz(-0.74218732) q[2];
sx q[2];
rz(1.0617219) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2877174) q[1];
sx q[1];
rz(-2.7734904) q[1];
sx q[1];
rz(2.6516857) q[1];
rz(-pi) q[2];
rz(-0.34890449) q[3];
sx q[3];
rz(-1.5530905) q[3];
sx q[3];
rz(-0.68428333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6172341) q[2];
sx q[2];
rz(-1.989216) q[2];
sx q[2];
rz(2.1481245) q[2];
rz(-2.1067545) q[3];
sx q[3];
rz(-2.5351758) q[3];
sx q[3];
rz(-0.16916999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33191037) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(1.5274973) q[0];
rz(2.2612259) q[1];
sx q[1];
rz(-2.7208734) q[1];
sx q[1];
rz(-0.31164935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0904734) q[0];
sx q[0];
rz(-1.3731806) q[0];
sx q[0];
rz(-1.86905) q[0];
rz(-1.6402354) q[2];
sx q[2];
rz(-0.23190325) q[2];
sx q[2];
rz(-2.7551485) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.27585718) q[1];
sx q[1];
rz(-1.9364662) q[1];
sx q[1];
rz(2.4128298) q[1];
rz(-pi) q[2];
rz(1.1470471) q[3];
sx q[3];
rz(-2.4422788) q[3];
sx q[3];
rz(-0.55845234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6722022) q[2];
sx q[2];
rz(-0.92350525) q[2];
sx q[2];
rz(1.7507318) q[2];
rz(-1.5270799) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(2.0670149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020318) q[0];
sx q[0];
rz(-1.1447516) q[0];
sx q[0];
rz(0.61258739) q[0];
rz(-0.9901498) q[1];
sx q[1];
rz(-1.9691111) q[1];
sx q[1];
rz(-1.490907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25967071) q[0];
sx q[0];
rz(-0.90357354) q[0];
sx q[0];
rz(-1.0259969) q[0];
rz(2.2305626) q[2];
sx q[2];
rz(-1.3675895) q[2];
sx q[2];
rz(3.0958297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0066291) q[1];
sx q[1];
rz(-1.4093168) q[1];
sx q[1];
rz(0.0063185255) q[1];
rz(-2.8464523) q[3];
sx q[3];
rz(-0.71986976) q[3];
sx q[3];
rz(-2.2146378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.426173) q[2];
sx q[2];
rz(-2.3121068) q[2];
sx q[2];
rz(-0.034991525) q[2];
rz(1.4404826) q[3];
sx q[3];
rz(-0.8927497) q[3];
sx q[3];
rz(0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4424425) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(1.9433446) q[1];
sx q[1];
rz(-1.6234963) q[1];
sx q[1];
rz(-2.1877098) q[1];
rz(-2.162355) q[2];
sx q[2];
rz(-2.1605347) q[2];
sx q[2];
rz(0.84183358) q[2];
rz(-0.71975868) q[3];
sx q[3];
rz(-1.2133141) q[3];
sx q[3];
rz(2.6837466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
