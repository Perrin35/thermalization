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
rz(-2.9094568) q[0];
sx q[0];
rz(-2.870626) q[0];
sx q[0];
rz(-1.969307) q[0];
rz(-1.6074033) q[1];
sx q[1];
rz(-1.489403) q[1];
sx q[1];
rz(2.593427) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0294103) q[0];
sx q[0];
rz(-1.8444841) q[0];
sx q[0];
rz(-3.0269483) q[0];
x q[1];
rz(2.1287969) q[2];
sx q[2];
rz(-2.0145973) q[2];
sx q[2];
rz(1.2324126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2506574) q[1];
sx q[1];
rz(-1.0708628) q[1];
sx q[1];
rz(-1.2229162) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1211419) q[3];
sx q[3];
rz(-1.5142875) q[3];
sx q[3];
rz(2.5380937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9008909) q[2];
sx q[2];
rz(-1.2071995) q[2];
sx q[2];
rz(1.5473676) q[2];
rz(-2.2089925) q[3];
sx q[3];
rz(-0.91351944) q[3];
sx q[3];
rz(2.3368321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.7050183) q[0];
sx q[0];
rz(-2.5387634) q[0];
sx q[0];
rz(0.21827179) q[0];
rz(-1.6328579) q[1];
sx q[1];
rz(-0.873133) q[1];
sx q[1];
rz(-2.0192718) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86271226) q[0];
sx q[0];
rz(-0.76970184) q[0];
sx q[0];
rz(0.93602009) q[0];
x q[1];
rz(-1.3706742) q[2];
sx q[2];
rz(-0.84017838) q[2];
sx q[2];
rz(-0.49337988) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52294815) q[1];
sx q[1];
rz(-2.1174333) q[1];
sx q[1];
rz(2.1039822) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25457766) q[3];
sx q[3];
rz(-1.2647332) q[3];
sx q[3];
rz(-1.7382415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0108769) q[2];
sx q[2];
rz(-1.1900095) q[2];
sx q[2];
rz(0.38401628) q[2];
rz(-0.555641) q[3];
sx q[3];
rz(-2.5540387) q[3];
sx q[3];
rz(-2.2491992) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4246849) q[0];
sx q[0];
rz(-0.45775828) q[0];
sx q[0];
rz(2.6440788) q[0];
rz(-2.6566907) q[1];
sx q[1];
rz(-1.5919911) q[1];
sx q[1];
rz(0.44152322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2427485) q[0];
sx q[0];
rz(-1.4024807) q[0];
sx q[0];
rz(0.90710137) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21887987) q[2];
sx q[2];
rz(-2.0538123) q[2];
sx q[2];
rz(0.85798664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.42558266) q[1];
sx q[1];
rz(-1.0908608) q[1];
sx q[1];
rz(-2.369057) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2763292) q[3];
sx q[3];
rz(-0.40832106) q[3];
sx q[3];
rz(0.43864076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8170844) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(-0.95477742) q[2];
rz(-3.1214516) q[3];
sx q[3];
rz(-2.3432799) q[3];
sx q[3];
rz(-0.32625833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.7666053) q[0];
sx q[0];
rz(-2.6191481) q[0];
sx q[0];
rz(0.0047542714) q[0];
rz(-0.44113723) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(-1.7316679) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4213181) q[0];
sx q[0];
rz(-1.7022812) q[0];
sx q[0];
rz(-1.907503) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27556117) q[2];
sx q[2];
rz(-2.2803734) q[2];
sx q[2];
rz(1.0312652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9121352) q[1];
sx q[1];
rz(-1.3452824) q[1];
sx q[1];
rz(2.0769801) q[1];
x q[2];
rz(-1.4291841) q[3];
sx q[3];
rz(-0.19053121) q[3];
sx q[3];
rz(-1.9857565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8498174) q[2];
sx q[2];
rz(-1.3128277) q[2];
sx q[2];
rz(2.4449091) q[2];
rz(-1.6126532) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(-0.67970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1729537) q[0];
sx q[0];
rz(-2.8346859) q[0];
sx q[0];
rz(-0.90721834) q[0];
rz(0.13547678) q[1];
sx q[1];
rz(-1.2716581) q[1];
sx q[1];
rz(-2.4438593) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0402148) q[0];
sx q[0];
rz(-1.4411949) q[0];
sx q[0];
rz(1.1762397) q[0];
rz(-pi) q[1];
rz(-1.5595857) q[2];
sx q[2];
rz(-0.19489842) q[2];
sx q[2];
rz(0.29981183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6271883) q[1];
sx q[1];
rz(-1.3980734) q[1];
sx q[1];
rz(-2.808213) q[1];
rz(0.99643965) q[3];
sx q[3];
rz(-0.054704156) q[3];
sx q[3];
rz(2.8131054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1873143) q[2];
sx q[2];
rz(-2.4772187) q[2];
sx q[2];
rz(-2.7962255) q[2];
rz(2.6464388) q[3];
sx q[3];
rz(-2.5457355) q[3];
sx q[3];
rz(-3.0143747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0536163) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(-1.9867058) q[0];
rz(-0.79479533) q[1];
sx q[1];
rz(-1.2966917) q[1];
sx q[1];
rz(0.48387873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6538972) q[0];
sx q[0];
rz(-1.1654108) q[0];
sx q[0];
rz(0.8542819) q[0];
rz(1.0740499) q[2];
sx q[2];
rz(-1.4306465) q[2];
sx q[2];
rz(1.4432009) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9660873) q[1];
sx q[1];
rz(-2.6350807) q[1];
sx q[1];
rz(-1.1706074) q[1];
x q[2];
rz(-1.9779786) q[3];
sx q[3];
rz(-1.0308427) q[3];
sx q[3];
rz(-2.0288717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7530219) q[2];
sx q[2];
rz(-1.2349671) q[2];
sx q[2];
rz(1.1080326) q[2];
rz(-1.5293416) q[3];
sx q[3];
rz(-0.11950167) q[3];
sx q[3];
rz(0.39335462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.954708) q[0];
sx q[0];
rz(-2.055838) q[0];
sx q[0];
rz(-2.4672274) q[0];
rz(-0.87633324) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(-0.032329917) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5139149) q[0];
sx q[0];
rz(-1.6925294) q[0];
sx q[0];
rz(-1.4343778) q[0];
rz(-2.9170119) q[2];
sx q[2];
rz(-1.4201284) q[2];
sx q[2];
rz(1.2008787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9738071) q[1];
sx q[1];
rz(-1.9085437) q[1];
sx q[1];
rz(0.54904292) q[1];
rz(-2.9495839) q[3];
sx q[3];
rz(-0.98658326) q[3];
sx q[3];
rz(-1.8580798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9548367) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(-2.1046861) q[2];
rz(-1.9131276) q[3];
sx q[3];
rz(-0.90932536) q[3];
sx q[3];
rz(-2.9508446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4036338) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(2.8926358) q[0];
rz(0.20949334) q[1];
sx q[1];
rz(-1.341235) q[1];
sx q[1];
rz(0.6627717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25191987) q[0];
sx q[0];
rz(-3.126069) q[0];
sx q[0];
rz(0.70308103) q[0];
rz(-1.0900108) q[2];
sx q[2];
rz(-2.2741246) q[2];
sx q[2];
rz(-2.472773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31964916) q[1];
sx q[1];
rz(-1.3256097) q[1];
sx q[1];
rz(-0.82464928) q[1];
rz(-0.39318496) q[3];
sx q[3];
rz(-2.4454456) q[3];
sx q[3];
rz(-2.7659211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9111619) q[2];
sx q[2];
rz(-1.5713567) q[2];
sx q[2];
rz(-0.39361185) q[2];
rz(-2.7502934) q[3];
sx q[3];
rz(-0.53520447) q[3];
sx q[3];
rz(-0.75214255) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1289718) q[0];
sx q[0];
rz(-0.18652815) q[0];
sx q[0];
rz(0.57998002) q[0];
rz(2.773556) q[1];
sx q[1];
rz(-2.3263704) q[1];
sx q[1];
rz(-3.0267402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57844394) q[0];
sx q[0];
rz(-0.73860335) q[0];
sx q[0];
rz(-1.369801) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0124577) q[2];
sx q[2];
rz(-1.7071305) q[2];
sx q[2];
rz(1.1669056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8424994) q[1];
sx q[1];
rz(-0.4173435) q[1];
sx q[1];
rz(-2.4853766) q[1];
rz(-pi) q[2];
rz(0.87646342) q[3];
sx q[3];
rz(-1.254458) q[3];
sx q[3];
rz(2.9522459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7241235) q[2];
sx q[2];
rz(-2.8026411) q[2];
sx q[2];
rz(0.69619703) q[2];
rz(-2.0255069) q[3];
sx q[3];
rz(-0.79056549) q[3];
sx q[3];
rz(2.4000786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29257947) q[0];
sx q[0];
rz(-0.96939033) q[0];
sx q[0];
rz(-2.868929) q[0];
rz(0.69372454) q[1];
sx q[1];
rz(-1.1169746) q[1];
sx q[1];
rz(-0.56347096) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.97992) q[0];
sx q[0];
rz(-1.663999) q[0];
sx q[0];
rz(0.037089238) q[0];
rz(-pi) q[1];
rz(-0.54948893) q[2];
sx q[2];
rz(-2.8656883) q[2];
sx q[2];
rz(2.9018108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.882) q[1];
sx q[1];
rz(-1.881885) q[1];
sx q[1];
rz(0.085373665) q[1];
rz(-pi) q[2];
rz(-0.93592398) q[3];
sx q[3];
rz(-0.96158035) q[3];
sx q[3];
rz(-3.0425231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1571265) q[2];
sx q[2];
rz(-2.0968585) q[2];
sx q[2];
rz(-2.5926479) q[2];
rz(2.1460311) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(-2.5115749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41848771) q[0];
sx q[0];
rz(-0.99011078) q[0];
sx q[0];
rz(-0.44575442) q[0];
rz(0.86722974) q[1];
sx q[1];
rz(-2.0384616) q[1];
sx q[1];
rz(-1.9834317) q[1];
rz(1.8015447) q[2];
sx q[2];
rz(-1.1044772) q[2];
sx q[2];
rz(0.92369631) q[2];
rz(-0.54789644) q[3];
sx q[3];
rz(-1.4461645) q[3];
sx q[3];
rz(-0.89492284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
