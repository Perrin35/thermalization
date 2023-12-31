OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(4.1242546) q[0];
sx q[0];
rz(10.186515) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(-0.74365562) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0775454) q[0];
sx q[0];
rz(-2.8923058) q[0];
sx q[0];
rz(-1.8403948) q[0];
rz(-pi) q[1];
rz(-2.488963) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(1.3053615) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8929157) q[1];
sx q[1];
rz(-0.49202737) q[1];
sx q[1];
rz(0.9534652) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2067912) q[3];
sx q[3];
rz(-1.6379106) q[3];
sx q[3];
rz(1.8458927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-3.1006295) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7777268) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(1.8198387) q[0];
x q[1];
rz(2.8616521) q[2];
sx q[2];
rz(-1.2739812) q[2];
sx q[2];
rz(-2.9505626) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37521024) q[1];
sx q[1];
rz(-0.76645215) q[1];
sx q[1];
rz(-0.68739989) q[1];
x q[2];
rz(-2.7621272) q[3];
sx q[3];
rz(-1.1960408) q[3];
sx q[3];
rz(-1.7750164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.162398) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(0.67726642) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30401858) q[0];
sx q[0];
rz(-1.4025592) q[0];
sx q[0];
rz(-1.6234682) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5696987) q[2];
sx q[2];
rz(-1.564431) q[2];
sx q[2];
rz(-0.15417834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2689506) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(0.93572576) q[1];
x q[2];
rz(-2.4694091) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-0.046860524) q[2];
rz(-2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-0.051483367) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.9690008) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(-1.4021224) q[0];
x q[1];
rz(-3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(1.4215353) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7395775) q[1];
sx q[1];
rz(-2.6699454) q[1];
sx q[1];
rz(-0.86400835) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7309534) q[3];
sx q[3];
rz(-0.43678108) q[3];
sx q[3];
rz(-2.3763451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.8160965) q[0];
rz(-2.5505113) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-0.65471929) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7514873) q[0];
sx q[0];
rz(-0.57919466) q[0];
sx q[0];
rz(-2.4672227) q[0];
x q[1];
rz(1.7582558) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(0.86134855) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66203413) q[1];
sx q[1];
rz(-2.2231632) q[1];
sx q[1];
rz(-1.0865092) q[1];
rz(-pi) q[2];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.4978131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59693903) q[0];
sx q[0];
rz(-2.4785846) q[0];
sx q[0];
rz(-1.4859499) q[0];
rz(-pi) q[1];
rz(2.1645855) q[2];
sx q[2];
rz(-1.0736246) q[2];
sx q[2];
rz(1.1360816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0648246) q[1];
sx q[1];
rz(-1.2829363) q[1];
sx q[1];
rz(-1.0899781) q[1];
rz(2.4160556) q[3];
sx q[3];
rz(-0.36645884) q[3];
sx q[3];
rz(3.0612502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(-2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(2.3506929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0349802) q[0];
sx q[0];
rz(-0.62823409) q[0];
sx q[0];
rz(-1.7321222) q[0];
rz(-pi) q[1];
rz(-2.3979264) q[2];
sx q[2];
rz(-1.8218092) q[2];
sx q[2];
rz(0.92168346) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8154527) q[1];
sx q[1];
rz(-0.37323144) q[1];
sx q[1];
rz(-1.7883854) q[1];
rz(-pi) q[2];
rz(-0.88364717) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(-1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2124436) q[0];
sx q[0];
rz(-1.7663167) q[0];
sx q[0];
rz(1.6385727) q[0];
rz(-pi) q[1];
rz(1.5723096) q[2];
sx q[2];
rz(-1.5655893) q[2];
sx q[2];
rz(-1.7092012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4916363) q[1];
sx q[1];
rz(-1.8815814) q[1];
sx q[1];
rz(-1.4375356) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59928008) q[3];
sx q[3];
rz(-1.6370956) q[3];
sx q[3];
rz(1.146572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-0.33561486) q[2];
rz(-0.32661682) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(2.1597247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40030038) q[0];
sx q[0];
rz(-1.9809082) q[0];
sx q[0];
rz(1.6923231) q[0];
x q[1];
rz(-1.4023196) q[2];
sx q[2];
rz(-2.5661764) q[2];
sx q[2];
rz(1.9926496) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2541483) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(1.385034) q[1];
rz(2.8544159) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(-1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(-2.1441933) q[2];
rz(-0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(0.51666623) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(2.8732252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8440588) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(-0.78155078) q[0];
rz(2.7029413) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(1.7739319) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8473709) q[1];
sx q[1];
rz(-2.1547744) q[1];
sx q[1];
rz(1.2482615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26465613) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(0.86268007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664809) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(1.8742758) q[2];
sx q[2];
rz(-1.0721285) q[2];
sx q[2];
rz(-0.68821651) q[2];
rz(0.56223829) q[3];
sx q[3];
rz(-2.7769965) q[3];
sx q[3];
rz(-2.1806352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
