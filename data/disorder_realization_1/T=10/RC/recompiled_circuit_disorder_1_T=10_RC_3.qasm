OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(1.948184) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0573187) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(2.1291332) q[0];
x q[1];
rz(1.869007) q[2];
sx q[2];
rz(-2.0993005) q[2];
sx q[2];
rz(2.3117711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8661583) q[1];
sx q[1];
rz(-0.94238102) q[1];
sx q[1];
rz(-0.98316146) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0317806) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(0.11085489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6821735) q[2];
sx q[2];
rz(-0.023962263) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(-1.4131644) q[3];
sx q[3];
rz(-1.1011522) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5380149) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(-1.0128101) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-2.467997) q[1];
sx q[1];
rz(-2.0181296) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22966188) q[0];
sx q[0];
rz(-2.6768885) q[0];
sx q[0];
rz(-1.6994516) q[0];
rz(-1.6337109) q[2];
sx q[2];
rz(-2.3496369) q[2];
sx q[2];
rz(-2.6346249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9065735) q[1];
sx q[1];
rz(-1.0636914) q[1];
sx q[1];
rz(-0.97096918) q[1];
rz(-0.65621891) q[3];
sx q[3];
rz(-1.1074293) q[3];
sx q[3];
rz(-0.046063395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.2751689) q[0];
rz(2.4480942) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(2.0085874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1904859) q[0];
sx q[0];
rz(-1.5674942) q[0];
sx q[0];
rz(1.7130873) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0622382) q[2];
sx q[2];
rz(-0.8478176) q[2];
sx q[2];
rz(1.6194956) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36724597) q[1];
sx q[1];
rz(-0.62421747) q[1];
sx q[1];
rz(1.063785) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5127701) q[3];
sx q[3];
rz(-1.0477133) q[3];
sx q[3];
rz(0.76997013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2514078) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.2934925) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(1.8815276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816198) q[0];
sx q[0];
rz(-0.078475229) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(0.89598957) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-3.0060351) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4281222) q[0];
sx q[0];
rz(-0.8849511) q[0];
sx q[0];
rz(-1.1629521) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1866456) q[2];
sx q[2];
rz(-1.640056) q[2];
sx q[2];
rz(0.70089507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23952661) q[1];
sx q[1];
rz(-0.24589989) q[1];
sx q[1];
rz(-0.86218254) q[1];
x q[2];
rz(-2.0849864) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(-1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(-2.2616852) q[2];
rz(0.044163477) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(0.28863171) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-1.0132382) q[0];
rz(3.0918616) q[1];
sx q[1];
rz(-0.91369349) q[1];
sx q[1];
rz(2.0577046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5729382) q[0];
sx q[0];
rz(-1.7540115) q[0];
sx q[0];
rz(-1.818548) q[0];
rz(-pi) q[1];
rz(2.8903557) q[2];
sx q[2];
rz(-1.305797) q[2];
sx q[2];
rz(-2.3064409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10478445) q[1];
sx q[1];
rz(-1.0462865) q[1];
sx q[1];
rz(3.0166237) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.089133457) q[3];
sx q[3];
rz(-2.6260178) q[3];
sx q[3];
rz(-2.6775132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-0.24442913) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(-2.6385245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2844834) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-0.094141468) q[0];
rz(2.969818) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0953513) q[0];
sx q[0];
rz(-1.6086676) q[0];
sx q[0];
rz(0.34237679) q[0];
x q[1];
rz(-1.8108098) q[2];
sx q[2];
rz(-0.95108205) q[2];
sx q[2];
rz(2.0680708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9858866) q[1];
sx q[1];
rz(-2.967318) q[1];
sx q[1];
rz(-1.3143015) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5395245) q[3];
sx q[3];
rz(-0.074723738) q[3];
sx q[3];
rz(2.5089695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.133693) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-0.80319476) q[2];
rz(1.1903654) q[3];
sx q[3];
rz(-1.9093711) q[3];
sx q[3];
rz(0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.068709277) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(-1.8849467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.885868) q[0];
sx q[0];
rz(-0.2547383) q[0];
sx q[0];
rz(1.7314918) q[0];
rz(-2.0357382) q[2];
sx q[2];
rz(-1.7591811) q[2];
sx q[2];
rz(-0.18907324) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9659121) q[1];
sx q[1];
rz(-1.6351846) q[1];
sx q[1];
rz(0.28190159) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7729633) q[3];
sx q[3];
rz(-0.7437403) q[3];
sx q[3];
rz(1.7891235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98823035) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(1.777565) q[2];
rz(0.91056943) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751188) q[0];
sx q[0];
rz(-2.5771038) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(0.072487436) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(-2.7546308) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21803074) q[0];
sx q[0];
rz(-0.98919981) q[0];
sx q[0];
rz(0.27115718) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53051051) q[2];
sx q[2];
rz(-2.1973917) q[2];
sx q[2];
rz(-2.8811787) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4657198) q[1];
sx q[1];
rz(-0.67589251) q[1];
sx q[1];
rz(3.0216316) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7510795) q[3];
sx q[3];
rz(-0.53005855) q[3];
sx q[3];
rz(-2.183225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(-0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.4322759) q[3];
sx q[3];
rz(1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14116645) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(-2.4138342) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5987199) q[0];
sx q[0];
rz(-0.96519404) q[0];
sx q[0];
rz(2.4332895) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6463037) q[2];
sx q[2];
rz(-2.1914748) q[2];
sx q[2];
rz(-0.6492614) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94433632) q[1];
sx q[1];
rz(-1.1068871) q[1];
sx q[1];
rz(1.0524366) q[1];
rz(-2.4008972) q[3];
sx q[3];
rz(-1.9441838) q[3];
sx q[3];
rz(-1.5965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14770517) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(1.0797427) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-0.22288999) q[1];
sx q[1];
rz(1.7396897) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3831543) q[0];
sx q[0];
rz(-2.8935195) q[0];
sx q[0];
rz(-2.1092578) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18677588) q[2];
sx q[2];
rz(-1.597607) q[2];
sx q[2];
rz(1.9865799) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.460443) q[1];
sx q[1];
rz(-0.76421684) q[1];
sx q[1];
rz(-2.0185673) q[1];
rz(-pi) q[2];
rz(2.8614282) q[3];
sx q[3];
rz(-0.6503085) q[3];
sx q[3];
rz(0.10661099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6587276) q[2];
sx q[2];
rz(-1.3663224) q[2];
sx q[2];
rz(1.520291) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(0.6974535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(2.2254754) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(0.79402906) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(2.0879073) q[3];
sx q[3];
rz(-1.5862982) q[3];
sx q[3];
rz(1.3261212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
