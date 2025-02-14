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
rz(2.0514367) q[0];
sx q[0];
rz(3.2618599) q[0];
sx q[0];
rz(9.5016228) q[0];
rz(2.285217) q[1];
sx q[1];
rz(5.2353274) q[1];
sx q[1];
rz(8.6040001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4030006) q[0];
sx q[0];
rz(-2.0401017) q[0];
sx q[0];
rz(-1.655939) q[0];
rz(-pi) q[1];
rz(-2.2256741) q[2];
sx q[2];
rz(-1.4400463) q[2];
sx q[2];
rz(-0.98721013) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1917853) q[1];
sx q[1];
rz(-0.92647431) q[1];
sx q[1];
rz(0.61442805) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54475875) q[3];
sx q[3];
rz(-1.7228456) q[3];
sx q[3];
rz(1.461538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0277752) q[2];
sx q[2];
rz(-0.57089388) q[2];
sx q[2];
rz(1.5296193) q[2];
rz(-0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(-2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73794) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(1.2427484) q[0];
rz(-2.2620762) q[1];
sx q[1];
rz(-2.1915235) q[1];
sx q[1];
rz(-1.3290728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97016818) q[0];
sx q[0];
rz(-2.1546083) q[0];
sx q[0];
rz(1.0310571) q[0];
rz(-pi) q[1];
rz(1.0813367) q[2];
sx q[2];
rz(-2.0149517) q[2];
sx q[2];
rz(1.7760488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0893007) q[1];
sx q[1];
rz(-1.1049905) q[1];
sx q[1];
rz(-2.0897081) q[1];
rz(0.19402282) q[3];
sx q[3];
rz(-2.4388154) q[3];
sx q[3];
rz(-1.8234092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8326412) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(1.857117) q[2];
rz(0.66713157) q[3];
sx q[3];
rz(-0.37225538) q[3];
sx q[3];
rz(2.4767806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.2892147) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(-1.941823) q[0];
rz(1.86357) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(-2.3866381) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4830355) q[0];
sx q[0];
rz(-2.6847976) q[0];
sx q[0];
rz(1.5798588) q[0];
x q[1];
rz(-0.16904449) q[2];
sx q[2];
rz(-0.72578207) q[2];
sx q[2];
rz(-0.58453416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9782516) q[1];
sx q[1];
rz(-0.74865197) q[1];
sx q[1];
rz(1.6221552) q[1];
rz(2.4944896) q[3];
sx q[3];
rz(-0.80296338) q[3];
sx q[3];
rz(2.3784007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7167012) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(1.971604) q[2];
rz(-2.8814426) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(-0.20497504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78885704) q[0];
sx q[0];
rz(-1.641888) q[0];
sx q[0];
rz(0.66057551) q[0];
rz(-0.13064101) q[1];
sx q[1];
rz(-0.6159997) q[1];
sx q[1];
rz(-2.7052243) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6141967) q[0];
sx q[0];
rz(-2.578311) q[0];
sx q[0];
rz(-2.0296996) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5549468) q[2];
sx q[2];
rz(-1.9244951) q[2];
sx q[2];
rz(-2.2207137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4727488) q[1];
sx q[1];
rz(-1.9276697) q[1];
sx q[1];
rz(-2.1118739) q[1];
x q[2];
rz(1.4722669) q[3];
sx q[3];
rz(-1.9103582) q[3];
sx q[3];
rz(1.3738541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.074162) q[2];
sx q[2];
rz(-2.6219411) q[2];
sx q[2];
rz(-2.0036428) q[2];
rz(-0.84749311) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.725089) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013019) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(2.9345595) q[0];
rz(1.8197458) q[1];
sx q[1];
rz(-1.1489979) q[1];
sx q[1];
rz(2.102899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.032271) q[0];
sx q[0];
rz(-2.3242852) q[0];
sx q[0];
rz(-0.82839806) q[0];
rz(0.62049753) q[2];
sx q[2];
rz(-1.7417241) q[2];
sx q[2];
rz(-0.88384274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6398479) q[1];
sx q[1];
rz(-1.5540489) q[1];
sx q[1];
rz(0.63629758) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4624427) q[3];
sx q[3];
rz(-1.8757678) q[3];
sx q[3];
rz(1.0602575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.98628226) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(2.5347064) q[2];
rz(0.33489975) q[3];
sx q[3];
rz(-1.5523942) q[3];
sx q[3];
rz(1.5437532) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63810054) q[0];
sx q[0];
rz(-1.7146716) q[0];
sx q[0];
rz(-0.18216369) q[0];
rz(-2.4496574) q[1];
sx q[1];
rz(-1.4357932) q[1];
sx q[1];
rz(-1.2729794) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6516354) q[0];
sx q[0];
rz(-3.0916442) q[0];
sx q[0];
rz(1.7836216) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7164549) q[2];
sx q[2];
rz(-1.7985059) q[2];
sx q[2];
rz(0.53066724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0611736) q[1];
sx q[1];
rz(-2.3234832) q[1];
sx q[1];
rz(-1.2123176) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9732861) q[3];
sx q[3];
rz(-2.7790439) q[3];
sx q[3];
rz(-2.1358638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89340297) q[2];
sx q[2];
rz(-2.4290163) q[2];
sx q[2];
rz(-0.87598652) q[2];
rz(2.857483) q[3];
sx q[3];
rz(-2.1967389) q[3];
sx q[3];
rz(0.42377728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1856325) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(2.60485) q[0];
rz(0.56450379) q[1];
sx q[1];
rz(-0.91465488) q[1];
sx q[1];
rz(-1.5435262) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48115981) q[0];
sx q[0];
rz(-2.8114722) q[0];
sx q[0];
rz(-0.31417234) q[0];
rz(0.41588628) q[2];
sx q[2];
rz(-1.3817817) q[2];
sx q[2];
rz(-3.031446) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98087039) q[1];
sx q[1];
rz(-2.0288784) q[1];
sx q[1];
rz(1.0700109) q[1];
rz(2.2143977) q[3];
sx q[3];
rz(-2.4193086) q[3];
sx q[3];
rz(2.3658906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87968612) q[2];
sx q[2];
rz(-1.8330816) q[2];
sx q[2];
rz(1.7730664) q[2];
rz(-0.3174828) q[3];
sx q[3];
rz(-1.8902794) q[3];
sx q[3];
rz(-2.4449463) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19318652) q[0];
sx q[0];
rz(-0.53006154) q[0];
sx q[0];
rz(-1.0411221) q[0];
rz(2.5913473) q[1];
sx q[1];
rz(-2.7698275) q[1];
sx q[1];
rz(-2.9464029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0597417) q[0];
sx q[0];
rz(-1.4343111) q[0];
sx q[0];
rz(1.6249715) q[0];
rz(-pi) q[1];
rz(2.0521021) q[2];
sx q[2];
rz(-1.6688108) q[2];
sx q[2];
rz(-2.5070407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4976552) q[1];
sx q[1];
rz(-1.8574839) q[1];
sx q[1];
rz(-3.1189256) q[1];
rz(-pi) q[2];
rz(0.82542586) q[3];
sx q[3];
rz(-1.9147827) q[3];
sx q[3];
rz(1.0937364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0027085) q[2];
sx q[2];
rz(-2.0033074) q[2];
sx q[2];
rz(1.2646593) q[2];
rz(-1.5539315) q[3];
sx q[3];
rz(-1.9374266) q[3];
sx q[3];
rz(-1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
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
rz(0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(1.6881855) q[1];
sx q[1];
rz(-1.6994349) q[1];
sx q[1];
rz(-2.2834987) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5018117) q[0];
sx q[0];
rz(-2.1970891) q[0];
sx q[0];
rz(0.39299742) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4420322) q[2];
sx q[2];
rz(-1.6673281) q[2];
sx q[2];
rz(-1.0573204) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8112642) q[1];
sx q[1];
rz(-1.1712046) q[1];
sx q[1];
rz(0.43021221) q[1];
rz(-pi) q[2];
rz(-1.7210049) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(-2.5774389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66576362) q[2];
sx q[2];
rz(-2.2593311) q[2];
sx q[2];
rz(-2.8909454) q[2];
rz(2.2541239) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(0.34621507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86824054) q[0];
sx q[0];
rz(-0.98889416) q[0];
sx q[0];
rz(2.2762779) q[0];
rz(2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(1.685198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0883852) q[0];
sx q[0];
rz(-1.6376978) q[0];
sx q[0];
rz(-1.7663203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9695332) q[2];
sx q[2];
rz(-2.2340748) q[2];
sx q[2];
rz(-3.1236529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2715312) q[1];
sx q[1];
rz(-1.9215596) q[1];
sx q[1];
rz(-0.49141541) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51442149) q[3];
sx q[3];
rz(-1.2687917) q[3];
sx q[3];
rz(-0.3946886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3978079) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(-1.3249116) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(2.3769784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.926173) q[0];
sx q[0];
rz(-1.8126491) q[0];
sx q[0];
rz(1.5480315) q[0];
rz(0.46319766) q[1];
sx q[1];
rz(-1.6013655) q[1];
sx q[1];
rz(-0.2688437) q[1];
rz(-1.2590969) q[2];
sx q[2];
rz(-2.1915956) q[2];
sx q[2];
rz(-1.6804463) q[2];
rz(-0.67288053) q[3];
sx q[3];
rz(-1.4757446) q[3];
sx q[3];
rz(-0.12991005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
