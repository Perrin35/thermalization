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
rz(-0.51195872) q[0];
sx q[0];
rz(6.6867642) q[0];
sx q[0];
rz(9.4288958) q[0];
rz(-2.4601958) q[1];
sx q[1];
rz(-1.8713142) q[1];
sx q[1];
rz(-1.4407925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82872921) q[0];
sx q[0];
rz(-1.848319) q[0];
sx q[0];
rz(2.9636431) q[0];
rz(-2.3902293) q[2];
sx q[2];
rz(-1.3168502) q[2];
sx q[2];
rz(-0.51515173) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5914247) q[1];
sx q[1];
rz(-1.9240849) q[1];
sx q[1];
rz(2.4387906) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1295092) q[3];
sx q[3];
rz(-1.3971431) q[3];
sx q[3];
rz(-2.1326667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9422841) q[2];
sx q[2];
rz(-2.0950623) q[2];
sx q[2];
rz(0.47041565) q[2];
rz(-2.2453902) q[3];
sx q[3];
rz(-1.4678518) q[3];
sx q[3];
rz(1.5508274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50297058) q[0];
sx q[0];
rz(-2.6574385) q[0];
sx q[0];
rz(-2.7726987) q[0];
rz(0.46345723) q[1];
sx q[1];
rz(-1.2838485) q[1];
sx q[1];
rz(2.8772112) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1702829) q[0];
sx q[0];
rz(-1.4173039) q[0];
sx q[0];
rz(0.92043368) q[0];
rz(2.8910341) q[2];
sx q[2];
rz(-1.7224285) q[2];
sx q[2];
rz(0.79094584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0707339) q[1];
sx q[1];
rz(-1.8740591) q[1];
sx q[1];
rz(1.8766848) q[1];
x q[2];
rz(-0.12118487) q[3];
sx q[3];
rz(-2.3555992) q[3];
sx q[3];
rz(1.3898897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5276864) q[2];
sx q[2];
rz(-2.9759585) q[2];
sx q[2];
rz(1.3045093) q[2];
rz(-2.8218609) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(-0.50362292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9673135) q[0];
sx q[0];
rz(-2.9496084) q[0];
sx q[0];
rz(2.018003) q[0];
rz(-0.38791052) q[1];
sx q[1];
rz(-1.265641) q[1];
sx q[1];
rz(-0.33043114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215805) q[0];
sx q[0];
rz(-2.8224484) q[0];
sx q[0];
rz(-2.1681251) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0082012) q[2];
sx q[2];
rz(-0.61163227) q[2];
sx q[2];
rz(0.69452121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.53239218) q[1];
sx q[1];
rz(-0.48621854) q[1];
sx q[1];
rz(2.4052252) q[1];
rz(-pi) q[2];
x q[2];
rz(2.167424) q[3];
sx q[3];
rz(-1.0181352) q[3];
sx q[3];
rz(2.139758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.419751) q[2];
sx q[2];
rz(-1.1417049) q[2];
sx q[2];
rz(1.1775796) q[2];
rz(1.8925331) q[3];
sx q[3];
rz(-1.306059) q[3];
sx q[3];
rz(3.1031928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053442001) q[0];
sx q[0];
rz(-1.7120687) q[0];
sx q[0];
rz(-0.088726774) q[0];
rz(-2.7170722) q[1];
sx q[1];
rz(-0.71966925) q[1];
sx q[1];
rz(-1.4094062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8363894) q[0];
sx q[0];
rz(-0.38581784) q[0];
sx q[0];
rz(-2.0842537) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1698977) q[2];
sx q[2];
rz(-1.9847775) q[2];
sx q[2];
rz(-1.5091015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.030050412) q[1];
sx q[1];
rz(-2.4218028) q[1];
sx q[1];
rz(-2.5538302) q[1];
rz(-2.0982101) q[3];
sx q[3];
rz(-2.5549915) q[3];
sx q[3];
rz(-2.9709904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7523664) q[2];
sx q[2];
rz(-2.0615292) q[2];
sx q[2];
rz(2.7092095) q[2];
rz(0.45807517) q[3];
sx q[3];
rz(-1.658344) q[3];
sx q[3];
rz(1.0079916) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78419375) q[0];
sx q[0];
rz(-2.9965239) q[0];
sx q[0];
rz(2.8302637) q[0];
rz(1.9642824) q[1];
sx q[1];
rz(-1.177634) q[1];
sx q[1];
rz(0.48431531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70894079) q[0];
sx q[0];
rz(-1.6659032) q[0];
sx q[0];
rz(-1.4115788) q[0];
x q[1];
rz(-1.823678) q[2];
sx q[2];
rz(-2.8399979) q[2];
sx q[2];
rz(-2.4689134) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.53953417) q[1];
sx q[1];
rz(-1.9440398) q[1];
sx q[1];
rz(-2.9887524) q[1];
x q[2];
rz(-1.442904) q[3];
sx q[3];
rz(-1.4223243) q[3];
sx q[3];
rz(0.080303513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8208661) q[2];
sx q[2];
rz(-1.4297239) q[2];
sx q[2];
rz(2.7166264) q[2];
rz(-3.037437) q[3];
sx q[3];
rz(-2.7724373) q[3];
sx q[3];
rz(2.4764376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686309) q[0];
sx q[0];
rz(-2.5017128) q[0];
sx q[0];
rz(-0.1567008) q[0];
rz(-2.8796097) q[1];
sx q[1];
rz(-1.1527088) q[1];
sx q[1];
rz(1.6798457) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77685415) q[0];
sx q[0];
rz(-0.97524736) q[0];
sx q[0];
rz(2.9315346) q[0];
rz(-2.2672923) q[2];
sx q[2];
rz(-1.7477711) q[2];
sx q[2];
rz(2.1483473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9380381) q[1];
sx q[1];
rz(-2.9239836) q[1];
sx q[1];
rz(0.59534351) q[1];
rz(-2.3912638) q[3];
sx q[3];
rz(-1.273386) q[3];
sx q[3];
rz(0.70575037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4713952) q[2];
sx q[2];
rz(-1.4963701) q[2];
sx q[2];
rz(-2.4246598) q[2];
rz(2.9251621) q[3];
sx q[3];
rz(-1.2856893) q[3];
sx q[3];
rz(1.1855679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
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
rz(2.7316498) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(1.8753847) q[0];
rz(0.75824291) q[1];
sx q[1];
rz(-1.9634602) q[1];
sx q[1];
rz(-1.0479124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1850994) q[0];
sx q[0];
rz(-0.66732115) q[0];
sx q[0];
rz(2.064608) q[0];
x q[1];
rz(0.4518544) q[2];
sx q[2];
rz(-1.0596709) q[2];
sx q[2];
rz(0.035511927) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3057237) q[1];
sx q[1];
rz(-0.36986884) q[1];
sx q[1];
rz(-1.3772411) q[1];
rz(-pi) q[2];
rz(1.8147477) q[3];
sx q[3];
rz(-1.1317963) q[3];
sx q[3];
rz(-1.148996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3332112) q[2];
sx q[2];
rz(-2.1344678) q[2];
sx q[2];
rz(-3.132931) q[2];
rz(-2.1374785) q[3];
sx q[3];
rz(-0.4929556) q[3];
sx q[3];
rz(-2.139411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.778331) q[0];
sx q[0];
rz(-2.9982428) q[0];
sx q[0];
rz(-2.1891201) q[0];
rz(-1.6078423) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(1.0362157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8795963) q[0];
sx q[0];
rz(-0.72234266) q[0];
sx q[0];
rz(-2.45558) q[0];
rz(-pi) q[1];
rz(0.63663738) q[2];
sx q[2];
rz(-1.271406) q[2];
sx q[2];
rz(-2.4518397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1324824) q[1];
sx q[1];
rz(-2.3308737) q[1];
sx q[1];
rz(-2.2894163) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1187068) q[3];
sx q[3];
rz(-2.2452659) q[3];
sx q[3];
rz(1.7728286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2967534) q[2];
sx q[2];
rz(-2.4256458) q[2];
sx q[2];
rz(-2.1412264) q[2];
rz(-1.2760466) q[3];
sx q[3];
rz(-1.6701271) q[3];
sx q[3];
rz(-2.0744417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.100383) q[0];
sx q[0];
rz(-0.90827933) q[0];
sx q[0];
rz(2.8874183) q[0];
rz(1.3379478) q[1];
sx q[1];
rz(-0.86054069) q[1];
sx q[1];
rz(-2.3635704) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36150071) q[0];
sx q[0];
rz(-0.92467148) q[0];
sx q[0];
rz(-1.161399) q[0];
rz(-pi) q[1];
rz(-2.4619589) q[2];
sx q[2];
rz(-0.25654116) q[2];
sx q[2];
rz(-2.7367965) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3624448) q[1];
sx q[1];
rz(-1.0082642) q[1];
sx q[1];
rz(-2.0570807) q[1];
rz(-pi) q[2];
rz(1.3555525) q[3];
sx q[3];
rz(-1.6674893) q[3];
sx q[3];
rz(-1.7413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.079387) q[2];
sx q[2];
rz(-1.5176682) q[2];
sx q[2];
rz(-0.13266955) q[2];
rz(1.1072055) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(1.6975105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4001813) q[0];
sx q[0];
rz(-1.9453229) q[0];
sx q[0];
rz(-0.74914002) q[0];
rz(-2.6465042) q[1];
sx q[1];
rz(-1.2352713) q[1];
sx q[1];
rz(1.7503768) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2783954) q[0];
sx q[0];
rz(-2.1957198) q[0];
sx q[0];
rz(1.004659) q[0];
rz(0.83413627) q[2];
sx q[2];
rz(-2.5747402) q[2];
sx q[2];
rz(-1.8994272) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5914982) q[1];
sx q[1];
rz(-1.6143394) q[1];
sx q[1];
rz(-2.7344879) q[1];
rz(-pi) q[2];
rz(2.0999317) q[3];
sx q[3];
rz(-1.1253998) q[3];
sx q[3];
rz(-2.9187849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5007925) q[2];
sx q[2];
rz(-1.4747138) q[2];
sx q[2];
rz(2.2643209) q[2];
rz(-2.6194465) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(1.6845901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40182879) q[0];
sx q[0];
rz(-0.54720989) q[0];
sx q[0];
rz(-1.3450958) q[0];
rz(1.2783891) q[1];
sx q[1];
rz(-2.8291193) q[1];
sx q[1];
rz(3.1183174) q[1];
rz(1.6017492) q[2];
sx q[2];
rz(-2.2014115) q[2];
sx q[2];
rz(1.9104107) q[2];
rz(-2.6865339) q[3];
sx q[3];
rz(-1.9225592) q[3];
sx q[3];
rz(-0.98071702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
