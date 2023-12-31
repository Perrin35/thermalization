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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418158) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(-3.0738897) q[0];
rz(-pi) q[1];
rz(-0.6526297) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(1.8362311) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2603862) q[1];
sx q[1];
rz(-1.8477866) q[1];
sx q[1];
rz(1.1587515) q[1];
x q[2];
rz(-1.7573962) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649439) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-2.9108677) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-3.1006295) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728535) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(0.98451891) q[0];
x q[1];
rz(1.2626921) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(-1.2958796) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7296655) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(2.5018442) q[1];
rz(-pi) q[2];
x q[2];
rz(2.326194) q[3];
sx q[3];
rz(-0.5268464) q[3];
sx q[3];
rz(-0.53838733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(2.6290821) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(0.67726642) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2579502) q[0];
sx q[0];
rz(-1.5188688) q[0];
sx q[0];
rz(-2.9731263) q[0];
rz(2.9708423) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(0.016575459) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2689506) q[1];
sx q[1];
rz(-0.22224717) q[1];
sx q[1];
rz(0.93572576) q[1];
rz(-pi) q[2];
rz(2.6636395) q[3];
sx q[3];
rz(-1.2198997) q[3];
sx q[3];
rz(1.4357476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(0.046860524) q[2];
rz(-0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2617944) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(1.7394702) q[0];
rz(-1.5863717) q[2];
sx q[2];
rz(-2.1767463) q[2];
sx q[2];
rz(1.7473999) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6231411) q[1];
sx q[1];
rz(-1.271283) q[1];
sx q[1];
rz(1.9407942) q[1];
rz(-pi) q[2];
rz(-1.138932) q[3];
sx q[3];
rz(-1.5032839) q[3];
sx q[3];
rz(2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-0.65471929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39010534) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(-2.4672227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7227313) q[2];
sx q[2];
rz(-1.7423811) q[2];
sx q[2];
rz(-0.78531839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4795585) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(-2.0550834) q[1];
x q[2];
rz(2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.6437795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2346674) q[0];
sx q[0];
rz(-1.6229796) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97700714) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-2.0055111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.076768) q[1];
sx q[1];
rz(-1.8586564) q[1];
sx q[1];
rz(2.0516146) q[1];
rz(-pi) q[2];
rz(-1.8201581) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(2.4621778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(-2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-2.3506929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40490155) q[0];
sx q[0];
rz(-1.4762523) q[0];
sx q[0];
rz(-2.1928284) q[0];
x q[1];
rz(-1.2355455) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(-2.7170979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.447688) q[1];
sx q[1];
rz(-1.649592) q[1];
sx q[1];
rz(1.2055956) q[1];
rz(2.0049719) q[3];
sx q[3];
rz(-0.73528157) q[3];
sx q[3];
rz(-2.6649464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(-0.25407243) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62846047) q[0];
sx q[0];
rz(-1.6372794) q[0];
sx q[0];
rz(-0.19595887) q[0];
x q[1];
rz(-1.569283) q[2];
sx q[2];
rz(-1.5760033) q[2];
sx q[2];
rz(1.7092012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1797807) q[1];
sx q[1];
rz(-1.4439549) q[1];
sx q[1];
rz(0.3133874) q[1];
rz(-pi) q[2];
rz(-0.11717511) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(-0.52091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(-0.32661682) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(2.1597247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69753416) q[0];
sx q[0];
rz(-0.42675787) q[0];
sx q[0];
rz(0.27192893) q[0];
x q[1];
rz(1.7392731) q[2];
sx q[2];
rz(-0.57541621) q[2];
sx q[2];
rz(1.1489431) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6135785) q[1];
sx q[1];
rz(-0.19436969) q[1];
sx q[1];
rz(1.8730875) q[1];
x q[2];
rz(0.28717678) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(-1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(-0.99739933) q[2];
rz(0.50619566) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(0.51666623) q[0];
rz(-0.38756469) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70290138) q[0];
sx q[0];
rz(-2.1500906) q[0];
sx q[0];
rz(-0.86433522) q[0];
x q[1];
rz(2.1860113) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(2.5126484) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25071535) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(0.4472181) q[1];
x q[2];
rz(1.1770505) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(-0.39214373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(-2.6396991) q[2];
sx q[2];
rz(-2.5645651) q[2];
sx q[2];
rz(-1.2679451) q[2];
rz(2.8292538) q[3];
sx q[3];
rz(-1.3795508) q[3];
sx q[3];
rz(3.0637904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
