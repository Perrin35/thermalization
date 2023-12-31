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
rz(2.397937) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640472) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(1.8403948) q[0];
rz(-pi) q[1];
x q[1];
rz(1.450199) q[2];
sx q[2];
rz(-0.92168027) q[2];
sx q[2];
rz(-0.33855864) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88120645) q[1];
sx q[1];
rz(-1.8477866) q[1];
sx q[1];
rz(-1.9828412) q[1];
rz(-1.7573962) q[3];
sx q[3];
rz(-2.7717234) q[3];
sx q[3];
rz(0.44934011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(2.4690348) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(0.040963106) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7777268) q[0];
sx q[0];
rz(-1.4085318) q[0];
sx q[0];
rz(-1.3217539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8789005) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(1.2958796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47536182) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(-2.1192141) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7621272) q[3];
sx q[3];
rz(-1.1960408) q[3];
sx q[3];
rz(1.3665762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(-2.6290821) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(2.6913753) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(2.4643262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8375741) q[0];
sx q[0];
rz(-1.4025592) q[0];
sx q[0];
rz(1.5181245) q[0];
x q[1];
rz(1.5718939) q[2];
sx q[2];
rz(-1.564431) q[2];
sx q[2];
rz(0.15417834) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87264204) q[1];
sx q[1];
rz(-0.22224717) q[1];
sx q[1];
rz(-0.93572576) q[1];
rz(-1.9618109) q[3];
sx q[3];
rz(-2.0174332) q[3];
sx q[3];
rz(0.041165813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.9690008) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6068891) q[0];
sx q[0];
rz(-0.67236116) q[0];
sx q[0];
rz(0.21557233) q[0];
rz(2.5355859) q[2];
sx q[2];
rz(-1.5835985) q[2];
sx q[2];
rz(2.9738604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.40201515) q[1];
sx q[1];
rz(-2.6699454) q[1];
sx q[1];
rz(-0.86400835) q[1];
rz(-pi) q[2];
rz(-1.138932) q[3];
sx q[3];
rz(-1.5032839) q[3];
sx q[3];
rz(-0.95089144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-0.48689294) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(0.025432767) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-2.4868734) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3712758) q[0];
sx q[0];
rz(-1.2219984) q[0];
sx q[0];
rz(-2.6692997) q[0];
rz(-pi) q[1];
rz(-1.3833369) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(-2.2802441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5419546) q[1];
sx q[1];
rz(-1.1918187) q[1];
sx q[1];
rz(-2.42948) q[1];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.9545385) q[3];
sx q[3];
rz(-1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(3.0736198) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2346674) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57967474) q[2];
sx q[2];
rz(-2.0849166) q[2];
sx q[2];
rz(-0.74619734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0648246) q[1];
sx q[1];
rz(-1.8586564) q[1];
sx q[1];
rz(1.0899781) q[1];
rz(-pi) q[2];
rz(1.8201581) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(-0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.41247535) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(2.3506929) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.233477) q[0];
sx q[0];
rz(-2.18963) q[0];
sx q[0];
rz(-0.11615642) q[0];
rz(-1.2355455) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(0.4244948) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0485718) q[1];
sx q[1];
rz(-1.934811) q[1];
sx q[1];
rz(0.084333468) q[1];
x q[2];
rz(0.36356504) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(-0.082106575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5131322) q[0];
sx q[0];
rz(-1.5043133) q[0];
sx q[0];
rz(2.9456338) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28283624) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(1.9920414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(-2.7493613) q[1];
rz(-pi) q[2];
rz(3.0244175) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(2.1597247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0197524) q[0];
sx q[0];
rz(-1.6822018) q[0];
sx q[0];
rz(2.7287672) q[0];
rz(-pi) q[1];
rz(-2.1397212) q[2];
sx q[2];
rz(-1.6621727) q[2];
sx q[2];
rz(2.5779974) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2541483) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(-1.7565586) q[1];
x q[2];
rz(-1.3217661) q[3];
sx q[3];
rz(-0.8753652) q[3];
sx q[3];
rz(-1.4679366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2173569) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(-2.754028) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-2.8732252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3049406) q[0];
sx q[0];
rz(-0.9965082) q[0];
sx q[0];
rz(-0.71026295) q[0];
rz(-0.95558138) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(-0.62894422) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8908773) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(2.6943745) q[1];
rz(-pi) q[2];
rz(0.99276944) q[3];
sx q[3];
rz(-1.7939995) q[3];
sx q[3];
rz(-2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(2.5777585) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(0.64175516) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(1.8019567) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(-1.8742758) q[2];
sx q[2];
rz(-2.0694642) q[2];
sx q[2];
rz(2.4533761) q[2];
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
