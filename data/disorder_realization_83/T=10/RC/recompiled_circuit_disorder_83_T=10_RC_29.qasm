OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(2.3798556) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(-0.74365562) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0640472) q[0];
sx q[0];
rz(-2.8923058) q[0];
sx q[0];
rz(1.8403948) q[0];
rz(-pi) q[1];
rz(1.6913937) q[2];
sx q[2];
rz(-2.2199124) q[2];
sx q[2];
rz(2.803034) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2603862) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(1.1587515) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2067912) q[3];
sx q[3];
rz(-1.5036821) q[3];
sx q[3];
rz(1.2957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7068229) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(-0.14262959) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728535) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(-2.1570737) q[0];
rz(-pi) q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(-2.1737614) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.47536182) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(-1.0223785) q[1];
x q[2];
rz(1.9713692) q[3];
sx q[3];
rz(-1.2188606) q[3];
sx q[3];
rz(0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(-1.846116) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(2.4643262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5325461) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(-0.30058582) q[0];
rz(2.9708423) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(-3.1250172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9159569) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(-3.008328) q[1];
rz(-2.6636395) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(-1.705845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(-2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.9690008) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2617944) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(-1.7394702) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6060068) q[2];
sx q[2];
rz(-1.5835985) q[2];
sx q[2];
rz(-2.9738604) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9753032) q[1];
sx q[1];
rz(-1.9235833) q[1];
sx q[1];
rz(2.8217479) q[1];
x q[2];
rz(1.138932) q[3];
sx q[3];
rz(-1.6383088) q[3];
sx q[3];
rz(2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(-1.8160965) q[0];
rz(-2.5505113) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(2.4868734) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39010534) q[0];
sx q[0];
rz(-0.57919466) q[0];
sx q[0];
rz(-0.67436995) q[0];
x q[1];
rz(2.7388217) q[2];
sx q[2];
rz(-2.6908814) q[2];
sx q[2];
rz(2.7223301) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4795585) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(2.0550834) q[1];
rz(2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(1.4978131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70443557) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(0.066083834) q[0];
x q[1];
rz(2.1645855) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-1.1360816) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99243473) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(1.0013594) q[1];
rz(-pi) q[2];
rz(1.8201581) q[3];
sx q[3];
rz(-1.8421838) q[3];
sx q[3];
rz(-0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(3.1141172) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-0.79089975) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.233477) q[0];
sx q[0];
rz(-0.95196264) q[0];
sx q[0];
rz(0.11615642) q[0];
rz(-0.74366624) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(0.92168346) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0485718) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(-3.0572592) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0049719) q[3];
sx q[3];
rz(-2.4063111) q[3];
sx q[3];
rz(-2.6649464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(-2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62846047) q[0];
sx q[0];
rz(-1.5043133) q[0];
sx q[0];
rz(-2.9456338) q[0];
rz(-pi) q[1];
rz(-2.8587564) q[2];
sx q[2];
rz(-3.1361702) q[2];
sx q[2];
rz(1.1495513) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64995631) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(-1.4375356) q[1];
x q[2];
rz(3.0244175) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-2.1597247) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0197524) q[0];
sx q[0];
rz(-1.6822018) q[0];
sx q[0];
rz(0.41282546) q[0];
x q[1];
rz(-2.1397212) q[2];
sx q[2];
rz(-1.6621727) q[2];
sx q[2];
rz(-0.56359529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2541483) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(-1.7565586) q[1];
rz(-pi) q[2];
rz(1.8198265) q[3];
sx q[3];
rz(-0.8753652) q[3];
sx q[3];
rz(-1.4679366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(-2.1441933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.6249264) q[0];
rz(-2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(2.8732252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.836652) q[0];
sx q[0];
rz(-0.9965082) q[0];
sx q[0];
rz(2.4313297) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1860113) q[2];
sx q[2];
rz(-2.4590543) q[2];
sx q[2];
rz(-0.62894422) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.25071535) q[1];
sx q[1];
rz(-0.6579352) q[1];
sx q[1];
rz(0.4472181) q[1];
rz(-0.99276944) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(0.85152599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(-1.1901723) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-1.2673169) q[2];
sx q[2];
rz(-1.0721285) q[2];
sx q[2];
rz(-0.68821651) q[2];
rz(-2.5793544) q[3];
sx q[3];
rz(-2.7769965) q[3];
sx q[3];
rz(-2.1806352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
