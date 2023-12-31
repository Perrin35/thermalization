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
rz(-2.1589307) q[0];
sx q[0];
rz(0.76173705) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(-0.74365562) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418158) q[0];
sx q[0];
rz(-1.3306949) q[0];
sx q[0];
rz(-0.067702985) q[0];
rz(-pi) q[1];
rz(2.9843569) q[2];
sx q[2];
rz(-2.4829645) q[2];
sx q[2];
rz(-0.14070357) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88120645) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(1.1587515) q[1];
rz(-1.7573962) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(-0.44934011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(3.1006295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728535) q[0];
sx q[0];
rz(-0.29631796) q[0];
sx q[0];
rz(2.1570737) q[0];
rz(-pi) q[1];
rz(1.2626921) q[2];
sx q[2];
rz(-1.8381881) q[2];
sx q[2];
rz(-1.2958796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6662308) q[1];
sx q[1];
rz(-2.1365709) q[1];
sx q[1];
rz(1.0223785) q[1];
rz(-pi) q[2];
rz(0.81539865) q[3];
sx q[3];
rz(-0.5268464) q[3];
sx q[3];
rz(0.53838733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(-2.6290821) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
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
rz(-2.9653774) q[0];
sx q[0];
rz(-2.8410068) q[0];
rz(-pi) q[1];
rz(-0.0063653221) q[2];
sx q[2];
rz(-1.5696988) q[2];
sx q[2];
rz(1.7249677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22563572) q[1];
sx q[1];
rz(-1.3924053) q[1];
sx q[1];
rz(-3.008328) q[1];
rz(-pi) q[2];
rz(1.9618109) q[3];
sx q[3];
rz(-1.1241594) q[3];
sx q[3];
rz(0.041165813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(-0.81165195) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(2.9714382) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.9690008) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(-1.7394702) q[0];
rz(-0.02247359) q[2];
sx q[2];
rz(-2.5354676) q[2];
sx q[2];
rz(1.4215353) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9753032) q[1];
sx q[1];
rz(-1.9235833) q[1];
sx q[1];
rz(2.8217479) q[1];
rz(-1.138932) q[3];
sx q[3];
rz(-1.6383088) q[3];
sx q[3];
rz(-2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(-1.9481109) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(-1.8160965) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7703169) q[0];
sx q[0];
rz(-1.9195942) q[0];
sx q[0];
rz(2.6692997) q[0];
x q[1];
rz(1.3833369) q[2];
sx q[2];
rz(-1.9831295) q[2];
sx q[2];
rz(0.86134855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5419546) q[1];
sx q[1];
rz(-1.9497739) q[1];
sx q[1];
rz(-2.42948) q[1];
rz(-pi) q[2];
rz(-2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(3.1298266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(0.58978224) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(1.6437795) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59693903) q[0];
sx q[0];
rz(-0.66300809) q[0];
sx q[0];
rz(1.6556428) q[0];
rz(-pi) q[1];
rz(-0.57967474) q[2];
sx q[2];
rz(-2.0849166) q[2];
sx q[2];
rz(-2.3953953) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1491579) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(1.0013594) q[1];
rz(-pi) q[2];
rz(-0.27961126) q[3];
sx q[3];
rz(-1.3307443) q[3];
sx q[3];
rz(2.1820501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(3.1141172) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291173) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-0.11418848) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(0.79089975) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9081156) q[0];
sx q[0];
rz(-0.95196264) q[0];
sx q[0];
rz(0.11615642) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9060471) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(0.4244948) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0485718) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(3.0572592) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1366208) q[3];
sx q[3];
rz(-0.73528157) q[3];
sx q[3];
rz(0.47664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(-2.8875202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2124436) q[0];
sx q[0];
rz(-1.375276) q[0];
sx q[0];
rz(-1.50302) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0052070219) q[2];
sx q[2];
rz(-1.569283) q[2];
sx q[2];
rz(-0.13841275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.904774) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(-2.7493613) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11717511) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(-2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(2.8059778) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(0.98186791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1218402) q[0];
sx q[0];
rz(-1.4593908) q[0];
sx q[0];
rz(2.7287672) q[0];
rz(-pi) q[1];
rz(1.4023196) q[2];
sx q[2];
rz(-2.5661764) q[2];
sx q[2];
rz(-1.9926496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5280142) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(1.2685052) q[1];
rz(2.8544159) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(-1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(2.1441933) q[2];
rz(0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
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
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70290138) q[0];
sx q[0];
rz(-2.1500906) q[0];
sx q[0];
rz(-0.86433522) q[0];
x q[1];
rz(0.98476121) q[2];
sx q[2];
rz(-1.1981989) q[2];
sx q[2];
rz(-0.44024703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4587935) q[1];
sx q[1];
rz(-1.3031712) q[1];
sx q[1];
rz(0.60857472) q[1];
x q[2];
rz(1.1770505) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(2.7494489) q[3];
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
rz(-0.95919132) q[3];
sx q[3];
rz(2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.8019567) q[1];
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
