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
rz(2.0775454) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(1.8403948) q[0];
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
rz(-pi) q[0];
x q[0];
rz(0.88120645) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(1.1587515) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3841964) q[3];
sx q[3];
rz(-2.7717234) q[3];
sx q[3];
rz(-2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(-3.1006295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7728535) q[0];
sx q[0];
rz(-0.29631796) q[0];
sx q[0];
rz(-2.1570737) q[0];
x q[1];
rz(-0.27994056) q[2];
sx q[2];
rz(-1.2739812) q[2];
sx q[2];
rz(-2.9505626) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4119271) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(2.5018442) q[1];
x q[2];
rz(-1.9713692) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(-2.7924201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9791947) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(2.6913753) q[0];
rz(1.846116) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(2.4643262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2579502) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(-0.16846637) q[0];
rz(-pi) q[1];
rz(1.5696987) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(-2.9874143) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87264204) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(-0.93572576) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6636395) q[3];
sx q[3];
rz(-1.2198997) q[3];
sx q[3];
rz(1.705845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(-2.3299407) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(0.17015447) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.1725918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5347036) q[0];
sx q[0];
rz(-0.67236116) q[0];
sx q[0];
rz(0.21557233) q[0];
x q[1];
rz(0.02247359) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(-1.7200574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6231411) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(1.2007984) q[1];
x q[2];
rz(1.4106393) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(-0.48689294) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(2.5505113) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(-0.65471929) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3712758) q[0];
sx q[0];
rz(-1.2219984) q[0];
sx q[0];
rz(2.6692997) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7388217) q[2];
sx q[2];
rz(-0.45071128) q[2];
sx q[2];
rz(2.7223301) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4795585) q[1];
sx q[1];
rz(-0.91842945) q[1];
sx q[1];
rz(2.0550834) q[1];
rz(-pi) q[2];
rz(1.6642903) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-0.5711242) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(-1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.4978131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59693903) q[0];
sx q[0];
rz(-2.4785846) q[0];
sx q[0];
rz(1.4859499) q[0];
rz(-pi) q[1];
rz(-2.3408893) q[2];
sx q[2];
rz(-2.3869037) q[2];
sx q[2];
rz(-2.961219) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99243473) q[1];
sx q[1];
rz(-0.5545534) q[1];
sx q[1];
rz(2.1402332) q[1];
x q[2];
rz(1.3214345) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(-0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291173) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(0.79089975) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0349802) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(1.4094704) q[0];
rz(-2.7795243) q[2];
sx q[2];
rz(-0.77713359) q[2];
sx q[2];
rz(-2.2287378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6939047) q[1];
sx q[1];
rz(-1.649592) q[1];
sx q[1];
rz(1.9359971) q[1];
rz(-pi) q[2];
rz(-0.36356504) q[3];
sx q[3];
rz(-0.91655469) q[3];
sx q[3];
rz(-0.082106575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-0.86501914) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62846047) q[0];
sx q[0];
rz(-1.6372794) q[0];
sx q[0];
rz(2.9456338) q[0];
rz(-pi) q[1];
rz(-1.5723096) q[2];
sx q[2];
rz(-1.5655893) q[2];
sx q[2];
rz(-1.4323915) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(-0.39223139) q[1];
rz(-pi) q[2];
rz(1.4905606) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(2.762592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(-2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(0.98186791) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1218402) q[0];
sx q[0];
rz(-1.6822018) q[0];
sx q[0];
rz(0.41282546) q[0];
rz(-3.0332546) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(-0.94891753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6135785) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(-1.2685052) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3217661) q[3];
sx q[3];
rz(-2.2662275) q[3];
sx q[3];
rz(1.4679366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(2.1441933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(0.034328073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8440588) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(-0.78155078) q[0];
x q[1];
rz(-0.95558138) q[2];
sx q[2];
rz(-2.4590543) q[2];
sx q[2];
rz(-2.5126484) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.25071535) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(2.6943745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9645421) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-2.5777585) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
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
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(0.50189353) q[2];
sx q[2];
rz(-2.5645651) q[2];
sx q[2];
rz(-1.2679451) q[2];
rz(2.5793544) q[3];
sx q[3];
rz(-0.36459618) q[3];
sx q[3];
rz(0.96095745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
