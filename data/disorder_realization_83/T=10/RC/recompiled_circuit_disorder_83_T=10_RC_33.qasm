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
rz(2.397937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24510358) q[0];
sx q[0];
rz(-1.5050383) q[0];
sx q[0];
rz(1.3301646) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15723575) q[2];
sx q[2];
rz(-0.6586282) q[2];
sx q[2];
rz(3.0008891) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.248677) q[1];
sx q[1];
rz(-2.6495653) q[1];
sx q[1];
rz(0.9534652) q[1];
rz(1.3841964) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(3.0649439) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(3.1006295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7777268) q[0];
sx q[0];
rz(-1.4085318) q[0];
sx q[0];
rz(1.8198387) q[0];
rz(-2.3054753) q[2];
sx q[2];
rz(-2.7364519) q[2];
sx q[2];
rz(-2.1737614) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47536182) q[1];
sx q[1];
rz(-1.0050217) q[1];
sx q[1];
rz(2.1192141) q[1];
rz(-pi) q[2];
rz(-2.7621272) q[3];
sx q[3];
rz(-1.1960408) q[3];
sx q[3];
rz(1.3665762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-2.5734148) q[2];
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
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-2.4643262) q[1];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.5718939) q[2];
sx q[2];
rz(-1.564431) q[2];
sx q[2];
rz(-0.15417834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3213768) q[1];
sx q[1];
rz(-1.4396588) q[1];
sx q[1];
rz(-1.3908435) q[1];
x q[2];
rz(1.1797818) q[3];
sx q[3];
rz(-1.1241594) q[3];
sx q[3];
rz(-0.041165813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(-0.051483367) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(1.4021224) q[0];
x q[1];
rz(1.5863717) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(-1.3941927) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6231411) q[1];
sx q[1];
rz(-1.271283) q[1];
sx q[1];
rz(-1.2007984) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.074313642) q[3];
sx q[3];
rz(-1.1399817) q[3];
sx q[3];
rz(-2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(-1.8160965) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(-0.65471929) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514873) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(-2.4672227) q[0];
rz(1.7582558) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(-2.2802441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3760738) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(0.54733025) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7563285) q[3];
sx q[3];
rz(-1.4841199) q[3];
sx q[3];
rz(-3.1298266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(-1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(1.4978131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90692524) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(-0.90953565) q[0];
rz(-pi) q[1];
rz(2.3408893) q[2];
sx q[2];
rz(-2.3869037) q[2];
sx q[2];
rz(-0.18037361) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1491579) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(2.1402332) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3214345) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(-0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(-2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(0.79089975) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7366911) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(2.1928284) q[0];
rz(-pi) q[1];
rz(2.7795243) q[2];
sx q[2];
rz(-2.3644591) q[2];
sx q[2];
rz(-2.2287378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6939047) q[1];
sx q[1];
rz(-1.649592) q[1];
sx q[1];
rz(-1.2055956) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88364717) q[3];
sx q[3];
rz(-1.8568608) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-0.86501914) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(-2.8875202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2124436) q[0];
sx q[0];
rz(-1.375276) q[0];
sx q[0];
rz(-1.6385727) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.569283) q[2];
sx q[2];
rz(-1.5760033) q[2];
sx q[2];
rz(-1.4323915) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64995631) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(1.7040571) q[1];
rz(3.0244175) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(2.8149758) q[3];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083387233) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-0.98186791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0197524) q[0];
sx q[0];
rz(-1.6822018) q[0];
sx q[0];
rz(-0.41282546) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0332546) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(2.1926751) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8874444) q[1];
sx q[1];
rz(-1.6283298) q[1];
sx q[1];
rz(-1.385034) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4307067) q[3];
sx q[3];
rz(-1.7611739) q[3];
sx q[3];
rz(0.058660942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2173569) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(2.6249264) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(-2.8732252) q[1];
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
rz(-pi) q[1];
rz(0.43865135) q[2];
sx q[2];
rz(-2.1119366) q[2];
sx q[2];
rz(-1.3676608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8908773) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(2.6943745) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26465613) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(2.2789126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(2.5777585) q[2];
rz(-1.1901723) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(1.8742758) q[2];
sx q[2];
rz(-1.0721285) q[2];
sx q[2];
rz(-0.68821651) q[2];
rz(-1.7715122) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];