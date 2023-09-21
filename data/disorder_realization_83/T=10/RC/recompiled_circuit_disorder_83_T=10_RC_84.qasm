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
rz(1.3418158) q[0];
sx q[0];
rz(-1.3306949) q[0];
sx q[0];
rz(0.067702985) q[0];
x q[1];
rz(1.6913937) q[2];
sx q[2];
rz(-2.2199124) q[2];
sx q[2];
rz(-0.33855864) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88120645) q[1];
sx q[1];
rz(-1.8477866) q[1];
sx q[1];
rz(-1.9828412) q[1];
rz(-pi) q[2];
rz(-1.3841964) q[3];
sx q[3];
rz(-2.7717234) q[3];
sx q[3];
rz(2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7068229) q[2];
sx q[2];
rz(-0.40199026) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(0.23072492) q[0];
rz(-1.8143066) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(-0.040963106) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9757257) q[0];
sx q[0];
rz(-1.8164993) q[0];
sx q[0];
rz(2.9742572) q[0];
x q[1];
rz(0.27994056) q[2];
sx q[2];
rz(-1.8676114) q[2];
sx q[2];
rz(-2.9505626) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7296655) q[1];
sx q[1];
rz(-2.0265059) q[1];
sx q[1];
rz(2.5018442) q[1];
x q[2];
rz(0.37946545) q[3];
sx q[3];
rz(-1.9455519) q[3];
sx q[3];
rz(-1.3665762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(2.5734148) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(0.67726642) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2579502) q[0];
sx q[0];
rz(-1.5188688) q[0];
sx q[0];
rz(0.16846637) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9708423) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(0.016575459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87264204) q[1];
sx q[1];
rz(-0.22224717) q[1];
sx q[1];
rz(-0.93572576) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1797818) q[3];
sx q[3];
rz(-1.1241594) q[3];
sx q[3];
rz(3.1004268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(0.046860524) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(1.9690008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5347036) q[0];
sx q[0];
rz(-0.67236116) q[0];
sx q[0];
rz(0.21557233) q[0];
rz(1.555221) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(-1.7473999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40201515) q[1];
sx q[1];
rz(-2.6699454) q[1];
sx q[1];
rz(-2.2775843) q[1];
rz(-pi) q[2];
rz(1.7309534) q[3];
sx q[3];
rz(-0.43678108) q[3];
sx q[3];
rz(-2.3763451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2525758) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.837773) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3712758) q[0];
sx q[0];
rz(-1.2219984) q[0];
sx q[0];
rz(-2.6692997) q[0];
rz(-0.41886139) q[2];
sx q[2];
rz(-1.7423811) q[2];
sx q[2];
rz(0.78531839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5419546) q[1];
sx q[1];
rz(-1.9497739) q[1];
sx q[1];
rz(0.71211262) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(-0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.8213182) q[1];
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
rz(-1.6229796) q[0];
sx q[0];
rz(0.90953565) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57967474) q[2];
sx q[2];
rz(-2.0849166) q[2];
sx q[2];
rz(-2.3953953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1491579) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(-2.1402332) q[1];
rz(2.8619814) q[3];
sx q[3];
rz(-1.8108484) q[3];
sx q[3];
rz(0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(3.0274042) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(-2.3506929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0349802) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(-1.4094704) q[0];
rz(0.36206836) q[2];
sx q[2];
rz(-0.77713359) q[2];
sx q[2];
rz(-2.2287378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0930209) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(0.084333468) q[1];
x q[2];
rz(0.36356504) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(3.0594861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(-2.8875202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92914903) q[0];
sx q[0];
rz(-1.375276) q[0];
sx q[0];
rz(-1.50302) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0052070219) q[2];
sx q[2];
rz(-1.5723096) q[2];
sx q[2];
rz(0.13841275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.904774) q[1];
sx q[1];
rz(-2.8042951) q[1];
sx q[1];
rz(-2.7493613) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6510321) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(0.3790006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-0.33561486) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-2.1597247) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69753416) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(0.27192893) q[0];
x q[1];
rz(0.1083381) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(-0.94891753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.3853449) q[1];
sx q[1];
rz(-0.058538392) q[1];
rz(-1.8198265) q[3];
sx q[3];
rz(-2.2662275) q[3];
sx q[3];
rz(1.673656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92423576) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(-0.99739933) q[2];
rz(-0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(3.1072646) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2845594) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(-2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(2.8732252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70290138) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(-2.2772574) q[0];
rz(-pi) q[1];
rz(-0.95558138) q[2];
sx q[2];
rz(-2.4590543) q[2];
sx q[2];
rz(-2.5126484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2942218) q[1];
sx q[1];
rz(-2.1547744) q[1];
sx q[1];
rz(-1.2482615) q[1];
rz(-pi) q[2];
rz(-2.8769365) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(-2.2789126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.8019567) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(2.6396991) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
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