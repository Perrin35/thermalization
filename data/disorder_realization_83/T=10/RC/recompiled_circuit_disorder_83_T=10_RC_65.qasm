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
rz(0.76454437) q[1];
sx q[1];
rz(-2.0643056) q[1];
sx q[1];
rz(0.74365562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3418158) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(-0.067702985) q[0];
rz(-2.9843569) q[2];
sx q[2];
rz(-2.4829645) q[2];
sx q[2];
rz(0.14070357) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88120645) q[1];
sx q[1];
rz(-1.8477866) q[1];
sx q[1];
rz(-1.1587515) q[1];
rz(-pi) q[2];
rz(1.7573962) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(0.44934011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(0.23072492) q[0];
rz(-1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(3.1006295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.165867) q[0];
sx q[0];
rz(-1.8164993) q[0];
sx q[0];
rz(2.9742572) q[0];
rz(-2.8616521) q[2];
sx q[2];
rz(-1.8676114) q[2];
sx q[2];
rz(-2.9505626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4119271) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(-0.63974849) q[1];
rz(-pi) q[2];
rz(1.9713692) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(2.7924201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(1.846116) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(0.67726642) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579502) q[0];
sx q[0];
rz(-1.5188688) q[0];
sx q[0];
rz(2.9731263) q[0];
rz(-pi) q[1];
rz(-2.9708423) q[2];
sx q[2];
rz(-0.0064592529) q[2];
sx q[2];
rz(-3.1250172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2689506) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(-0.93572576) q[1];
x q[2];
rz(-2.6636395) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(1.4357476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-0.051483367) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576142) q[0];
sx q[0];
rz(-1.4371705) q[0];
sx q[0];
rz(2.4806116) q[0];
rz(-pi) q[1];
rz(2.5355859) q[2];
sx q[2];
rz(-1.5579941) q[2];
sx q[2];
rz(0.16773227) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40201515) q[1];
sx q[1];
rz(-0.47164729) q[1];
sx q[1];
rz(-0.86400835) q[1];
rz(-0.074313642) q[3];
sx q[3];
rz(-2.001611) q[3];
sx q[3];
rz(2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-0.48689294) q[2];
rz(-1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30381969) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.8160965) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(2.4868734) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7692208) q[0];
sx q[0];
rz(-1.1290316) q[0];
sx q[0];
rz(-1.9584993) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40277092) q[2];
sx q[2];
rz(-0.45071128) q[2];
sx q[2];
rz(2.7223301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59963804) q[1];
sx q[1];
rz(-1.1918187) q[1];
sx q[1];
rz(2.42948) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(-1.3202745) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70443557) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(0.066083834) q[0];
x q[1];
rz(-0.80070337) q[2];
sx q[2];
rz(-0.75468894) q[2];
sx q[2];
rz(0.18037361) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7945929) q[1];
sx q[1];
rz(-1.1113249) q[1];
sx q[1];
rz(-0.32230349) q[1];
rz(-pi) q[2];
rz(-1.8201581) q[3];
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
rz(-1.4279782) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(-0.79089975) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1066125) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(1.7321222) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9060471) q[2];
sx q[2];
rz(-2.2860048) q[2];
sx q[2];
rz(-2.7170979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6939047) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(-1.9359971) q[1];
rz(-1.1366208) q[3];
sx q[3];
rz(-0.73528157) q[3];
sx q[3];
rz(-2.6649464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(-2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi/2) q[2];
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
x q[1];
rz(0.28283624) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(-1.1495513) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2368187) q[1];
sx q[1];
rz(-0.33729759) q[1];
sx q[1];
rz(2.7493613) q[1];
rz(3.0244175) q[3];
sx q[3];
rz(-2.5391038) q[3];
sx q[3];
rz(-2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-2.1597247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40030038) q[0];
sx q[0];
rz(-1.1606845) q[0];
sx q[0];
rz(1.4492695) q[0];
rz(3.0332546) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(0.94891753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6135785) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(-1.8730875) q[1];
rz(1.3217661) q[3];
sx q[3];
rz(-0.8753652) q[3];
sx q[3];
rz(-1.673656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2845594) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(0.51666623) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(-2.8732252) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4386913) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(0.86433522) q[0];
x q[1];
rz(-0.95558138) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(-0.62894422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8473709) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(1.8933312) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99276944) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(-0.85152599) q[3];
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
rz(-1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(2.6396991) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
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
