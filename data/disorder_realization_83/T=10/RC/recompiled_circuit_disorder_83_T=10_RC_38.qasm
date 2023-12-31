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
rz(2.397937) q[1];
rz(-pi/2) q[2];
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
x q[1];
rz(-2.9843569) q[2];
sx q[2];
rz(-0.6586282) q[2];
sx q[2];
rz(-0.14070357) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.570959) q[1];
sx q[1];
rz(-1.966241) q[1];
sx q[1];
rz(-2.8407437) q[1];
x q[2];
rz(-1.9348014) q[3];
sx q[3];
rz(-1.6379106) q[3];
sx q[3];
rz(-1.2957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-2.4690348) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-0.67115152) q[1];
sx q[1];
rz(0.040963106) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3687392) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(-0.98451891) q[0];
x q[1];
rz(2.3054753) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(0.9678313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37521024) q[1];
sx q[1];
rz(-0.76645215) q[1];
sx q[1];
rz(2.4541928) q[1];
rz(-pi) q[2];
rz(0.37946545) q[3];
sx q[3];
rz(-1.9455519) q[3];
sx q[3];
rz(1.7750164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9791947) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(2.6290821) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(1.846116) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-2.4643262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5325461) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(-0.30058582) q[0];
rz(-pi) q[1];
rz(2.9708423) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(-3.1250172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22563572) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(3.008328) q[1];
x q[2];
rz(2.4694091) q[3];
sx q[3];
rz(-2.5568092) q[3];
sx q[3];
rz(-2.4206846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(-2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.9690008) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8797982) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(1.7394702) q[0];
x q[1];
rz(-1.555221) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(1.7473999) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6231411) q[1];
sx q[1];
rz(-1.271283) q[1];
sx q[1];
rz(1.2007984) q[1];
rz(0.074313642) q[3];
sx q[3];
rz(-2.001611) q[3];
sx q[3];
rz(-2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(2.5505113) q[1];
sx q[1];
rz(-0.68060827) q[1];
sx q[1];
rz(0.65471929) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39010534) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(-0.67436995) q[0];
x q[1];
rz(2.7388217) q[2];
sx q[2];
rz(-0.45071128) q[2];
sx q[2];
rz(0.41926256) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.66203413) q[1];
sx q[1];
rz(-2.2231632) q[1];
sx q[1];
rz(1.0865092) q[1];
x q[2];
rz(-1.6642903) q[3];
sx q[3];
rz(-1.9545385) q[3];
sx q[3];
rz(-1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(-2.5518104) q[3];
sx q[3];
rz(-0.46001205) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59693903) q[0];
sx q[0];
rz(-0.66300809) q[0];
sx q[0];
rz(-1.4859499) q[0];
x q[1];
rz(-0.97700714) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-1.1360816) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99243473) q[1];
sx q[1];
rz(-2.5870393) q[1];
sx q[1];
rz(-1.0013594) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8619814) q[3];
sx q[3];
rz(-1.3307443) q[3];
sx q[3];
rz(0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(3.1141172) q[2];
rz(-0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(2.3506929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9081156) q[0];
sx q[0];
rz(-0.95196264) q[0];
sx q[0];
rz(-3.0254362) q[0];
rz(-0.36206836) q[2];
sx q[2];
rz(-2.3644591) q[2];
sx q[2];
rz(0.91285489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0930209) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(-0.084333468) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2579455) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(-1.7162232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(0.67251742) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(-2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5131322) q[0];
sx q[0];
rz(-1.6372794) q[0];
sx q[0];
rz(0.19595887) q[0];
x q[1];
rz(0.28283624) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(-1.1495513) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1797807) q[1];
sx q[1];
rz(-1.6976377) q[1];
sx q[1];
rz(-2.8282053) q[1];
rz(-pi) q[2];
rz(-1.6510321) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(-0.3790006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(2.1597247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69753416) q[0];
sx q[0];
rz(-0.42675787) q[0];
sx q[0];
rz(2.8696637) q[0];
rz(-pi) q[1];
rz(0.1083381) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(2.1926751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5280142) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(-1.8730875) q[1];
x q[2];
rz(0.28717678) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-0.034328073) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(2.6249264) q[0];
rz(-2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(2.8732252) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4386913) q[0];
sx q[0];
rz(-2.1500906) q[0];
sx q[0];
rz(-0.86433522) q[0];
rz(-2.1568314) q[2];
sx q[2];
rz(-1.9433937) q[2];
sx q[2];
rz(0.44024703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8908773) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(0.4472181) q[1];
x q[2];
rz(1.1770505) q[3];
sx q[3];
rz(-0.6150133) q[3];
sx q[3];
rz(-2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.664809) q[0];
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
rz(0.31233882) q[3];
sx q[3];
rz(-1.7620419) q[3];
sx q[3];
rz(-0.077802303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
