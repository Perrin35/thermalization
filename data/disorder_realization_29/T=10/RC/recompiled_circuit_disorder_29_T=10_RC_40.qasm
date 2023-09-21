OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(-0.32615647) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0155976) q[0];
sx q[0];
rz(-1.3552357) q[0];
sx q[0];
rz(-0.21685812) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36122303) q[2];
sx q[2];
rz(-2.5075956) q[2];
sx q[2];
rz(-2.0915973) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3781158) q[1];
sx q[1];
rz(-0.33707481) q[1];
sx q[1];
rz(0.76428767) q[1];
rz(-0.013734038) q[3];
sx q[3];
rz(-2.3547958) q[3];
sx q[3];
rz(2.9247583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(0.88511434) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(2.2639993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7859902) q[0];
sx q[0];
rz(-0.88760932) q[0];
sx q[0];
rz(-1.3366633) q[0];
rz(2.8009731) q[2];
sx q[2];
rz(-0.38481958) q[2];
sx q[2];
rz(0.22573267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21786015) q[1];
sx q[1];
rz(-1.4638312) q[1];
sx q[1];
rz(-0.033543368) q[1];
rz(-pi) q[2];
rz(1.9545994) q[3];
sx q[3];
rz(-1.2591397) q[3];
sx q[3];
rz(2.4812738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39891222) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(1.1304643) q[2];
rz(1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.6931504) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(-0.6668123) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(3.0677632) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4742972) q[0];
sx q[0];
rz(-1.4814261) q[0];
sx q[0];
rz(-1.6325634) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29088144) q[2];
sx q[2];
rz(-1.2419495) q[2];
sx q[2];
rz(1.8287303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0970043) q[1];
sx q[1];
rz(-1.6267596) q[1];
sx q[1];
rz(1.3929277) q[1];
rz(-pi) q[2];
x q[2];
rz(2.79014) q[3];
sx q[3];
rz(-1.4593399) q[3];
sx q[3];
rz(-0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6510216) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(2.0329287) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(-2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(-1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80285145) q[0];
sx q[0];
rz(-2.2995298) q[0];
sx q[0];
rz(1.6943504) q[0];
rz(0.32707733) q[2];
sx q[2];
rz(-1.2750669) q[2];
sx q[2];
rz(-2.0138274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1632299) q[1];
sx q[1];
rz(-1.1095188) q[1];
sx q[1];
rz(2.7510838) q[1];
x q[2];
rz(-1.4401682) q[3];
sx q[3];
rz(-1.6960973) q[3];
sx q[3];
rz(-0.13392042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58549762) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(2.0193224) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595903) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(-0.37297747) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81985695) q[0];
sx q[0];
rz(-0.61367354) q[0];
sx q[0];
rz(-0.86921285) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8226132) q[2];
sx q[2];
rz(-1.8251112) q[2];
sx q[2];
rz(-0.11676678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59606325) q[1];
sx q[1];
rz(-1.5585594) q[1];
sx q[1];
rz(-1.6227325) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1042215) q[3];
sx q[3];
rz(-1.08764) q[3];
sx q[3];
rz(-2.0869568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0405154) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-0.80580795) q[2];
rz(2.6082883) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(-2.1300952) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.3202753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0549714) q[0];
sx q[0];
rz(-1.939659) q[0];
sx q[0];
rz(-0.47199179) q[0];
x q[1];
rz(2.8684902) q[2];
sx q[2];
rz(-0.58758508) q[2];
sx q[2];
rz(-2.7762129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7399866) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(-0.92958881) q[1];
rz(-pi) q[2];
rz(2.2269507) q[3];
sx q[3];
rz(-1.802889) q[3];
sx q[3];
rz(-2.7616012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-2.5406204) q[2];
rz(-0.48505923) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(1.4453567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(0.55554187) q[0];
rz(3.1069966) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.3909891) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44156528) q[0];
sx q[0];
rz(-1.0374984) q[0];
sx q[0];
rz(1.9576859) q[0];
rz(-pi) q[1];
rz(-2.7752635) q[2];
sx q[2];
rz(-1.3372256) q[2];
sx q[2];
rz(-0.69586588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3607189) q[1];
sx q[1];
rz(-0.90457537) q[1];
sx q[1];
rz(-2.9800376) q[1];
rz(-pi) q[2];
rz(3.0656747) q[3];
sx q[3];
rz(-0.77709353) q[3];
sx q[3];
rz(1.0958375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3283219) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.288712) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(-0.95343268) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(-1.4377726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8673082) q[0];
sx q[0];
rz(-1.058393) q[0];
sx q[0];
rz(-2.9318277) q[0];
rz(2.4687924) q[2];
sx q[2];
rz(-2.9252508) q[2];
sx q[2];
rz(-2.1397482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2807323) q[1];
sx q[1];
rz(-0.8287462) q[1];
sx q[1];
rz(2.4651205) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.612547) q[3];
sx q[3];
rz(-2.6937727) q[3];
sx q[3];
rz(0.24275045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4961204) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(-2.9628741) q[2];
rz(-0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(-0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(2.0478915) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(2.0827983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537542) q[0];
sx q[0];
rz(-2.5231579) q[0];
sx q[0];
rz(-1.0879602) q[0];
rz(-pi) q[1];
rz(-1.0979707) q[2];
sx q[2];
rz(-1.2038004) q[2];
sx q[2];
rz(2.8851913) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.345929) q[1];
sx q[1];
rz(-1.1946031) q[1];
sx q[1];
rz(-0.13161195) q[1];
rz(2.6685647) q[3];
sx q[3];
rz(-1.2086476) q[3];
sx q[3];
rz(0.87272296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.4808562) q[2];
rz(2.7311834) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(-2.9836392) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(2.8503382) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390022) q[0];
sx q[0];
rz(-1.0800835) q[0];
sx q[0];
rz(-1.9309994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45551331) q[2];
sx q[2];
rz(-1.3580772) q[2];
sx q[2];
rz(2.6754975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1283975) q[1];
sx q[1];
rz(-2.5880475) q[1];
sx q[1];
rz(-3.058606) q[1];
rz(1.7581975) q[3];
sx q[3];
rz(-0.43863505) q[3];
sx q[3];
rz(0.24500971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46135205) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(1.1414026) q[2];
rz(-1.6067778) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-2.4070807) q[2];
sx q[2];
rz(-2.8399158) q[2];
sx q[2];
rz(0.5010571) q[2];
rz(-0.42320078) q[3];
sx q[3];
rz(-1.3907708) q[3];
sx q[3];
rz(1.6650865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];