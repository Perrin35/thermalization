OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6095088) q[0];
sx q[0];
rz(-2.3713787) q[0];
sx q[0];
rz(-2.8340333) q[0];
rz(-pi) q[1];
rz(-2.5277532) q[2];
sx q[2];
rz(-1.5547353) q[2];
sx q[2];
rz(-1.2889372) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7588501) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(1.43169) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59431608) q[3];
sx q[3];
rz(-1.4694957) q[3];
sx q[3];
rz(3.1241824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(1.356396) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(1.9869841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(-2.4029762) q[0];
rz(-pi) q[1];
rz(0.39312675) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(-1.4002422) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6007538) q[1];
sx q[1];
rz(-0.80889091) q[1];
sx q[1];
rz(-1.6879338) q[1];
rz(1.8335908) q[3];
sx q[3];
rz(-1.7495219) q[3];
sx q[3];
rz(2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(-0.88341218) q[2];
rz(-0.47131053) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(2.3538891) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8283591) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(1.6261684) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(1.0916969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1356782) q[0];
sx q[0];
rz(-2.1164829) q[0];
sx q[0];
rz(-2.2360327) q[0];
rz(2.0902363) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(1.5067593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.162902) q[1];
sx q[1];
rz(-2.2802417) q[1];
sx q[1];
rz(0.49180007) q[1];
rz(-pi) q[2];
rz(1.6352429) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(-2.8116022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83051935) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(0.31035796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.719602) q[0];
sx q[0];
rz(-1.70277) q[0];
sx q[0];
rz(2.7601526) q[0];
rz(-1.4448302) q[2];
sx q[2];
rz(-2.2519886) q[2];
sx q[2];
rz(-1.7484401) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.031361) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(-1.5178174) q[1];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.4471874) q[3];
sx q[3];
rz(-1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0115396) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(-1.0774353) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(1.594054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-2.8919343) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72083005) q[0];
sx q[0];
rz(-1.8846858) q[0];
sx q[0];
rz(-1.3563966) q[0];
rz(-pi) q[1];
rz(2.444961) q[2];
sx q[2];
rz(-1.4259035) q[2];
sx q[2];
rz(0.87755132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40904564) q[1];
sx q[1];
rz(-2.2605719) q[1];
sx q[1];
rz(-1.9802666) q[1];
x q[2];
rz(0.16964511) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(2.5559705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8683118) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(0.67374054) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.3195066) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34981397) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(-0.42516431) q[1];
sx q[1];
rz(-0.95562569) q[1];
sx q[1];
rz(1.649883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0817889) q[0];
sx q[0];
rz(-2.7002618) q[0];
sx q[0];
rz(2.9102737) q[0];
rz(-2.4620373) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(-pi) q[2];
rz(-1.8317354) q[3];
sx q[3];
rz(-0.3762227) q[3];
sx q[3];
rz(2.6747243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-3.0498665) q[2];
rz(-0.84364676) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-3.1076028) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2720374) q[0];
sx q[0];
rz(-0.79622686) q[0];
sx q[0];
rz(-1.7982593) q[0];
rz(0.88862822) q[2];
sx q[2];
rz(-2.3805328) q[2];
sx q[2];
rz(1.467848) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36044869) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(-3.0767246) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61492413) q[3];
sx q[3];
rz(-1.691754) q[3];
sx q[3];
rz(1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5380481) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-2.9984737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8975163) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(1.4454909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60364265) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(-3.135878) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1922853) q[2];
sx q[2];
rz(-2.0645803) q[2];
sx q[2];
rz(0.91044237) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6519424) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(-2.5850992) q[1];
rz(-pi) q[2];
rz(-0.66283488) q[3];
sx q[3];
rz(-1.0532866) q[3];
sx q[3];
rz(-2.585632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80609926) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(-0.30977419) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99684925) q[0];
sx q[0];
rz(-2.092917) q[0];
sx q[0];
rz(1.4247308) q[0];
rz(1.6236213) q[2];
sx q[2];
rz(-0.17715684) q[2];
sx q[2];
rz(-1.0747386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27948353) q[1];
sx q[1];
rz(-0.21394193) q[1];
sx q[1];
rz(1.2941542) q[1];
rz(-pi) q[2];
rz(0.36848948) q[3];
sx q[3];
rz(-2.512305) q[3];
sx q[3];
rz(-3.1143509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.9343728) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(-0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(-1.4996128) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8382032) q[0];
sx q[0];
rz(-1.8614385) q[0];
sx q[0];
rz(0.10586664) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76403107) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(-1.3438091) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3824532) q[1];
sx q[1];
rz(-1.1083974) q[1];
sx q[1];
rz(-0.23666246) q[1];
x q[2];
rz(-1.3979982) q[3];
sx q[3];
rz(-1.5683335) q[3];
sx q[3];
rz(-0.60009225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88400921) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(-0.36869129) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35836999) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(3.070667) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-1.0735738) q[2];
sx q[2];
rz(-1.9855269) q[2];
sx q[2];
rz(-1.2863458) q[2];
rz(0.55340135) q[3];
sx q[3];
rz(-2.3407866) q[3];
sx q[3];
rz(1.8368807) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
