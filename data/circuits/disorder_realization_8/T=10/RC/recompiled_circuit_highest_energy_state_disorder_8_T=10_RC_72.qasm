OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.8166135) q[0];
sx q[0];
rz(-0.54083523) q[0];
sx q[0];
rz(1.1582561) q[0];
rz(0.090016063) q[1];
sx q[1];
rz(-2.6675192) q[1];
sx q[1];
rz(2.3064244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1285889) q[0];
sx q[0];
rz(-2.3478201) q[0];
sx q[0];
rz(2.6966266) q[0];
rz(2.7273601) q[2];
sx q[2];
rz(-2.4576839) q[2];
sx q[2];
rz(1.4642844) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0660634) q[1];
sx q[1];
rz(-1.3261686) q[1];
sx q[1];
rz(-0.3800769) q[1];
rz(-pi) q[2];
rz(-2.4832544) q[3];
sx q[3];
rz(-1.1000203) q[3];
sx q[3];
rz(1.8792764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0917255) q[2];
sx q[2];
rz(-0.89605248) q[2];
sx q[2];
rz(0.33207616) q[2];
rz(2.8927228) q[3];
sx q[3];
rz(-1.9460461) q[3];
sx q[3];
rz(-2.2539049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2700972) q[0];
sx q[0];
rz(-2.8506554) q[0];
sx q[0];
rz(-2.6089597) q[0];
rz(-0.10781413) q[1];
sx q[1];
rz(-2.0262521) q[1];
sx q[1];
rz(-0.32049387) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6924393) q[0];
sx q[0];
rz(-0.65523132) q[0];
sx q[0];
rz(1.4661319) q[0];
rz(-2.2079289) q[2];
sx q[2];
rz(-0.28544989) q[2];
sx q[2];
rz(-0.97412005) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35831603) q[1];
sx q[1];
rz(-0.22138295) q[1];
sx q[1];
rz(2.3795147) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32933195) q[3];
sx q[3];
rz(-1.7747702) q[3];
sx q[3];
rz(2.4402809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8104441) q[2];
sx q[2];
rz(-2.836477) q[2];
sx q[2];
rz(-1.6395052) q[2];
rz(-2.8034927) q[3];
sx q[3];
rz(-2.2564042) q[3];
sx q[3];
rz(-2.6147208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51760393) q[0];
sx q[0];
rz(-1.232134) q[0];
sx q[0];
rz(2.7841618) q[0];
rz(-2.7492211) q[1];
sx q[1];
rz(-0.79634276) q[1];
sx q[1];
rz(1.2145112) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75823513) q[0];
sx q[0];
rz(-1.5534983) q[0];
sx q[0];
rz(-1.7123459) q[0];
rz(-pi) q[1];
rz(1.7262801) q[2];
sx q[2];
rz(-2.1124055) q[2];
sx q[2];
rz(2.3692148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7347265) q[1];
sx q[1];
rz(-0.38741437) q[1];
sx q[1];
rz(-1.563036) q[1];
x q[2];
rz(-2.1585095) q[3];
sx q[3];
rz(-1.6020892) q[3];
sx q[3];
rz(0.33558095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7963205) q[2];
sx q[2];
rz(-0.97769633) q[2];
sx q[2];
rz(0.39682445) q[2];
rz(1.8848298) q[3];
sx q[3];
rz(-1.7233012) q[3];
sx q[3];
rz(0.61029148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9339555) q[0];
sx q[0];
rz(-1.3960681) q[0];
sx q[0];
rz(-1.2285832) q[0];
rz(-1.2127016) q[1];
sx q[1];
rz(-2.1804501) q[1];
sx q[1];
rz(-0.62087762) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29781326) q[0];
sx q[0];
rz(-0.50619805) q[0];
sx q[0];
rz(-1.9202581) q[0];
rz(-pi) q[1];
rz(0.75350301) q[2];
sx q[2];
rz(-0.22543365) q[2];
sx q[2];
rz(1.6627251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.634269) q[1];
sx q[1];
rz(-2.7719927) q[1];
sx q[1];
rz(-2.3564767) q[1];
x q[2];
rz(-1.6895164) q[3];
sx q[3];
rz(-1.190589) q[3];
sx q[3];
rz(0.79352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2584201) q[2];
sx q[2];
rz(-1.1383388) q[2];
sx q[2];
rz(2.9849198) q[2];
rz(-2.2251718) q[3];
sx q[3];
rz(-1.0333002) q[3];
sx q[3];
rz(2.5887183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85585344) q[0];
sx q[0];
rz(-1.1612949) q[0];
sx q[0];
rz(-0.39749417) q[0];
rz(-3.1187348) q[1];
sx q[1];
rz(-2.6409179) q[1];
sx q[1];
rz(0.674725) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8431339) q[0];
sx q[0];
rz(-1.3359016) q[0];
sx q[0];
rz(-2.9265334) q[0];
rz(0.12395383) q[2];
sx q[2];
rz(-1.6646241) q[2];
sx q[2];
rz(1.4012865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51395638) q[1];
sx q[1];
rz(-1.7412724) q[1];
sx q[1];
rz(-1.7254616) q[1];
rz(-pi) q[2];
rz(0.27366112) q[3];
sx q[3];
rz(-1.869264) q[3];
sx q[3];
rz(0.65200114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61458331) q[2];
sx q[2];
rz(-0.60636568) q[2];
sx q[2];
rz(2.442339) q[2];
rz(-1.3462542) q[3];
sx q[3];
rz(-2.5739539) q[3];
sx q[3];
rz(0.7114555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4166477) q[0];
sx q[0];
rz(-0.41794932) q[0];
sx q[0];
rz(0.39837343) q[0];
rz(-2.1669972) q[1];
sx q[1];
rz(-1.3037126) q[1];
sx q[1];
rz(-1.4385361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2319039) q[0];
sx q[0];
rz(-1.6758907) q[0];
sx q[0];
rz(-2.8382481) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6104524) q[2];
sx q[2];
rz(-0.25601124) q[2];
sx q[2];
rz(-2.2289952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8827445) q[1];
sx q[1];
rz(-2.2387894) q[1];
sx q[1];
rz(-2.8440688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8277728) q[3];
sx q[3];
rz(-2.0517618) q[3];
sx q[3];
rz(-2.0723267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.67821104) q[2];
sx q[2];
rz(-1.3621829) q[2];
sx q[2];
rz(1.8360651) q[2];
rz(1.3156923) q[3];
sx q[3];
rz(-1.7782327) q[3];
sx q[3];
rz(2.2731884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797794) q[0];
sx q[0];
rz(-0.4158622) q[0];
sx q[0];
rz(-1.4965936) q[0];
rz(-2.1553701) q[1];
sx q[1];
rz(-1.58135) q[1];
sx q[1];
rz(-2.644002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80579306) q[0];
sx q[0];
rz(-1.3168653) q[0];
sx q[0];
rz(3.0428314) q[0];
x q[1];
rz(2.1946218) q[2];
sx q[2];
rz(-1.1188036) q[2];
sx q[2];
rz(-0.78540451) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10863241) q[1];
sx q[1];
rz(-1.4147621) q[1];
sx q[1];
rz(-1.1724657) q[1];
rz(-pi) q[2];
rz(1.6646181) q[3];
sx q[3];
rz(-1.2872496) q[3];
sx q[3];
rz(2.97992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8387973) q[2];
sx q[2];
rz(-1.5668198) q[2];
sx q[2];
rz(-2.8509169) q[2];
rz(0.21026462) q[3];
sx q[3];
rz(-0.99431521) q[3];
sx q[3];
rz(2.1073585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656089) q[0];
sx q[0];
rz(-1.1701595) q[0];
sx q[0];
rz(1.1822816) q[0];
rz(-2.9871509) q[1];
sx q[1];
rz(-1.7223822) q[1];
sx q[1];
rz(-0.90528893) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8715031) q[0];
sx q[0];
rz(-2.8686028) q[0];
sx q[0];
rz(1.6156625) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8087093) q[2];
sx q[2];
rz(-1.1255956) q[2];
sx q[2];
rz(-0.92724909) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0876079) q[1];
sx q[1];
rz(-1.4248669) q[1];
sx q[1];
rz(2.3838359) q[1];
rz(-1.9117113) q[3];
sx q[3];
rz(-1.3671759) q[3];
sx q[3];
rz(3.1217755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0514544) q[2];
sx q[2];
rz(-1.5608414) q[2];
sx q[2];
rz(2.1913989) q[2];
rz(3.0692302) q[3];
sx q[3];
rz(-1.3772929) q[3];
sx q[3];
rz(-2.4322521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41910928) q[0];
sx q[0];
rz(-2.1868732) q[0];
sx q[0];
rz(-1.0954274) q[0];
rz(-2.0504045) q[1];
sx q[1];
rz(-0.90091101) q[1];
sx q[1];
rz(-0.99517623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5916717) q[0];
sx q[0];
rz(-0.66837817) q[0];
sx q[0];
rz(0.10381283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2500317) q[2];
sx q[2];
rz(-2.7195647) q[2];
sx q[2];
rz(-1.0996549) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6217664) q[1];
sx q[1];
rz(-1.6556181) q[1];
sx q[1];
rz(1.6165401) q[1];
x q[2];
rz(2.9561958) q[3];
sx q[3];
rz(-1.3747921) q[3];
sx q[3];
rz(0.40962266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5114078) q[2];
sx q[2];
rz(-0.46773043) q[2];
sx q[2];
rz(-2.4616145) q[2];
rz(-1.6857111) q[3];
sx q[3];
rz(-2.2472491) q[3];
sx q[3];
rz(1.9623914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42530123) q[0];
sx q[0];
rz(-0.93451262) q[0];
sx q[0];
rz(2.5262078) q[0];
rz(1.8062493) q[1];
sx q[1];
rz(-1.5348624) q[1];
sx q[1];
rz(-1.7937484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3894661) q[0];
sx q[0];
rz(-1.9258317) q[0];
sx q[0];
rz(1.5219164) q[0];
rz(-pi) q[1];
rz(-0.87839076) q[2];
sx q[2];
rz(-2.2579402) q[2];
sx q[2];
rz(-0.85604446) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19619689) q[1];
sx q[1];
rz(-2.7295503) q[1];
sx q[1];
rz(2.7315188) q[1];
rz(-pi) q[2];
rz(-0.12935454) q[3];
sx q[3];
rz(-2.3812897) q[3];
sx q[3];
rz(3.1227675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8936257) q[2];
sx q[2];
rz(-1.3004356) q[2];
sx q[2];
rz(2.628053) q[2];
rz(2.9945471) q[3];
sx q[3];
rz(-2.5976318) q[3];
sx q[3];
rz(1.2646382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2697987) q[0];
sx q[0];
rz(-0.88903058) q[0];
sx q[0];
rz(-0.24118184) q[0];
rz(1.8409894) q[1];
sx q[1];
rz(-1.069297) q[1];
sx q[1];
rz(-0.020513608) q[1];
rz(-1.945449) q[2];
sx q[2];
rz(-1.7032663) q[2];
sx q[2];
rz(1.8864529) q[2];
rz(-0.5823395) q[3];
sx q[3];
rz(-2.4262541) q[3];
sx q[3];
rz(-2.2179009) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
