OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.46848133) q[0];
sx q[0];
rz(3.0335479) q[0];
sx q[0];
rz(11.57668) q[0];
rz(-1.5762848) q[1];
sx q[1];
rz(-1.58374) q[1];
sx q[1];
rz(-1.6852112) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1030758) q[0];
sx q[0];
rz(-1.4355735) q[0];
sx q[0];
rz(-1.200202) q[0];
x q[1];
rz(0.26881071) q[2];
sx q[2];
rz(-0.13829002) q[2];
sx q[2];
rz(-2.2594096) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.072468199) q[1];
sx q[1];
rz(-2.1241423) q[1];
sx q[1];
rz(-2.8072559) q[1];
rz(-pi) q[2];
rz(3.0979324) q[3];
sx q[3];
rz(-2.147445) q[3];
sx q[3];
rz(-0.94983627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9911389) q[2];
sx q[2];
rz(-3.1296215) q[2];
sx q[2];
rz(-1.0452622) q[2];
rz(2.1446877) q[3];
sx q[3];
rz(-0.0054587047) q[3];
sx q[3];
rz(-1.8256942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.588722) q[0];
sx q[0];
rz(-1.9010239) q[0];
sx q[0];
rz(-1.7887822) q[0];
rz(3.1006587) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(-1.5997684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38986015) q[0];
sx q[0];
rz(-0.19987389) q[0];
sx q[0];
rz(-0.78011192) q[0];
rz(1.5920609) q[2];
sx q[2];
rz(-1.5870811) q[2];
sx q[2];
rz(-2.2546538) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1499298) q[1];
sx q[1];
rz(-2.5444391) q[1];
sx q[1];
rz(-1.6258904) q[1];
x q[2];
rz(1.0066693) q[3];
sx q[3];
rz(-2.4897277) q[3];
sx q[3];
rz(-1.3019555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.505595) q[2];
sx q[2];
rz(-0.032568585) q[2];
sx q[2];
rz(0.51844281) q[2];
rz(2.017766) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(0.43784416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951303) q[0];
sx q[0];
rz(-0.050124425) q[0];
sx q[0];
rz(-1.6019524) q[0];
rz(-0.72499544) q[1];
sx q[1];
rz(-3.1097737) q[1];
sx q[1];
rz(-0.67951387) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5121661) q[0];
sx q[0];
rz(-2.0039194) q[0];
sx q[0];
rz(0.64164759) q[0];
x q[1];
rz(0.21444397) q[2];
sx q[2];
rz(-1.5753928) q[2];
sx q[2];
rz(1.4397804) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0371548) q[1];
sx q[1];
rz(-1.3111115) q[1];
sx q[1];
rz(-1.2472769) q[1];
rz(-pi) q[2];
rz(2.1587426) q[3];
sx q[3];
rz(-1.790739) q[3];
sx q[3];
rz(1.6538805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2346377) q[2];
sx q[2];
rz(-0.24181557) q[2];
sx q[2];
rz(-0.50092906) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-3.1152476) q[3];
sx q[3];
rz(-0.19550368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2723715) q[0];
sx q[0];
rz(-0.048373241) q[0];
sx q[0];
rz(0.79917556) q[0];
rz(0.29798206) q[1];
sx q[1];
rz(-0.25845343) q[1];
sx q[1];
rz(-2.2395649) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69660115) q[0];
sx q[0];
rz(-1.4334911) q[0];
sx q[0];
rz(-0.75604646) q[0];
rz(2.8935585) q[2];
sx q[2];
rz(-2.5861861) q[2];
sx q[2];
rz(0.75128864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9060256) q[1];
sx q[1];
rz(-0.14793554) q[1];
sx q[1];
rz(-1.6355913) q[1];
rz(-pi) q[2];
rz(2.1893327) q[3];
sx q[3];
rz(-0.094933184) q[3];
sx q[3];
rz(-0.003530276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4681089) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(1.6710949) q[2];
rz(-2.9521613) q[3];
sx q[3];
rz(-3.093284) q[3];
sx q[3];
rz(2.3725177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22537941) q[0];
sx q[0];
rz(-0.12752859) q[0];
sx q[0];
rz(-2.7498229) q[0];
rz(2.0562992) q[1];
sx q[1];
rz(-3.1317874) q[1];
sx q[1];
rz(0.36878961) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1802496) q[0];
sx q[0];
rz(-1.9260487) q[0];
sx q[0];
rz(0.60549462) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1541751) q[2];
sx q[2];
rz(-1.81798) q[2];
sx q[2];
rz(-0.32933035) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5318067) q[1];
sx q[1];
rz(-1.5641677) q[1];
sx q[1];
rz(3.1408491) q[1];
rz(-pi) q[2];
rz(-0.27009876) q[3];
sx q[3];
rz(-2.4135347) q[3];
sx q[3];
rz(-0.61157507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58086777) q[2];
sx q[2];
rz(-3.0256425) q[2];
sx q[2];
rz(1.7725393) q[2];
rz(0.12618682) q[3];
sx q[3];
rz(-2.7044665) q[3];
sx q[3];
rz(1.060846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5996025) q[0];
sx q[0];
rz(-0.76102155) q[0];
sx q[0];
rz(-1.5498932) q[0];
rz(-0.97269336) q[1];
sx q[1];
rz(-2.9284271) q[1];
sx q[1];
rz(1.0290283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1465107) q[0];
sx q[0];
rz(-1.2702939) q[0];
sx q[0];
rz(2.5142558) q[0];
rz(-pi) q[1];
rz(-0.53035276) q[2];
sx q[2];
rz(-2.8674922) q[2];
sx q[2];
rz(-2.0567187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29063836) q[1];
sx q[1];
rz(-0.096157638) q[1];
sx q[1];
rz(1.5753217) q[1];
rz(-pi) q[2];
rz(-0.036110445) q[3];
sx q[3];
rz(-1.7381769) q[3];
sx q[3];
rz(0.13425628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35919967) q[2];
sx q[2];
rz(-1.426037) q[2];
sx q[2];
rz(-0.4314118) q[2];
rz(2.1443478) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(-2.5844311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3417086) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(-2.0836015) q[0];
rz(-0.82408389) q[1];
sx q[1];
rz(-1.7015142e-05) q[1];
sx q[1];
rz(-0.81738671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5800487) q[0];
sx q[0];
rz(-1.9684682) q[0];
sx q[0];
rz(2.9133767) q[0];
rz(-0.36211966) q[2];
sx q[2];
rz(-0.019500168) q[2];
sx q[2];
rz(1.3355314) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83113545) q[1];
sx q[1];
rz(-1.5316446) q[1];
sx q[1];
rz(-1.3292946) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1002186) q[3];
sx q[3];
rz(-2.2499871) q[3];
sx q[3];
rz(2.8619932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99007964) q[2];
sx q[2];
rz(-1.9033056) q[2];
sx q[2];
rz(1.5129169) q[2];
rz(0.40142909) q[3];
sx q[3];
rz(-3.1135058) q[3];
sx q[3];
rz(2.1872971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7734739) q[0];
sx q[0];
rz(-3.0768657) q[0];
sx q[0];
rz(-1.7777959) q[0];
rz(-3.0996481) q[1];
sx q[1];
rz(-3.0105803) q[1];
sx q[1];
rz(-1.0036453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59597385) q[0];
sx q[0];
rz(-0.64028469) q[0];
sx q[0];
rz(2.9856634) q[0];
rz(-pi) q[1];
rz(1.9081536) q[2];
sx q[2];
rz(-1.5631096) q[2];
sx q[2];
rz(-0.53042646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.072362547) q[1];
sx q[1];
rz(-1.8025928) q[1];
sx q[1];
rz(-2.9332471) q[1];
rz(-pi) q[2];
rz(-0.3866638) q[3];
sx q[3];
rz(-2.4485976) q[3];
sx q[3];
rz(2.3946144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7046788) q[2];
sx q[2];
rz(-0.059160058) q[2];
sx q[2];
rz(1.3979647) q[2];
rz(-0.49020234) q[3];
sx q[3];
rz(-3.0981045) q[3];
sx q[3];
rz(1.9628261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040084664) q[0];
sx q[0];
rz(-0.12797102) q[0];
sx q[0];
rz(-0.21139938) q[0];
rz(2.0231817) q[1];
sx q[1];
rz(-0.011592955) q[1];
sx q[1];
rz(1.5107907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64563066) q[0];
sx q[0];
rz(-1.2617153) q[0];
sx q[0];
rz(-1.8739971) q[0];
rz(-pi) q[1];
rz(1.1338992) q[2];
sx q[2];
rz(-0.99183509) q[2];
sx q[2];
rz(1.4329662) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.69204077) q[1];
sx q[1];
rz(-0.15555432) q[1];
sx q[1];
rz(3.0306308) q[1];
rz(-pi) q[2];
rz(-0.16638593) q[3];
sx q[3];
rz(-1.4825487) q[3];
sx q[3];
rz(-1.4180753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7819034) q[2];
sx q[2];
rz(-0.017904559) q[2];
sx q[2];
rz(-2.8622799) q[2];
rz(0.8638047) q[3];
sx q[3];
rz(-0.0044718663) q[3];
sx q[3];
rz(2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2281247) q[0];
sx q[0];
rz(-1.1044015) q[0];
sx q[0];
rz(-1.8260691) q[0];
rz(-2.3663991) q[1];
sx q[1];
rz(-0.23861353) q[1];
sx q[1];
rz(-1.7528037) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7167146) q[0];
sx q[0];
rz(-1.8535063) q[0];
sx q[0];
rz(1.7260051) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49225537) q[2];
sx q[2];
rz(-2.6894719) q[2];
sx q[2];
rz(-1.7017572) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99936472) q[1];
sx q[1];
rz(-1.570163) q[1];
sx q[1];
rz(3.1406979) q[1];
rz(-pi) q[2];
rz(1.989886) q[3];
sx q[3];
rz(-1.1858318) q[3];
sx q[3];
rz(0.46483913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76643252) q[2];
sx q[2];
rz(-0.015724026) q[2];
sx q[2];
rz(-1.4520221) q[2];
rz(2.8826513) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(1.9848721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6257085) q[0];
sx q[0];
rz(-0.71737552) q[0];
sx q[0];
rz(1.4364568) q[0];
rz(1.4868078) q[1];
sx q[1];
rz(-0.27324067) q[1];
sx q[1];
rz(-2.9437093) q[1];
rz(-2.941359) q[2];
sx q[2];
rz(-1.5641134) q[2];
sx q[2];
rz(-1.3431298) q[2];
rz(-3.1155125) q[3];
sx q[3];
rz(-1.9175538) q[3];
sx q[3];
rz(-3.1160311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
