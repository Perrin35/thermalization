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
rz(1.4577515) q[0];
sx q[0];
rz(-3.0384851) q[0];
sx q[0];
rz(-0.74080324) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(3.8771602) q[1];
sx q[1];
rz(9.1280042) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73707092) q[0];
sx q[0];
rz(-1.9849791) q[0];
sx q[0];
rz(-0.41612723) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8662017) q[2];
sx q[2];
rz(-2.4045334) q[2];
sx q[2];
rz(0.136497) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7267043) q[1];
sx q[1];
rz(-1.0203531) q[1];
sx q[1];
rz(-1.5617983) q[1];
x q[2];
rz(-1.9363251) q[3];
sx q[3];
rz(-0.46896471) q[3];
sx q[3];
rz(-0.381857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2510117) q[2];
sx q[2];
rz(-0.0958395) q[2];
sx q[2];
rz(1.4061692) q[2];
rz(-2.3928394) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(2.9296618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9894079) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(-2.8023791) q[0];
rz(-2.5665414) q[1];
sx q[1];
rz(-1.1851858) q[1];
sx q[1];
rz(-2.7599879) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71020102) q[0];
sx q[0];
rz(-1.9862735) q[0];
sx q[0];
rz(0.25449591) q[0];
x q[1];
rz(2.223338) q[2];
sx q[2];
rz(-0.41792575) q[2];
sx q[2];
rz(-1.1993091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2927865) q[1];
sx q[1];
rz(-1.5804351) q[1];
sx q[1];
rz(1.9697492) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80581237) q[3];
sx q[3];
rz(-1.1639813) q[3];
sx q[3];
rz(-0.47570566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7316458) q[2];
sx q[2];
rz(-0.68845981) q[2];
sx q[2];
rz(-0.48639578) q[2];
rz(0.54388034) q[3];
sx q[3];
rz(-1.0433652) q[3];
sx q[3];
rz(-0.16415183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51089066) q[0];
sx q[0];
rz(-0.4158026) q[0];
sx q[0];
rz(-2.1422332) q[0];
rz(0.73127812) q[1];
sx q[1];
rz(-2.6731698) q[1];
sx q[1];
rz(0.17459248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098487735) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(-1.5763603) q[0];
rz(-0.83589696) q[2];
sx q[2];
rz(-0.23698254) q[2];
sx q[2];
rz(1.8191172) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.120003) q[1];
sx q[1];
rz(-1.6232511) q[1];
sx q[1];
rz(-1.4823556) q[1];
rz(-2.0295983) q[3];
sx q[3];
rz(-1.3743041) q[3];
sx q[3];
rz(-2.8216763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38873765) q[2];
sx q[2];
rz(-0.85192215) q[2];
sx q[2];
rz(-1.8641776) q[2];
rz(-2.3169005) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(-0.1674913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.326062) q[0];
sx q[0];
rz(-0.46634665) q[0];
sx q[0];
rz(-1.9933568) q[0];
rz(0.3321906) q[1];
sx q[1];
rz(-2.5416608) q[1];
sx q[1];
rz(-2.0567599) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6713632) q[0];
sx q[0];
rz(-1.7106904) q[0];
sx q[0];
rz(-1.8663919) q[0];
rz(-pi) q[1];
rz(2.8021028) q[2];
sx q[2];
rz(-2.2953646) q[2];
sx q[2];
rz(2.9629415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7521811) q[1];
sx q[1];
rz(-1.6232921) q[1];
sx q[1];
rz(0.73867259) q[1];
x q[2];
rz(-2.6779019) q[3];
sx q[3];
rz(-1.9640018) q[3];
sx q[3];
rz(-0.56305199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3499902) q[2];
sx q[2];
rz(-2.377066) q[2];
sx q[2];
rz(0.0066268607) q[2];
rz(-0.22348063) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(0.82911432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70581907) q[0];
sx q[0];
rz(-2.3888102) q[0];
sx q[0];
rz(-1.1821795) q[0];
rz(-1.301282) q[1];
sx q[1];
rz(-0.92295206) q[1];
sx q[1];
rz(-0.15929793) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.015291) q[0];
sx q[0];
rz(-1.8582004) q[0];
sx q[0];
rz(-1.5319583) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2550687) q[2];
sx q[2];
rz(-1.7509394) q[2];
sx q[2];
rz(-2.2314928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9827798) q[1];
sx q[1];
rz(-2.4430664) q[1];
sx q[1];
rz(2.34312) q[1];
rz(-pi) q[2];
rz(-2.716655) q[3];
sx q[3];
rz(-2.1930088) q[3];
sx q[3];
rz(0.10729161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2948239) q[2];
sx q[2];
rz(-0.20192768) q[2];
sx q[2];
rz(-1.798604) q[2];
rz(-2.790847) q[3];
sx q[3];
rz(-2.0028159) q[3];
sx q[3];
rz(0.32575592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2419234) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(-2.9747466) q[0];
rz(0.22653656) q[1];
sx q[1];
rz(-1.3275361) q[1];
sx q[1];
rz(2.66364) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27875459) q[0];
sx q[0];
rz(-2.4669381) q[0];
sx q[0];
rz(2.1341896) q[0];
rz(1.208838) q[2];
sx q[2];
rz(-2.2565292) q[2];
sx q[2];
rz(0.63206965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2061658) q[1];
sx q[1];
rz(-0.7784673) q[1];
sx q[1];
rz(0.73164083) q[1];
x q[2];
rz(0.43059723) q[3];
sx q[3];
rz(-1.9816508) q[3];
sx q[3];
rz(-1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6106674) q[2];
sx q[2];
rz(-1.3767367) q[2];
sx q[2];
rz(-1.2923856) q[2];
rz(2.2409706) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(1.5952516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47386277) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(-1.5245755) q[0];
rz(1.698311) q[1];
sx q[1];
rz(-0.89569211) q[1];
sx q[1];
rz(0.5009833) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058280073) q[0];
sx q[0];
rz(-0.73375637) q[0];
sx q[0];
rz(-1.0007798) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6691241) q[2];
sx q[2];
rz(-0.61206619) q[2];
sx q[2];
rz(1.2304359) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9588291) q[1];
sx q[1];
rz(-1.0729959) q[1];
sx q[1];
rz(-2.2341841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4464842) q[3];
sx q[3];
rz(-1.2642702) q[3];
sx q[3];
rz(0.55858726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78918004) q[2];
sx q[2];
rz(-3.0173306) q[2];
sx q[2];
rz(-0.46335709) q[2];
rz(-0.067666791) q[3];
sx q[3];
rz(-1.8384408) q[3];
sx q[3];
rz(-0.1629924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25245923) q[0];
sx q[0];
rz(-1.2754138) q[0];
sx q[0];
rz(0.5109936) q[0];
rz(-0.64396089) q[1];
sx q[1];
rz(-2.0168346) q[1];
sx q[1];
rz(-2.8535829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32569186) q[0];
sx q[0];
rz(-1.3790695) q[0];
sx q[0];
rz(1.6098538) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2637469) q[2];
sx q[2];
rz(-1.0509509) q[2];
sx q[2];
rz(2.9016888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5254613) q[1];
sx q[1];
rz(-1.0122293) q[1];
sx q[1];
rz(1.8957047) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6828641) q[3];
sx q[3];
rz(-1.1945416) q[3];
sx q[3];
rz(-2.3628836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1945232) q[2];
sx q[2];
rz(-1.9788519) q[2];
sx q[2];
rz(-0.21491773) q[2];
rz(-1.8042709) q[3];
sx q[3];
rz(-2.603172) q[3];
sx q[3];
rz(-2.9609093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8682078) q[0];
sx q[0];
rz(-0.93623638) q[0];
sx q[0];
rz(2.0599763) q[0];
rz(-2.1514905) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(2.8066011) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8222447) q[0];
sx q[0];
rz(-1.8935504) q[0];
sx q[0];
rz(3.1060973) q[0];
rz(-pi) q[1];
x q[1];
rz(1.556646) q[2];
sx q[2];
rz(-2.7196537) q[2];
sx q[2];
rz(2.5463605) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8364063) q[1];
sx q[1];
rz(-1.4855774) q[1];
sx q[1];
rz(-3.042074) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5906628) q[3];
sx q[3];
rz(-2.7405313) q[3];
sx q[3];
rz(-0.10894081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0830903) q[2];
sx q[2];
rz(-0.52512705) q[2];
sx q[2];
rz(-2.7527909) q[2];
rz(0.36924103) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(0.84454876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19287547) q[0];
sx q[0];
rz(-2.188864) q[0];
sx q[0];
rz(-2.4396851) q[0];
rz(1.6096055) q[1];
sx q[1];
rz(-2.46789) q[1];
sx q[1];
rz(-0.38175499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4629342) q[0];
sx q[0];
rz(-1.6525499) q[0];
sx q[0];
rz(-1.6340748) q[0];
rz(-pi) q[1];
rz(-1.4456621) q[2];
sx q[2];
rz(-1.1627253) q[2];
sx q[2];
rz(3.0681075) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4873969) q[1];
sx q[1];
rz(-1.3302478) q[1];
sx q[1];
rz(-1.3122561) q[1];
x q[2];
rz(0.66017229) q[3];
sx q[3];
rz(-0.91554612) q[3];
sx q[3];
rz(-1.8546263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6751042) q[2];
sx q[2];
rz(-1.0305104) q[2];
sx q[2];
rz(2.4934736) q[2];
rz(-2.134792) q[3];
sx q[3];
rz(-1.1454134) q[3];
sx q[3];
rz(-3.0692611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5299912) q[0];
sx q[0];
rz(-1.7399104) q[0];
sx q[0];
rz(-2.4763784) q[0];
rz(-2.8499659) q[1];
sx q[1];
rz(-1.3529774) q[1];
sx q[1];
rz(-1.1381961) q[1];
rz(1.2164581) q[2];
sx q[2];
rz(-2.0632498) q[2];
sx q[2];
rz(-1.9634631) q[2];
rz(-1.2286366) q[3];
sx q[3];
rz(-2.6989326) q[3];
sx q[3];
rz(-1.1748675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
