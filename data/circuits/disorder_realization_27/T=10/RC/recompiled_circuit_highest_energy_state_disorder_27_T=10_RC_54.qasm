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
rz(-0.24212317) q[0];
sx q[0];
rz(-0.72532931) q[0];
sx q[0];
rz(0.90670937) q[0];
rz(2.6746305) q[1];
sx q[1];
rz(-0.90277201) q[1];
sx q[1];
rz(0.24388193) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8236602) q[0];
sx q[0];
rz(-2.3711304) q[0];
sx q[0];
rz(-2.5567358) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93893013) q[2];
sx q[2];
rz(-0.9900695) q[2];
sx q[2];
rz(1.1685358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0816312) q[1];
sx q[1];
rz(-1.3880127) q[1];
sx q[1];
rz(-1.6441397) q[1];
rz(-pi) q[2];
rz(-0.21190819) q[3];
sx q[3];
rz(-0.82900199) q[3];
sx q[3];
rz(-2.172989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4563518) q[2];
sx q[2];
rz(-2.6754003) q[2];
sx q[2];
rz(2.1166128) q[2];
rz(0.13925615) q[3];
sx q[3];
rz(-1.1956297) q[3];
sx q[3];
rz(-2.9774184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76992947) q[0];
sx q[0];
rz(-2.0491845) q[0];
sx q[0];
rz(-1.2051693) q[0];
rz(-1.4681762) q[1];
sx q[1];
rz(-1.877715) q[1];
sx q[1];
rz(1.7211627) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98200765) q[0];
sx q[0];
rz(-2.3075342) q[0];
sx q[0];
rz(-2.9247051) q[0];
x q[1];
rz(-1.7718138) q[2];
sx q[2];
rz(-2.8716785) q[2];
sx q[2];
rz(-0.52598276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.08733315) q[1];
sx q[1];
rz(-1.2337605) q[1];
sx q[1];
rz(2.9059306) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74052625) q[3];
sx q[3];
rz(-1.9364001) q[3];
sx q[3];
rz(-1.4359563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.18985192) q[2];
sx q[2];
rz(-2.9517089) q[2];
sx q[2];
rz(-2.3404549) q[2];
rz(0.43870157) q[3];
sx q[3];
rz(-1.8935685) q[3];
sx q[3];
rz(2.0285105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37757117) q[0];
sx q[0];
rz(-2.8442597) q[0];
sx q[0];
rz(1.1391033) q[0];
rz(2.8746919) q[1];
sx q[1];
rz(-1.8903939) q[1];
sx q[1];
rz(-2.7762754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6613839) q[0];
sx q[0];
rz(-1.669642) q[0];
sx q[0];
rz(0.5151202) q[0];
rz(-pi) q[1];
rz(2.5796125) q[2];
sx q[2];
rz(-2.1393074) q[2];
sx q[2];
rz(1.3696483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85161419) q[1];
sx q[1];
rz(-1.9383926) q[1];
sx q[1];
rz(-3.0795891) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85478304) q[3];
sx q[3];
rz(-2.0826591) q[3];
sx q[3];
rz(-1.3123684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2306564) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(-0.12348565) q[2];
rz(2.2423045) q[3];
sx q[3];
rz(-1.4495918) q[3];
sx q[3];
rz(1.9062769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7742291) q[0];
sx q[0];
rz(-1.4643359) q[0];
sx q[0];
rz(-0.88957077) q[0];
rz(-0.74596897) q[1];
sx q[1];
rz(-1.4624701) q[1];
sx q[1];
rz(-1.3847345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69834405) q[0];
sx q[0];
rz(-0.067509934) q[0];
sx q[0];
rz(-1.0288074) q[0];
x q[1];
rz(-2.8207702) q[2];
sx q[2];
rz(-2.4745382) q[2];
sx q[2];
rz(2.6503682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8306668) q[1];
sx q[1];
rz(-1.3841024) q[1];
sx q[1];
rz(0.6142389) q[1];
rz(-0.05352002) q[3];
sx q[3];
rz(-2.60703) q[3];
sx q[3];
rz(3.0381615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3978465) q[2];
sx q[2];
rz(-0.83796871) q[2];
sx q[2];
rz(2.1261334) q[2];
rz(-0.022653496) q[3];
sx q[3];
rz(-1.3709143) q[3];
sx q[3];
rz(3.0034351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37462336) q[0];
sx q[0];
rz(-1.8648819) q[0];
sx q[0];
rz(-0.20481566) q[0];
rz(-2.9264033) q[1];
sx q[1];
rz(-2.9209825) q[1];
sx q[1];
rz(-2.2664216) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78007215) q[0];
sx q[0];
rz(-2.8018537) q[0];
sx q[0];
rz(2.0302515) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4007447) q[2];
sx q[2];
rz(-2.399246) q[2];
sx q[2];
rz(-0.60630732) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2332747) q[1];
sx q[1];
rz(-1.868778) q[1];
sx q[1];
rz(-1.3402935) q[1];
x q[2];
rz(-2.2017893) q[3];
sx q[3];
rz(-1.4550337) q[3];
sx q[3];
rz(0.94652688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63336664) q[2];
sx q[2];
rz(-2.2534011) q[2];
sx q[2];
rz(1.7249031) q[2];
rz(-2.21777) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(-1.1317322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0210339) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-2.293204) q[0];
rz(1.5658762) q[1];
sx q[1];
rz(-1.3278241) q[1];
sx q[1];
rz(2.1956992) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1740055) q[0];
sx q[0];
rz(-1.8552836) q[0];
sx q[0];
rz(2.4422285) q[0];
rz(-0.20393238) q[2];
sx q[2];
rz(-0.82101594) q[2];
sx q[2];
rz(1.7199788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1887363) q[1];
sx q[1];
rz(-0.99299255) q[1];
sx q[1];
rz(2.0835447) q[1];
rz(0.35870798) q[3];
sx q[3];
rz(-0.51878319) q[3];
sx q[3];
rz(2.0174873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1762323) q[2];
sx q[2];
rz(-1.9051899) q[2];
sx q[2];
rz(-1.7410834) q[2];
rz(1.6598353) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(-0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7401212) q[0];
sx q[0];
rz(-2.6417702) q[0];
sx q[0];
rz(-1.2116145) q[0];
rz(2.6648193) q[1];
sx q[1];
rz(-1.8649273) q[1];
sx q[1];
rz(1.635294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20902987) q[0];
sx q[0];
rz(-1.8537921) q[0];
sx q[0];
rz(-2.1542633) q[0];
rz(-pi) q[1];
rz(-1.2477307) q[2];
sx q[2];
rz(-1.808721) q[2];
sx q[2];
rz(-1.0987154) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7249178) q[1];
sx q[1];
rz(-0.70958269) q[1];
sx q[1];
rz(2.9516944) q[1];
rz(1.8651857) q[3];
sx q[3];
rz(-2.0352425) q[3];
sx q[3];
rz(0.73976952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8022884) q[2];
sx q[2];
rz(-1.6243287) q[2];
sx q[2];
rz(-0.18292546) q[2];
rz(2.4798992) q[3];
sx q[3];
rz(-1.9293834) q[3];
sx q[3];
rz(2.4367512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45109192) q[0];
sx q[0];
rz(-1.6766312) q[0];
sx q[0];
rz(2.9943384) q[0];
rz(-2.9946949) q[1];
sx q[1];
rz(-2.2203827) q[1];
sx q[1];
rz(-1.9619092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3578276) q[0];
sx q[0];
rz(-1.9067147) q[0];
sx q[0];
rz(0.442105) q[0];
x q[1];
rz(-1.0919763) q[2];
sx q[2];
rz(-1.9648332) q[2];
sx q[2];
rz(-2.0872598) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.98683607) q[1];
sx q[1];
rz(-2.1636226) q[1];
sx q[1];
rz(-1.703771) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0249765) q[3];
sx q[3];
rz(-0.75748235) q[3];
sx q[3];
rz(2.1471786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94072056) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(-2.3737523) q[2];
rz(-2.614894) q[3];
sx q[3];
rz(-2.0898762) q[3];
sx q[3];
rz(-2.5832978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4621157) q[0];
sx q[0];
rz(-0.18590346) q[0];
sx q[0];
rz(2.0515077) q[0];
rz(1.7915626) q[1];
sx q[1];
rz(-1.4005902) q[1];
sx q[1];
rz(-0.42492351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23608828) q[0];
sx q[0];
rz(-0.47471913) q[0];
sx q[0];
rz(-1.4464054) q[0];
x q[1];
rz(0.64984729) q[2];
sx q[2];
rz(-1.51553) q[2];
sx q[2];
rz(-1.856232) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38389123) q[1];
sx q[1];
rz(-0.5846068) q[1];
sx q[1];
rz(-0.65442185) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5492776) q[3];
sx q[3];
rz(-0.80192536) q[3];
sx q[3];
rz(1.3682113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2316042) q[2];
sx q[2];
rz(-0.51077545) q[2];
sx q[2];
rz(-2.6166022) q[2];
rz(-1.7867583) q[3];
sx q[3];
rz(-0.89559186) q[3];
sx q[3];
rz(-1.7790599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43596426) q[0];
sx q[0];
rz(-0.66052496) q[0];
sx q[0];
rz(-3.0349162) q[0];
rz(-2.0872929) q[1];
sx q[1];
rz(-1.432447) q[1];
sx q[1];
rz(-1.7342825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46074902) q[0];
sx q[0];
rz(-1.1011194) q[0];
sx q[0];
rz(0.90523714) q[0];
x q[1];
rz(2.4881192) q[2];
sx q[2];
rz(-1.1879041) q[2];
sx q[2];
rz(-3.058397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6165054) q[1];
sx q[1];
rz(-2.8697851) q[1];
sx q[1];
rz(2.1735783) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94707113) q[3];
sx q[3];
rz(-0.92309381) q[3];
sx q[3];
rz(-0.050681678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6537031) q[2];
sx q[2];
rz(-2.6370878) q[2];
sx q[2];
rz(-2.8686236) q[2];
rz(-2.2702787) q[3];
sx q[3];
rz(-1.6815691) q[3];
sx q[3];
rz(0.13450809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3684261) q[0];
sx q[0];
rz(-1.9482524) q[0];
sx q[0];
rz(0.85309749) q[0];
rz(0.38289616) q[1];
sx q[1];
rz(-1.2877512) q[1];
sx q[1];
rz(3.126694) q[1];
rz(-2.0128925) q[2];
sx q[2];
rz(-2.0876035) q[2];
sx q[2];
rz(0.58135396) q[2];
rz(-0.79991585) q[3];
sx q[3];
rz(-1.6408995) q[3];
sx q[3];
rz(-1.2560918) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
