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
rz(1.2122756) q[0];
sx q[0];
rz(4.2733856) q[0];
sx q[0];
rz(10.797664) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(-1.7665266) q[1];
sx q[1];
rz(2.3553203) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7776913) q[0];
sx q[0];
rz(-0.4253952) q[0];
sx q[0];
rz(1.0546272) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2271871) q[2];
sx q[2];
rz(-1.1160276) q[2];
sx q[2];
rz(0.17780534) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3226124) q[1];
sx q[1];
rz(-0.29085813) q[1];
sx q[1];
rz(1.4664654) q[1];
rz(-pi) q[2];
rz(-2.1699675) q[3];
sx q[3];
rz(-1.2279945) q[3];
sx q[3];
rz(2.9620217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0218411) q[2];
sx q[2];
rz(-1.7181052) q[2];
sx q[2];
rz(1.8915668) q[2];
rz(-0.73578468) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(1.5031987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4982872) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(0.071320891) q[0];
rz(-1.1530863) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(0.52871314) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.816975) q[0];
sx q[0];
rz(-1.2551089) q[0];
sx q[0];
rz(-1.0782918) q[0];
rz(-pi) q[1];
rz(-2.2430482) q[2];
sx q[2];
rz(-0.29783422) q[2];
sx q[2];
rz(-0.98984776) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9829927) q[1];
sx q[1];
rz(-2.0865177) q[1];
sx q[1];
rz(2.2186568) q[1];
rz(1.9857446) q[3];
sx q[3];
rz(-1.7768716) q[3];
sx q[3];
rz(1.8725439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66136393) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(0.028623494) q[2];
rz(-2.7875767) q[3];
sx q[3];
rz(-2.386644) q[3];
sx q[3];
rz(0.55907512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4383168) q[0];
sx q[0];
rz(-2.1385312) q[0];
sx q[0];
rz(-0.89135528) q[0];
rz(-0.50093961) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(-3.0827789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0552154) q[0];
sx q[0];
rz(-3.0432467) q[0];
sx q[0];
rz(0.50758657) q[0];
rz(-0.28532736) q[2];
sx q[2];
rz(-2.0801736) q[2];
sx q[2];
rz(2.2338108) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0173024) q[1];
sx q[1];
rz(-2.2180386) q[1];
sx q[1];
rz(0.94396293) q[1];
x q[2];
rz(0.18043442) q[3];
sx q[3];
rz(-1.5431607) q[3];
sx q[3];
rz(0.35865006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3567051) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(0.2365665) q[2];
rz(-0.92075721) q[3];
sx q[3];
rz(-1.0509138) q[3];
sx q[3];
rz(1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22223602) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(-0.87523571) q[0];
rz(1.5454166) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(-0.045305591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169679) q[0];
sx q[0];
rz(-0.12147203) q[0];
sx q[0];
rz(-1.0018906) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95618177) q[2];
sx q[2];
rz(-0.76757694) q[2];
sx q[2];
rz(2.0337348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6958624) q[1];
sx q[1];
rz(-0.67477422) q[1];
sx q[1];
rz(1.4314907) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8604467) q[3];
sx q[3];
rz(-2.5389606) q[3];
sx q[3];
rz(1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8594325) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(1.947594) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(0.94849006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10876656) q[0];
sx q[0];
rz(-2.325401) q[0];
sx q[0];
rz(2.4591675) q[0];
rz(1.2884864) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.3312181) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2829943) q[0];
sx q[0];
rz(-0.72688519) q[0];
sx q[0];
rz(-1.5461393) q[0];
rz(-2.140215) q[2];
sx q[2];
rz(-0.45854202) q[2];
sx q[2];
rz(1.4024782) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.694931) q[1];
sx q[1];
rz(-0.38685054) q[1];
sx q[1];
rz(-0.37987654) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0430452) q[3];
sx q[3];
rz(-1.1748136) q[3];
sx q[3];
rz(2.2137303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1149301) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(0.66693532) q[2];
rz(0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(-2.8426898) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894593) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(-0.70372787) q[0];
rz(2.9723889) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(-0.15883787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7089091) q[0];
sx q[0];
rz(-1.7129717) q[0];
sx q[0];
rz(-0.12169038) q[0];
rz(-pi) q[1];
rz(2.4689697) q[2];
sx q[2];
rz(-2.9458698) q[2];
sx q[2];
rz(-1.4162404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1353246) q[1];
sx q[1];
rz(-2.6416831) q[1];
sx q[1];
rz(2.8630524) q[1];
rz(-pi) q[2];
rz(-2.8377526) q[3];
sx q[3];
rz(-0.66038495) q[3];
sx q[3];
rz(2.8840898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3682897) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(2.2840195) q[2];
rz(1.1693906) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17539772) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(2.4229557) q[0];
rz(-0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(-0.66863543) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.74092) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(-0.4192062) q[0];
rz(-0.35923194) q[2];
sx q[2];
rz(-2.1678535) q[2];
sx q[2];
rz(-3.0199241) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23827502) q[1];
sx q[1];
rz(-1.2438477) q[1];
sx q[1];
rz(1.1823149) q[1];
rz(-3.1304006) q[3];
sx q[3];
rz(-2.5103223) q[3];
sx q[3];
rz(0.61345631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64688993) q[2];
sx q[2];
rz(-1.1855482) q[2];
sx q[2];
rz(2.1051712) q[2];
rz(-2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(-2.3488267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(-0.72599894) q[0];
rz(1.1622102) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(1.1964218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0045753) q[0];
sx q[0];
rz(-1.5423623) q[0];
sx q[0];
rz(1.7445376) q[0];
rz(-pi) q[1];
rz(1.4216719) q[2];
sx q[2];
rz(-2.8310378) q[2];
sx q[2];
rz(-0.20937411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70564044) q[1];
sx q[1];
rz(-1.1060103) q[1];
sx q[1];
rz(-2.57354) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2231908) q[3];
sx q[3];
rz(-2.2112527) q[3];
sx q[3];
rz(2.9909822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90765816) q[2];
sx q[2];
rz(-1.5479167) q[2];
sx q[2];
rz(0.93351239) q[2];
rz(2.524611) q[3];
sx q[3];
rz(-1.0022997) q[3];
sx q[3];
rz(-2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1652949) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(3.0883375) q[0];
rz(-2.8540197) q[1];
sx q[1];
rz(-2.2122999) q[1];
sx q[1];
rz(-0.25921777) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1179781) q[0];
sx q[0];
rz(-0.058246944) q[0];
sx q[0];
rz(1.9969166) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0087148) q[2];
sx q[2];
rz(-1.8009552) q[2];
sx q[2];
rz(-1.8425737) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51356835) q[1];
sx q[1];
rz(-1.5092938) q[1];
sx q[1];
rz(-0.81925591) q[1];
rz(0.15786981) q[3];
sx q[3];
rz(-0.77220687) q[3];
sx q[3];
rz(-1.3881573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5549434) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-2.9602236) q[2];
rz(-0.3950611) q[3];
sx q[3];
rz(-0.22160465) q[3];
sx q[3];
rz(-0.980353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3996537) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(-2.7594866) q[0];
rz(2.8451879) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(1.2688676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5345602) q[0];
sx q[0];
rz(-1.1447313) q[0];
sx q[0];
rz(2.0924278) q[0];
rz(-pi) q[1];
rz(1.4051132) q[2];
sx q[2];
rz(-0.68000092) q[2];
sx q[2];
rz(0.14762893) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2603392) q[1];
sx q[1];
rz(-1.2311544) q[1];
sx q[1];
rz(-0.022625523) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4488698) q[3];
sx q[3];
rz(-1.1200179) q[3];
sx q[3];
rz(-1.4565695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84539139) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(-2.7517547) q[2];
rz(-1.8440394) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98904499) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(0.46089725) q[2];
sx q[2];
rz(-2.1620038) q[2];
sx q[2];
rz(2.9481859) q[2];
rz(1.9080347) q[3];
sx q[3];
rz(-0.81546406) q[3];
sx q[3];
rz(-2.245695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
