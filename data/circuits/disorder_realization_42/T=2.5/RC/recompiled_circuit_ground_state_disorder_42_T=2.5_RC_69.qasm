OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3990134) q[0];
sx q[0];
rz(-1.1630381) q[0];
sx q[0];
rz(2.0030339) q[0];
rz(2.5055655) q[1];
sx q[1];
rz(-2.6316402) q[1];
sx q[1];
rz(2.5165601) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.807686) q[0];
sx q[0];
rz(-1.429883) q[0];
sx q[0];
rz(1.5265092) q[0];
rz(-1.5364067) q[2];
sx q[2];
rz(-2.6011438) q[2];
sx q[2];
rz(-2.2480223) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33746943) q[1];
sx q[1];
rz(-1.1589246) q[1];
sx q[1];
rz(-2.1571674) q[1];
rz(0.19842438) q[3];
sx q[3];
rz(-2.3188667) q[3];
sx q[3];
rz(-2.2446475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92224377) q[2];
sx q[2];
rz(-0.85180989) q[2];
sx q[2];
rz(2.7395978) q[2];
rz(2.3123907) q[3];
sx q[3];
rz(-1.821937) q[3];
sx q[3];
rz(1.7792262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7336693) q[0];
sx q[0];
rz(-0.55183721) q[0];
sx q[0];
rz(-2.0358987) q[0];
rz(0.78877527) q[1];
sx q[1];
rz(-2.5194247) q[1];
sx q[1];
rz(2.6355991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5010586) q[0];
sx q[0];
rz(-1.6468684) q[0];
sx q[0];
rz(-2.9058558) q[0];
rz(-2.3680192) q[2];
sx q[2];
rz(-2.1974652) q[2];
sx q[2];
rz(-0.0093912436) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2179394) q[1];
sx q[1];
rz(-2.7193177) q[1];
sx q[1];
rz(2.1417066) q[1];
x q[2];
rz(-0.45428975) q[3];
sx q[3];
rz(-1.807782) q[3];
sx q[3];
rz(-1.911357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2371696) q[2];
sx q[2];
rz(-1.0639031) q[2];
sx q[2];
rz(-0.87192956) q[2];
rz(2.0875841) q[3];
sx q[3];
rz(-1.0796615) q[3];
sx q[3];
rz(-1.7006629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0525518) q[0];
sx q[0];
rz(-2.4854923) q[0];
sx q[0];
rz(-1.9144527) q[0];
rz(-0.50651208) q[1];
sx q[1];
rz(-2.3498693) q[1];
sx q[1];
rz(-0.92481771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4445405) q[0];
sx q[0];
rz(-0.18029515) q[0];
sx q[0];
rz(-0.96526115) q[0];
x q[1];
rz(-1.501899) q[2];
sx q[2];
rz(-2.0259948) q[2];
sx q[2];
rz(-2.8719939) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9021685) q[1];
sx q[1];
rz(-1.268506) q[1];
sx q[1];
rz(0.92642529) q[1];
x q[2];
rz(-0.4731725) q[3];
sx q[3];
rz(-0.71801502) q[3];
sx q[3];
rz(1.1188467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0792803) q[2];
sx q[2];
rz(-1.5996876) q[2];
sx q[2];
rz(-0.49528948) q[2];
rz(1.7700206) q[3];
sx q[3];
rz(-2.3973231) q[3];
sx q[3];
rz(1.3220968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1247193) q[0];
sx q[0];
rz(-1.9833516) q[0];
sx q[0];
rz(2.25368) q[0];
rz(1.7601695) q[1];
sx q[1];
rz(-1.0180749) q[1];
sx q[1];
rz(-2.6240614) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4221188) q[0];
sx q[0];
rz(-1.5567915) q[0];
sx q[0];
rz(-1.3097358) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4231269) q[2];
sx q[2];
rz(-1.9794165) q[2];
sx q[2];
rz(1.1312616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0486517) q[1];
sx q[1];
rz(-3.042964) q[1];
sx q[1];
rz(-1.4709378) q[1];
rz(-pi) q[2];
rz(-2.4740287) q[3];
sx q[3];
rz(-2.5680411) q[3];
sx q[3];
rz(1.9907601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.069328221) q[2];
sx q[2];
rz(-1.1257409) q[2];
sx q[2];
rz(-2.1605055) q[2];
rz(0.4782933) q[3];
sx q[3];
rz(-2.7893453) q[3];
sx q[3];
rz(0.52811629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9454055) q[0];
sx q[0];
rz(-1.4237175) q[0];
sx q[0];
rz(0.032935306) q[0];
rz(-1.5252652) q[1];
sx q[1];
rz(-2.567629) q[1];
sx q[1];
rz(-0.93564916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0230377) q[0];
sx q[0];
rz(-1.584125) q[0];
sx q[0];
rz(2.6715735) q[0];
rz(0.79765908) q[2];
sx q[2];
rz(-1.2204542) q[2];
sx q[2];
rz(1.2454741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0914382) q[1];
sx q[1];
rz(-2.3622516) q[1];
sx q[1];
rz(-3.0528751) q[1];
rz(1.3282975) q[3];
sx q[3];
rz(-2.1418013) q[3];
sx q[3];
rz(-1.5712488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66964883) q[2];
sx q[2];
rz(-2.7574597) q[2];
sx q[2];
rz(-1.6432537) q[2];
rz(-2.5390427) q[3];
sx q[3];
rz(-2.3157401) q[3];
sx q[3];
rz(-1.92441) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5947386) q[0];
sx q[0];
rz(-2.8475519) q[0];
sx q[0];
rz(-3.102741) q[0];
rz(0.36000577) q[1];
sx q[1];
rz(-1.729676) q[1];
sx q[1];
rz(-0.14883277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94939828) q[0];
sx q[0];
rz(-1.697505) q[0];
sx q[0];
rz(-1.4437946) q[0];
rz(-pi) q[1];
rz(-2.932933) q[2];
sx q[2];
rz(-2.5620915) q[2];
sx q[2];
rz(2.8148871) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56797261) q[1];
sx q[1];
rz(-0.96155969) q[1];
sx q[1];
rz(-2.1813758) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5555796) q[3];
sx q[3];
rz(-1.3162321) q[3];
sx q[3];
rz(2.3814122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0570602) q[2];
sx q[2];
rz(-0.82814211) q[2];
sx q[2];
rz(0.35489902) q[2];
rz(2.0830294) q[3];
sx q[3];
rz(-0.94616977) q[3];
sx q[3];
rz(0.53275776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31838378) q[0];
sx q[0];
rz(-1.7192778) q[0];
sx q[0];
rz(2.3329155) q[0];
rz(3.0478364) q[1];
sx q[1];
rz(-0.84051991) q[1];
sx q[1];
rz(-2.7309928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79390271) q[0];
sx q[0];
rz(-1.9436033) q[0];
sx q[0];
rz(-1.1469141) q[0];
rz(-pi) q[1];
rz(1.3099395) q[2];
sx q[2];
rz(-0.93376389) q[2];
sx q[2];
rz(2.2230679) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0714662) q[1];
sx q[1];
rz(-1.1845329) q[1];
sx q[1];
rz(1.1995308) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15488829) q[3];
sx q[3];
rz(-2.4714176) q[3];
sx q[3];
rz(2.6763499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14313301) q[2];
sx q[2];
rz(-0.3930347) q[2];
sx q[2];
rz(-0.03579363) q[2];
rz(-1.834747) q[3];
sx q[3];
rz(-1.9948317) q[3];
sx q[3];
rz(0.80865639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3150113) q[0];
sx q[0];
rz(-2.1666574) q[0];
sx q[0];
rz(-0.63631979) q[0];
rz(2.4782205) q[1];
sx q[1];
rz(-1.1095108) q[1];
sx q[1];
rz(2.5136307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7610891) q[0];
sx q[0];
rz(-0.081898339) q[0];
sx q[0];
rz(-2.7837672) q[0];
rz(0.2070932) q[2];
sx q[2];
rz(-0.40488714) q[2];
sx q[2];
rz(1.5804497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.99599528) q[1];
sx q[1];
rz(-0.72628262) q[1];
sx q[1];
rz(-0.65205814) q[1];
x q[2];
rz(0.18744882) q[3];
sx q[3];
rz(-1.6050242) q[3];
sx q[3];
rz(-1.6899861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5593354) q[2];
sx q[2];
rz(-0.11056837) q[2];
sx q[2];
rz(-1.4165233) q[2];
rz(-2.0420117) q[3];
sx q[3];
rz(-2.1018335) q[3];
sx q[3];
rz(-2.288868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13903862) q[0];
sx q[0];
rz(-2.250493) q[0];
sx q[0];
rz(2.52676) q[0];
rz(-0.50827208) q[1];
sx q[1];
rz(-0.95242396) q[1];
sx q[1];
rz(0.27613786) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50867453) q[0];
sx q[0];
rz(-2.6946697) q[0];
sx q[0];
rz(1.8742976) q[0];
x q[1];
rz(-0.15900826) q[2];
sx q[2];
rz(-1.0463042) q[2];
sx q[2];
rz(0.25689143) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4003168) q[1];
sx q[1];
rz(-1.1172677) q[1];
sx q[1];
rz(-0.67429115) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7998447) q[3];
sx q[3];
rz(-2.2015988) q[3];
sx q[3];
rz(-2.3113826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95163661) q[2];
sx q[2];
rz(-1.0389682) q[2];
sx q[2];
rz(0.098380066) q[2];
rz(-1.9010057) q[3];
sx q[3];
rz(-2.0799347) q[3];
sx q[3];
rz(-0.043717535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8561309) q[0];
sx q[0];
rz(-1.0482482) q[0];
sx q[0];
rz(-2.7833126) q[0];
rz(-1.0048535) q[1];
sx q[1];
rz(-0.88468164) q[1];
sx q[1];
rz(-0.34214941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44624871) q[0];
sx q[0];
rz(-0.98606743) q[0];
sx q[0];
rz(1.8323932) q[0];
rz(2.1120295) q[2];
sx q[2];
rz(-2.4839253) q[2];
sx q[2];
rz(2.5785411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0562065) q[1];
sx q[1];
rz(-1.6780938) q[1];
sx q[1];
rz(-2.5683968) q[1];
x q[2];
rz(1.9072745) q[3];
sx q[3];
rz(-1.668317) q[3];
sx q[3];
rz(0.3995801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3442012) q[2];
sx q[2];
rz(-2.6579393) q[2];
sx q[2];
rz(0.3717711) q[2];
rz(2.9504635) q[3];
sx q[3];
rz(-1.976795) q[3];
sx q[3];
rz(-0.42310664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38901781) q[0];
sx q[0];
rz(-0.49674635) q[0];
sx q[0];
rz(2.4865271) q[0];
rz(-2.7816506) q[1];
sx q[1];
rz(-1.5725726) q[1];
sx q[1];
rz(-1.6307065) q[1];
rz(3.0790569) q[2];
sx q[2];
rz(-1.2690496) q[2];
sx q[2];
rz(-3.0435111) q[2];
rz(-3.1265518) q[3];
sx q[3];
rz(-0.81285601) q[3];
sx q[3];
rz(2.6668919) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
