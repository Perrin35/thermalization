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
rz(-0.51195872) q[0];
sx q[0];
rz(-2.7380138) q[0];
sx q[0];
rz(3.1374748) q[0];
rz(0.68139684) q[1];
sx q[1];
rz(-1.2702785) q[1];
sx q[1];
rz(1.4407925) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4096298) q[0];
sx q[0];
rz(-2.8131631) q[0];
sx q[0];
rz(-1.0148125) q[0];
rz(-pi) q[1];
rz(2.7782562) q[2];
sx q[2];
rz(-2.3564995) q[2];
sx q[2];
rz(1.3183644) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8358678) q[1];
sx q[1];
rz(-0.9192217) q[1];
sx q[1];
rz(-1.1206133) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5020237) q[3];
sx q[3];
rz(-2.9675238) q[3];
sx q[3];
rz(-0.93910142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.19930856) q[2];
sx q[2];
rz(-1.0465304) q[2];
sx q[2];
rz(-2.671177) q[2];
rz(2.2453902) q[3];
sx q[3];
rz(-1.4678518) q[3];
sx q[3];
rz(1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6386221) q[0];
sx q[0];
rz(-2.6574385) q[0];
sx q[0];
rz(2.7726987) q[0];
rz(-2.6781354) q[1];
sx q[1];
rz(-1.8577441) q[1];
sx q[1];
rz(0.2643815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6579012) q[0];
sx q[0];
rz(-0.92934791) q[0];
sx q[0];
rz(0.19199706) q[0];
rz(-1.7272378) q[2];
sx q[2];
rz(-1.3231734) q[2];
sx q[2];
rz(2.3231017) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0708587) q[1];
sx q[1];
rz(-1.8740591) q[1];
sx q[1];
rz(1.2649078) q[1];
x q[2];
rz(-0.78231298) q[3];
sx q[3];
rz(-1.6564329) q[3];
sx q[3];
rz(-0.2667564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5276864) q[2];
sx q[2];
rz(-2.9759585) q[2];
sx q[2];
rz(-1.8370834) q[2];
rz(-2.8218609) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(2.6379697) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17427915) q[0];
sx q[0];
rz(-0.19198424) q[0];
sx q[0];
rz(-1.1235896) q[0];
rz(0.38791052) q[1];
sx q[1];
rz(-1.8759517) q[1];
sx q[1];
rz(2.8111615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200121) q[0];
sx q[0];
rz(-0.3191443) q[0];
sx q[0];
rz(-0.97346755) q[0];
rz(3.0082012) q[2];
sx q[2];
rz(-2.5299604) q[2];
sx q[2];
rz(-0.69452121) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6092005) q[1];
sx q[1];
rz(-0.48621854) q[1];
sx q[1];
rz(2.4052252) q[1];
rz(0.73883812) q[3];
sx q[3];
rz(-2.3519302) q[3];
sx q[3];
rz(-1.9146321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.419751) q[2];
sx q[2];
rz(-1.1417049) q[2];
sx q[2];
rz(-1.9640131) q[2];
rz(1.8925331) q[3];
sx q[3];
rz(-1.8355337) q[3];
sx q[3];
rz(0.038399847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881507) q[0];
sx q[0];
rz(-1.7120687) q[0];
sx q[0];
rz(-0.088726774) q[0];
rz(2.7170722) q[1];
sx q[1];
rz(-0.71966925) q[1];
sx q[1];
rz(-1.7321865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30520327) q[0];
sx q[0];
rz(-0.38581784) q[0];
sx q[0];
rz(2.0842537) q[0];
rz(-pi) q[1];
x q[1];
rz(2.652651) q[2];
sx q[2];
rz(-2.1132815) q[2];
sx q[2];
rz(-2.9352293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3864435) q[1];
sx q[1];
rz(-2.1514822) q[1];
sx q[1];
rz(-2.0232852) q[1];
rz(-pi) q[2];
rz(-2.0982101) q[3];
sx q[3];
rz(-0.58660117) q[3];
sx q[3];
rz(2.9709904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38922629) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(0.43238315) q[2];
rz(-0.45807517) q[3];
sx q[3];
rz(-1.4832486) q[3];
sx q[3];
rz(-2.1336011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78419375) q[0];
sx q[0];
rz(-0.14506871) q[0];
sx q[0];
rz(2.8302637) q[0];
rz(1.9642824) q[1];
sx q[1];
rz(-1.177634) q[1];
sx q[1];
rz(0.48431531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70894079) q[0];
sx q[0];
rz(-1.4756894) q[0];
sx q[0];
rz(-1.7300138) q[0];
rz(1.3179146) q[2];
sx q[2];
rz(-0.30159471) q[2];
sx q[2];
rz(2.4689134) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97515177) q[1];
sx q[1];
rz(-1.7130392) q[1];
sx q[1];
rz(1.1935545) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9919162) q[3];
sx q[3];
rz(-1.6972741) q[3];
sx q[3];
rz(-1.6320796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3207265) q[2];
sx q[2];
rz(-1.7118688) q[2];
sx q[2];
rz(-2.7166264) q[2];
rz(3.037437) q[3];
sx q[3];
rz(-2.7724373) q[3];
sx q[3];
rz(-2.4764376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686309) q[0];
sx q[0];
rz(-0.63987982) q[0];
sx q[0];
rz(-0.1567008) q[0];
rz(-2.8796097) q[1];
sx q[1];
rz(-1.1527088) q[1];
sx q[1];
rz(-1.461747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77685415) q[0];
sx q[0];
rz(-2.1663453) q[0];
sx q[0];
rz(-2.9315346) q[0];
rz(1.842672) q[2];
sx q[2];
rz(-2.4266234) q[2];
sx q[2];
rz(-0.78503099) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20355454) q[1];
sx q[1];
rz(-2.9239836) q[1];
sx q[1];
rz(-0.59534351) q[1];
x q[2];
rz(-0.42243345) q[3];
sx q[3];
rz(-0.79630536) q[3];
sx q[3];
rz(0.56047201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4713952) q[2];
sx q[2];
rz(-1.6452226) q[2];
sx q[2];
rz(0.71693286) q[2];
rz(0.21643058) q[3];
sx q[3];
rz(-1.2856893) q[3];
sx q[3];
rz(1.9560248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7316498) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(1.2662079) q[0];
rz(-0.75824291) q[1];
sx q[1];
rz(-1.9634602) q[1];
sx q[1];
rz(-2.0936802) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58430481) q[0];
sx q[0];
rz(-2.1471229) q[0];
sx q[0];
rz(0.35741522) q[0];
x q[1];
rz(-2.2323147) q[2];
sx q[2];
rz(-2.4729441) q[2];
sx q[2];
rz(0.74597154) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0985154) q[1];
sx q[1];
rz(-1.9334353) q[1];
sx q[1];
rz(0.074438268) q[1];
rz(-2.6664958) q[3];
sx q[3];
rz(-0.49834278) q[3];
sx q[3];
rz(-1.4628177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3332112) q[2];
sx q[2];
rz(-2.1344678) q[2];
sx q[2];
rz(-3.132931) q[2];
rz(-1.0041142) q[3];
sx q[3];
rz(-2.6486371) q[3];
sx q[3];
rz(1.0021817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36326161) q[0];
sx q[0];
rz(-0.14334981) q[0];
sx q[0];
rz(-2.1891201) q[0];
rz(-1.5337503) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(2.105377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0909611) q[0];
sx q[0];
rz(-1.0337752) q[0];
sx q[0];
rz(-1.0616674) q[0];
x q[1];
rz(-0.63663738) q[2];
sx q[2];
rz(-1.8701866) q[2];
sx q[2];
rz(-2.4518397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0091102) q[1];
sx q[1];
rz(-0.81071893) q[1];
sx q[1];
rz(2.2894163) q[1];
rz(-pi) q[2];
rz(1.1187068) q[3];
sx q[3];
rz(-2.2452659) q[3];
sx q[3];
rz(1.3687641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8448392) q[2];
sx q[2];
rz(-2.4256458) q[2];
sx q[2];
rz(1.0003662) q[2];
rz(1.2760466) q[3];
sx q[3];
rz(-1.6701271) q[3];
sx q[3];
rz(2.0744417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.100383) q[0];
sx q[0];
rz(-2.2333133) q[0];
sx q[0];
rz(-0.25417438) q[0];
rz(-1.8036448) q[1];
sx q[1];
rz(-2.281052) q[1];
sx q[1];
rz(-0.77802229) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8786273) q[0];
sx q[0];
rz(-2.3926982) q[0];
sx q[0];
rz(-0.48567943) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7341975) q[2];
sx q[2];
rz(-1.7694574) q[2];
sx q[2];
rz(-0.29115788) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48315037) q[1];
sx q[1];
rz(-1.9772234) q[1];
sx q[1];
rz(-0.61950018) q[1];
rz(-1.7860402) q[3];
sx q[3];
rz(-1.4741033) q[3];
sx q[3];
rz(-1.4002864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.079387) q[2];
sx q[2];
rz(-1.5176682) q[2];
sx q[2];
rz(-3.0089231) q[2];
rz(2.0343871) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(1.4440822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4001813) q[0];
sx q[0];
rz(-1.1962698) q[0];
sx q[0];
rz(-0.74914002) q[0];
rz(2.6465042) q[1];
sx q[1];
rz(-1.9063213) q[1];
sx q[1];
rz(-1.3912158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063569) q[0];
sx q[0];
rz(-2.0208404) q[0];
sx q[0];
rz(-2.4343632) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83413627) q[2];
sx q[2];
rz(-2.5747402) q[2];
sx q[2];
rz(-1.2421654) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0616845) q[1];
sx q[1];
rz(-0.40929738) q[1];
sx q[1];
rz(-0.10959919) q[1];
x q[2];
rz(-2.6364348) q[3];
sx q[3];
rz(-1.0978175) q[3];
sx q[3];
rz(1.5947756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64080015) q[2];
sx q[2];
rz(-1.6668789) q[2];
sx q[2];
rz(-0.87727171) q[2];
rz(2.6194465) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(1.4570025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40182879) q[0];
sx q[0];
rz(-0.54720989) q[0];
sx q[0];
rz(-1.3450958) q[0];
rz(1.8632035) q[1];
sx q[1];
rz(-0.31247333) q[1];
sx q[1];
rz(-0.023275274) q[1];
rz(0.042365778) q[2];
sx q[2];
rz(-0.63127098) q[2];
sx q[2];
rz(1.9628738) q[2];
rz(2.6865339) q[3];
sx q[3];
rz(-1.2190335) q[3];
sx q[3];
rz(2.1608756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
