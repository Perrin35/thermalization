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
rz(1.3923378) q[0];
sx q[0];
rz(-2.7633986) q[0];
sx q[0];
rz(-0.35051546) q[0];
rz(0.89138436) q[1];
sx q[1];
rz(3.9337629) q[1];
sx q[1];
rz(13.73929) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.354686) q[0];
sx q[0];
rz(-1.3808492) q[0];
sx q[0];
rz(2.9482916) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1831398) q[2];
sx q[2];
rz(-0.26587379) q[2];
sx q[2];
rz(2.079351) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4306372) q[1];
sx q[1];
rz(-1.733755) q[1];
sx q[1];
rz(0.32117543) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7426026) q[3];
sx q[3];
rz(-1.1049924) q[3];
sx q[3];
rz(-1.4143296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9709836) q[2];
sx q[2];
rz(-1.8656518) q[2];
sx q[2];
rz(-0.21273908) q[2];
rz(-0.085266026) q[3];
sx q[3];
rz(-1.8906967) q[3];
sx q[3];
rz(1.8703478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79412115) q[0];
sx q[0];
rz(-2.329282) q[0];
sx q[0];
rz(2.0982657) q[0];
rz(0.13283816) q[1];
sx q[1];
rz(-2.9265407) q[1];
sx q[1];
rz(-1.9827693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3029859) q[0];
sx q[0];
rz(-1.9667454) q[0];
sx q[0];
rz(1.7274169) q[0];
rz(-pi) q[1];
rz(-0.01387502) q[2];
sx q[2];
rz(-1.4094947) q[2];
sx q[2];
rz(-2.814295) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1082127) q[1];
sx q[1];
rz(-1.8361109) q[1];
sx q[1];
rz(2.6379774) q[1];
x q[2];
rz(-0.93166931) q[3];
sx q[3];
rz(-0.51995819) q[3];
sx q[3];
rz(0.66321841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9449238) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(-0.72466737) q[2];
rz(-0.66257462) q[3];
sx q[3];
rz(-2.1734838) q[3];
sx q[3];
rz(-0.29964963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.9579983) q[0];
sx q[0];
rz(-1.2677001) q[0];
sx q[0];
rz(-0.9147574) q[0];
rz(-2.5547408) q[1];
sx q[1];
rz(-1.7210759) q[1];
sx q[1];
rz(2.9920726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0256074) q[0];
sx q[0];
rz(-1.7113026) q[0];
sx q[0];
rz(2.0951659) q[0];
x q[1];
rz(-2.3517866) q[2];
sx q[2];
rz(-2.6457204) q[2];
sx q[2];
rz(-2.7117031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4853277) q[1];
sx q[1];
rz(-2.1187003) q[1];
sx q[1];
rz(2.5119971) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87823509) q[3];
sx q[3];
rz(-2.1562169) q[3];
sx q[3];
rz(0.4957605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2933423) q[2];
sx q[2];
rz(-1.205227) q[2];
sx q[2];
rz(-0.80043522) q[2];
rz(-0.46659255) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(-2.485937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34224299) q[0];
sx q[0];
rz(-1.290134) q[0];
sx q[0];
rz(2.2829862) q[0];
rz(-0.41060064) q[1];
sx q[1];
rz(-2.938439) q[1];
sx q[1];
rz(2.1536749) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2938474) q[0];
sx q[0];
rz(-1.721707) q[0];
sx q[0];
rz(1.9035089) q[0];
rz(-0.12148492) q[2];
sx q[2];
rz(-1.0787691) q[2];
sx q[2];
rz(-0.58771261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6776442) q[1];
sx q[1];
rz(-1.8419187) q[1];
sx q[1];
rz(-1.8810924) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8186422) q[3];
sx q[3];
rz(-1.3846372) q[3];
sx q[3];
rz(1.8017839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9504488) q[2];
sx q[2];
rz(-1.7473651) q[2];
sx q[2];
rz(2.6960755) q[2];
rz(-0.71349239) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(0.92002404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.54818654) q[0];
sx q[0];
rz(-1.0819409) q[0];
sx q[0];
rz(1.5418381) q[0];
rz(1.8707188) q[1];
sx q[1];
rz(-1.4022695) q[1];
sx q[1];
rz(-0.52938968) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7177757) q[0];
sx q[0];
rz(-2.7090906) q[0];
sx q[0];
rz(2.1958817) q[0];
rz(-1.7725792) q[2];
sx q[2];
rz(-2.1786961) q[2];
sx q[2];
rz(2.4948134) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31570753) q[1];
sx q[1];
rz(-1.8216019) q[1];
sx q[1];
rz(-0.52108379) q[1];
rz(-pi) q[2];
rz(1.2182203) q[3];
sx q[3];
rz(-1.7099172) q[3];
sx q[3];
rz(-1.9162045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9993837) q[2];
sx q[2];
rz(-1.9500407) q[2];
sx q[2];
rz(-1.3231529) q[2];
rz(0.38149825) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(2.748446) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2926272) q[0];
sx q[0];
rz(-2.3830074) q[0];
sx q[0];
rz(1.0211771) q[0];
rz(-1.609833) q[1];
sx q[1];
rz(-0.94667089) q[1];
sx q[1];
rz(-2.3381332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96413104) q[0];
sx q[0];
rz(-1.877907) q[0];
sx q[0];
rz(1.403128) q[0];
x q[1];
rz(-2.7293936) q[2];
sx q[2];
rz(-1.4934818) q[2];
sx q[2];
rz(-0.22082034) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66113055) q[1];
sx q[1];
rz(-1.4509038) q[1];
sx q[1];
rz(0.57173034) q[1];
rz(-1.6709953) q[3];
sx q[3];
rz(-2.1952944) q[3];
sx q[3];
rz(1.9621094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44567406) q[2];
sx q[2];
rz(-2.5542104) q[2];
sx q[2];
rz(2.7109801) q[2];
rz(0.63498354) q[3];
sx q[3];
rz(-3.1228784) q[3];
sx q[3];
rz(1.2682605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1465313) q[0];
sx q[0];
rz(-2.596031) q[0];
sx q[0];
rz(-1.5978285) q[0];
rz(2.457288) q[1];
sx q[1];
rz(-1.4477718) q[1];
sx q[1];
rz(0.79944557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22173026) q[0];
sx q[0];
rz(-1.5654025) q[0];
sx q[0];
rz(3.0649867) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2861023) q[2];
sx q[2];
rz(-1.6281307) q[2];
sx q[2];
rz(2.3279026) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3189074) q[1];
sx q[1];
rz(-0.46055183) q[1];
sx q[1];
rz(3.1189633) q[1];
rz(-1.0897899) q[3];
sx q[3];
rz(-0.79339992) q[3];
sx q[3];
rz(0.97360669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6058309) q[2];
sx q[2];
rz(-2.3039218) q[2];
sx q[2];
rz(2.3568995) q[2];
rz(3.0253518) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(-1.2522662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5410974) q[0];
sx q[0];
rz(-1.2299812) q[0];
sx q[0];
rz(-1.0294718) q[0];
rz(3.0211499) q[1];
sx q[1];
rz(-1.8414626) q[1];
sx q[1];
rz(0.57074237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3523173) q[0];
sx q[0];
rz(-1.1162151) q[0];
sx q[0];
rz(-0.42735703) q[0];
x q[1];
rz(-3.013754) q[2];
sx q[2];
rz(-0.72261506) q[2];
sx q[2];
rz(1.3764718) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51008049) q[1];
sx q[1];
rz(-1.2464379) q[1];
sx q[1];
rz(-2.0381169) q[1];
x q[2];
rz(1.9071155) q[3];
sx q[3];
rz(-1.3988931) q[3];
sx q[3];
rz(0.055489648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.268078) q[2];
sx q[2];
rz(-0.87468481) q[2];
sx q[2];
rz(-0.52379215) q[2];
rz(1.9889471) q[3];
sx q[3];
rz(-2.5609784) q[3];
sx q[3];
rz(1.4279648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48083392) q[0];
sx q[0];
rz(-1.2013712) q[0];
sx q[0];
rz(-2.7177366) q[0];
rz(2.2747874) q[1];
sx q[1];
rz(-1.024217) q[1];
sx q[1];
rz(0.20326916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.568726) q[0];
sx q[0];
rz(-0.6558658) q[0];
sx q[0];
rz(0.076506581) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5166984) q[2];
sx q[2];
rz(-1.2612169) q[2];
sx q[2];
rz(2.5961909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8228559) q[1];
sx q[1];
rz(-1.4278432) q[1];
sx q[1];
rz(-0.28495423) q[1];
rz(-0.8747845) q[3];
sx q[3];
rz(-2.751707) q[3];
sx q[3];
rz(-2.3830151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.086143494) q[2];
sx q[2];
rz(-1.2148427) q[2];
sx q[2];
rz(-2.8622368) q[2];
rz(1.8481567) q[3];
sx q[3];
rz(-1.8297628) q[3];
sx q[3];
rz(1.9259341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16348895) q[0];
sx q[0];
rz(-0.80253974) q[0];
sx q[0];
rz(-1.4168903) q[0];
rz(-2.2143927) q[1];
sx q[1];
rz(-1.5798774) q[1];
sx q[1];
rz(1.7318116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15148189) q[0];
sx q[0];
rz(-1.6674433) q[0];
sx q[0];
rz(-2.2352909) q[0];
x q[1];
rz(2.0538834) q[2];
sx q[2];
rz(-1.8894686) q[2];
sx q[2];
rz(-2.0234194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5621582) q[1];
sx q[1];
rz(-0.82565281) q[1];
sx q[1];
rz(2.5316115) q[1];
rz(-1.7396443) q[3];
sx q[3];
rz(-2.8956684) q[3];
sx q[3];
rz(-2.3000474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37567821) q[2];
sx q[2];
rz(-1.8199074) q[2];
sx q[2];
rz(0.970617) q[2];
rz(2.4961903) q[3];
sx q[3];
rz(-2.161945) q[3];
sx q[3];
rz(-1.0055044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29393016) q[0];
sx q[0];
rz(-1.6456589) q[0];
sx q[0];
rz(1.8346067) q[0];
rz(-2.7597799) q[1];
sx q[1];
rz(-1.7996856) q[1];
sx q[1];
rz(-2.3736384) q[1];
rz(2.3293498) q[2];
sx q[2];
rz(-1.9612189) q[2];
sx q[2];
rz(-1.8283755) q[2];
rz(-2.0313203) q[3];
sx q[3];
rz(-1.9721748) q[3];
sx q[3];
rz(1.0655793) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
