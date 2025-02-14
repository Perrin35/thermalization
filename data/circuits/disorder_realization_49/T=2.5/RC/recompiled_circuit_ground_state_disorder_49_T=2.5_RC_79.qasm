OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1385652) q[0];
sx q[0];
rz(-0.87831098) q[0];
sx q[0];
rz(-0.83100975) q[0];
rz(2.5660958) q[1];
sx q[1];
rz(-0.73605186) q[1];
sx q[1];
rz(2.4017258) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66798009) q[0];
sx q[0];
rz(-1.6428609) q[0];
sx q[0];
rz(-2.9524809) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0027356) q[2];
sx q[2];
rz(-1.4502) q[2];
sx q[2];
rz(-1.9386148) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4783096) q[1];
sx q[1];
rz(-1.7170719) q[1];
sx q[1];
rz(2.0520794) q[1];
x q[2];
rz(2.7571477) q[3];
sx q[3];
rz(-1.3923402) q[3];
sx q[3];
rz(0.78027356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2089219) q[2];
sx q[2];
rz(-2.966556) q[2];
sx q[2];
rz(2.8026061) q[2];
rz(-0.50421667) q[3];
sx q[3];
rz(-2.1239069) q[3];
sx q[3];
rz(-2.0174111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9826688) q[0];
sx q[0];
rz(-2.447154) q[0];
sx q[0];
rz(0.50952953) q[0];
rz(-1.5024028) q[1];
sx q[1];
rz(-2.8645611) q[1];
sx q[1];
rz(0.94430077) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4128542) q[0];
sx q[0];
rz(-1.3711689) q[0];
sx q[0];
rz(2.9908604) q[0];
x q[1];
rz(-0.048205094) q[2];
sx q[2];
rz(-1.4314326) q[2];
sx q[2];
rz(2.4037619) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2118476) q[1];
sx q[1];
rz(-2.3073688) q[1];
sx q[1];
rz(-2.3843308) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89280309) q[3];
sx q[3];
rz(-1.2965186) q[3];
sx q[3];
rz(2.1712042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3462191) q[2];
sx q[2];
rz(-1.5495164) q[2];
sx q[2];
rz(-2.6056371) q[2];
rz(2.0761944) q[3];
sx q[3];
rz(-0.25748101) q[3];
sx q[3];
rz(-2.4901701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1123493) q[0];
sx q[0];
rz(-1.9922682) q[0];
sx q[0];
rz(-2.405622) q[0];
rz(0.53572267) q[1];
sx q[1];
rz(-1.2993206) q[1];
sx q[1];
rz(-0.89964286) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7539983) q[0];
sx q[0];
rz(-1.7090461) q[0];
sx q[0];
rz(1.686245) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2399142) q[2];
sx q[2];
rz(-2.5106259) q[2];
sx q[2];
rz(-0.081693782) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9316072) q[1];
sx q[1];
rz(-1.6495541) q[1];
sx q[1];
rz(1.0754536) q[1];
x q[2];
rz(-0.24921649) q[3];
sx q[3];
rz(-1.936541) q[3];
sx q[3];
rz(2.5178227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6003517) q[2];
sx q[2];
rz(-1.1949801) q[2];
sx q[2];
rz(2.3386653) q[2];
rz(-0.7488572) q[3];
sx q[3];
rz(-1.3978981) q[3];
sx q[3];
rz(2.9822541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51711851) q[0];
sx q[0];
rz(-1.7298537) q[0];
sx q[0];
rz(0.13036048) q[0];
rz(-1.2654001) q[1];
sx q[1];
rz(-1.1771026) q[1];
sx q[1];
rz(0.1098384) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43183655) q[0];
sx q[0];
rz(-2.0227814) q[0];
sx q[0];
rz(-2.7664037) q[0];
x q[1];
rz(-1.2939343) q[2];
sx q[2];
rz(-2.1199391) q[2];
sx q[2];
rz(2.800966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3567645) q[1];
sx q[1];
rz(-2.3127529) q[1];
sx q[1];
rz(-2.1020562) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5073152) q[3];
sx q[3];
rz(-1.1697672) q[3];
sx q[3];
rz(2.4672535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.639223) q[2];
sx q[2];
rz(-1.9618192) q[2];
sx q[2];
rz(-2.7184674) q[2];
rz(2.1028178) q[3];
sx q[3];
rz(-2.7602502) q[3];
sx q[3];
rz(1.0264621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.0005242) q[0];
sx q[0];
rz(-1.8924014) q[0];
sx q[0];
rz(2.8619859) q[0];
rz(0.33310834) q[1];
sx q[1];
rz(-2.6393642) q[1];
sx q[1];
rz(2.8531029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97388291) q[0];
sx q[0];
rz(-1.6339193) q[0];
sx q[0];
rz(1.0991251) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1568858) q[2];
sx q[2];
rz(-2.8749646) q[2];
sx q[2];
rz(-1.3150584) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7964871) q[1];
sx q[1];
rz(-3.0287841) q[1];
sx q[1];
rz(-2.2276001) q[1];
x q[2];
rz(-1.9506504) q[3];
sx q[3];
rz(-0.65841802) q[3];
sx q[3];
rz(1.0059716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5491817) q[2];
sx q[2];
rz(-1.9224527) q[2];
sx q[2];
rz(2.9782817) q[2];
rz(-1.1641938) q[3];
sx q[3];
rz(-2.4653698) q[3];
sx q[3];
rz(-1.7827079) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0613681) q[0];
sx q[0];
rz(-3.1243262) q[0];
sx q[0];
rz(1.9153216) q[0];
rz(-1.3310883) q[1];
sx q[1];
rz(-2.2115579) q[1];
sx q[1];
rz(-0.61940449) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3400522) q[0];
sx q[0];
rz(-2.0503133) q[0];
sx q[0];
rz(2.1229486) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.370843) q[2];
sx q[2];
rz(-1.3226349) q[2];
sx q[2];
rz(-1.2436109) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88404578) q[1];
sx q[1];
rz(-2.9830898) q[1];
sx q[1];
rz(-1.201215) q[1];
x q[2];
rz(-1.9483637) q[3];
sx q[3];
rz(-1.7488297) q[3];
sx q[3];
rz(-1.5291027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2254534) q[2];
sx q[2];
rz(-1.4832387) q[2];
sx q[2];
rz(0.43530604) q[2];
rz(-0.12886038) q[3];
sx q[3];
rz(-1.83056) q[3];
sx q[3];
rz(0.35663566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72047609) q[0];
sx q[0];
rz(-0.55835503) q[0];
sx q[0];
rz(-0.55225736) q[0];
rz(-2.5241959) q[1];
sx q[1];
rz(-1.9300902) q[1];
sx q[1];
rz(-1.6389821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.033806697) q[0];
sx q[0];
rz(-1.5272041) q[0];
sx q[0];
rz(-1.1826452) q[0];
x q[1];
rz(-1.8322629) q[2];
sx q[2];
rz(-0.51957031) q[2];
sx q[2];
rz(2.7594523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.43728033) q[1];
sx q[1];
rz(-1.4131792) q[1];
sx q[1];
rz(-1.0606517) q[1];
rz(0.076528744) q[3];
sx q[3];
rz(-2.8815443) q[3];
sx q[3];
rz(-1.6116774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5611394) q[2];
sx q[2];
rz(-0.20623198) q[2];
sx q[2];
rz(2.0533766) q[2];
rz(-0.020180833) q[3];
sx q[3];
rz(-1.2729278) q[3];
sx q[3];
rz(1.5816241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7839171) q[0];
sx q[0];
rz(-2.7516784) q[0];
sx q[0];
rz(-2.1141323) q[0];
rz(-2.0383535) q[1];
sx q[1];
rz(-1.5733066) q[1];
sx q[1];
rz(-1.2148414) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2526557) q[0];
sx q[0];
rz(-1.2644469) q[0];
sx q[0];
rz(-0.10821786) q[0];
rz(-0.80963366) q[2];
sx q[2];
rz(-0.61695951) q[2];
sx q[2];
rz(3.0669341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.852927) q[1];
sx q[1];
rz(-0.92354362) q[1];
sx q[1];
rz(-2.701328) q[1];
rz(-pi) q[2];
rz(-1.070511) q[3];
sx q[3];
rz(-2.1246582) q[3];
sx q[3];
rz(1.1244233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2601629) q[2];
sx q[2];
rz(-1.9378928) q[2];
sx q[2];
rz(-1.0677968) q[2];
rz(-2.2504375) q[3];
sx q[3];
rz(-2.6813337) q[3];
sx q[3];
rz(-2.0021745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.7282309) q[0];
sx q[0];
rz(-2.3869393) q[0];
sx q[0];
rz(-0.2555787) q[0];
rz(-0.74343395) q[1];
sx q[1];
rz(-1.5836704) q[1];
sx q[1];
rz(1.5240634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179068) q[0];
sx q[0];
rz(-2.0605378) q[0];
sx q[0];
rz(2.6262002) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7268098) q[2];
sx q[2];
rz(-0.66294248) q[2];
sx q[2];
rz(-2.489733) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2801622) q[1];
sx q[1];
rz(-1.4373684) q[1];
sx q[1];
rz(0.89806865) q[1];
rz(2.8761707) q[3];
sx q[3];
rz(-2.1947104) q[3];
sx q[3];
rz(1.3359631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1114379) q[2];
sx q[2];
rz(-1.43575) q[2];
sx q[2];
rz(-1.6072404) q[2];
rz(-2.0700908) q[3];
sx q[3];
rz(-1.6000308) q[3];
sx q[3];
rz(-1.967954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2509505) q[0];
sx q[0];
rz(-0.28674704) q[0];
sx q[0];
rz(-0.51837921) q[0];
rz(0.80825835) q[1];
sx q[1];
rz(-1.4603442) q[1];
sx q[1];
rz(1.6917276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5980412) q[0];
sx q[0];
rz(-1.8596453) q[0];
sx q[0];
rz(3.1310215) q[0];
rz(-0.53624714) q[2];
sx q[2];
rz(-1.9384346) q[2];
sx q[2];
rz(-0.73208955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.58123484) q[1];
sx q[1];
rz(-2.5311845) q[1];
sx q[1];
rz(-2.6950652) q[1];
rz(0.57761044) q[3];
sx q[3];
rz(-2.0311714) q[3];
sx q[3];
rz(-2.1652997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0290252) q[2];
sx q[2];
rz(-1.2351278) q[2];
sx q[2];
rz(0.9355363) q[2];
rz(-1.4714636) q[3];
sx q[3];
rz(-1.4751438) q[3];
sx q[3];
rz(-1.0940301) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6846631) q[0];
sx q[0];
rz(-2.5419432) q[0];
sx q[0];
rz(2.3210617) q[0];
rz(-0.44878557) q[1];
sx q[1];
rz(-1.1460591) q[1];
sx q[1];
rz(-1.9312327) q[1];
rz(-1.2032897) q[2];
sx q[2];
rz(-1.6151645) q[2];
sx q[2];
rz(0.9565959) q[2];
rz(-1.9121691) q[3];
sx q[3];
rz(-2.1044272) q[3];
sx q[3];
rz(2.9499346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
