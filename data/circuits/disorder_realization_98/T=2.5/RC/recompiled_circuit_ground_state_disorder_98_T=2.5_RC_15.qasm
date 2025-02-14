OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(-0.21259354) q[0];
sx q[0];
rz(0.73417443) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(-0.34595481) q[1];
sx q[1];
rz(3.1168361) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7914331) q[0];
sx q[0];
rz(-1.6140249) q[0];
sx q[0];
rz(0.062383609) q[0];
x q[1];
rz(-2.0169746) q[2];
sx q[2];
rz(-0.25814787) q[2];
sx q[2];
rz(3.019697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7721482) q[1];
sx q[1];
rz(-1.8696864) q[1];
sx q[1];
rz(-1.5159831) q[1];
rz(-pi) q[2];
rz(1.3583899) q[3];
sx q[3];
rz(-1.5161361) q[3];
sx q[3];
rz(0.89059356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2572702) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(-2.8327668) q[2];
rz(-2.3158) q[3];
sx q[3];
rz(-1.0079931) q[3];
sx q[3];
rz(0.76396137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17077133) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(1.0276851) q[0];
rz(-2.3254501) q[1];
sx q[1];
rz(-2.0232537) q[1];
sx q[1];
rz(0.81248409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4554169) q[0];
sx q[0];
rz(-0.38639613) q[0];
sx q[0];
rz(-1.1375582) q[0];
rz(-pi) q[1];
rz(1.5565926) q[2];
sx q[2];
rz(-1.1974632) q[2];
sx q[2];
rz(2.4844784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1907249) q[1];
sx q[1];
rz(-2.0800033) q[1];
sx q[1];
rz(2.9835761) q[1];
rz(-pi) q[2];
rz(-0.79576335) q[3];
sx q[3];
rz(-0.47945346) q[3];
sx q[3];
rz(-2.457452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5673148) q[2];
sx q[2];
rz(-1.652176) q[2];
sx q[2];
rz(2.2226649) q[2];
rz(-2.610176) q[3];
sx q[3];
rz(-2.0960977) q[3];
sx q[3];
rz(-0.04118583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.4384005) q[0];
sx q[0];
rz(-2.0719318) q[0];
sx q[0];
rz(-2.0042787) q[0];
rz(-1.3905585) q[1];
sx q[1];
rz(-0.5568234) q[1];
sx q[1];
rz(0.73296076) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61444401) q[0];
sx q[0];
rz(-1.9839459) q[0];
sx q[0];
rz(-1.2507417) q[0];
x q[1];
rz(2.5733092) q[2];
sx q[2];
rz(-2.097762) q[2];
sx q[2];
rz(-2.2672841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3324686) q[1];
sx q[1];
rz(-2.4753503) q[1];
sx q[1];
rz(-0.82102832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5571345) q[3];
sx q[3];
rz(-0.44876305) q[3];
sx q[3];
rz(-0.45873935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0040032337) q[2];
sx q[2];
rz(-1.4121476) q[2];
sx q[2];
rz(-1.0388733) q[2];
rz(-2.1393356) q[3];
sx q[3];
rz(-1.301731) q[3];
sx q[3];
rz(2.8514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9349979) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(2.2344053) q[0];
rz(-2.9838003) q[1];
sx q[1];
rz(-2.1267499) q[1];
sx q[1];
rz(-1.8603604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84128252) q[0];
sx q[0];
rz(-2.1774946) q[0];
sx q[0];
rz(-0.9228031) q[0];
x q[1];
rz(2.0220535) q[2];
sx q[2];
rz(-0.76821487) q[2];
sx q[2];
rz(0.31203416) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7042027) q[1];
sx q[1];
rz(-1.9949732) q[1];
sx q[1];
rz(0.68628879) q[1];
x q[2];
rz(-0.82366039) q[3];
sx q[3];
rz(-2.3170217) q[3];
sx q[3];
rz(2.0093105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24718757) q[2];
sx q[2];
rz(-0.93623585) q[2];
sx q[2];
rz(1.8518764) q[2];
rz(1.7973409) q[3];
sx q[3];
rz(-1.4569837) q[3];
sx q[3];
rz(0.70014203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808788) q[0];
sx q[0];
rz(-2.6432156) q[0];
sx q[0];
rz(-2.1233249) q[0];
rz(-2.0385888) q[1];
sx q[1];
rz(-0.7119199) q[1];
sx q[1];
rz(-1.8011372) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0606275) q[0];
sx q[0];
rz(-0.61765352) q[0];
sx q[0];
rz(-2.4897442) q[0];
rz(-pi) q[1];
rz(-1.6685772) q[2];
sx q[2];
rz(-1.1927529) q[2];
sx q[2];
rz(-0.60286544) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.75234883) q[1];
sx q[1];
rz(-1.6977786) q[1];
sx q[1];
rz(-1.701322) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0304673) q[3];
sx q[3];
rz(-1.6202462) q[3];
sx q[3];
rz(1.7211308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5420142) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(-1.0020024) q[2];
rz(-0.69139785) q[3];
sx q[3];
rz(-1.9611497) q[3];
sx q[3];
rz(1.6208167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-1.6028676) q[0];
sx q[0];
rz(-2.0724917) q[0];
sx q[0];
rz(-0.62527239) q[0];
rz(1.9484776) q[1];
sx q[1];
rz(-0.68980491) q[1];
sx q[1];
rz(-2.5742721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7517487) q[0];
sx q[0];
rz(-1.3951256) q[0];
sx q[0];
rz(2.3542464) q[0];
rz(-0.82639931) q[2];
sx q[2];
rz(-1.5652839) q[2];
sx q[2];
rz(2.7720087) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4299791) q[1];
sx q[1];
rz(-0.95250722) q[1];
sx q[1];
rz(0.031876335) q[1];
rz(-pi) q[2];
rz(2.3472957) q[3];
sx q[3];
rz(-0.24883379) q[3];
sx q[3];
rz(0.90422309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83711964) q[2];
sx q[2];
rz(-1.8581055) q[2];
sx q[2];
rz(-1.5817969) q[2];
rz(-0.61257735) q[3];
sx q[3];
rz(-1.9809096) q[3];
sx q[3];
rz(2.9858203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0953858) q[0];
sx q[0];
rz(-1.950773) q[0];
sx q[0];
rz(1.1268536) q[0];
rz(1.8719748) q[1];
sx q[1];
rz(-1.1879299) q[1];
sx q[1];
rz(-2.6205305) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714037) q[0];
sx q[0];
rz(-0.91556433) q[0];
sx q[0];
rz(-2.6960424) q[0];
x q[1];
rz(-1.2418141) q[2];
sx q[2];
rz(-1.6018724) q[2];
sx q[2];
rz(-0.46361332) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4849629) q[1];
sx q[1];
rz(-0.99492517) q[1];
sx q[1];
rz(0.63159512) q[1];
rz(-pi) q[2];
rz(-2.7016958) q[3];
sx q[3];
rz(-0.58618136) q[3];
sx q[3];
rz(-2.2324454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1178939) q[2];
sx q[2];
rz(-1.3641027) q[2];
sx q[2];
rz(-1.914631) q[2];
rz(-1.4620694) q[3];
sx q[3];
rz(-1.5173802) q[3];
sx q[3];
rz(0.44155651) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11314497) q[0];
sx q[0];
rz(-2.1118836) q[0];
sx q[0];
rz(2.5944769) q[0];
rz(-0.088134915) q[1];
sx q[1];
rz(-1.3476177) q[1];
sx q[1];
rz(-0.42253447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28332253) q[0];
sx q[0];
rz(-1.2898603) q[0];
sx q[0];
rz(-3.1188117) q[0];
rz(-0.59178517) q[2];
sx q[2];
rz(-0.6066423) q[2];
sx q[2];
rz(-3.0559412) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7393028) q[1];
sx q[1];
rz(-1.6268432) q[1];
sx q[1];
rz(1.6399553) q[1];
x q[2];
rz(1.3874153) q[3];
sx q[3];
rz(-1.7088582) q[3];
sx q[3];
rz(-2.3473397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2253458) q[2];
sx q[2];
rz(-1.6371181) q[2];
sx q[2];
rz(-0.28727356) q[2];
rz(-0.54287994) q[3];
sx q[3];
rz(-0.77016872) q[3];
sx q[3];
rz(-0.058102593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420237) q[0];
sx q[0];
rz(-0.68411198) q[0];
sx q[0];
rz(2.3642819) q[0];
rz(-0.9264535) q[1];
sx q[1];
rz(-1.9994206) q[1];
sx q[1];
rz(-1.3949589) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.004121) q[0];
sx q[0];
rz(-0.51252675) q[0];
sx q[0];
rz(0.88051535) q[0];
x q[1];
rz(0.35819004) q[2];
sx q[2];
rz(-0.6300504) q[2];
sx q[2];
rz(-1.9262528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7662188) q[1];
sx q[1];
rz(-1.2412984) q[1];
sx q[1];
rz(-1.0887515) q[1];
rz(-pi) q[2];
rz(0.80628245) q[3];
sx q[3];
rz(-0.7976992) q[3];
sx q[3];
rz(-2.9258941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7159783) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(-1.0905637) q[2];
rz(-0.96450949) q[3];
sx q[3];
rz(-1.0114074) q[3];
sx q[3];
rz(-1.2151659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57508093) q[0];
sx q[0];
rz(-2.8031741) q[0];
sx q[0];
rz(-1.5100719) q[0];
rz(-3.1072726) q[1];
sx q[1];
rz(-1.3395373) q[1];
sx q[1];
rz(-0.43201772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2902108) q[0];
sx q[0];
rz(-0.54199666) q[0];
sx q[0];
rz(1.4567514) q[0];
rz(-pi) q[1];
rz(1.2197184) q[2];
sx q[2];
rz(-1.6475999) q[2];
sx q[2];
rz(-1.6258282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4544983) q[1];
sx q[1];
rz(-1.7155572) q[1];
sx q[1];
rz(1.6705546) q[1];
rz(1.179078) q[3];
sx q[3];
rz(-2.5234902) q[3];
sx q[3];
rz(0.51392344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4711275) q[2];
sx q[2];
rz(-1.1521143) q[2];
sx q[2];
rz(-1.0738037) q[2];
rz(0.18887575) q[3];
sx q[3];
rz(-2.5329068) q[3];
sx q[3];
rz(-2.7591738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0976681) q[0];
sx q[0];
rz(-0.32079874) q[0];
sx q[0];
rz(-0.70377845) q[0];
rz(1.7882998) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(-2.399171) q[2];
sx q[2];
rz(-2.0423642) q[2];
sx q[2];
rz(-2.5004417) q[2];
rz(0.13808098) q[3];
sx q[3];
rz(-0.46317536) q[3];
sx q[3];
rz(-2.4899766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
