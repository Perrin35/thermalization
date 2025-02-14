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
rz(-2.4401378) q[0];
sx q[0];
rz(-0.48818809) q[0];
sx q[0];
rz(0.39946431) q[0];
rz(-1.0723298) q[1];
sx q[1];
rz(-2.2430099) q[1];
sx q[1];
rz(2.8314765) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43543571) q[0];
sx q[0];
rz(-2.6118267) q[0];
sx q[0];
rz(-0.72491531) q[0];
x q[1];
rz(0.691857) q[2];
sx q[2];
rz(-0.82946316) q[2];
sx q[2];
rz(-2.8645309) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9609144) q[1];
sx q[1];
rz(-0.49094683) q[1];
sx q[1];
rz(0.21638201) q[1];
x q[2];
rz(-1.8350527) q[3];
sx q[3];
rz(-2.744361) q[3];
sx q[3];
rz(-2.5108932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8862137) q[2];
sx q[2];
rz(-1.9510521) q[2];
sx q[2];
rz(2.1204685) q[2];
rz(1.4897664) q[3];
sx q[3];
rz(-1.8332053) q[3];
sx q[3];
rz(-2.3511353) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0806231) q[0];
sx q[0];
rz(-1.0401833) q[0];
sx q[0];
rz(-2.8676497) q[0];
rz(-0.25320369) q[1];
sx q[1];
rz(-0.72414032) q[1];
sx q[1];
rz(-0.23987548) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.138809) q[0];
sx q[0];
rz(-1.0595508) q[0];
sx q[0];
rz(-0.89936043) q[0];
x q[1];
rz(-0.48440964) q[2];
sx q[2];
rz(-2.5229215) q[2];
sx q[2];
rz(1.2728438) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4644313) q[1];
sx q[1];
rz(-1.7993235) q[1];
sx q[1];
rz(-2.7941969) q[1];
rz(1.3680458) q[3];
sx q[3];
rz(-0.79907954) q[3];
sx q[3];
rz(-0.54982027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46951744) q[2];
sx q[2];
rz(-1.3152452) q[2];
sx q[2];
rz(-1.7496656) q[2];
rz(-2.4972656) q[3];
sx q[3];
rz(-1.9985806) q[3];
sx q[3];
rz(-2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7740087) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(2.9665663) q[0];
rz(1.4352098) q[1];
sx q[1];
rz(-1.8670466) q[1];
sx q[1];
rz(-0.79644901) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0284517) q[0];
sx q[0];
rz(-2.7243834) q[0];
sx q[0];
rz(0.88930486) q[0];
rz(-pi) q[1];
rz(2.5416614) q[2];
sx q[2];
rz(-1.33432) q[2];
sx q[2];
rz(1.8452871) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63523102) q[1];
sx q[1];
rz(-2.0331675) q[1];
sx q[1];
rz(1.1504786) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9655592) q[3];
sx q[3];
rz(-2.6280132) q[3];
sx q[3];
rz(1.5638417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1548057) q[2];
sx q[2];
rz(-0.044175819) q[2];
sx q[2];
rz(-3.0955443) q[2];
rz(2.0264528) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(-2.0980289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.1107156) q[0];
sx q[0];
rz(-0.59030384) q[0];
sx q[0];
rz(2.9764771) q[0];
rz(-1.1394399) q[1];
sx q[1];
rz(-0.93430263) q[1];
sx q[1];
rz(-1.6868235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44581977) q[0];
sx q[0];
rz(-0.26560387) q[0];
sx q[0];
rz(0.37047343) q[0];
rz(-pi) q[1];
rz(-1.7795935) q[2];
sx q[2];
rz(-1.1399802) q[2];
sx q[2];
rz(0.58961678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36697179) q[1];
sx q[1];
rz(-2.4200388) q[1];
sx q[1];
rz(-1.5040892) q[1];
rz(2.1490578) q[3];
sx q[3];
rz(-0.49697486) q[3];
sx q[3];
rz(-0.24028905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4399061) q[2];
sx q[2];
rz(-2.2540698) q[2];
sx q[2];
rz(-2.9627724) q[2];
rz(1.5170826) q[3];
sx q[3];
rz(-2.8330018) q[3];
sx q[3];
rz(1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4894067) q[0];
sx q[0];
rz(-1.711373) q[0];
sx q[0];
rz(1.7328316) q[0];
rz(-2.5694555) q[1];
sx q[1];
rz(-1.9913543) q[1];
sx q[1];
rz(-2.9170654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1869769) q[0];
sx q[0];
rz(-0.77565926) q[0];
sx q[0];
rz(-2.5903754) q[0];
rz(-2.3064086) q[2];
sx q[2];
rz(-1.3408644) q[2];
sx q[2];
rz(1.060677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46426168) q[1];
sx q[1];
rz(-0.867093) q[1];
sx q[1];
rz(-1.000885) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8711617) q[3];
sx q[3];
rz(-1.1704418) q[3];
sx q[3];
rz(2.3107236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3276334) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(2.853788) q[2];
rz(2.6156901) q[3];
sx q[3];
rz(-1.1651499) q[3];
sx q[3];
rz(0.56720081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7228058) q[0];
sx q[0];
rz(-2.4961508) q[0];
sx q[0];
rz(-3.1268815) q[0];
rz(0.22444935) q[1];
sx q[1];
rz(-1.4597471) q[1];
sx q[1];
rz(-0.86123484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447198) q[0];
sx q[0];
rz(-1.5278491) q[0];
sx q[0];
rz(1.5687902) q[0];
rz(-pi) q[1];
rz(3.1003816) q[2];
sx q[2];
rz(-1.905426) q[2];
sx q[2];
rz(-1.24805) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3726681) q[1];
sx q[1];
rz(-2.3834991) q[1];
sx q[1];
rz(-2.7102986) q[1];
x q[2];
rz(-1.3545996) q[3];
sx q[3];
rz(-2.4830468) q[3];
sx q[3];
rz(-1.6525846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33846778) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(-3.0306446) q[2];
rz(-2.0153996) q[3];
sx q[3];
rz(-1.6323099) q[3];
sx q[3];
rz(0.66669983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768537) q[0];
sx q[0];
rz(-0.56346098) q[0];
sx q[0];
rz(1.2891084) q[0];
rz(1.8076757) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(-1.8905554) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7374501) q[0];
sx q[0];
rz(-2.2206428) q[0];
sx q[0];
rz(-2.823816) q[0];
rz(-pi) q[1];
rz(0.47908135) q[2];
sx q[2];
rz(-0.63036171) q[2];
sx q[2];
rz(-2.1674974) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6896162) q[1];
sx q[1];
rz(-0.9202846) q[1];
sx q[1];
rz(-0.30779771) q[1];
rz(2.9032307) q[3];
sx q[3];
rz(-1.2796448) q[3];
sx q[3];
rz(-2.1057768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0879849) q[2];
sx q[2];
rz(-1.8156275) q[2];
sx q[2];
rz(-1.717022) q[2];
rz(0.21814957) q[3];
sx q[3];
rz(-1.679922) q[3];
sx q[3];
rz(0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3394534) q[0];
sx q[0];
rz(-2.5886783) q[0];
sx q[0];
rz(-0.46808991) q[0];
rz(2.3098436) q[1];
sx q[1];
rz(-2.2305198) q[1];
sx q[1];
rz(-0.98240486) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8881852) q[0];
sx q[0];
rz(-1.5063852) q[0];
sx q[0];
rz(-3.0712434) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8562819) q[2];
sx q[2];
rz(-0.64183706) q[2];
sx q[2];
rz(-1.8571719) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.959572) q[1];
sx q[1];
rz(-0.38117711) q[1];
sx q[1];
rz(-2.2439754) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.053388552) q[3];
sx q[3];
rz(-1.3636737) q[3];
sx q[3];
rz(-2.6762251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19782797) q[2];
sx q[2];
rz(-0.48913893) q[2];
sx q[2];
rz(-0.93734199) q[2];
rz(1.4165953) q[3];
sx q[3];
rz(-1.2920486) q[3];
sx q[3];
rz(0.14911252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8775403) q[0];
sx q[0];
rz(-2.6314681) q[0];
sx q[0];
rz(-2.383655) q[0];
rz(-2.0856048) q[1];
sx q[1];
rz(-2.2426558) q[1];
sx q[1];
rz(-0.24076375) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4419542) q[0];
sx q[0];
rz(-2.505971) q[0];
sx q[0];
rz(-2.5997396) q[0];
rz(-2.3182436) q[2];
sx q[2];
rz(-1.8728796) q[2];
sx q[2];
rz(-0.27034098) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4283838) q[1];
sx q[1];
rz(-1.2782103) q[1];
sx q[1];
rz(-0.87134113) q[1];
x q[2];
rz(0.82262294) q[3];
sx q[3];
rz(-0.47815505) q[3];
sx q[3];
rz(1.218623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15464887) q[2];
sx q[2];
rz(-1.1661531) q[2];
sx q[2];
rz(1.9722975) q[2];
rz(0.5430921) q[3];
sx q[3];
rz(-1.2767867) q[3];
sx q[3];
rz(-2.5663466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52694046) q[0];
sx q[0];
rz(-1.696538) q[0];
sx q[0];
rz(0.27969435) q[0];
rz(2.0277297) q[1];
sx q[1];
rz(-2.0852456) q[1];
sx q[1];
rz(0.14811555) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046127) q[0];
sx q[0];
rz(-1.7610504) q[0];
sx q[0];
rz(-2.4664761) q[0];
rz(-1.4101548) q[2];
sx q[2];
rz(-0.88692188) q[2];
sx q[2];
rz(1.6647673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1821226) q[1];
sx q[1];
rz(-2.39678) q[1];
sx q[1];
rz(-0.710979) q[1];
x q[2];
rz(2.4579416) q[3];
sx q[3];
rz(-1.7507075) q[3];
sx q[3];
rz(2.266117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0327882) q[2];
sx q[2];
rz(-2.2970436) q[2];
sx q[2];
rz(-2.8631701) q[2];
rz(-1.6935211) q[3];
sx q[3];
rz(-1.672594) q[3];
sx q[3];
rz(1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0447822) q[0];
sx q[0];
rz(-2.3533996) q[0];
sx q[0];
rz(2.0871373) q[0];
rz(1.0983559) q[1];
sx q[1];
rz(-1.337468) q[1];
sx q[1];
rz(-2.1877847) q[1];
rz(2.5658812) q[2];
sx q[2];
rz(-1.7127219) q[2];
sx q[2];
rz(-0.59319152) q[2];
rz(1.7058115) q[3];
sx q[3];
rz(-1.6208456) q[3];
sx q[3];
rz(2.7593233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
