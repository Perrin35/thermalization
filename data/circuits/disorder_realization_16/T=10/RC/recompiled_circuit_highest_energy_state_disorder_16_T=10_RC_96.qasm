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
rz(-1.1666522) q[0];
sx q[0];
rz(5.8402099) q[0];
sx q[0];
rz(9.0415434) q[0];
rz(-1.2869599) q[1];
sx q[1];
rz(1.9998963) q[1];
sx q[1];
rz(13.286446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.274668) q[0];
sx q[0];
rz(-2.5733181) q[0];
sx q[0];
rz(2.3734762) q[0];
rz(-2.949657) q[2];
sx q[2];
rz(-1.3752116) q[2];
sx q[2];
rz(1.5337613) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0648918) q[1];
sx q[1];
rz(-1.2982839) q[1];
sx q[1];
rz(0.38857675) q[1];
rz(-pi) q[2];
rz(-0.96020697) q[3];
sx q[3];
rz(-1.7060602) q[3];
sx q[3];
rz(1.8253411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2893452) q[2];
sx q[2];
rz(-0.87163681) q[2];
sx q[2];
rz(-0.80725011) q[2];
rz(-1.6926258) q[3];
sx q[3];
rz(-2.2117386) q[3];
sx q[3];
rz(0.71678954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221231) q[0];
sx q[0];
rz(-1.8417646) q[0];
sx q[0];
rz(-1.8362554) q[0];
rz(2.3985825) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(1.8441127) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3083619) q[0];
sx q[0];
rz(-1.1740645) q[0];
sx q[0];
rz(0.19773592) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2114685) q[2];
sx q[2];
rz(-0.67652297) q[2];
sx q[2];
rz(-0.97201024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49238047) q[1];
sx q[1];
rz(-0.67427507) q[1];
sx q[1];
rz(-1.0843305) q[1];
x q[2];
rz(0.30001202) q[3];
sx q[3];
rz(-0.37676806) q[3];
sx q[3];
rz(-3.026791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7617191) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(-1.0404111) q[2];
rz(-2.8271293) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(-1.8837121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8085025) q[0];
sx q[0];
rz(-0.98121488) q[0];
sx q[0];
rz(-1.3370978) q[0];
rz(0.62726504) q[1];
sx q[1];
rz(-1.7341055) q[1];
sx q[1];
rz(-1.6076535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3298137) q[0];
sx q[0];
rz(-0.80429351) q[0];
sx q[0];
rz(-0.95383184) q[0];
x q[1];
rz(1.2731601) q[2];
sx q[2];
rz(-1.9797851) q[2];
sx q[2];
rz(-1.9252418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4348373) q[1];
sx q[1];
rz(-1.1011775) q[1];
sx q[1];
rz(-1.3634741) q[1];
x q[2];
rz(-0.19948761) q[3];
sx q[3];
rz(-1.9428245) q[3];
sx q[3];
rz(0.92687273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73438054) q[2];
sx q[2];
rz(-0.89443365) q[2];
sx q[2];
rz(-1.9341932) q[2];
rz(-2.1842128) q[3];
sx q[3];
rz(-1.2355685) q[3];
sx q[3];
rz(-2.4522219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6974238) q[0];
sx q[0];
rz(-2.1902695) q[0];
sx q[0];
rz(-3.0498258) q[0];
rz(0.16608876) q[1];
sx q[1];
rz(-1.2867915) q[1];
sx q[1];
rz(-1.0413569) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0386376) q[0];
sx q[0];
rz(-2.345076) q[0];
sx q[0];
rz(0.060925555) q[0];
x q[1];
rz(-2.741124) q[2];
sx q[2];
rz(-2.31268) q[2];
sx q[2];
rz(-2.9603855) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2745121) q[1];
sx q[1];
rz(-1.5364944) q[1];
sx q[1];
rz(-1.097789) q[1];
x q[2];
rz(1.6575069) q[3];
sx q[3];
rz(-1.711688) q[3];
sx q[3];
rz(2.2839462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2287717) q[2];
sx q[2];
rz(-0.81835881) q[2];
sx q[2];
rz(1.5405601) q[2];
rz(-1.0176963) q[3];
sx q[3];
rz(-1.8120268) q[3];
sx q[3];
rz(1.2525919) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8215264) q[0];
sx q[0];
rz(-1.7463266) q[0];
sx q[0];
rz(-0.23859247) q[0];
rz(-2.353031) q[1];
sx q[1];
rz(-0.32222727) q[1];
sx q[1];
rz(2.2408748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470578) q[0];
sx q[0];
rz(-1.3119895) q[0];
sx q[0];
rz(2.6954891) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3447755) q[2];
sx q[2];
rz(-0.39835762) q[2];
sx q[2];
rz(-2.4657754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5594192) q[1];
sx q[1];
rz(-1.4430077) q[1];
sx q[1];
rz(-2.726597) q[1];
rz(-pi) q[2];
rz(-2.4191148) q[3];
sx q[3];
rz(-2.1584508) q[3];
sx q[3];
rz(-2.8624299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0476524) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(1.368604) q[2];
rz(2.7641344) q[3];
sx q[3];
rz(-2.3182175) q[3];
sx q[3];
rz(1.1641758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0193598) q[0];
sx q[0];
rz(-0.43411175) q[0];
sx q[0];
rz(-0.22943108) q[0];
rz(2.5550628) q[1];
sx q[1];
rz(-2.6507381) q[1];
sx q[1];
rz(-1.4803001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9113113) q[0];
sx q[0];
rz(-2.3846941) q[0];
sx q[0];
rz(1.043399) q[0];
rz(1.3103799) q[2];
sx q[2];
rz(-1.0170817) q[2];
sx q[2];
rz(1.9268203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28162128) q[1];
sx q[1];
rz(-2.9088417) q[1];
sx q[1];
rz(0.23334664) q[1];
x q[2];
rz(-0.89648439) q[3];
sx q[3];
rz(-1.4004882) q[3];
sx q[3];
rz(-1.299082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4218939) q[2];
sx q[2];
rz(-2.5514558) q[2];
sx q[2];
rz(-1.3479007) q[2];
rz(0.19139309) q[3];
sx q[3];
rz(-2.8189711) q[3];
sx q[3];
rz(1.9026683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557068) q[0];
sx q[0];
rz(-1.6275591) q[0];
sx q[0];
rz(0.53644449) q[0];
rz(-3.0570807) q[1];
sx q[1];
rz(-2.5630496) q[1];
sx q[1];
rz(1.1892148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0482423) q[0];
sx q[0];
rz(-2.1879401) q[0];
sx q[0];
rz(0.18927745) q[0];
x q[1];
rz(1.9574478) q[2];
sx q[2];
rz(-2.2154567) q[2];
sx q[2];
rz(1.2520916) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7658733) q[1];
sx q[1];
rz(-0.60085591) q[1];
sx q[1];
rz(-2.145776) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0310845) q[3];
sx q[3];
rz(-1.7841859) q[3];
sx q[3];
rz(-1.8387356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1875923) q[2];
sx q[2];
rz(-1.2728929) q[2];
sx q[2];
rz(3.0583337) q[2];
rz(-1.8795053) q[3];
sx q[3];
rz(-0.83674651) q[3];
sx q[3];
rz(-2.4038103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5543723) q[0];
sx q[0];
rz(-1.502259) q[0];
sx q[0];
rz(2.3578405) q[0];
rz(1.4844249) q[1];
sx q[1];
rz(-0.52394167) q[1];
sx q[1];
rz(0.2934244) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3104738) q[0];
sx q[0];
rz(-1.2037414) q[0];
sx q[0];
rz(1.0793174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9181554) q[2];
sx q[2];
rz(-2.2268837) q[2];
sx q[2];
rz(-2.3272397) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1841976) q[1];
sx q[1];
rz(-2.7009472) q[1];
sx q[1];
rz(1.22681) q[1];
x q[2];
rz(2.226172) q[3];
sx q[3];
rz(-1.7599981) q[3];
sx q[3];
rz(2.9331741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7744814) q[2];
sx q[2];
rz(-1.7256871) q[2];
sx q[2];
rz(-1.492738) q[2];
rz(1.9914918) q[3];
sx q[3];
rz(-2.8425358) q[3];
sx q[3];
rz(-0.90112954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25006008) q[0];
sx q[0];
rz(-2.1136668) q[0];
sx q[0];
rz(-0.0082536396) q[0];
rz(0.0025657733) q[1];
sx q[1];
rz(-0.51104128) q[1];
sx q[1];
rz(-1.0645688) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9998523) q[0];
sx q[0];
rz(-2.2069262) q[0];
sx q[0];
rz(2.0881235) q[0];
rz(3.1334115) q[2];
sx q[2];
rz(-1.3213833) q[2];
sx q[2];
rz(-1.8298263) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9656532) q[1];
sx q[1];
rz(-2.4373131) q[1];
sx q[1];
rz(0.85051103) q[1];
rz(-pi) q[2];
rz(-1.4225619) q[3];
sx q[3];
rz(-0.98662107) q[3];
sx q[3];
rz(2.119182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.061444) q[2];
sx q[2];
rz(-1.4465569) q[2];
sx q[2];
rz(0.20305571) q[2];
rz(-1.2114581) q[3];
sx q[3];
rz(-2.0720427) q[3];
sx q[3];
rz(-3.0058461) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900629) q[0];
sx q[0];
rz(-1.1030581) q[0];
sx q[0];
rz(-0.67361012) q[0];
rz(2.8296962) q[1];
sx q[1];
rz(-1.2047647) q[1];
sx q[1];
rz(-1.9685251) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6098227) q[0];
sx q[0];
rz(-1.2780452) q[0];
sx q[0];
rz(-1.0171153) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3373702) q[2];
sx q[2];
rz(-0.82841043) q[2];
sx q[2];
rz(0.5898005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7562779) q[1];
sx q[1];
rz(-0.42149252) q[1];
sx q[1];
rz(1.384686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7133663) q[3];
sx q[3];
rz(-0.90598327) q[3];
sx q[3];
rz(-1.6237824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6761026) q[2];
sx q[2];
rz(-2.5238621) q[2];
sx q[2];
rz(-0.22107302) q[2];
rz(-2.5681833) q[3];
sx q[3];
rz(-2.2677877) q[3];
sx q[3];
rz(2.0360086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40972805) q[0];
sx q[0];
rz(-1.2936214) q[0];
sx q[0];
rz(-3.1142942) q[0];
rz(-1.9505386) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(-1.9223735) q[2];
sx q[2];
rz(-0.58813358) q[2];
sx q[2];
rz(1.0434601) q[2];
rz(2.7577362) q[3];
sx q[3];
rz(-1.0105269) q[3];
sx q[3];
rz(-0.15906048) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
