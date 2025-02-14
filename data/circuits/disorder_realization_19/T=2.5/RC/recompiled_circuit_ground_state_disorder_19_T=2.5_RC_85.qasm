OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.98915339) q[0];
sx q[0];
rz(-1.5664772) q[0];
sx q[0];
rz(-1.1251261) q[0];
rz(-2.1195124) q[1];
sx q[1];
rz(-2.4740969) q[1];
sx q[1];
rz(-0.8134841) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73972244) q[0];
sx q[0];
rz(-0.39416322) q[0];
sx q[0];
rz(-1.1959082) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9314693) q[2];
sx q[2];
rz(-1.1481592) q[2];
sx q[2];
rz(0.3061184) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58373815) q[1];
sx q[1];
rz(-1.1314794) q[1];
sx q[1];
rz(-1.735461) q[1];
rz(-pi) q[2];
rz(2.0735246) q[3];
sx q[3];
rz(-2.9554071) q[3];
sx q[3];
rz(-1.26621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4890613) q[2];
sx q[2];
rz(-1.4514613) q[2];
sx q[2];
rz(0.92864621) q[2];
rz(1.5422025) q[3];
sx q[3];
rz(-1.8079115) q[3];
sx q[3];
rz(-1.9699875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9376675) q[0];
sx q[0];
rz(-1.7610022) q[0];
sx q[0];
rz(-0.12705886) q[0];
rz(0.98310414) q[1];
sx q[1];
rz(-1.3652722) q[1];
sx q[1];
rz(-2.3703221) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8255804) q[0];
sx q[0];
rz(-0.86992747) q[0];
sx q[0];
rz(0.46937816) q[0];
rz(-pi) q[1];
x q[1];
rz(0.060734435) q[2];
sx q[2];
rz(-1.2337451) q[2];
sx q[2];
rz(3.0063546) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3202618) q[1];
sx q[1];
rz(-1.2148569) q[1];
sx q[1];
rz(1.9411646) q[1];
rz(-pi) q[2];
rz(-1.0326321) q[3];
sx q[3];
rz(-2.0797044) q[3];
sx q[3];
rz(0.63652388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9942921) q[2];
sx q[2];
rz(-1.2039801) q[2];
sx q[2];
rz(0.0083943923) q[2];
rz(-0.66347915) q[3];
sx q[3];
rz(-1.9050262) q[3];
sx q[3];
rz(2.8765163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10107772) q[0];
sx q[0];
rz(-2.3150257) q[0];
sx q[0];
rz(-0.44152942) q[0];
rz(0.98835522) q[1];
sx q[1];
rz(-1.1343196) q[1];
sx q[1];
rz(0.13557869) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5268742) q[0];
sx q[0];
rz(-1.5784987) q[0];
sx q[0];
rz(-3.1253424) q[0];
x q[1];
rz(-0.71696059) q[2];
sx q[2];
rz(-1.2768942) q[2];
sx q[2];
rz(-0.64168054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5244658) q[1];
sx q[1];
rz(-1.9248665) q[1];
sx q[1];
rz(-2.0984142) q[1];
rz(-2.9445678) q[3];
sx q[3];
rz(-0.84614119) q[3];
sx q[3];
rz(-1.0850832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93379891) q[2];
sx q[2];
rz(-1.7094882) q[2];
sx q[2];
rz(-3.1033893) q[2];
rz(0.52538747) q[3];
sx q[3];
rz(-2.5140258) q[3];
sx q[3];
rz(2.4562522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4811089) q[0];
sx q[0];
rz(-2.3479192) q[0];
sx q[0];
rz(-0.61087459) q[0];
rz(1.3350217) q[1];
sx q[1];
rz(-1.3742615) q[1];
sx q[1];
rz(0.23922051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67475457) q[0];
sx q[0];
rz(-1.309899) q[0];
sx q[0];
rz(-0.99715085) q[0];
x q[1];
rz(0.68065721) q[2];
sx q[2];
rz(-0.58154642) q[2];
sx q[2];
rz(-2.8342063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4799616) q[1];
sx q[1];
rz(-2.0212272) q[1];
sx q[1];
rz(-1.2110787) q[1];
rz(-pi) q[2];
rz(-2.753965) q[3];
sx q[3];
rz(-2.648733) q[3];
sx q[3];
rz(-0.29407497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7188344) q[2];
sx q[2];
rz(-1.6966635) q[2];
sx q[2];
rz(-0.0017496721) q[2];
rz(-2.8478029) q[3];
sx q[3];
rz(-1.9521451) q[3];
sx q[3];
rz(-0.26688117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9002429) q[0];
sx q[0];
rz(-1.3768063) q[0];
sx q[0];
rz(-0.84306651) q[0];
rz(-2.8040366) q[1];
sx q[1];
rz(-1.2755716) q[1];
sx q[1];
rz(1.6036124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8385411) q[0];
sx q[0];
rz(-1.3796796) q[0];
sx q[0];
rz(0.79940807) q[0];
x q[1];
rz(2.9672253) q[2];
sx q[2];
rz(-1.5182759) q[2];
sx q[2];
rz(-0.32584056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5953491) q[1];
sx q[1];
rz(-1.3145295) q[1];
sx q[1];
rz(-1.7728189) q[1];
x q[2];
rz(2.7790623) q[3];
sx q[3];
rz(-2.9177319) q[3];
sx q[3];
rz(-0.26765841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7097077) q[2];
sx q[2];
rz(-2.0543435) q[2];
sx q[2];
rz(1.0464) q[2];
rz(-2.8228068) q[3];
sx q[3];
rz(-0.71109486) q[3];
sx q[3];
rz(-0.54767245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043561291) q[0];
sx q[0];
rz(-1.8836319) q[0];
sx q[0];
rz(2.4936254) q[0];
rz(1.3994392) q[1];
sx q[1];
rz(-1.6439227) q[1];
sx q[1];
rz(1.8125777) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8662939) q[0];
sx q[0];
rz(-1.6801375) q[0];
sx q[0];
rz(-1.7924395) q[0];
rz(-pi) q[1];
rz(-0.29000303) q[2];
sx q[2];
rz(-1.776374) q[2];
sx q[2];
rz(-0.21943352) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4635515) q[1];
sx q[1];
rz(-1.8811418) q[1];
sx q[1];
rz(2.2928659) q[1];
rz(-pi) q[2];
rz(-0.51402199) q[3];
sx q[3];
rz(-1.2482621) q[3];
sx q[3];
rz(-0.0016101282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38137388) q[2];
sx q[2];
rz(-1.9332644) q[2];
sx q[2];
rz(-1.5927429) q[2];
rz(0.069843944) q[3];
sx q[3];
rz(-1.8826238) q[3];
sx q[3];
rz(-0.86863345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182564) q[0];
sx q[0];
rz(-1.9255487) q[0];
sx q[0];
rz(-2.1997531) q[0];
rz(2.9815004) q[1];
sx q[1];
rz(-1.4804877) q[1];
sx q[1];
rz(0.19217415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7640141) q[0];
sx q[0];
rz(-1.6068216) q[0];
sx q[0];
rz(1.7948304) q[0];
x q[1];
rz(1.501154) q[2];
sx q[2];
rz(-2.6187839) q[2];
sx q[2];
rz(-1.4774587) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6814179) q[1];
sx q[1];
rz(-2.3263086) q[1];
sx q[1];
rz(2.6667206) q[1];
x q[2];
rz(-1.0732157) q[3];
sx q[3];
rz(-1.5514152) q[3];
sx q[3];
rz(1.3598702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75752246) q[2];
sx q[2];
rz(-1.2382058) q[2];
sx q[2];
rz(-2.7080217) q[2];
rz(-1.5844257) q[3];
sx q[3];
rz(-1.6470563) q[3];
sx q[3];
rz(2.7092194) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6381391) q[0];
sx q[0];
rz(-1.764955) q[0];
sx q[0];
rz(1.3442511) q[0];
rz(1.9380219) q[1];
sx q[1];
rz(-1.1886339) q[1];
sx q[1];
rz(-1.3628091) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575131) q[0];
sx q[0];
rz(-1.5410893) q[0];
sx q[0];
rz(-0.021935181) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8261189) q[2];
sx q[2];
rz(-1.6996133) q[2];
sx q[2];
rz(-1.4908189) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2143933) q[1];
sx q[1];
rz(-1.3424557) q[1];
sx q[1];
rz(1.4424999) q[1];
rz(-2.8239488) q[3];
sx q[3];
rz(-2.2322725) q[3];
sx q[3];
rz(0.7484439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56889304) q[2];
sx q[2];
rz(-0.77304825) q[2];
sx q[2];
rz(0.9838689) q[2];
rz(-0.24108663) q[3];
sx q[3];
rz(-0.9328931) q[3];
sx q[3];
rz(0.69055313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6702061) q[0];
sx q[0];
rz(-2.3455878) q[0];
sx q[0];
rz(0.27467003) q[0];
rz(1.341691) q[1];
sx q[1];
rz(-0.62961737) q[1];
sx q[1];
rz(-2.0955657) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957841) q[0];
sx q[0];
rz(-1.2465451) q[0];
sx q[0];
rz(-0.62600531) q[0];
rz(-pi) q[1];
rz(1.3730896) q[2];
sx q[2];
rz(-1.4914163) q[2];
sx q[2];
rz(-2.6153836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99524263) q[1];
sx q[1];
rz(-2.1533794) q[1];
sx q[1];
rz(1.4928774) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.041954354) q[3];
sx q[3];
rz(-1.0681515) q[3];
sx q[3];
rz(1.6142538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2785953) q[2];
sx q[2];
rz(-1.3906761) q[2];
sx q[2];
rz(-0.39946237) q[2];
rz(2.5755889) q[3];
sx q[3];
rz(-0.96650201) q[3];
sx q[3];
rz(2.32617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28867662) q[0];
sx q[0];
rz(-1.7221907) q[0];
sx q[0];
rz(-0.5823108) q[0];
rz(0.18320006) q[1];
sx q[1];
rz(-0.87545005) q[1];
sx q[1];
rz(0.80642548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9574653) q[0];
sx q[0];
rz(-0.13617198) q[0];
sx q[0];
rz(-2.2720112) q[0];
x q[1];
rz(-0.66566531) q[2];
sx q[2];
rz(-2.7241926) q[2];
sx q[2];
rz(1.1102499) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2576221) q[1];
sx q[1];
rz(-0.7695573) q[1];
sx q[1];
rz(-0.29410024) q[1];
x q[2];
rz(-1.459645) q[3];
sx q[3];
rz(-1.2254997) q[3];
sx q[3];
rz(0.87792859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.23239423) q[2];
sx q[2];
rz(-2.4269203) q[2];
sx q[2];
rz(-0.8052899) q[2];
rz(-2.9946839) q[3];
sx q[3];
rz(-0.78607905) q[3];
sx q[3];
rz(-2.4033191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2522226) q[0];
sx q[0];
rz(-1.3852373) q[0];
sx q[0];
rz(-1.2607384) q[0];
rz(1.5994785) q[1];
sx q[1];
rz(-1.5972932) q[1];
sx q[1];
rz(1.6398026) q[1];
rz(-3.0270544) q[2];
sx q[2];
rz(-2.2907612) q[2];
sx q[2];
rz(-1.9775122) q[2];
rz(1.9143424) q[3];
sx q[3];
rz(-2.187163) q[3];
sx q[3];
rz(1.3507395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
