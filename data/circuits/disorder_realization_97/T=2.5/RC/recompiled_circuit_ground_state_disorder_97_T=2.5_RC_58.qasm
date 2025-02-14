OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.41807732) q[0];
sx q[0];
rz(3.1273841) q[0];
sx q[0];
rz(9.4655884) q[0];
rz(3.1006676) q[1];
sx q[1];
rz(-2.0857781) q[1];
sx q[1];
rz(-2.5850886) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4427821) q[0];
sx q[0];
rz(-2.7469146) q[0];
sx q[0];
rz(-2.0641607) q[0];
rz(2.7950613) q[2];
sx q[2];
rz(-2.541447) q[2];
sx q[2];
rz(2.481386) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79136363) q[1];
sx q[1];
rz(-1.204317) q[1];
sx q[1];
rz(0.54270737) q[1];
rz(-pi) q[2];
x q[2];
rz(0.09381525) q[3];
sx q[3];
rz(-1.3348186) q[3];
sx q[3];
rz(0.79506563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9875662) q[2];
sx q[2];
rz(-1.7511) q[2];
sx q[2];
rz(-1.7517368) q[2];
rz(-3.0912345) q[3];
sx q[3];
rz(-1.4202838) q[3];
sx q[3];
rz(-1.095298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9339226) q[0];
sx q[0];
rz(-1.4140797) q[0];
sx q[0];
rz(-3.0811144) q[0];
rz(0.98840797) q[1];
sx q[1];
rz(-1.6840839) q[1];
sx q[1];
rz(-2.9284533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7558862) q[0];
sx q[0];
rz(-1.1637628) q[0];
sx q[0];
rz(-0.72188053) q[0];
rz(0.11535203) q[2];
sx q[2];
rz(-1.7322503) q[2];
sx q[2];
rz(0.20736883) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7139529) q[1];
sx q[1];
rz(-0.25498495) q[1];
sx q[1];
rz(-0.1509576) q[1];
x q[2];
rz(2.0673236) q[3];
sx q[3];
rz(-1.4098216) q[3];
sx q[3];
rz(-2.9825008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.1468143) q[2];
sx q[2];
rz(-1.549336) q[2];
sx q[2];
rz(-3.1352622) q[2];
rz(-1.4396884) q[3];
sx q[3];
rz(-0.78841811) q[3];
sx q[3];
rz(0.44062135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45563662) q[0];
sx q[0];
rz(-2.7485479) q[0];
sx q[0];
rz(1.8259557) q[0];
rz(-1.3872967) q[1];
sx q[1];
rz(-0.60391128) q[1];
sx q[1];
rz(3.0953298) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14213698) q[0];
sx q[0];
rz(-1.6007413) q[0];
sx q[0];
rz(-1.5697451) q[0];
x q[1];
rz(1.9579024) q[2];
sx q[2];
rz(-2.0069242) q[2];
sx q[2];
rz(-1.8226464) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97291482) q[1];
sx q[1];
rz(-0.70446223) q[1];
sx q[1];
rz(0.55120123) q[1];
rz(1.336724) q[3];
sx q[3];
rz(-2.5137874) q[3];
sx q[3];
rz(2.1677728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5578654) q[2];
sx q[2];
rz(-1.1745619) q[2];
sx q[2];
rz(2.2306856) q[2];
rz(-1.2298498) q[3];
sx q[3];
rz(-2.3807447) q[3];
sx q[3];
rz(2.7228444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8987027) q[0];
sx q[0];
rz(-1.7430584) q[0];
sx q[0];
rz(-2.718495) q[0];
rz(0.93370071) q[1];
sx q[1];
rz(-1.2010778) q[1];
sx q[1];
rz(0.43412128) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86193) q[0];
sx q[0];
rz(-1.2796784) q[0];
sx q[0];
rz(1.5345957) q[0];
x q[1];
rz(-1.7291315) q[2];
sx q[2];
rz(-0.91671093) q[2];
sx q[2];
rz(0.20928644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4695376) q[1];
sx q[1];
rz(-2.5651882) q[1];
sx q[1];
rz(0.23643929) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0634123) q[3];
sx q[3];
rz(-0.23325167) q[3];
sx q[3];
rz(0.95913974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9551129) q[2];
sx q[2];
rz(-0.13804144) q[2];
sx q[2];
rz(-1.0008) q[2];
rz(-0.75438386) q[3];
sx q[3];
rz(-2.3013134) q[3];
sx q[3];
rz(1.8050571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5767625) q[0];
sx q[0];
rz(-0.63693988) q[0];
sx q[0];
rz(1.7621967) q[0];
rz(0.55343902) q[1];
sx q[1];
rz(-1.6575927) q[1];
sx q[1];
rz(-0.71054202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64122771) q[0];
sx q[0];
rz(-0.7754179) q[0];
sx q[0];
rz(-1.3432608) q[0];
rz(-pi) q[1];
rz(0.24181409) q[2];
sx q[2];
rz(-1.6266357) q[2];
sx q[2];
rz(-0.18130359) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2031999) q[1];
sx q[1];
rz(-1.6050996) q[1];
sx q[1];
rz(-1.3399948) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8666302) q[3];
sx q[3];
rz(-2.5841613) q[3];
sx q[3];
rz(-2.9334588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1826627) q[2];
sx q[2];
rz(-1.6561597) q[2];
sx q[2];
rz(3.0976234) q[2];
rz(2.8751657) q[3];
sx q[3];
rz(-0.1259585) q[3];
sx q[3];
rz(0.21336475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63785315) q[0];
sx q[0];
rz(-1.6197236) q[0];
sx q[0];
rz(2.7730686) q[0];
rz(-0.35399327) q[1];
sx q[1];
rz(-1.1292255) q[1];
sx q[1];
rz(2.0781719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2475488) q[0];
sx q[0];
rz(-0.14078939) q[0];
sx q[0];
rz(-1.6611499) q[0];
rz(-pi) q[1];
rz(-3.0292461) q[2];
sx q[2];
rz(-1.8469166) q[2];
sx q[2];
rz(0.39673765) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8700347) q[1];
sx q[1];
rz(-1.9816625) q[1];
sx q[1];
rz(0.4510057) q[1];
rz(-pi) q[2];
rz(-2.3928349) q[3];
sx q[3];
rz(-1.4331265) q[3];
sx q[3];
rz(-2.2434821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4949558) q[2];
sx q[2];
rz(-1.8821303) q[2];
sx q[2];
rz(1.4259526) q[2];
rz(-0.87092733) q[3];
sx q[3];
rz(-2.0723074) q[3];
sx q[3];
rz(-0.23437414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17212412) q[0];
sx q[0];
rz(-0.5640465) q[0];
sx q[0];
rz(-3.1308351) q[0];
rz(2.5116008) q[1];
sx q[1];
rz(-2.0959581) q[1];
sx q[1];
rz(0.90525544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2242611) q[0];
sx q[0];
rz(-2.1035827) q[0];
sx q[0];
rz(-0.075418069) q[0];
rz(-1.0192746) q[2];
sx q[2];
rz(-2.5576303) q[2];
sx q[2];
rz(2.0083186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.9295302) q[1];
sx q[1];
rz(-1.7133949) q[1];
sx q[1];
rz(2.6349956) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91628381) q[3];
sx q[3];
rz(-1.9184417) q[3];
sx q[3];
rz(0.8485444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4524727) q[2];
sx q[2];
rz(-1.2484739) q[2];
sx q[2];
rz(-3.0242237) q[2];
rz(-0.72702414) q[3];
sx q[3];
rz(-0.89594642) q[3];
sx q[3];
rz(2.0015543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3177976) q[0];
sx q[0];
rz(-2.0841086) q[0];
sx q[0];
rz(2.8463038) q[0];
rz(-1.3914039) q[1];
sx q[1];
rz(-0.65381217) q[1];
sx q[1];
rz(-0.48694912) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21509296) q[0];
sx q[0];
rz(-1.7679259) q[0];
sx q[0];
rz(-1.6736945) q[0];
rz(2.0854961) q[2];
sx q[2];
rz(-1.2916512) q[2];
sx q[2];
rz(-2.5893958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0978969) q[1];
sx q[1];
rz(-2.0123031) q[1];
sx q[1];
rz(-2.312078) q[1];
x q[2];
rz(-1.7510909) q[3];
sx q[3];
rz(-1.3077985) q[3];
sx q[3];
rz(-0.20151787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.22750172) q[2];
sx q[2];
rz(-2.5469683) q[2];
sx q[2];
rz(1.5458858) q[2];
rz(-2.5325736) q[3];
sx q[3];
rz(-1.1142497) q[3];
sx q[3];
rz(-2.473089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.2862227) q[0];
sx q[0];
rz(-0.90317059) q[0];
sx q[0];
rz(-1.5089996) q[0];
rz(1.1434932) q[1];
sx q[1];
rz(-1.1610463) q[1];
sx q[1];
rz(-0.11018363) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78980434) q[0];
sx q[0];
rz(-0.98362255) q[0];
sx q[0];
rz(-1.4465062) q[0];
rz(-1.0855488) q[2];
sx q[2];
rz(-0.51328906) q[2];
sx q[2];
rz(-1.1977149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5563123) q[1];
sx q[1];
rz(-1.9686799) q[1];
sx q[1];
rz(-0.75979397) q[1];
rz(-0.0024282495) q[3];
sx q[3];
rz(-2.6442827) q[3];
sx q[3];
rz(-2.8893472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1630359) q[2];
sx q[2];
rz(-2.0608993) q[2];
sx q[2];
rz(-0.85868305) q[2];
rz(-1.5052172) q[3];
sx q[3];
rz(-1.7965763) q[3];
sx q[3];
rz(1.4428562) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34864146) q[0];
sx q[0];
rz(-1.071278) q[0];
sx q[0];
rz(-2.4559378) q[0];
rz(2.3220956) q[1];
sx q[1];
rz(-1.6623431) q[1];
sx q[1];
rz(-2.562838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1534655) q[0];
sx q[0];
rz(-2.1442226) q[0];
sx q[0];
rz(-0.98116409) q[0];
rz(1.9482934) q[2];
sx q[2];
rz(-0.60040081) q[2];
sx q[2];
rz(2.3410554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.600665) q[1];
sx q[1];
rz(-1.828421) q[1];
sx q[1];
rz(-2.7476946) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10921918) q[3];
sx q[3];
rz(-1.5223323) q[3];
sx q[3];
rz(0.27849659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3026352) q[2];
sx q[2];
rz(-2.9557648) q[2];
sx q[2];
rz(-0.1798943) q[2];
rz(2.5707865) q[3];
sx q[3];
rz(-2.4195318) q[3];
sx q[3];
rz(-0.65176898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30052341) q[0];
sx q[0];
rz(-2.4916334) q[0];
sx q[0];
rz(-1.4148225) q[0];
rz(-0.48814804) q[1];
sx q[1];
rz(-0.5562677) q[1];
sx q[1];
rz(1.5604896) q[1];
rz(-0.93964259) q[2];
sx q[2];
rz(-1.8710536) q[2];
sx q[2];
rz(-1.6761129) q[2];
rz(-2.6861784) q[3];
sx q[3];
rz(-2.5546163) q[3];
sx q[3];
rz(3.0163812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
