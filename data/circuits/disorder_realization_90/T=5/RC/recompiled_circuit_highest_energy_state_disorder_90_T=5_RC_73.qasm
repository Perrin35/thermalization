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
rz(-0.65547216) q[0];
sx q[0];
rz(-2.1894426) q[0];
sx q[0];
rz(0.093753554) q[0];
rz(1.3682415) q[1];
sx q[1];
rz(4.7028766) q[1];
sx q[1];
rz(8.75755) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8225559) q[0];
sx q[0];
rz(-1.0780436) q[0];
sx q[0];
rz(-2.7125554) q[0];
rz(-pi) q[1];
rz(-2.480443) q[2];
sx q[2];
rz(-2.8183658) q[2];
sx q[2];
rz(-1.2502535) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2161795) q[1];
sx q[1];
rz(-2.254018) q[1];
sx q[1];
rz(-2.776078) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5205403) q[3];
sx q[3];
rz(-0.9844616) q[3];
sx q[3];
rz(-1.037327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5253456) q[2];
sx q[2];
rz(-1.6051925) q[2];
sx q[2];
rz(3.1180535) q[2];
rz(-0.25447887) q[3];
sx q[3];
rz(-1.9813709) q[3];
sx q[3];
rz(2.3833073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1940521) q[0];
sx q[0];
rz(-1.3964615) q[0];
sx q[0];
rz(0.61554712) q[0];
rz(-2.5415892) q[1];
sx q[1];
rz(-1.871385) q[1];
sx q[1];
rz(-0.34872762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48744142) q[0];
sx q[0];
rz(-0.92986996) q[0];
sx q[0];
rz(1.4750255) q[0];
x q[1];
rz(1.9127513) q[2];
sx q[2];
rz(-1.2970902) q[2];
sx q[2];
rz(0.051816377) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8932027) q[1];
sx q[1];
rz(-2.5856308) q[1];
sx q[1];
rz(-0.979579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.058699) q[3];
sx q[3];
rz(-2.5827105) q[3];
sx q[3];
rz(-0.36039613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6500924) q[2];
sx q[2];
rz(-1.3050175) q[2];
sx q[2];
rz(0.65518641) q[2];
rz(2.9789467) q[3];
sx q[3];
rz(-0.9442257) q[3];
sx q[3];
rz(2.9338525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96838897) q[0];
sx q[0];
rz(-2.0252616) q[0];
sx q[0];
rz(3.047347) q[0];
rz(1.8415797) q[1];
sx q[1];
rz(-0.899122) q[1];
sx q[1];
rz(1.947044) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.822552) q[0];
sx q[0];
rz(-1.8739432) q[0];
sx q[0];
rz(-0.52355741) q[0];
x q[1];
rz(1.2066417) q[2];
sx q[2];
rz(-1.6192163) q[2];
sx q[2];
rz(-0.044755699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49385168) q[1];
sx q[1];
rz(-1.6147227) q[1];
sx q[1];
rz(-1.8274587) q[1];
x q[2];
rz(-1.1949092) q[3];
sx q[3];
rz(-2.2975911) q[3];
sx q[3];
rz(0.90968859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82075787) q[2];
sx q[2];
rz(-1.4471549) q[2];
sx q[2];
rz(0.40795946) q[2];
rz(-1.7631433) q[3];
sx q[3];
rz(-2.2429376) q[3];
sx q[3];
rz(-0.11817008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57678643) q[0];
sx q[0];
rz(-1.3917568) q[0];
sx q[0];
rz(2.3396662) q[0];
rz(2.7469514) q[1];
sx q[1];
rz(-2.092974) q[1];
sx q[1];
rz(0.035042979) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8405404) q[0];
sx q[0];
rz(-2.4581215) q[0];
sx q[0];
rz(0.95998623) q[0];
rz(0.50582992) q[2];
sx q[2];
rz(-0.37908812) q[2];
sx q[2];
rz(1.0729147) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1815419) q[1];
sx q[1];
rz(-1.6142577) q[1];
sx q[1];
rz(3.0994013) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7862042) q[3];
sx q[3];
rz(-1.3878763) q[3];
sx q[3];
rz(-1.2661883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0617712) q[2];
sx q[2];
rz(-1.062919) q[2];
sx q[2];
rz(1.4351832) q[2];
rz(-1.7363413) q[3];
sx q[3];
rz(-2.0515714) q[3];
sx q[3];
rz(2.0315571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46471304) q[0];
sx q[0];
rz(-1.9929303) q[0];
sx q[0];
rz(-1.1418463) q[0];
rz(-0.51072085) q[1];
sx q[1];
rz(-1.9074214) q[1];
sx q[1];
rz(1.4265192) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620112) q[0];
sx q[0];
rz(-1.0252165) q[0];
sx q[0];
rz(2.1929492) q[0];
rz(-pi) q[1];
rz(2.9147706) q[2];
sx q[2];
rz(-1.2661627) q[2];
sx q[2];
rz(0.90639988) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.74523679) q[1];
sx q[1];
rz(-2.5380683) q[1];
sx q[1];
rz(1.3717432) q[1];
rz(-pi) q[2];
rz(2.9742091) q[3];
sx q[3];
rz(-1.9731083) q[3];
sx q[3];
rz(-1.2237807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4003754) q[2];
sx q[2];
rz(-1.3072689) q[2];
sx q[2];
rz(0.2571787) q[2];
rz(-2.2687965) q[3];
sx q[3];
rz(-1.8200487) q[3];
sx q[3];
rz(1.1375554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207101) q[0];
sx q[0];
rz(-0.45419422) q[0];
sx q[0];
rz(-1.0783476) q[0];
rz(0.9067761) q[1];
sx q[1];
rz(-2.4557143) q[1];
sx q[1];
rz(-0.041898601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8804733) q[0];
sx q[0];
rz(-2.6676237) q[0];
sx q[0];
rz(1.4920477) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0187835) q[2];
sx q[2];
rz(-0.90722668) q[2];
sx q[2];
rz(-1.0994436) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4387233) q[1];
sx q[1];
rz(-1.8449191) q[1];
sx q[1];
rz(-1.1945748) q[1];
rz(-pi) q[2];
rz(-2.3056729) q[3];
sx q[3];
rz(-1.984388) q[3];
sx q[3];
rz(-1.6294162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.28356734) q[2];
sx q[2];
rz(-2.4391386) q[2];
sx q[2];
rz(0.89135998) q[2];
rz(-2.4274965) q[3];
sx q[3];
rz(-0.73914206) q[3];
sx q[3];
rz(1.063063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63356584) q[0];
sx q[0];
rz(-1.2440246) q[0];
sx q[0];
rz(-0.6749534) q[0];
rz(-1.630111) q[1];
sx q[1];
rz(-0.62050301) q[1];
sx q[1];
rz(-2.564548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2733611) q[0];
sx q[0];
rz(-1.8087862) q[0];
sx q[0];
rz(-1.0704051) q[0];
x q[1];
rz(2.9222832) q[2];
sx q[2];
rz(-2.6990934) q[2];
sx q[2];
rz(2.2377917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8220209) q[1];
sx q[1];
rz(-0.66950018) q[1];
sx q[1];
rz(-1.4386402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87793276) q[3];
sx q[3];
rz(-2.0725277) q[3];
sx q[3];
rz(2.2449765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6881037) q[2];
sx q[2];
rz(-2.0192396) q[2];
sx q[2];
rz(2.2171059) q[2];
rz(1.2201355) q[3];
sx q[3];
rz(-0.28885463) q[3];
sx q[3];
rz(-0.060062241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39153758) q[0];
sx q[0];
rz(-0.67111641) q[0];
sx q[0];
rz(1.3439939) q[0];
rz(-1.2229819) q[1];
sx q[1];
rz(-1.5593301) q[1];
sx q[1];
rz(1.5498243) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6236728) q[0];
sx q[0];
rz(-2.8123283) q[0];
sx q[0];
rz(-0.98401208) q[0];
rz(0.28002589) q[2];
sx q[2];
rz(-2.7529) q[2];
sx q[2];
rz(-2.1469146) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2651974) q[1];
sx q[1];
rz(-1.9361456) q[1];
sx q[1];
rz(-2.4202034) q[1];
rz(2.377691) q[3];
sx q[3];
rz(-1.6391067) q[3];
sx q[3];
rz(-2.1602283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0595155) q[2];
sx q[2];
rz(-2.1166708) q[2];
sx q[2];
rz(-2.9151741) q[2];
rz(-3.1021127) q[3];
sx q[3];
rz(-1.6624781) q[3];
sx q[3];
rz(-1.1994908) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74462849) q[0];
sx q[0];
rz(-1.9292984) q[0];
sx q[0];
rz(-0.82897559) q[0];
rz(3.0204311) q[1];
sx q[1];
rz(-2.1794901) q[1];
sx q[1];
rz(1.6896348) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0110553) q[0];
sx q[0];
rz(-2.4086039) q[0];
sx q[0];
rz(-0.20662465) q[0];
rz(2.0214861) q[2];
sx q[2];
rz(-2.3419437) q[2];
sx q[2];
rz(-0.19935184) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0740856) q[1];
sx q[1];
rz(-1.1065919) q[1];
sx q[1];
rz(-0.84949018) q[1];
rz(-pi) q[2];
rz(2.1249346) q[3];
sx q[3];
rz(-1.9151558) q[3];
sx q[3];
rz(-2.240802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5539603) q[2];
sx q[2];
rz(-0.82531896) q[2];
sx q[2];
rz(0.86453214) q[2];
rz(2.4899321) q[3];
sx q[3];
rz(-1.0385907) q[3];
sx q[3];
rz(0.017875044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.637735) q[0];
sx q[0];
rz(-0.88541579) q[0];
sx q[0];
rz(0.46022415) q[0];
rz(-1.3453311) q[1];
sx q[1];
rz(-2.5049152) q[1];
sx q[1];
rz(1.771079) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.329418) q[0];
sx q[0];
rz(-1.8746601) q[0];
sx q[0];
rz(1.162975) q[0];
rz(0.6376736) q[2];
sx q[2];
rz(-1.3402001) q[2];
sx q[2];
rz(-1.892145) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1980236) q[1];
sx q[1];
rz(-2.7881099) q[1];
sx q[1];
rz(-0.82228735) q[1];
x q[2];
rz(0.053456177) q[3];
sx q[3];
rz(-1.7776907) q[3];
sx q[3];
rz(-0.89689909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3687849) q[2];
sx q[2];
rz(-1.3489172) q[2];
sx q[2];
rz(-2.9985912) q[2];
rz(-0.16608206) q[3];
sx q[3];
rz(-2.4487285) q[3];
sx q[3];
rz(-0.79296976) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3731257) q[0];
sx q[0];
rz(-1.4906727) q[0];
sx q[0];
rz(1.3937108) q[0];
rz(1.985818) q[1];
sx q[1];
rz(-2.0230237) q[1];
sx q[1];
rz(0.3442234) q[1];
rz(2.106582) q[2];
sx q[2];
rz(-0.9101609) q[2];
sx q[2];
rz(2.6725651) q[2];
rz(-2.0438016) q[3];
sx q[3];
rz(-0.65752959) q[3];
sx q[3];
rz(0.65617954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
