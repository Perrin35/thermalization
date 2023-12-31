OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(-2.4175329) q[0];
sx q[0];
rz(1.6568503) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(0.88820052) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9125497) q[0];
sx q[0];
rz(-2.9648844) q[0];
sx q[0];
rz(-2.8852709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1076803) q[2];
sx q[2];
rz(-1.7134588) q[2];
sx q[2];
rz(-2.2906274) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29477316) q[1];
sx q[1];
rz(-1.5702489) q[1];
sx q[1];
rz(-0.16695395) q[1];
rz(-pi) q[2];
rz(-3.1021318) q[3];
sx q[3];
rz(-0.69898116) q[3];
sx q[3];
rz(2.3208502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-1.1179914) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(-2.4867687) q[0];
rz(-1.2163935) q[1];
sx q[1];
rz(-1.1647859) q[1];
sx q[1];
rz(0.12589802) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038371041) q[0];
sx q[0];
rz(-1.9403337) q[0];
sx q[0];
rz(-2.080337) q[0];
x q[1];
rz(-0.89181487) q[2];
sx q[2];
rz(-2.0795155) q[2];
sx q[2];
rz(-1.5985135) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.241908) q[1];
sx q[1];
rz(-0.13026127) q[1];
sx q[1];
rz(-1.8036519) q[1];
x q[2];
rz(-1.5829854) q[3];
sx q[3];
rz(-2.4180531) q[3];
sx q[3];
rz(0.42750588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(1.4240501) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(3.0575867) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(-1.9817339) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0214329) q[0];
sx q[0];
rz(-0.5172356) q[0];
sx q[0];
rz(-1.73818) q[0];
rz(-2.1986387) q[2];
sx q[2];
rz(-1.1477594) q[2];
sx q[2];
rz(-1.1689651) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4574036) q[1];
sx q[1];
rz(-1.1154798) q[1];
sx q[1];
rz(-1.1475569) q[1];
rz(1.0252762) q[3];
sx q[3];
rz(-1.5226411) q[3];
sx q[3];
rz(0.84850509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(0.1082871) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(-2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.165034) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(2.4404793) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(-0.28465095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4581504) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(0.083806888) q[0];
x q[1];
rz(1.2830164) q[2];
sx q[2];
rz(-1.7680401) q[2];
sx q[2];
rz(1.4892727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.22420158) q[1];
sx q[1];
rz(-0.14897878) q[1];
sx q[1];
rz(-2.1267183) q[1];
rz(-0.986226) q[3];
sx q[3];
rz(-2.3826736) q[3];
sx q[3];
rz(-2.330247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.231679) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(0.56751928) q[2];
rz(0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(-1.8150785) q[0];
rz(-1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-0.93793905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468231) q[0];
sx q[0];
rz(-1.4786647) q[0];
sx q[0];
rz(3.1196306) q[0];
rz(2.6229834) q[2];
sx q[2];
rz(-2.0687639) q[2];
sx q[2];
rz(-1.2155611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8541969) q[1];
sx q[1];
rz(-1.5970267) q[1];
sx q[1];
rz(1.6051952) q[1];
rz(-pi) q[2];
rz(-0.35477562) q[3];
sx q[3];
rz(-1.6924904) q[3];
sx q[3];
rz(2.0718758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(2.1002634) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(1.3180102) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(-2.5851137) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0628478) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(1.3486805) q[0];
rz(-pi) q[1];
rz(0.85530497) q[2];
sx q[2];
rz(-1.7024634) q[2];
sx q[2];
rz(-1.5065187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0007799) q[1];
sx q[1];
rz(-1.3197348) q[1];
sx q[1];
rz(-2.954133) q[1];
rz(-pi) q[2];
rz(-0.052080215) q[3];
sx q[3];
rz(-2.3651603) q[3];
sx q[3];
rz(-1.4610425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(-1.5103643) q[3];
sx q[3];
rz(-1.3204201) q[3];
sx q[3];
rz(-1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625793) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(2.2242916) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(0.50061217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0450889) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(1.7486497) q[0];
rz(0.31502864) q[2];
sx q[2];
rz(-2.0909967) q[2];
sx q[2];
rz(0.5859642) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2523633) q[1];
sx q[1];
rz(-1.6912618) q[1];
sx q[1];
rz(2.7389588) q[1];
rz(-pi) q[2];
rz(2.385419) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(2.365436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-2.4659757) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(-0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.6759466) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4439125) q[0];
sx q[0];
rz(-1.4364103) q[0];
sx q[0];
rz(2.0356376) q[0];
x q[1];
rz(1.6797811) q[2];
sx q[2];
rz(-1.3853067) q[2];
sx q[2];
rz(-2.2456004) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2579736) q[1];
sx q[1];
rz(-2.5273364) q[1];
sx q[1];
rz(-2.1780464) q[1];
rz(2.1727209) q[3];
sx q[3];
rz(-1.0191917) q[3];
sx q[3];
rz(1.5495891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(1.6652997) q[2];
rz(-2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(-0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575689) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(0.73721686) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(0.92528701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46426526) q[0];
sx q[0];
rz(-2.263501) q[0];
sx q[0];
rz(-0.6662743) q[0];
rz(1.4383573) q[2];
sx q[2];
rz(-2.2297511) q[2];
sx q[2];
rz(-3.1140285) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8268938) q[1];
sx q[1];
rz(-1.7539382) q[1];
sx q[1];
rz(-1.7032743) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0989283) q[3];
sx q[3];
rz(-2.3657551) q[3];
sx q[3];
rz(0.27234205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3725738) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(-2.1378689) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.2596624) q[0];
rz(2.8758077) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(-2.3840747) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.950338) q[0];
sx q[0];
rz(-1.7636824) q[0];
sx q[0];
rz(-3.0015776) q[0];
rz(-pi) q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.7014116) q[2];
sx q[2];
rz(-0.57934258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.65749189) q[1];
sx q[1];
rz(-1.3490281) q[1];
sx q[1];
rz(-2.0055254) q[1];
x q[2];
rz(-3.1086139) q[3];
sx q[3];
rz(-2.4313201) q[3];
sx q[3];
rz(0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-0.63888597) q[3];
sx q[3];
rz(-2.7975438) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(2.5095148) q[2];
sx q[2];
rz(-2.6161604) q[2];
sx q[2];
rz(-2.5642774) q[2];
rz(0.32140857) q[3];
sx q[3];
rz(-1.0151498) q[3];
sx q[3];
rz(-2.2890454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
