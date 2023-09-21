OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(1.6341524) q[0];
rz(-1.545067) q[1];
sx q[1];
rz(-2.5453321) q[1];
sx q[1];
rz(2.526386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98696729) q[0];
sx q[0];
rz(-2.0320503) q[0];
sx q[0];
rz(-2.2485562) q[0];
x q[1];
rz(2.7841714) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-2.0703966) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2088036) q[1];
sx q[1];
rz(-2.1425769) q[1];
sx q[1];
rz(-2.6245481) q[1];
x q[2];
rz(2.6184222) q[3];
sx q[3];
rz(-2.7912931) q[3];
sx q[3];
rz(1.3594128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(1.7017378) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.171339) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-0.030348226) q[0];
rz(-3.0753823) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(-1.5240086) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46166566) q[0];
sx q[0];
rz(-1.7704417) q[0];
sx q[0];
rz(0.0017077831) q[0];
rz(-pi) q[1];
rz(-1.1163887) q[2];
sx q[2];
rz(-1.5523567) q[2];
sx q[2];
rz(1.4510643) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1387353) q[1];
sx q[1];
rz(-2.5668975) q[1];
sx q[1];
rz(1.6981305) q[1];
x q[2];
rz(1.7582943) q[3];
sx q[3];
rz(-1.4783995) q[3];
sx q[3];
rz(2.932991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7559738) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(1.9937817) q[2];
rz(1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148934) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-0.31578627) q[0];
rz(-2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2825851) q[0];
sx q[0];
rz(-2.1827645) q[0];
sx q[0];
rz(-0.99143272) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71528541) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(2.8298024) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38110106) q[1];
sx q[1];
rz(-0.64844202) q[1];
sx q[1];
rz(1.2566503) q[1];
rz(-pi) q[2];
rz(-1.2139981) q[3];
sx q[3];
rz(-1.1988415) q[3];
sx q[3];
rz(1.7372204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(-0.034051731) q[2];
rz(-3.1241336) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(-2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(1.5267641) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(-0.70708752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0079460572) q[0];
sx q[0];
rz(-2.5292853) q[0];
sx q[0];
rz(0.72511073) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44666501) q[2];
sx q[2];
rz(-1.457505) q[2];
sx q[2];
rz(-0.54892533) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45338079) q[1];
sx q[1];
rz(-0.66426316) q[1];
sx q[1];
rz(-1.7654256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2624192) q[3];
sx q[3];
rz(-1.7763419) q[3];
sx q[3];
rz(1.2984848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(0.46009955) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(2.8869693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4655612) q[0];
sx q[0];
rz(-1.6499632) q[0];
sx q[0];
rz(-3.1278419) q[0];
x q[1];
rz(-2.1483634) q[2];
sx q[2];
rz(-1.5830056) q[2];
sx q[2];
rz(-2.4075367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.188272) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(0.377368) q[1];
x q[2];
rz(-0.43318627) q[3];
sx q[3];
rz(-2.232589) q[3];
sx q[3];
rz(-1.8720686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.022481) q[2];
sx q[2];
rz(-2.9412061) q[2];
sx q[2];
rz(-1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(1.0401789) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-2.9350231) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20139192) q[0];
sx q[0];
rz(-2.3513146) q[0];
sx q[0];
rz(1.2339562) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7467696) q[2];
sx q[2];
rz(-2.0645112) q[2];
sx q[2];
rz(-1.0600818) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5094604) q[1];
sx q[1];
rz(-1.3941947) q[1];
sx q[1];
rz(0.86061865) q[1];
x q[2];
rz(-1.1849095) q[3];
sx q[3];
rz(-0.69291249) q[3];
sx q[3];
rz(3.0544359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(2.3287866) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(-0.16684428) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(-1.8136224) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0921811) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(-1.1294424) q[0];
rz(0.94857256) q[2];
sx q[2];
rz(-1.9554536) q[2];
sx q[2];
rz(2.4228061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.59203397) q[1];
sx q[1];
rz(-0.70211071) q[1];
sx q[1];
rz(1.83105) q[1];
x q[2];
rz(-2.6286969) q[3];
sx q[3];
rz(-1.1937871) q[3];
sx q[3];
rz(-2.5412154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61838377) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(1.7112188) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946063) q[0];
sx q[0];
rz(-1.9817579) q[0];
sx q[0];
rz(-0.18187901) q[0];
rz(0.47422844) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(2.1906733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.788818) q[0];
sx q[0];
rz(-2.7633861) q[0];
sx q[0];
rz(-0.87120716) q[0];
rz(1.3457001) q[2];
sx q[2];
rz(-0.27268073) q[2];
sx q[2];
rz(2.8358012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2434395) q[1];
sx q[1];
rz(-2.9228518) q[1];
sx q[1];
rz(-2.8380413) q[1];
rz(-2.3202053) q[3];
sx q[3];
rz(-0.52762023) q[3];
sx q[3];
rz(1.9963095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(-2.5935796) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(-2.5089335) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178745) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(-0.65972796) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83803672) q[0];
sx q[0];
rz(-2.3065901) q[0];
sx q[0];
rz(-0.8291709) q[0];
rz(0.78166878) q[2];
sx q[2];
rz(-1.726892) q[2];
sx q[2];
rz(-0.26542703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1253423) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(2.7793105) q[1];
rz(2.1094443) q[3];
sx q[3];
rz(-1.6197455) q[3];
sx q[3];
rz(0.34070542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62347162) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-0.74111104) q[0];
sx q[0];
rz(-1.3056668) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(1.0356888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10782345) q[0];
sx q[0];
rz(-0.53350893) q[0];
sx q[0];
rz(0.67722042) q[0];
rz(-pi) q[1];
rz(-0.91949384) q[2];
sx q[2];
rz(-0.80923015) q[2];
sx q[2];
rz(2.2724255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90696883) q[1];
sx q[1];
rz(-1.3303489) q[1];
sx q[1];
rz(2.8333227) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74929897) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(2.628905) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37968996) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.5079386) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(2.7307636) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(-0.012398331) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
