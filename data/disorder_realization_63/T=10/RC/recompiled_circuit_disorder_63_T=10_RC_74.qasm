OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(-0.35342616) q[0];
sx q[0];
rz(1.0647635) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(-0.53607166) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9809496) q[0];
sx q[0];
rz(-2.0619218) q[0];
sx q[0];
rz(0.26382291) q[0];
rz(-0.067692368) q[2];
sx q[2];
rz(-2.2798685) q[2];
sx q[2];
rz(2.6807705) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2546665) q[1];
sx q[1];
rz(-2.1517422) q[1];
sx q[1];
rz(-1.7960153) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6081798) q[3];
sx q[3];
rz(-0.80848137) q[3];
sx q[3];
rz(1.4097139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(-0.70835152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3268711) q[0];
sx q[0];
rz(-2.673809) q[0];
sx q[0];
rz(2.4844556) q[0];
rz(0.96946851) q[2];
sx q[2];
rz(-1.6919486) q[2];
sx q[2];
rz(-3.0960992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(-0.80177387) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80531081) q[3];
sx q[3];
rz(-2.478963) q[3];
sx q[3];
rz(-0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-2.8498245) q[2];
rz(-0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3290688) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(-1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-0.57166878) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0288491) q[0];
sx q[0];
rz(-2.1035806) q[0];
sx q[0];
rz(1.4677731) q[0];
rz(0.35927202) q[2];
sx q[2];
rz(-1.9810988) q[2];
sx q[2];
rz(-1.8112195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0464222) q[1];
sx q[1];
rz(-1.4383525) q[1];
sx q[1];
rz(-2.4080647) q[1];
rz(-pi) q[2];
x q[2];
rz(1.709135) q[3];
sx q[3];
rz(-2.0190911) q[3];
sx q[3];
rz(1.5642779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9812575) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(-0.70880115) q[2];
rz(2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70808327) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(0.79950142) q[0];
rz(0.049830534) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(-2.9611011) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7364396) q[0];
sx q[0];
rz(-2.7442051) q[0];
sx q[0];
rz(-2.9050164) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7724491) q[2];
sx q[2];
rz(-1.7338848) q[2];
sx q[2];
rz(-0.8085608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6199477) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(-0.32089969) q[1];
x q[2];
rz(1.0749712) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(1.6881975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.9767438) q[2];
sx q[2];
rz(-1.583741) q[2];
rz(1.7662988) q[3];
sx q[3];
rz(-1.8245274) q[3];
sx q[3];
rz(-2.2089675) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76064008) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(-0.98651648) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.2170075) q[1];
sx q[1];
rz(2.8900237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35206301) q[0];
sx q[0];
rz(-1.3955727) q[0];
sx q[0];
rz(0.87900889) q[0];
rz(-2.6816363) q[2];
sx q[2];
rz(-1.6486247) q[2];
sx q[2];
rz(-0.83173448) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0136452) q[1];
sx q[1];
rz(-2.9442759) q[1];
sx q[1];
rz(3.0017588) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98925791) q[3];
sx q[3];
rz(-2.3495418) q[3];
sx q[3];
rz(-2.0056412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.95785561) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-0.95932275) q[3];
sx q[3];
rz(2.5103536) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3866766) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(1.9592346) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.052554) q[0];
sx q[0];
rz(-2.7900643) q[0];
sx q[0];
rz(-1.0765443) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82627798) q[2];
sx q[2];
rz(-1.7457361) q[2];
sx q[2];
rz(0.28184055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1054354) q[1];
sx q[1];
rz(-1.565355) q[1];
sx q[1];
rz(-2.2219707) q[1];
rz(1.844448) q[3];
sx q[3];
rz(-1.2840464) q[3];
sx q[3];
rz(-1.2424038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6958375) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(-1.5234891) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80124827) q[0];
sx q[0];
rz(-0.72754117) q[0];
sx q[0];
rz(1.9968541) q[0];
rz(-pi) q[1];
rz(-2.0231831) q[2];
sx q[2];
rz(-2.2707319) q[2];
sx q[2];
rz(1.5024904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2318864) q[1];
sx q[1];
rz(-1.6344374) q[1];
sx q[1];
rz(1.015889) q[1];
x q[2];
rz(-2.7695914) q[3];
sx q[3];
rz(-1.5270546) q[3];
sx q[3];
rz(1.352076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-0.19006426) q[2];
rz(2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(-0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0311325) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-0.53623143) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3535239) q[0];
sx q[0];
rz(-0.81809645) q[0];
sx q[0];
rz(1.8188994) q[0];
x q[1];
rz(1.0529222) q[2];
sx q[2];
rz(-1.1793943) q[2];
sx q[2];
rz(-2.1640167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4211385) q[1];
sx q[1];
rz(-0.23067833) q[1];
sx q[1];
rz(-2.8645664) q[1];
rz(-1.6471345) q[3];
sx q[3];
rz(-1.9266955) q[3];
sx q[3];
rz(-0.59698856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(-2.9013157) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(-0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48052045) q[0];
sx q[0];
rz(-1.6854223) q[0];
sx q[0];
rz(-1.4474439) q[0];
x q[1];
rz(1.4264832) q[2];
sx q[2];
rz(-1.2735954) q[2];
sx q[2];
rz(0.65629634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7749412) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(1.8371546) q[1];
rz(-pi) q[2];
rz(2.3064763) q[3];
sx q[3];
rz(-1.1856688) q[3];
sx q[3];
rz(-1.9851828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(1.4871037) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47637981) q[0];
sx q[0];
rz(-1.3276275) q[0];
sx q[0];
rz(2.5961844) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0087183) q[2];
sx q[2];
rz(-2.808411) q[2];
sx q[2];
rz(-2.9269232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6259389) q[1];
sx q[1];
rz(-1.3733555) q[1];
sx q[1];
rz(2.291912) q[1];
x q[2];
rz(2.8285633) q[3];
sx q[3];
rz(-2.1306681) q[3];
sx q[3];
rz(-2.5525639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(-0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(-1.2623513) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1508355) q[0];
sx q[0];
rz(-2.3322454) q[0];
sx q[0];
rz(-1.1582751) q[0];
rz(1.5785718) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(-0.14317748) q[2];
sx q[2];
rz(-1.8191874) q[2];
sx q[2];
rz(-3.0468265) q[2];
rz(-1.8347918) q[3];
sx q[3];
rz(-1.055607) q[3];
sx q[3];
rz(-1.9029688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];