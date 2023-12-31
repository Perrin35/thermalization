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
rz(1.5965257) q[1];
sx q[1];
rz(2.5453321) q[1];
sx q[1];
rz(8.8095713) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98696729) q[0];
sx q[0];
rz(-2.0320503) q[0];
sx q[0];
rz(-2.2485562) q[0];
rz(-0.35742128) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(1.071196) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7429744) q[1];
sx q[1];
rz(-2.3906039) q[1];
sx q[1];
rz(2.2258334) q[1];
rz(-pi) q[2];
x q[2];
rz(1.390236) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(0.80871049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.6117233) q[2];
sx q[2];
rz(-2.8033076) q[2];
rz(-1.4398549) q[3];
sx q[3];
rz(-0.9153291) q[3];
sx q[3];
rz(-2.2556944) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(3.1112444) q[0];
rz(0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(1.5240086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.679927) q[0];
sx q[0];
rz(-1.7704417) q[0];
sx q[0];
rz(-0.0017077831) q[0];
rz(-pi) q[1];
rz(-1.1163887) q[2];
sx q[2];
rz(-1.5523567) q[2];
sx q[2];
rz(-1.6905284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6807032) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(-2.1417888) q[1];
rz(-pi) q[2];
rz(-3.0475572) q[3];
sx q[3];
rz(-1.3841076) q[3];
sx q[3];
rz(-1.3446913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38561884) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(1.1478109) q[2];
rz(-1.8148445) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148934) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(-0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4323498) q[0];
sx q[0];
rz(-2.3253257) q[0];
sx q[0];
rz(2.4791251) q[0];
rz(-pi) q[1];
rz(-0.88148586) q[2];
sx q[2];
rz(-0.93886095) q[2];
sx q[2];
rz(-1.880868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6985059) q[1];
sx q[1];
rz(-1.7585187) q[1];
sx q[1];
rz(-0.94633533) q[1];
rz(2.7470845) q[3];
sx q[3];
rz(-1.9022226) q[3];
sx q[3];
rz(-0.031771914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(2.0461369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5144192) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(-2.4345051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3084761) q[0];
sx q[0];
rz(-2.0154698) q[0];
sx q[0];
rz(2.0067257) q[0];
x q[1];
rz(0.25755067) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(2.3515153) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45338079) q[1];
sx q[1];
rz(-2.4773295) q[1];
sx q[1];
rz(1.7654256) q[1];
rz(-pi) q[2];
rz(-2.9261758) q[3];
sx q[3];
rz(-1.2691174) q[3];
sx q[3];
rz(0.33723436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84918555) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(0.24266711) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(-1.3943577) q[0];
rz(-2.0460515) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(-2.8869693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1063227) q[0];
sx q[0];
rz(-1.584504) q[0];
sx q[0];
rz(1.4916219) q[0];
rz(-pi) q[1];
x q[1];
rz(0.014572797) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(2.2968959) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6756145) q[1];
sx q[1];
rz(-1.3295768) q[1];
sx q[1];
rz(-0.66959186) q[1];
x q[2];
rz(2.7084064) q[3];
sx q[3];
rz(-0.90900366) q[3];
sx q[3];
rz(1.8720686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.3767892) q[2];
rz(1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45143932) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(-1.0401789) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(0.20656955) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6109989) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(-2.3321652) q[0];
rz(-pi) q[1];
rz(1.0429522) q[2];
sx q[2];
rz(-1.2252508) q[2];
sx q[2];
rz(-2.8258459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0530015) q[1];
sx q[1];
rz(-2.2676761) q[1];
sx q[1];
rz(-2.9104396) q[1];
rz(-pi) q[2];
rz(-0.30287403) q[3];
sx q[3];
rz(-2.2040963) q[3];
sx q[3];
rz(0.57297046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.841659) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(2.858813) q[2];
rz(-2.3287866) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(2.9747484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2151826) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.8136224) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049411557) q[0];
sx q[0];
rz(-1.4565399) q[0];
sx q[0];
rz(1.1294424) q[0];
x q[1];
rz(0.46220025) q[2];
sx q[2];
rz(-2.1415347) q[2];
sx q[2];
rz(0.58909033) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5495587) q[1];
sx q[1];
rz(-2.4394819) q[1];
sx q[1];
rz(-1.3105427) q[1];
rz(-pi) q[2];
rz(1.1442723) q[3];
sx q[3];
rz(-2.0445619) q[3];
sx q[3];
rz(-1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(-1.7810812) q[2];
rz(-1.7112188) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(-0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(2.6673642) q[1];
sx q[1];
rz(-1.0206181) q[1];
sx q[1];
rz(-2.1906733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2598341) q[0];
sx q[0];
rz(-1.330733) q[0];
sx q[0];
rz(-1.2756707) q[0];
x q[1];
rz(-1.7958926) q[2];
sx q[2];
rz(-0.27268073) q[2];
sx q[2];
rz(2.8358012) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37590313) q[1];
sx q[1];
rz(-1.5058869) q[1];
sx q[1];
rz(-0.20903559) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82138737) q[3];
sx q[3];
rz(-2.6139724) q[3];
sx q[3];
rz(1.1452831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2237079) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(3.0656832) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.8173822) q[3];
sx q[3];
rz(-0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(0.41704047) q[1];
sx q[1];
rz(-1.7224256) q[1];
sx q[1];
rz(2.4818647) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8575681) q[0];
sx q[0];
rz(-1.0463456) q[0];
sx q[0];
rz(-0.88733034) q[0];
rz(-0.78166878) q[2];
sx q[2];
rz(-1.4147007) q[2];
sx q[2];
rz(2.8761656) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0162504) q[1];
sx q[1];
rz(-0.94031912) q[1];
sx q[1];
rz(-0.36228212) q[1];
x q[2];
rz(-2.1094443) q[3];
sx q[3];
rz(-1.6197455) q[3];
sx q[3];
rz(2.8008872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(-2.1155604) q[2];
rz(2.9927411) q[3];
sx q[3];
rz(-2.1089349) q[3];
sx q[3];
rz(-3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(-1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-2.1059039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64338387) q[0];
sx q[0];
rz(-1.9783101) q[0];
sx q[0];
rz(-1.9252752) q[0];
rz(-2.2220988) q[2];
sx q[2];
rz(-2.3323625) q[2];
sx q[2];
rz(-0.86916718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2346238) q[1];
sx q[1];
rz(-1.8112438) q[1];
sx q[1];
rz(0.30827) q[1];
rz(2.8519248) q[3];
sx q[3];
rz(-2.3710459) q[3];
sx q[3];
rz(-0.25911261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.1516085) q[2];
rz(-0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(2.7307636) q[2];
sx q[2];
rz(-1.5317691) q[2];
sx q[2];
rz(-0.017824235) q[2];
rz(3.1291943) q[3];
sx q[3];
rz(-0.61840246) q[3];
sx q[3];
rz(1.7411504) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
