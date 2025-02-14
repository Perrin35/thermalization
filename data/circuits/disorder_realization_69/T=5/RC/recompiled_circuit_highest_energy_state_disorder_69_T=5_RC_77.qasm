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
rz(0.34613553) q[0];
sx q[0];
rz(-1.5032285) q[0];
sx q[0];
rz(0.7315973) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(2.7076758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2303378) q[0];
sx q[0];
rz(-0.83090913) q[0];
sx q[0];
rz(-3.0940542) q[0];
rz(-pi) q[1];
rz(-2.5662759) q[2];
sx q[2];
rz(-1.0566718) q[2];
sx q[2];
rz(0.34808394) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0632532) q[1];
sx q[1];
rz(-1.9615478) q[1];
sx q[1];
rz(0.88576742) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1103575) q[3];
sx q[3];
rz(-1.7618351) q[3];
sx q[3];
rz(1.7241378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71243858) q[2];
sx q[2];
rz(-1.466208) q[2];
sx q[2];
rz(-1.0038556) q[2];
rz(-0.53576523) q[3];
sx q[3];
rz(-0.86447132) q[3];
sx q[3];
rz(-0.0011477688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67343229) q[0];
sx q[0];
rz(-2.4591481) q[0];
sx q[0];
rz(-2.4144507) q[0];
rz(-0.06761059) q[1];
sx q[1];
rz(-1.7865684) q[1];
sx q[1];
rz(1.1236069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1557187) q[0];
sx q[0];
rz(-2.576568) q[0];
sx q[0];
rz(-0.54604097) q[0];
x q[1];
rz(-3.0762663) q[2];
sx q[2];
rz(-1.6931173) q[2];
sx q[2];
rz(2.9768012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3930291) q[1];
sx q[1];
rz(-1.8504228) q[1];
sx q[1];
rz(-3.0381027) q[1];
rz(-1.0960078) q[3];
sx q[3];
rz(-1.9663621) q[3];
sx q[3];
rz(2.6717348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8565389) q[2];
sx q[2];
rz(-1.5679789) q[2];
sx q[2];
rz(-2.6598568) q[2];
rz(-0.5528062) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(-3.0520181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8168617) q[0];
sx q[0];
rz(-0.4751927) q[0];
sx q[0];
rz(-3.0480296) q[0];
rz(1.3803253) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(0.63724744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475315) q[0];
sx q[0];
rz(-1.1406745) q[0];
sx q[0];
rz(-0.57292666) q[0];
rz(-pi) q[1];
rz(-2.5334444) q[2];
sx q[2];
rz(-1.7575193) q[2];
sx q[2];
rz(3.1352459) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67575021) q[1];
sx q[1];
rz(-1.9681962) q[1];
sx q[1];
rz(3.1134362) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1430238) q[3];
sx q[3];
rz(-0.29353729) q[3];
sx q[3];
rz(-2.0484201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9351585) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(2.4578102) q[2];
rz(2.911496) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(-0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0919331) q[0];
sx q[0];
rz(-0.50685087) q[0];
sx q[0];
rz(1.4403213) q[0];
rz(-0.47538844) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(-0.58116523) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064539486) q[0];
sx q[0];
rz(-2.4145899) q[0];
sx q[0];
rz(1.6995656) q[0];
rz(-pi) q[1];
rz(-1.1240684) q[2];
sx q[2];
rz(-2.6368195) q[2];
sx q[2];
rz(-3.1044132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95422259) q[1];
sx q[1];
rz(-2.0285602) q[1];
sx q[1];
rz(1.9083896) q[1];
rz(-pi) q[2];
rz(-1.4122333) q[3];
sx q[3];
rz(-1.0845636) q[3];
sx q[3];
rz(2.9610046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0541957) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(-2.6206214) q[2];
rz(-1.4969131) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5533376) q[0];
sx q[0];
rz(-2.8659358) q[0];
sx q[0];
rz(1.1302554) q[0];
rz(2.0626119) q[1];
sx q[1];
rz(-1.1993473) q[1];
sx q[1];
rz(-2.030453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95100194) q[0];
sx q[0];
rz(-1.7365121) q[0];
sx q[0];
rz(1.1072876) q[0];
rz(-pi) q[1];
rz(-0.52972858) q[2];
sx q[2];
rz(-1.6090983) q[2];
sx q[2];
rz(-1.2438347) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4823522) q[1];
sx q[1];
rz(-0.9352881) q[1];
sx q[1];
rz(-0.19584943) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2258198) q[3];
sx q[3];
rz(-2.0213599) q[3];
sx q[3];
rz(-0.612606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7378716) q[2];
sx q[2];
rz(-2.626494) q[2];
sx q[2];
rz(-2.1007288) q[2];
rz(-0.8199842) q[3];
sx q[3];
rz(-1.9085725) q[3];
sx q[3];
rz(-0.6238873) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89989221) q[0];
sx q[0];
rz(-0.59881678) q[0];
sx q[0];
rz(-2.4851121) q[0];
rz(1.3735324) q[1];
sx q[1];
rz(-2.3479925) q[1];
sx q[1];
rz(2.0097282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053018) q[0];
sx q[0];
rz(-1.7680918) q[0];
sx q[0];
rz(-0.70081607) q[0];
x q[1];
rz(0.67645881) q[2];
sx q[2];
rz(-2.0148811) q[2];
sx q[2];
rz(-1.9537587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1161728) q[1];
sx q[1];
rz(-2.4301404) q[1];
sx q[1];
rz(-1.5283714) q[1];
rz(-pi) q[2];
rz(-3.1127315) q[3];
sx q[3];
rz(-0.15046826) q[3];
sx q[3];
rz(-0.34666016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.446283) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(-2.9122635) q[3];
sx q[3];
rz(-2.1181483) q[3];
sx q[3];
rz(-2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0685773) q[0];
sx q[0];
rz(-1.1443161) q[0];
sx q[0];
rz(1.2314433) q[0];
rz(0.07864174) q[1];
sx q[1];
rz(-1.4631203) q[1];
sx q[1];
rz(2.9873649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336727) q[0];
sx q[0];
rz(-1.3929954) q[0];
sx q[0];
rz(-2.7975797) q[0];
rz(-pi) q[1];
rz(1.9811822) q[2];
sx q[2];
rz(-1.2912116) q[2];
sx q[2];
rz(-1.2321763) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7423305) q[1];
sx q[1];
rz(-0.62472099) q[1];
sx q[1];
rz(0.66466753) q[1];
rz(-pi) q[2];
rz(2.6811872) q[3];
sx q[3];
rz(-0.45259991) q[3];
sx q[3];
rz(-0.79595882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33588931) q[2];
sx q[2];
rz(-1.615639) q[2];
sx q[2];
rz(1.6501144) q[2];
rz(1.7283745) q[3];
sx q[3];
rz(-1.3946984) q[3];
sx q[3];
rz(-1.570805) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042628057) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(-1.4549103) q[0];
rz(0.97575724) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(-1.3386493) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0473235) q[0];
sx q[0];
rz(-1.1828831) q[0];
sx q[0];
rz(2.3045425) q[0];
x q[1];
rz(2.8076914) q[2];
sx q[2];
rz(-1.0907764) q[2];
sx q[2];
rz(1.0261818) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8725135) q[1];
sx q[1];
rz(-0.95035205) q[1];
sx q[1];
rz(-0.2618813) q[1];
rz(-pi) q[2];
rz(-2.678431) q[3];
sx q[3];
rz(-0.71802959) q[3];
sx q[3];
rz(-2.5754186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6300388) q[2];
sx q[2];
rz(-1.0147107) q[2];
sx q[2];
rz(0.83941984) q[2];
rz(-1.2342341) q[3];
sx q[3];
rz(-1.0890361) q[3];
sx q[3];
rz(-0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79427528) q[0];
sx q[0];
rz(-1.3823771) q[0];
sx q[0];
rz(2.7224702) q[0];
rz(-1.6436815) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(2.7083414) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8041291) q[0];
sx q[0];
rz(-0.47352284) q[0];
sx q[0];
rz(1.3086523) q[0];
x q[1];
rz(-0.017288329) q[2];
sx q[2];
rz(-2.4845036) q[2];
sx q[2];
rz(1.8985871) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9135838) q[1];
sx q[1];
rz(-0.59803793) q[1];
sx q[1];
rz(1.0526471) q[1];
rz(-pi) q[2];
rz(-2.8178704) q[3];
sx q[3];
rz(-2.5383213) q[3];
sx q[3];
rz(-2.3761185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7822632) q[2];
sx q[2];
rz(-1.1507582) q[2];
sx q[2];
rz(2.9456054) q[2];
rz(-1.1228784) q[3];
sx q[3];
rz(-2.4298318) q[3];
sx q[3];
rz(-1.6897197) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6561683) q[0];
sx q[0];
rz(-0.66101414) q[0];
sx q[0];
rz(2.0458903) q[0];
rz(0.51086673) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(2.9313472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71111403) q[0];
sx q[0];
rz(-2.1415882) q[0];
sx q[0];
rz(-2.9080703) q[0];
rz(-1.4214244) q[2];
sx q[2];
rz(-0.66022849) q[2];
sx q[2];
rz(0.82158839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8556274) q[1];
sx q[1];
rz(-1.1447765) q[1];
sx q[1];
rz(0.53526874) q[1];
rz(-pi) q[2];
rz(-0.56783592) q[3];
sx q[3];
rz(-0.68151268) q[3];
sx q[3];
rz(2.5208254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0290587) q[2];
sx q[2];
rz(-1.6363279) q[2];
sx q[2];
rz(1.0851592) q[2];
rz(0.30760136) q[3];
sx q[3];
rz(-2.8234973) q[3];
sx q[3];
rz(-0.65348452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45162421) q[0];
sx q[0];
rz(-1.987048) q[0];
sx q[0];
rz(-1.2183627) q[0];
rz(0.65144173) q[1];
sx q[1];
rz(-2.1919498) q[1];
sx q[1];
rz(0.776074) q[1];
rz(0.87985676) q[2];
sx q[2];
rz(-1.5825888) q[2];
sx q[2];
rz(1.3144944) q[2];
rz(-0.50363398) q[3];
sx q[3];
rz(-1.0009652) q[3];
sx q[3];
rz(2.9611369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
