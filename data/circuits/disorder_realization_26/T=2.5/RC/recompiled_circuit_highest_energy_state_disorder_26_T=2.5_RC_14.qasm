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
rz(0.8910203) q[0];
sx q[0];
rz(-1.2863337) q[0];
sx q[0];
rz(0.88598716) q[0];
rz(-2.6180144) q[1];
sx q[1];
rz(-0.55966592) q[1];
sx q[1];
rz(1.3203415) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.490075) q[0];
sx q[0];
rz(-1.0371672) q[0];
sx q[0];
rz(-0.20488157) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2314531) q[2];
sx q[2];
rz(-2.3187713) q[2];
sx q[2];
rz(1.0642327) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1883386) q[1];
sx q[1];
rz(-1.0085114) q[1];
sx q[1];
rz(-0.072093318) q[1];
x q[2];
rz(2.3137847) q[3];
sx q[3];
rz(-0.95658619) q[3];
sx q[3];
rz(0.18060623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3628799) q[2];
sx q[2];
rz(-1.820463) q[2];
sx q[2];
rz(-2.2134181) q[2];
rz(1.5556395) q[3];
sx q[3];
rz(-2.0944244) q[3];
sx q[3];
rz(1.9716523) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8856119) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(-2.0804491) q[0];
rz(1.9288918) q[1];
sx q[1];
rz(-0.78293982) q[1];
sx q[1];
rz(-1.8276259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43976682) q[0];
sx q[0];
rz(-1.9226388) q[0];
sx q[0];
rz(-0.74811305) q[0];
rz(-2.1232067) q[2];
sx q[2];
rz(-2.4168042) q[2];
sx q[2];
rz(-1.9375436) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3281341) q[1];
sx q[1];
rz(-2.5812006) q[1];
sx q[1];
rz(1.2499534) q[1];
rz(-pi) q[2];
rz(0.98269083) q[3];
sx q[3];
rz(-0.81702166) q[3];
sx q[3];
rz(2.6213561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5709915) q[2];
sx q[2];
rz(-2.1676895) q[2];
sx q[2];
rz(-2.2368597) q[2];
rz(1.5500801) q[3];
sx q[3];
rz(-1.855987) q[3];
sx q[3];
rz(-1.8486283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631113) q[0];
sx q[0];
rz(-0.08992973) q[0];
sx q[0];
rz(-0.91823804) q[0];
rz(-0.69084424) q[1];
sx q[1];
rz(-1.3965239) q[1];
sx q[1];
rz(0.7965368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7920216) q[0];
sx q[0];
rz(-0.22788985) q[0];
sx q[0];
rz(1.9877276) q[0];
rz(-pi) q[1];
rz(-1.5112707) q[2];
sx q[2];
rz(-1.5280485) q[2];
sx q[2];
rz(1.8190365) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5480876) q[1];
sx q[1];
rz(-1.4626164) q[1];
sx q[1];
rz(-0.87203474) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91751601) q[3];
sx q[3];
rz(-1.2923354) q[3];
sx q[3];
rz(2.3115932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.010926509) q[2];
sx q[2];
rz(-2.1911759) q[2];
sx q[2];
rz(1.8939023) q[2];
rz(2.583336) q[3];
sx q[3];
rz(-1.4562675) q[3];
sx q[3];
rz(3.1202417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70304865) q[0];
sx q[0];
rz(-3.092364) q[0];
sx q[0];
rz(2.4141648) q[0];
rz(0.1618596) q[1];
sx q[1];
rz(-1.1715803) q[1];
sx q[1];
rz(-1.1579827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25140554) q[0];
sx q[0];
rz(-1.5047538) q[0];
sx q[0];
rz(1.6679881) q[0];
x q[1];
rz(-1.4314992) q[2];
sx q[2];
rz(-1.5698395) q[2];
sx q[2];
rz(-0.19569163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6111485) q[1];
sx q[1];
rz(-1.4475249) q[1];
sx q[1];
rz(-1.7632496) q[1];
x q[2];
rz(2.6154989) q[3];
sx q[3];
rz(-1.4988314) q[3];
sx q[3];
rz(1.1854894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8689279) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(-1.7499917) q[2];
rz(0.23022716) q[3];
sx q[3];
rz(-1.9672491) q[3];
sx q[3];
rz(-2.9496884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4153083) q[0];
sx q[0];
rz(-3.0400161) q[0];
sx q[0];
rz(-2.0368077) q[0];
rz(3.0305908) q[1];
sx q[1];
rz(-2.0260725) q[1];
sx q[1];
rz(-1.312779) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60393023) q[0];
sx q[0];
rz(-0.27034187) q[0];
sx q[0];
rz(-2.2946847) q[0];
rz(-pi) q[1];
rz(0.13938015) q[2];
sx q[2];
rz(-2.5733893) q[2];
sx q[2];
rz(-2.2903493) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1066619) q[1];
sx q[1];
rz(-2.6900869) q[1];
sx q[1];
rz(-1.0978218) q[1];
rz(-pi) q[2];
rz(0.41683414) q[3];
sx q[3];
rz(-2.2456944) q[3];
sx q[3];
rz(-1.0639497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.50292) q[2];
sx q[2];
rz(-2.0827677) q[2];
sx q[2];
rz(-1.0271094) q[2];
rz(2.1549759) q[3];
sx q[3];
rz(-2.8592181) q[3];
sx q[3];
rz(1.7178887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8091549) q[0];
sx q[0];
rz(-1.9092535) q[0];
sx q[0];
rz(-0.80068457) q[0];
rz(0.77146161) q[1];
sx q[1];
rz(-2.0636676) q[1];
sx q[1];
rz(1.7652184) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05089914) q[0];
sx q[0];
rz(-2.2138023) q[0];
sx q[0];
rz(1.0674632) q[0];
rz(-pi) q[1];
rz(-1.6402836) q[2];
sx q[2];
rz(-0.46705267) q[2];
sx q[2];
rz(-2.3882773) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4478057) q[1];
sx q[1];
rz(-0.70560938) q[1];
sx q[1];
rz(3.1211057) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82628754) q[3];
sx q[3];
rz(-1.0184231) q[3];
sx q[3];
rz(-2.1635522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8419522) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(3.1237777) q[2];
rz(0.31339112) q[3];
sx q[3];
rz(-2.119795) q[3];
sx q[3];
rz(-0.71948403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7844836) q[0];
sx q[0];
rz(-1.5672368) q[0];
sx q[0];
rz(-2.9743279) q[0];
rz(0.74370614) q[1];
sx q[1];
rz(-2.2552762) q[1];
sx q[1];
rz(-0.13626616) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1377099) q[0];
sx q[0];
rz(-1.2943177) q[0];
sx q[0];
rz(-2.2599758) q[0];
rz(0.17301128) q[2];
sx q[2];
rz(-1.8699416) q[2];
sx q[2];
rz(1.1227705) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72754656) q[1];
sx q[1];
rz(-2.7767477) q[1];
sx q[1];
rz(-2.9529497) q[1];
rz(-1.5149917) q[3];
sx q[3];
rz(-1.2254834) q[3];
sx q[3];
rz(-1.2962411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83884376) q[2];
sx q[2];
rz(-0.15041298) q[2];
sx q[2];
rz(1.0813084) q[2];
rz(0.14351621) q[3];
sx q[3];
rz(-1.84294) q[3];
sx q[3];
rz(1.4474086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3926587) q[0];
sx q[0];
rz(-1.8086139) q[0];
sx q[0];
rz(-0.64312154) q[0];
rz(-0.92203036) q[1];
sx q[1];
rz(-0.26607251) q[1];
sx q[1];
rz(-1.6304852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47925835) q[0];
sx q[0];
rz(-1.3915724) q[0];
sx q[0];
rz(-1.916612) q[0];
rz(-pi) q[1];
rz(1.950932) q[2];
sx q[2];
rz(-2.8415749) q[2];
sx q[2];
rz(2.0643108) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.585941) q[1];
sx q[1];
rz(-1.5787447) q[1];
sx q[1];
rz(-1.7674602) q[1];
rz(-pi) q[2];
rz(1.891252) q[3];
sx q[3];
rz(-2.0600187) q[3];
sx q[3];
rz(1.7699522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.94613218) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(-1.8411609) q[2];
rz(-0.91222936) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(-0.31360489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4881956) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(0.61757863) q[0];
rz(3.0600582) q[1];
sx q[1];
rz(-0.50059861) q[1];
sx q[1];
rz(-0.68731442) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.222037) q[0];
sx q[0];
rz(-1.8074146) q[0];
sx q[0];
rz(0.048892269) q[0];
x q[1];
rz(1.6825283) q[2];
sx q[2];
rz(-1.8364753) q[2];
sx q[2];
rz(-0.39171644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3442823) q[1];
sx q[1];
rz(-0.51206368) q[1];
sx q[1];
rz(-2.5355829) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6507023) q[3];
sx q[3];
rz(-0.54466313) q[3];
sx q[3];
rz(-2.2063257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11757892) q[2];
sx q[2];
rz(-1.964183) q[2];
sx q[2];
rz(-0.071361072) q[2];
rz(-1.9060382) q[3];
sx q[3];
rz(-2.8046298) q[3];
sx q[3];
rz(2.5753042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5679034) q[0];
sx q[0];
rz(-0.25147831) q[0];
sx q[0];
rz(0.13791826) q[0];
rz(0.57016405) q[1];
sx q[1];
rz(-1.0823931) q[1];
sx q[1];
rz(-3.1261442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041067657) q[0];
sx q[0];
rz(-0.78162949) q[0];
sx q[0];
rz(2.0360721) q[0];
x q[1];
rz(-0.49039109) q[2];
sx q[2];
rz(-2.3182959) q[2];
sx q[2];
rz(-3.1237912) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8739104) q[1];
sx q[1];
rz(-1.4060258) q[1];
sx q[1];
rz(2.6450637) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.055850765) q[3];
sx q[3];
rz(-1.7540364) q[3];
sx q[3];
rz(2.2574772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7323759) q[2];
sx q[2];
rz(-2.3541383) q[2];
sx q[2];
rz(2.709205) q[2];
rz(-2.2435097) q[3];
sx q[3];
rz(-2.2663074) q[3];
sx q[3];
rz(1.5842277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6569923) q[0];
sx q[0];
rz(-2.0130172) q[0];
sx q[0];
rz(0.76242557) q[0];
rz(2.5905329) q[1];
sx q[1];
rz(-1.3403799) q[1];
sx q[1];
rz(-0.47245477) q[1];
rz(-0.36334262) q[2];
sx q[2];
rz(-1.6132334) q[2];
sx q[2];
rz(1.2237534) q[2];
rz(-2.6656125) q[3];
sx q[3];
rz(-2.130571) q[3];
sx q[3];
rz(-0.68031753) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
