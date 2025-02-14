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
rz(-2.6391368) q[0];
sx q[0];
rz(-1.1075736) q[0];
sx q[0];
rz(1.7662319) q[0];
rz(-3.4604685) q[1];
sx q[1];
rz(4.7825216) q[1];
sx q[1];
rz(8.6934269) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36424822) q[0];
sx q[0];
rz(-2.4349182) q[0];
sx q[0];
rz(-1.0663435) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1362304) q[2];
sx q[2];
rz(-2.4796133) q[2];
sx q[2];
rz(2.299451) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.22905502) q[1];
sx q[1];
rz(-2.3542535) q[1];
sx q[1];
rz(-3.0773735) q[1];
x q[2];
rz(-2.9792487) q[3];
sx q[3];
rz(-1.4777017) q[3];
sx q[3];
rz(-2.7422991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.88087624) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-0.94240776) q[2];
rz(-0.85855329) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-0.058636531) q[0];
rz(0.05642852) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(1.1172392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7579382) q[0];
sx q[0];
rz(-0.87791872) q[0];
sx q[0];
rz(2.6206159) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2710167) q[2];
sx q[2];
rz(-1.6938871) q[2];
sx q[2];
rz(1.276615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2721418) q[1];
sx q[1];
rz(-1.488416) q[1];
sx q[1];
rz(-1.8100965) q[1];
rz(-1.4499217) q[3];
sx q[3];
rz(-1.0981005) q[3];
sx q[3];
rz(2.8877088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(2.7351725) q[2];
rz(0.31600076) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86421788) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.3360485) q[0];
rz(-0.27851963) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(2.1727402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1135546) q[0];
sx q[0];
rz(-2.4062706) q[0];
sx q[0];
rz(2.5707629) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0090221) q[2];
sx q[2];
rz(-0.86350212) q[2];
sx q[2];
rz(-0.85894859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75894657) q[1];
sx q[1];
rz(-0.29932705) q[1];
sx q[1];
rz(-1.7141729) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3125456) q[3];
sx q[3];
rz(-2.2792993) q[3];
sx q[3];
rz(1.4227305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1442673) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(-2.9664795) q[2];
rz(1.3052321) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.85873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(-2.112222) q[0];
rz(-0.81946212) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(-2.3447461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681172) q[0];
sx q[0];
rz(-1.1321018) q[0];
sx q[0];
rz(1.3583769) q[0];
x q[1];
rz(-1.0428456) q[2];
sx q[2];
rz(-1.5656131) q[2];
sx q[2];
rz(-2.5739365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.12415826) q[1];
sx q[1];
rz(-2.7982111) q[1];
sx q[1];
rz(-2.2209973) q[1];
rz(0.66241499) q[3];
sx q[3];
rz(-1.9180505) q[3];
sx q[3];
rz(-1.2353237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0113819) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(-1.9127362) q[2];
rz(0.43904385) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48290408) q[0];
sx q[0];
rz(-1.5331601) q[0];
sx q[0];
rz(-2.3332692) q[0];
rz(-1.7841548) q[1];
sx q[1];
rz(-0.789855) q[1];
sx q[1];
rz(-0.70560169) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8238753) q[0];
sx q[0];
rz(-0.30254811) q[0];
sx q[0];
rz(-0.11647336) q[0];
rz(2.91053) q[2];
sx q[2];
rz(-2.2602644) q[2];
sx q[2];
rz(-1.5901515) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29632211) q[1];
sx q[1];
rz(-2.8883838) q[1];
sx q[1];
rz(-1.5601539) q[1];
rz(-pi) q[2];
rz(-2.3427395) q[3];
sx q[3];
rz(-1.0980513) q[3];
sx q[3];
rz(2.1780739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48050532) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(-2.535848) q[2];
rz(3.0211966) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(-2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26009387) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(-3.1392198) q[0];
rz(-0.35898769) q[1];
sx q[1];
rz(-1.2009883) q[1];
sx q[1];
rz(2.733309) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6946163) q[0];
sx q[0];
rz(-0.13152619) q[0];
sx q[0];
rz(1.294554) q[0];
rz(-pi) q[1];
rz(2.177816) q[2];
sx q[2];
rz(-1.3492108) q[2];
sx q[2];
rz(2.2100248) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1241345) q[1];
sx q[1];
rz(-1.5038905) q[1];
sx q[1];
rz(0.23048877) q[1];
rz(-pi) q[2];
rz(-2.5465132) q[3];
sx q[3];
rz(-2.236955) q[3];
sx q[3];
rz(2.4196168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6458873) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(-3.0118946) q[2];
rz(2.7221223) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(-0.81370846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3892589) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(-2.4580521) q[0];
rz(-0.93841249) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(1.9155115) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968564) q[0];
sx q[0];
rz(-2.4312907) q[0];
sx q[0];
rz(-0.5819178) q[0];
x q[1];
rz(-2.8713538) q[2];
sx q[2];
rz(-0.61737379) q[2];
sx q[2];
rz(2.5021105) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.214512) q[1];
sx q[1];
rz(-2.0190372) q[1];
sx q[1];
rz(-0.31772504) q[1];
rz(-2.5874373) q[3];
sx q[3];
rz(-2.1423577) q[3];
sx q[3];
rz(0.038427834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.73416) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(0.70719353) q[2];
rz(-1.4736942) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68607512) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(2.2177875) q[0];
rz(-1.0651945) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(1.1346029) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89399983) q[0];
sx q[0];
rz(-1.1957279) q[0];
sx q[0];
rz(3.0984237) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.258111) q[2];
sx q[2];
rz(-2.0728014) q[2];
sx q[2];
rz(-2.4740681) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2345786) q[1];
sx q[1];
rz(-1.5173787) q[1];
sx q[1];
rz(-1.0397041) q[1];
rz(-pi) q[2];
rz(0.52319877) q[3];
sx q[3];
rz(-0.6776132) q[3];
sx q[3];
rz(-2.1385156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.0040466641) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-2.498632) q[2];
rz(-0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.74649015) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-0.80628959) q[0];
rz(0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-0.26201216) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1073955) q[0];
sx q[0];
rz(-0.55343628) q[0];
sx q[0];
rz(1.8897927) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4100894) q[2];
sx q[2];
rz(-2.2101058) q[2];
sx q[2];
rz(-2.8455545) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6422185) q[1];
sx q[1];
rz(-2.0009577) q[1];
sx q[1];
rz(-0.22689928) q[1];
rz(-pi) q[2];
rz(2.8132961) q[3];
sx q[3];
rz(-1.0919763) q[3];
sx q[3];
rz(1.6678068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7664465) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(2.2625066) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(-2.5118929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64163369) q[0];
sx q[0];
rz(-2.2370506) q[0];
sx q[0];
rz(-1.0532265) q[0];
rz(-pi) q[1];
rz(2.1037322) q[2];
sx q[2];
rz(-1.3383972) q[2];
sx q[2];
rz(0.17017066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15805298) q[1];
sx q[1];
rz(-1/(6*pi)) q[1];
sx q[1];
rz(-2.1237462) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0694169) q[3];
sx q[3];
rz(-2.5540941) q[3];
sx q[3];
rz(1.0621493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(-0.68010124) q[2];
rz(-0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730561) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(0.14280351) q[1];
sx q[1];
rz(-1.0719943) q[1];
sx q[1];
rz(0.73678585) q[1];
rz(0.057351107) q[2];
sx q[2];
rz(-1.8125368) q[2];
sx q[2];
rz(-2.9619458) q[2];
rz(-1.5132001) q[3];
sx q[3];
rz(-2.5313898) q[3];
sx q[3];
rz(1.678086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
