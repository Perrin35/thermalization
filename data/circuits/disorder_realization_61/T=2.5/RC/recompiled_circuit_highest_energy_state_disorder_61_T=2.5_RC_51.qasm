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
rz(0.50245589) q[0];
sx q[0];
rz(-2.034019) q[0];
sx q[0];
rz(-1.7662319) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(0.73135102) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99217447) q[0];
sx q[0];
rz(-0.96620027) q[0];
sx q[0];
rz(0.3913619) q[0];
x q[1];
rz(-1.574975) q[2];
sx q[2];
rz(-0.90882817) q[2];
sx q[2];
rz(-2.3062492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3870942) q[1];
sx q[1];
rz(-1.6162786) q[1];
sx q[1];
rz(-2.3552853) q[1];
rz(-1.4764686) q[3];
sx q[3];
rz(-1.7324311) q[3];
sx q[3];
rz(1.9548655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88087624) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(0.94240776) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-2.5311354) q[3];
sx q[3];
rz(-2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(3.079916) q[0];
sx q[0];
rz(-0.55773568) q[0];
sx q[0];
rz(-0.058636531) q[0];
rz(-0.05642852) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-2.0243534) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53848828) q[0];
sx q[0];
rz(-1.177801) q[0];
sx q[0];
rz(2.3343143) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8705759) q[2];
sx q[2];
rz(-1.4477056) q[2];
sx q[2];
rz(1.276615) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5176425) q[1];
sx q[1];
rz(-0.25282598) q[1];
sx q[1];
rz(-1.9060017) q[1];
rz(-pi) q[2];
rz(2.9100203) q[3];
sx q[3];
rz(-2.6548214) q[3];
sx q[3];
rz(0.51460217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77659082) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(2.7351725) q[2];
rz(-0.31600076) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(0.78280226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2773748) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.8055441) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(-0.96885243) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1287945) q[0];
sx q[0];
rz(-1.199882) q[0];
sx q[0];
rz(-0.65058915) q[0];
rz(-0.79040852) q[2];
sx q[2];
rz(-1.9876754) q[2];
sx q[2];
rz(-2.0415155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4668256) q[1];
sx q[1];
rz(-1.6129426) q[1];
sx q[1];
rz(-1.8672297) q[1];
rz(0.28964759) q[3];
sx q[3];
rz(-0.74637369) q[3];
sx q[3];
rz(-1.0371072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1442673) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(-0.17511314) q[2];
rz(1.8363606) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(2.112222) q[0];
rz(2.3221305) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(0.79684657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7356257) q[0];
sx q[0];
rz(-1.7628364) q[0];
sx q[0];
rz(0.44749949) q[0];
x q[1];
rz(0.0060002319) q[2];
sx q[2];
rz(-2.0987392) q[2];
sx q[2];
rz(-1.0061629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5863492) q[1];
sx q[1];
rz(-1.8420911) q[1];
sx q[1];
rz(-2.9284412) q[1];
rz(2.6096973) q[3];
sx q[3];
rz(-0.73557702) q[3];
sx q[3];
rz(3.0656505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13021079) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(-1.2288564) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(0.80832344) q[0];
rz(1.7841548) q[1];
sx q[1];
rz(-0.789855) q[1];
sx q[1];
rz(-2.435991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7019136) q[0];
sx q[0];
rz(-1.2703623) q[0];
sx q[0];
rz(-1.6070532) q[0];
x q[1];
rz(-2.2735177) q[2];
sx q[2];
rz(-1.7484312) q[2];
sx q[2];
rz(0.16791074) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8562634) q[1];
sx q[1];
rz(-1.3176022) q[1];
sx q[1];
rz(-0.002753792) q[1];
rz(-pi) q[2];
rz(0.6198778) q[3];
sx q[3];
rz(-0.90074632) q[3];
sx q[3];
rz(1.0244964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6610873) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(-0.60574469) q[2];
rz(-0.12039603) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(0.81407434) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26009387) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(0.0023728097) q[0];
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
rz(-1.8470386) q[0];
rz(-pi) q[1];
rz(1.9469366) q[2];
sx q[2];
rz(-0.6413817) q[2];
sx q[2];
rz(0.94571404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.830891) q[1];
sx q[1];
rz(-2.9017555) q[1];
sx q[1];
rz(-0.28530052) q[1];
rz(0.95127912) q[3];
sx q[3];
rz(-2.2797717) q[3];
sx q[3];
rz(-3.0326642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4957054) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(3.0118946) q[2];
rz(2.7221223) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(0.81370846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3892589) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(2.4580521) q[0];
rz(-0.93841249) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(-1.2260812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3049604) q[0];
sx q[0];
rz(-1.2042521) q[0];
sx q[0];
rz(2.5185597) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3834882) q[2];
sx q[2];
rz(-0.97895998) q[2];
sx q[2];
rz(2.1746152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.576215) q[1];
sx q[1];
rz(-0.54311308) q[1];
sx q[1];
rz(0.99467038) q[1];
x q[2];
rz(0.92323233) q[3];
sx q[3];
rz(-1.1123163) q[3];
sx q[3];
rz(-1.9322559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.73416) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(0.70719353) q[2];
rz(-1.6678984) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(-1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555175) q[0];
sx q[0];
rz(-0.34611836) q[0];
sx q[0];
rz(-0.92380512) q[0];
rz(-2.0763981) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(-1.1346029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69261937) q[0];
sx q[0];
rz(-1.6109626) q[0];
sx q[0];
rz(1.9461826) q[0];
x q[1];
rz(0.8834817) q[2];
sx q[2];
rz(-1.0687912) q[2];
sx q[2];
rz(-0.66752455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3871613) q[1];
sx q[1];
rz(-0.5335156) q[1];
sx q[1];
rz(-1.4656161) q[1];
rz(1.1884965) q[3];
sx q[3];
rz(-0.99671066) q[3];
sx q[3];
rz(0.36575438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.137546) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-0.64296067) q[2];
rz(-0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(2.3353031) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-2.8795805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952103) q[0];
sx q[0];
rz(-1.4052009) q[0];
sx q[0];
rz(1.0402334) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4960619) q[2];
sx q[2];
rz(-1.4420266) q[2];
sx q[2];
rz(-1.9632531) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6422185) q[1];
sx q[1];
rz(-1.1406349) q[1];
sx q[1];
rz(0.22689928) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.069182) q[3];
sx q[3];
rz(-1.2805689) q[3];
sx q[3];
rz(0.25267664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(-0.14287512) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(-2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751462) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(-0.65650702) q[0];
rz(0.87908602) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(-0.62969977) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.499959) q[0];
sx q[0];
rz(-0.90454209) q[0];
sx q[0];
rz(1.0532265) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26816396) q[2];
sx q[2];
rz(-1.0536453) q[2];
sx q[2];
rz(-1.6059396) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39552626) q[1];
sx q[1];
rz(-1.5256572) q[1];
sx q[1];
rz(-0.027881134) q[1];
rz(-pi) q[2];
rz(-1.042243) q[3];
sx q[3];
rz(-1.3011329) q[3];
sx q[3];
rz(2.2049513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(2.4614914) q[2];
rz(2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(-2.0228001) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685365) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(2.9987891) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-3.0842415) q[2];
sx q[2];
rz(-1.8125368) q[2];
sx q[2];
rz(-2.9619458) q[2];
rz(2.18022) q[3];
sx q[3];
rz(-1.6037887) q[3];
sx q[3];
rz(-2.9870839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
