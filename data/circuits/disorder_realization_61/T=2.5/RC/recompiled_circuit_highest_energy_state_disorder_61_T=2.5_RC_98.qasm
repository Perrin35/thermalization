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
rz(4.2491663) q[0];
sx q[0];
rz(7.6585461) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(-2.4102416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99217447) q[0];
sx q[0];
rz(-2.1753924) q[0];
sx q[0];
rz(-0.3913619) q[0];
x q[1];
rz(0.0053623036) q[2];
sx q[2];
rz(-0.66197936) q[2];
sx q[2];
rz(-2.299451) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0034084) q[1];
sx q[1];
rz(-2.3560682) q[1];
sx q[1];
rz(-1.6351321) q[1];
x q[2];
rz(0.16234397) q[3];
sx q[3];
rz(-1.4777017) q[3];
sx q[3];
rz(0.39929354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.88087624) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(2.1991849) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(0.058636531) q[0];
rz(-0.05642852) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-2.0243534) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53848828) q[0];
sx q[0];
rz(-1.9637917) q[0];
sx q[0];
rz(-0.80727838) q[0];
rz(-pi) q[1];
rz(1.9675206) q[2];
sx q[2];
rz(-0.32336059) q[2];
sx q[2];
rz(0.67229313) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27858116) q[1];
sx q[1];
rz(-1.3323235) q[1];
sx q[1];
rz(-3.0568074) q[1];
rz(-pi) q[2];
rz(-1.6916709) q[3];
sx q[3];
rz(-2.0434922) q[3];
sx q[3];
rz(2.8877088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77659082) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(-2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86421788) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.3360485) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(-0.96885243) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012798158) q[0];
sx q[0];
rz(-1.199882) q[0];
sx q[0];
rz(0.65058915) q[0];
rz(-pi) q[1];
rz(2.1325705) q[2];
sx q[2];
rz(-0.86350212) q[2];
sx q[2];
rz(-2.2826441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4668256) q[1];
sx q[1];
rz(-1.6129426) q[1];
sx q[1];
rz(1.8672297) q[1];
rz(1.8290471) q[3];
sx q[3];
rz(-2.2792993) q[3];
sx q[3];
rz(1.7188622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1442673) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(0.17511314) q[2];
rz(1.8363606) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.99821943) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(1.0293707) q[0];
rz(-2.3221305) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(-0.79684657) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5433079) q[0];
sx q[0];
rz(-0.48439041) q[0];
sx q[0];
rz(-0.422307) q[0];
rz(-pi) q[1];
rz(0.0060002319) q[2];
sx q[2];
rz(-2.0987392) q[2];
sx q[2];
rz(2.1354298) q[2];
rz(-pi) q[3];
x q[3];
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
rz(-pi) q[2];
x q[2];
rz(0.53189532) q[3];
sx q[3];
rz(-2.4060156) q[3];
sx q[3];
rz(3.0656505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0113819) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(1.2288564) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(-1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586886) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(-0.80832344) q[0];
rz(-1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(2.435991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8238753) q[0];
sx q[0];
rz(-0.30254811) q[0];
sx q[0];
rz(-3.0251193) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23106261) q[2];
sx q[2];
rz(-0.88132826) q[2];
sx q[2];
rz(-1.5901515) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29632211) q[1];
sx q[1];
rz(-0.25320881) q[1];
sx q[1];
rz(1.5814387) q[1];
x q[2];
rz(2.5217149) q[3];
sx q[3];
rz(-2.2408463) q[3];
sx q[3];
rz(1.0244964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48050532) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(-2.535848) q[2];
rz(3.0211966) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(0.81407434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8814988) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(0.0023728097) q[0];
rz(-0.35898769) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(-2.733309) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4160894) q[0];
sx q[0];
rz(-1.6973087) q[0];
sx q[0];
rz(0.036065412) q[0];
rz(-pi) q[1];
rz(-1.194656) q[2];
sx q[2];
rz(-2.500211) q[2];
sx q[2];
rz(-0.94571404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.017458113) q[1];
sx q[1];
rz(-1.5038905) q[1];
sx q[1];
rz(0.23048877) q[1];
rz(-2.5465132) q[3];
sx q[3];
rz(-0.90463764) q[3];
sx q[3];
rz(0.72197589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4957054) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(-3.0118946) q[2];
rz(-2.7221223) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(-0.81370846) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7523338) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(0.68354052) q[0];
rz(-2.2031802) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(1.2260812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6594567) q[0];
sx q[0];
rz(-0.99471751) q[0];
sx q[0];
rz(-1.129219) q[0];
rz(-0.27023882) q[2];
sx q[2];
rz(-0.61737379) q[2];
sx q[2];
rz(0.63948217) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49785532) q[1];
sx q[1];
rz(-1.2853936) q[1];
sx q[1];
rz(-1.1021815) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88507401) q[3];
sx q[3];
rz(-0.7739017) q[3];
sx q[3];
rz(0.89088551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40743264) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(-2.4343991) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-2.4652822) q[3];
sx q[3];
rz(-1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68607512) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(0.92380512) q[0];
rz(-2.0763981) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(-2.0069897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4489733) q[0];
sx q[0];
rz(-1.53063) q[0];
sx q[0];
rz(1.19541) q[0];
x q[1];
rz(-2.258111) q[2];
sx q[2];
rz(-2.0728014) q[2];
sx q[2];
rz(0.66752455) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75443134) q[1];
sx q[1];
rz(-0.5335156) q[1];
sx q[1];
rz(-1.4656161) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60875684) q[3];
sx q[3];
rz(-1.8894102) q[3];
sx q[3];
rz(-0.99005885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.137546) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(-0.64296067) q[2];
rz(-0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(-1.0954866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-0.80628959) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-0.6178304) q[1];
sx q[1];
rz(-0.26201216) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4777139) q[0];
sx q[0];
rz(-2.0933525) q[0];
sx q[0];
rz(-2.9502003) q[0];
x q[1];
rz(0.21199439) q[2];
sx q[2];
rz(-0.65644453) q[2];
sx q[2];
rz(-2.580263) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6422185) q[1];
sx q[1];
rz(-2.0009577) q[1];
sx q[1];
rz(-0.22689928) q[1];
rz(-pi) q[2];
rz(1.069182) q[3];
sx q[3];
rz(-1.8610238) q[3];
sx q[3];
rz(0.25267664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(-0.90397942) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(0.8814632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.7664465) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(2.4850856) q[0];
rz(0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(0.62969977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.499959) q[0];
sx q[0];
rz(-2.2370506) q[0];
sx q[0];
rz(-2.0883661) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8734287) q[2];
sx q[2];
rz(-1.0536453) q[2];
sx q[2];
rz(-1.535653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7460664) q[1];
sx q[1];
rz(-1.6159355) q[1];
sx q[1];
rz(-3.1137115) q[1];
rz(1.0694169) q[3];
sx q[3];
rz(-0.58749858) q[3];
sx q[3];
rz(-1.0621493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(-2.4614914) q[2];
rz(-2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730561) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(2.9987891) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-1.7992146) q[2];
sx q[2];
rz(-2.8932718) q[2];
sx q[2];
rz(0.415033) q[2];
rz(1.6283926) q[3];
sx q[3];
rz(-2.5313898) q[3];
sx q[3];
rz(1.678086) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
