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
rz(-1.3753608) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(-2.4102416) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80901805) q[0];
sx q[0];
rz(-1.2515731) q[0];
sx q[0];
rz(2.2126424) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4796203) q[2];
sx q[2];
rz(-1.5740924) q[2];
sx q[2];
rz(2.4087083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.13818422) q[1];
sx q[1];
rz(-2.3560682) q[1];
sx q[1];
rz(-1.6351321) q[1];
rz(-2.6177789) q[3];
sx q[3];
rz(-2.9546545) q[3];
sx q[3];
rz(-2.4863249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88087624) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(2.1991849) q[2];
rz(-2.2830394) q[3];
sx q[3];
rz(-2.5311354) q[3];
sx q[3];
rz(-0.71949351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.079916) q[0];
sx q[0];
rz(-0.55773568) q[0];
sx q[0];
rz(3.0829561) q[0];
rz(3.0851641) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(1.1172392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3836544) q[0];
sx q[0];
rz(-2.2636739) q[0];
sx q[0];
rz(-2.6206159) q[0];
rz(-pi) q[1];
rz(-1.9675206) q[2];
sx q[2];
rz(-0.32336059) q[2];
sx q[2];
rz(-0.67229313) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.27858116) q[1];
sx q[1];
rz(-1.8092691) q[1];
sx q[1];
rz(0.084785297) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6916709) q[3];
sx q[3];
rz(-2.0434922) q[3];
sx q[3];
rz(0.25388381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3650018) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(0.4064202) q[2];
rz(-0.31600076) q[3];
sx q[3];
rz(-1.6812811) q[3];
sx q[3];
rz(-0.78280226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2773748) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(1.3360485) q[0];
rz(0.27851963) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(0.96885243) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8271885) q[0];
sx q[0];
rz(-2.1705856) q[0];
sx q[0];
rz(1.1161854) q[0];
x q[1];
rz(-0.79040852) q[2];
sx q[2];
rz(-1.9876754) q[2];
sx q[2];
rz(1.1000772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75894657) q[1];
sx q[1];
rz(-0.29932705) q[1];
sx q[1];
rz(-1.4274197) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28964759) q[3];
sx q[3];
rz(-0.74637369) q[3];
sx q[3];
rz(1.0371072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1442673) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(0.17511314) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-2.1988018) q[3];
sx q[3];
rz(1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433732) q[0];
sx q[0];
rz(-2.0046736) q[0];
sx q[0];
rz(-1.0293707) q[0];
rz(0.81946212) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(2.3447461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5982847) q[0];
sx q[0];
rz(-0.48439041) q[0];
sx q[0];
rz(-2.7192857) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0428456) q[2];
sx q[2];
rz(-1.5656131) q[2];
sx q[2];
rz(-0.56765616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0681035) q[1];
sx q[1];
rz(-1.3655546) q[1];
sx q[1];
rz(-1.2935335) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53189532) q[3];
sx q[3];
rz(-0.73557702) q[3];
sx q[3];
rz(3.0656505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13021079) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(1.9127362) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(2.3332692) q[0];
rz(1.7841548) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(2.435991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1418512) q[0];
sx q[0];
rz(-1.6054285) q[0];
sx q[0];
rz(-0.3006199) q[0];
x q[1];
rz(-0.86807495) q[2];
sx q[2];
rz(-1.7484312) q[2];
sx q[2];
rz(-0.16791074) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8452705) q[1];
sx q[1];
rz(-2.8883838) q[1];
sx q[1];
rz(-1.5814387) q[1];
rz(-pi) q[2];
rz(0.93813809) q[3];
sx q[3];
rz(-2.2626503) q[3];
sx q[3];
rz(2.9718338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48050532) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(2.535848) q[2];
rz(0.12039603) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8814988) q[0];
sx q[0];
rz(-0.81387481) q[0];
sx q[0];
rz(3.1392198) q[0];
rz(2.782605) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(0.40828362) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9914382) q[0];
sx q[0];
rz(-1.5350193) q[0];
sx q[0];
rz(1.4442025) q[0];
rz(-pi) q[1];
rz(1.9469366) q[2];
sx q[2];
rz(-0.6413817) q[2];
sx q[2];
rz(0.94571404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3107016) q[1];
sx q[1];
rz(-2.9017555) q[1];
sx q[1];
rz(0.28530052) q[1];
rz(-pi) q[2];
rz(0.95127912) q[3];
sx q[3];
rz(-0.86182098) q[3];
sx q[3];
rz(-0.10892841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4957054) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(0.1296981) q[2];
rz(-2.7221223) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(-2.3278842) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3892589) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(-2.4580521) q[0];
rz(-2.2031802) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(-1.9155115) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3049604) q[0];
sx q[0];
rz(-1.2042521) q[0];
sx q[0];
rz(2.5185597) q[0];
rz(-pi) q[1];
rz(0.60003321) q[2];
sx q[2];
rz(-1.4156315) q[2];
sx q[2];
rz(0.70916353) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9270806) q[1];
sx q[1];
rz(-2.0190372) q[1];
sx q[1];
rz(-0.31772504) q[1];
x q[2];
rz(0.88507401) q[3];
sx q[3];
rz(-2.367691) q[3];
sx q[3];
rz(0.89088551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40743264) q[2];
sx q[2];
rz(-1.4340883) q[2];
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
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
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
rz(-1.1895475) q[1];
sx q[1];
rz(-1.1346029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69261937) q[0];
sx q[0];
rz(-1.6109626) q[0];
sx q[0];
rz(1.9461826) q[0];
x q[1];
rz(2.258111) q[2];
sx q[2];
rz(-1.0687912) q[2];
sx q[2];
rz(0.66752455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5091619) q[1];
sx q[1];
rz(-2.1010509) q[1];
sx q[1];
rz(3.0796618) q[1];
rz(-pi) q[2];
rz(0.60875684) q[3];
sx q[3];
rz(-1.2521825) q[3];
sx q[3];
rz(-0.99005885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.7377487) q[2];
sx q[2];
rz(0.64296067) q[2];
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
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3951025) q[0];
sx q[0];
rz(-0.13245067) q[0];
sx q[0];
rz(-0.80628959) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(-2.8795805) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6638787) q[0];
sx q[0];
rz(-1.0482402) q[0];
sx q[0];
rz(0.19139238) q[0];
rz(1.7315032) q[2];
sx q[2];
rz(-2.2101058) q[2];
sx q[2];
rz(0.29603816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9741874) q[1];
sx q[1];
rz(-1.7767118) q[1];
sx q[1];
rz(2.0108825) q[1];
rz(2.0724107) q[3];
sx q[3];
rz(-1.2805689) q[3];
sx q[3];
rz(-2.888916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-0.8814632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(2.3751462) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(0.62969977) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7555576) q[0];
sx q[0];
rz(-0.81869253) q[0];
sx q[0];
rz(0.56171239) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.26816396) q[2];
sx q[2];
rz(-2.0879474) q[2];
sx q[2];
rz(1.535653) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15805298) q[1];
sx q[1];
rz(-3.088542) q[1];
sx q[1];
rz(2.1237462) q[1];
rz(-pi) q[2];
rz(-0.30977003) q[3];
sx q[3];
rz(-1.0632205) q[3];
sx q[3];
rz(-2.6617756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(-2.4614914) q[2];
rz(-0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730561) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(1.8129195) q[2];
sx q[2];
rz(-1.5151146) q[2];
sx q[2];
rz(1.7641868) q[2];
rz(0.040228514) q[3];
sx q[3];
rz(-0.9617525) q[3];
sx q[3];
rz(1.7483275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
