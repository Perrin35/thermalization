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
rz(0.73135102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3325746) q[0];
sx q[0];
rz(-1.8900196) q[0];
sx q[0];
rz(-2.2126424) q[0];
rz(-pi) q[1];
rz(3.1362304) q[2];
sx q[2];
rz(-2.4796133) q[2];
sx q[2];
rz(0.84214166) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9125376) q[1];
sx q[1];
rz(-2.3542535) q[1];
sx q[1];
rz(-0.064219193) q[1];
x q[2];
rz(1.4764686) q[3];
sx q[3];
rz(-1.7324311) q[3];
sx q[3];
rz(1.1867271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88087624) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(-0.94240776) q[2];
rz(0.85855329) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(-2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079916) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-0.058636531) q[0];
rz(-0.05642852) q[1];
sx q[1];
rz(-1.3899048) q[1];
sx q[1];
rz(-1.1172392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6031044) q[0];
sx q[0];
rz(-1.9637917) q[0];
sx q[0];
rz(2.3343143) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1740721) q[2];
sx q[2];
rz(-2.8182321) q[2];
sx q[2];
rz(0.67229313) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.27858116) q[1];
sx q[1];
rz(-1.3323235) q[1];
sx q[1];
rz(0.084785297) q[1];
x q[2];
rz(-2.6659219) q[3];
sx q[3];
rz(-1.4632309) q[3];
sx q[3];
rz(1.8799262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3650018) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(2.7351725) q[2];
rz(2.8255919) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(-2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86421788) q[0];
sx q[0];
rz(-0.45806956) q[0];
sx q[0];
rz(-1.3360485) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(2.1727402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1287945) q[0];
sx q[0];
rz(-1.199882) q[0];
sx q[0];
rz(-0.65058915) q[0];
x q[1];
rz(0.79040852) q[2];
sx q[2];
rz(-1.1539173) q[2];
sx q[2];
rz(-2.0415155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67476701) q[1];
sx q[1];
rz(-1.5286501) q[1];
sx q[1];
rz(1.274363) q[1];
rz(1.3125456) q[3];
sx q[3];
rz(-2.2792993) q[3];
sx q[3];
rz(-1.7188622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1442673) q[2];
sx q[2];
rz(-0.48339016) q[2];
sx q[2];
rz(2.9664795) q[2];
rz(1.3052321) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(-1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1433732) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(-1.0293707) q[0];
rz(-0.81946212) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(0.79684657) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681172) q[0];
sx q[0];
rz(-2.0094909) q[0];
sx q[0];
rz(1.3583769) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1355924) q[2];
sx q[2];
rz(-1.0428535) q[2];
sx q[2];
rz(-2.1354298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12415826) q[1];
sx q[1];
rz(-0.34338152) q[1];
sx q[1];
rz(0.92059532) q[1];
rz(-2.6096973) q[3];
sx q[3];
rz(-0.73557702) q[3];
sx q[3];
rz(-3.0656505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0113819) q[2];
sx q[2];
rz(-3.104976) q[2];
sx q[2];
rz(1.9127362) q[2];
rz(2.7025488) q[3];
sx q[3];
rz(-1.565758) q[3];
sx q[3];
rz(-1.2307833) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586886) q[0];
sx q[0];
rz(-1.5331601) q[0];
sx q[0];
rz(0.80832344) q[0];
rz(1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(0.70560169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1418512) q[0];
sx q[0];
rz(-1.5361642) q[0];
sx q[0];
rz(0.3006199) q[0];
rz(-pi) q[1];
rz(2.2735177) q[2];
sx q[2];
rz(-1.3931615) q[2];
sx q[2];
rz(0.16791074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29632211) q[1];
sx q[1];
rz(-2.8883838) q[1];
sx q[1];
rz(1.5814387) q[1];
rz(-pi) q[2];
rz(-0.79885317) q[3];
sx q[3];
rz(-1.0980513) q[3];
sx q[3];
rz(-2.1780739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6610873) q[2];
sx q[2];
rz(-1.2837807) q[2];
sx q[2];
rz(0.60574469) q[2];
rz(0.12039603) q[3];
sx q[3];
rz(-1.8431289) q[3];
sx q[3];
rz(2.3275183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8814988) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(-0.0023728097) q[0];
rz(-0.35898769) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(-2.733309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9914382) q[0];
sx q[0];
rz(-1.5350193) q[0];
sx q[0];
rz(-1.6973901) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96377669) q[2];
sx q[2];
rz(-1.7923819) q[2];
sx q[2];
rz(-2.2100248) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.830891) q[1];
sx q[1];
rz(-2.9017555) q[1];
sx q[1];
rz(-2.8562921) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95127912) q[3];
sx q[3];
rz(-0.86182098) q[3];
sx q[3];
rz(0.10892841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6458873) q[2];
sx q[2];
rz(-2.3453823) q[2];
sx q[2];
rz(3.0118946) q[2];
rz(-0.4194704) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(-2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7523338) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(0.68354052) q[0];
rz(2.2031802) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(1.2260812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48213595) q[0];
sx q[0];
rz(-0.99471751) q[0];
sx q[0];
rz(2.0123737) q[0];
x q[1];
rz(0.60003321) q[2];
sx q[2];
rz(-1.7259611) q[2];
sx q[2];
rz(-0.70916353) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56537762) q[1];
sx q[1];
rz(-0.54311308) q[1];
sx q[1];
rz(0.99467038) q[1];
x q[2];
rz(-0.55415537) q[3];
sx q[3];
rz(-2.1423577) q[3];
sx q[3];
rz(3.1031648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.73416) q[2];
sx q[2];
rz(-1.4340883) q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68607512) q[0];
sx q[0];
rz(-0.34611836) q[0];
sx q[0];
rz(2.2177875) q[0];
rz(2.0763981) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(-1.1346029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89399983) q[0];
sx q[0];
rz(-1.1957279) q[0];
sx q[0];
rz(-3.0984237) q[0];
x q[1];
rz(0.61750268) q[2];
sx q[2];
rz(-0.98101014) q[2];
sx q[2];
rz(-1.8621572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5091619) q[1];
sx q[1];
rz(-2.1010509) q[1];
sx q[1];
rz(-3.0796618) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52319877) q[3];
sx q[3];
rz(-0.6776132) q[3];
sx q[3];
rz(2.1385156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0040466641) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(0.64296067) q[2];
rz(0.17063394) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(1.0954866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.3951025) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(2.3353031) q[0];
rz(-3.131648) q[1];
sx q[1];
rz(-0.6178304) q[1];
sx q[1];
rz(0.26201216) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0341971) q[0];
sx q[0];
rz(-2.5881564) q[0];
sx q[0];
rz(1.8897927) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9295983) q[2];
sx q[2];
rz(-0.65644453) q[2];
sx q[2];
rz(-2.580263) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9741874) q[1];
sx q[1];
rz(-1.7767118) q[1];
sx q[1];
rz(-2.0108825) q[1];
rz(1.069182) q[3];
sx q[3];
rz(-1.8610238) q[3];
sx q[3];
rz(-2.888916) q[3];
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
rz(-0.90397942) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(2.2601295) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7664465) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(0.65650702) q[0];
rz(0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(0.62969977) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5507767) q[0];
sx q[0];
rz(-1.1712946) q[0];
sx q[0];
rz(-2.4062064) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0378605) q[2];
sx q[2];
rz(-1.3383972) q[2];
sx q[2];
rz(0.17017066) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9650642) q[1];
sx q[1];
rz(-1.5429436) q[1];
sx q[1];
rz(-1.5256397) q[1];
rz(-pi) q[2];
rz(-1.0694169) q[3];
sx q[3];
rz(-2.5540941) q[3];
sx q[3];
rz(2.0794433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5188344) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(2.4614914) q[2];
rz(-2.425219) q[3];
sx q[3];
rz(-3.0383737) q[3];
sx q[3];
rz(2.0228001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685365) q[0];
sx q[0];
rz(-1.2328883) q[0];
sx q[0];
rz(0.20986025) q[0];
rz(0.14280351) q[1];
sx q[1];
rz(-1.0719943) q[1];
sx q[1];
rz(0.73678585) q[1];
rz(-0.057351107) q[2];
sx q[2];
rz(-1.3290559) q[2];
sx q[2];
rz(0.17964687) q[2];
rz(1.5132001) q[3];
sx q[3];
rz(-0.61020281) q[3];
sx q[3];
rz(-1.4635066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
