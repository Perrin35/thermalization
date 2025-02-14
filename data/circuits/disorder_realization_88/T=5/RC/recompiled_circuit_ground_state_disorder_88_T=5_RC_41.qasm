OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4360566) q[0];
sx q[0];
rz(-0.85890618) q[0];
sx q[0];
rz(0.98281759) q[0];
rz(1.1801899) q[1];
sx q[1];
rz(-2.4153914) q[1];
sx q[1];
rz(-0.89816165) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8456062) q[0];
sx q[0];
rz(-1.4665415) q[0];
sx q[0];
rz(2.8235675) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9299081) q[2];
sx q[2];
rz(-2.0835597) q[2];
sx q[2];
rz(2.7732234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4136988) q[1];
sx q[1];
rz(-2.4174712) q[1];
sx q[1];
rz(0.65774386) q[1];
rz(-pi) q[2];
rz(-2.5407422) q[3];
sx q[3];
rz(-2.966189) q[3];
sx q[3];
rz(3.1186171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0183705) q[2];
sx q[2];
rz(-1.3192588) q[2];
sx q[2];
rz(2.3940864) q[2];
rz(-1.8788762) q[3];
sx q[3];
rz(-1.6044173) q[3];
sx q[3];
rz(0.023716299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30775192) q[0];
sx q[0];
rz(-3.096014) q[0];
sx q[0];
rz(1.0652834) q[0];
rz(0.65159687) q[1];
sx q[1];
rz(-0.35531303) q[1];
sx q[1];
rz(-1.5078872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20017524) q[0];
sx q[0];
rz(-1.8719409) q[0];
sx q[0];
rz(0.99267545) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47453158) q[2];
sx q[2];
rz(-1.2515278) q[2];
sx q[2];
rz(-2.266707) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44017556) q[1];
sx q[1];
rz(-1.9721627) q[1];
sx q[1];
rz(2.077515) q[1];
x q[2];
rz(-1.5298858) q[3];
sx q[3];
rz(-2.5299151) q[3];
sx q[3];
rz(3.1168695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6642586) q[2];
sx q[2];
rz(-1.4724255) q[2];
sx q[2];
rz(3.0490457) q[2];
rz(0.44068286) q[3];
sx q[3];
rz(-2.7350072) q[3];
sx q[3];
rz(1.1455166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4785081) q[0];
sx q[0];
rz(-0.19190754) q[0];
sx q[0];
rz(-1.9051911) q[0];
rz(-0.23794404) q[1];
sx q[1];
rz(-1.4711719) q[1];
sx q[1];
rz(2.1740289) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2500754) q[0];
sx q[0];
rz(-2.2359497) q[0];
sx q[0];
rz(-0.7789156) q[0];
x q[1];
rz(-1.2252878) q[2];
sx q[2];
rz(-1.0566525) q[2];
sx q[2];
rz(1.9793069) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5284651) q[1];
sx q[1];
rz(-1.8741762) q[1];
sx q[1];
rz(-0.19618285) q[1];
x q[2];
rz(3.0839663) q[3];
sx q[3];
rz(-1.1579683) q[3];
sx q[3];
rz(2.0183764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6588916) q[2];
sx q[2];
rz(-2.7693558) q[2];
sx q[2];
rz(-2.2230395) q[2];
rz(-0.77361584) q[3];
sx q[3];
rz(-0.88671237) q[3];
sx q[3];
rz(-0.95019597) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2660265) q[0];
sx q[0];
rz(-0.73930621) q[0];
sx q[0];
rz(0.29228041) q[0];
rz(-2.267011) q[1];
sx q[1];
rz(-2.5129109) q[1];
sx q[1];
rz(-0.94327092) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1312201) q[0];
sx q[0];
rz(-1.2772075) q[0];
sx q[0];
rz(2.2468727) q[0];
rz(-pi) q[1];
rz(2.9246632) q[2];
sx q[2];
rz(-1.1116905) q[2];
sx q[2];
rz(0.87068671) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29674655) q[1];
sx q[1];
rz(-1.2592717) q[1];
sx q[1];
rz(2.2492337) q[1];
x q[2];
rz(1.1499829) q[3];
sx q[3];
rz(-0.61987034) q[3];
sx q[3];
rz(3.0734143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5248519) q[2];
sx q[2];
rz(-1.1139694) q[2];
sx q[2];
rz(-0.88391602) q[2];
rz(2.0058477) q[3];
sx q[3];
rz(-0.92848778) q[3];
sx q[3];
rz(-0.69042027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8166872) q[0];
sx q[0];
rz(-3.0563834) q[0];
sx q[0];
rz(-2.9499522) q[0];
rz(-1.8457671) q[1];
sx q[1];
rz(-1.192966) q[1];
sx q[1];
rz(0.4506909) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2795581) q[0];
sx q[0];
rz(-0.64230317) q[0];
sx q[0];
rz(-1.4181644) q[0];
rz(-pi) q[1];
rz(-0.27235548) q[2];
sx q[2];
rz(-1.5030454) q[2];
sx q[2];
rz(2.1960047) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1289336) q[1];
sx q[1];
rz(-1.6427354) q[1];
sx q[1];
rz(2.1484445) q[1];
x q[2];
rz(-0.21416625) q[3];
sx q[3];
rz(-0.29931047) q[3];
sx q[3];
rz(-2.4714937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6777163) q[2];
sx q[2];
rz(-1.8666942) q[2];
sx q[2];
rz(2.1837168) q[2];
rz(0.052834474) q[3];
sx q[3];
rz(-2.5962679) q[3];
sx q[3];
rz(2.7132645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4868454) q[0];
sx q[0];
rz(-2.0992278) q[0];
sx q[0];
rz(0.47252193) q[0];
rz(-2.3782702) q[1];
sx q[1];
rz(-2.4210052) q[1];
sx q[1];
rz(-2.8514013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4241353) q[0];
sx q[0];
rz(-0.73834921) q[0];
sx q[0];
rz(-2.8108869) q[0];
rz(-pi) q[1];
rz(0.31748469) q[2];
sx q[2];
rz(-2.0986663) q[2];
sx q[2];
rz(1.1880529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.063263254) q[1];
sx q[1];
rz(-1.4681879) q[1];
sx q[1];
rz(0.66184499) q[1];
x q[2];
rz(-0.17866023) q[3];
sx q[3];
rz(-1.3746486) q[3];
sx q[3];
rz(2.2771108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7758238) q[2];
sx q[2];
rz(-0.86161986) q[2];
sx q[2];
rz(-0.59930581) q[2];
rz(1.722909) q[3];
sx q[3];
rz(-1.320188) q[3];
sx q[3];
rz(2.1883709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391984) q[0];
sx q[0];
rz(-1.9559487) q[0];
sx q[0];
rz(-1.3108569) q[0];
rz(1.5015548) q[1];
sx q[1];
rz(-2.7412667) q[1];
sx q[1];
rz(0.60360533) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5854322) q[0];
sx q[0];
rz(-1.7045665) q[0];
sx q[0];
rz(-1.9942392) q[0];
rz(-pi) q[1];
rz(1.7896176) q[2];
sx q[2];
rz(-1.0299152) q[2];
sx q[2];
rz(-2.4336124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.43638602) q[1];
sx q[1];
rz(-1.5028186) q[1];
sx q[1];
rz(-1.8118565) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9333618) q[3];
sx q[3];
rz(-2.3891267) q[3];
sx q[3];
rz(-0.83112874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32460585) q[2];
sx q[2];
rz(-2.2548455) q[2];
sx q[2];
rz(-0.19503197) q[2];
rz(2.4845691) q[3];
sx q[3];
rz(-1.7200836) q[3];
sx q[3];
rz(3.0278964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3312155) q[0];
sx q[0];
rz(-1.7278142) q[0];
sx q[0];
rz(0.072176607) q[0];
rz(-0.78308925) q[1];
sx q[1];
rz(-2.5091645) q[1];
sx q[1];
rz(-2.2236688) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24195053) q[0];
sx q[0];
rz(-0.52142087) q[0];
sx q[0];
rz(2.7160286) q[0];
rz(-pi) q[1];
rz(-0.36717271) q[2];
sx q[2];
rz(-0.65362602) q[2];
sx q[2];
rz(-0.12993654) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9787679) q[1];
sx q[1];
rz(-0.37436327) q[1];
sx q[1];
rz(-2.4434213) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25592978) q[3];
sx q[3];
rz(-1.157981) q[3];
sx q[3];
rz(0.61724058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5607295) q[2];
sx q[2];
rz(-1.2141289) q[2];
sx q[2];
rz(1.5011935) q[2];
rz(0.8864657) q[3];
sx q[3];
rz(-1.804988) q[3];
sx q[3];
rz(0.10543536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5378872) q[0];
sx q[0];
rz(-0.34905809) q[0];
sx q[0];
rz(-2.6339997) q[0];
rz(-0.72788584) q[1];
sx q[1];
rz(-1.4898841) q[1];
sx q[1];
rz(1.1061888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3547886) q[0];
sx q[0];
rz(-1.3193325) q[0];
sx q[0];
rz(1.2261816) q[0];
rz(1.395325) q[2];
sx q[2];
rz(-1.9242491) q[2];
sx q[2];
rz(0.33909106) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4708097) q[1];
sx q[1];
rz(-1.9341262) q[1];
sx q[1];
rz(-2.7881507) q[1];
rz(2.15739) q[3];
sx q[3];
rz(-1.1569258) q[3];
sx q[3];
rz(-2.4936287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20405208) q[2];
sx q[2];
rz(-2.0378518) q[2];
sx q[2];
rz(-0.44720116) q[2];
rz(-1.4372545) q[3];
sx q[3];
rz(-0.81366003) q[3];
sx q[3];
rz(-1.6566431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5337885) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(-2.6265889) q[0];
rz(1.1732514) q[1];
sx q[1];
rz(-1.2898022) q[1];
sx q[1];
rz(2.692093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5052348) q[0];
sx q[0];
rz(-0.32669386) q[0];
sx q[0];
rz(-1.6411843) q[0];
rz(-pi) q[1];
rz(0.88926824) q[2];
sx q[2];
rz(-0.73822108) q[2];
sx q[2];
rz(0.29730931) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.892657) q[1];
sx q[1];
rz(-0.79467862) q[1];
sx q[1];
rz(0.67130295) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2902255) q[3];
sx q[3];
rz(-1.3847305) q[3];
sx q[3];
rz(2.1894941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0939193) q[2];
sx q[2];
rz(-2.8670222) q[2];
sx q[2];
rz(0.74335113) q[2];
rz(-0.36176935) q[3];
sx q[3];
rz(-1.0310562) q[3];
sx q[3];
rz(-0.37775347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7985228) q[0];
sx q[0];
rz(-1.1953851) q[0];
sx q[0];
rz(0.60085798) q[0];
rz(2.9877904) q[1];
sx q[1];
rz(-1.3187131) q[1];
sx q[1];
rz(-2.9153894) q[1];
rz(-0.13300399) q[2];
sx q[2];
rz(-1.2659182) q[2];
sx q[2];
rz(1.5230303) q[2];
rz(-2.5220925) q[3];
sx q[3];
rz(-1.2589069) q[3];
sx q[3];
rz(-0.3862602) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
