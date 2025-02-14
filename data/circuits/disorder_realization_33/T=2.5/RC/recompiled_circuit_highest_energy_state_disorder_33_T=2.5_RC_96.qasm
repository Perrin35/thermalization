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
rz(0.1641195) q[0];
sx q[0];
rz(-2.1174705) q[0];
sx q[0];
rz(0.48409387) q[0];
rz(-0.76264277) q[1];
sx q[1];
rz(-1.5781382) q[1];
sx q[1];
rz(0.054952316) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8380175) q[0];
sx q[0];
rz(-0.9009046) q[0];
sx q[0];
rz(3.1017257) q[0];
x q[1];
rz(1.3447433) q[2];
sx q[2];
rz(-2.885427) q[2];
sx q[2];
rz(-2.0992172) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72803264) q[1];
sx q[1];
rz(-1.7031341) q[1];
sx q[1];
rz(1.669426) q[1];
x q[2];
rz(-2.8047885) q[3];
sx q[3];
rz(-0.97645491) q[3];
sx q[3];
rz(-1.6860733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.310828) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(-2.2667387) q[3];
sx q[3];
rz(-2.0503876) q[3];
sx q[3];
rz(2.8680475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65983588) q[0];
sx q[0];
rz(-0.6898703) q[0];
sx q[0];
rz(-0.39875317) q[0];
rz(-0.2440456) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(-2.0358548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92764054) q[0];
sx q[0];
rz(-2.0114779) q[0];
sx q[0];
rz(-1.8233612) q[0];
x q[1];
rz(2.2244338) q[2];
sx q[2];
rz(-2.0431402) q[2];
sx q[2];
rz(-1.2231959) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9404906) q[1];
sx q[1];
rz(-0.0041714287) q[1];
sx q[1];
rz(-0.68912403) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8375754) q[3];
sx q[3];
rz(-1.4543797) q[3];
sx q[3];
rz(2.4698225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.70046052) q[2];
sx q[2];
rz(-1.9596142) q[2];
sx q[2];
rz(-1.0750809) q[2];
rz(-2.4454146) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8657846) q[0];
sx q[0];
rz(-1.8901261) q[0];
sx q[0];
rz(-2.9248917) q[0];
rz(-1.7614583) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(2.5221141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19720896) q[0];
sx q[0];
rz(-1.5329307) q[0];
sx q[0];
rz(-3.1119431) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2080433) q[2];
sx q[2];
rz(-0.13158509) q[2];
sx q[2];
rz(-3.0818444) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8113831) q[1];
sx q[1];
rz(-2.7683677) q[1];
sx q[1];
rz(0.95287816) q[1];
rz(-2.1501174) q[3];
sx q[3];
rz(-2.1001171) q[3];
sx q[3];
rz(2.4059699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2565903) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(0.090506434) q[2];
rz(-0.45804405) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(-1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8103545) q[0];
sx q[0];
rz(-0.88307035) q[0];
sx q[0];
rz(-0.93165398) q[0];
rz(0.16904198) q[1];
sx q[1];
rz(-0.5642429) q[1];
sx q[1];
rz(-2.1021252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1202034) q[0];
sx q[0];
rz(-1.3303555) q[0];
sx q[0];
rz(-0.73657764) q[0];
rz(-pi) q[1];
rz(-2.2957689) q[2];
sx q[2];
rz(-2.8081144) q[2];
sx q[2];
rz(0.84104702) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2478777) q[1];
sx q[1];
rz(-1.3637241) q[1];
sx q[1];
rz(0.061092579) q[1];
rz(-pi) q[2];
rz(0.28431953) q[3];
sx q[3];
rz(-1.7756724) q[3];
sx q[3];
rz(-2.7903872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6197551) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(2.6452765) q[2];
rz(0.79658341) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(2.1697178) q[0];
rz(0.12174363) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(1.8121388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055453528) q[0];
sx q[0];
rz(-0.3124961) q[0];
sx q[0];
rz(-1.7715095) q[0];
rz(-pi) q[1];
rz(2.7260323) q[2];
sx q[2];
rz(-1.1141127) q[2];
sx q[2];
rz(1.6650852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.59898224) q[1];
sx q[1];
rz(-2.3282781) q[1];
sx q[1];
rz(2.5573362) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6869808) q[3];
sx q[3];
rz(-0.079616485) q[3];
sx q[3];
rz(1.5821567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10151265) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(-2.9902048) q[2];
rz(0.76465145) q[3];
sx q[3];
rz(-2.305856) q[3];
sx q[3];
rz(1.5966655) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695892) q[0];
sx q[0];
rz(-1.66865) q[0];
sx q[0];
rz(-1.0714916) q[0];
rz(0.53120652) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(-2.5154617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3817978) q[0];
sx q[0];
rz(-0.017649895) q[0];
sx q[0];
rz(1.2977029) q[0];
x q[1];
rz(-2.5444736) q[2];
sx q[2];
rz(-1.7749987) q[2];
sx q[2];
rz(1.1115505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5007361) q[1];
sx q[1];
rz(-1.8759512) q[1];
sx q[1];
rz(2.5049097) q[1];
x q[2];
rz(-2.0141861) q[3];
sx q[3];
rz(-1.1133685) q[3];
sx q[3];
rz(-1.2113038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7834187) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(-2.1273071) q[2];
rz(-1.8861534) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(1.0534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805873) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(2.2210806) q[0];
rz(0.11407425) q[1];
sx q[1];
rz(-1.736085) q[1];
sx q[1];
rz(1.1309518) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3944954) q[0];
sx q[0];
rz(-0.97309006) q[0];
sx q[0];
rz(1.5639485) q[0];
x q[1];
rz(-0.80318309) q[2];
sx q[2];
rz(-2.4121662) q[2];
sx q[2];
rz(-2.8620697) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9662094) q[1];
sx q[1];
rz(-1.7031324) q[1];
sx q[1];
rz(2.2369409) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.00131) q[3];
sx q[3];
rz(-0.52061235) q[3];
sx q[3];
rz(-1.8050107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-0.64212126) q[2];
sx q[2];
rz(2.7327909) q[2];
rz(0.18320228) q[3];
sx q[3];
rz(-1.5902404) q[3];
sx q[3];
rz(-2.0028152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278397) q[0];
sx q[0];
rz(-0.04700679) q[0];
sx q[0];
rz(-0.29079944) q[0];
rz(0.94789061) q[1];
sx q[1];
rz(-0.46698505) q[1];
sx q[1];
rz(1.9416169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24058293) q[0];
sx q[0];
rz(-1.9685317) q[0];
sx q[0];
rz(-0.57698864) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4418996) q[2];
sx q[2];
rz(-2.202575) q[2];
sx q[2];
rz(-0.80365411) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.34741286) q[1];
sx q[1];
rz(-2.844226) q[1];
sx q[1];
rz(1.7296687) q[1];
x q[2];
rz(0.10082106) q[3];
sx q[3];
rz(-1.3885048) q[3];
sx q[3];
rz(1.3829447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36134186) q[2];
sx q[2];
rz(-1.6929408) q[2];
sx q[2];
rz(-1.348749) q[2];
rz(1.6382943) q[3];
sx q[3];
rz(-0.88204757) q[3];
sx q[3];
rz(-0.1851113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(1.1143484) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(1.0741023) q[0];
rz(2.7117924) q[1];
sx q[1];
rz(-1.0440412) q[1];
sx q[1];
rz(0.46357402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060424711) q[0];
sx q[0];
rz(-1.3827818) q[0];
sx q[0];
rz(-3.0449165) q[0];
rz(-pi) q[1];
rz(-2.8797382) q[2];
sx q[2];
rz(-1.7657649) q[2];
sx q[2];
rz(-1.5443813) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.025109865) q[1];
sx q[1];
rz(-1.1949202) q[1];
sx q[1];
rz(2.2459338) q[1];
rz(-2.5829599) q[3];
sx q[3];
rz(-1.6470634) q[3];
sx q[3];
rz(-2.7879451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.27434906) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(-2.3331433) q[2];
rz(-0.48458734) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(-2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3633858) q[0];
sx q[0];
rz(-2.6860542) q[0];
sx q[0];
rz(-1.1325915) q[0];
rz(3.1032108) q[1];
sx q[1];
rz(-1.3382341) q[1];
sx q[1];
rz(-2.8448232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7420121) q[0];
sx q[0];
rz(-2.098408) q[0];
sx q[0];
rz(0.3216089) q[0];
rz(1.1266668) q[2];
sx q[2];
rz(-1.9812366) q[2];
sx q[2];
rz(-0.12516147) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1047639) q[1];
sx q[1];
rz(-1.6460168) q[1];
sx q[1];
rz(-2.9192218) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88480437) q[3];
sx q[3];
rz(-1.2594885) q[3];
sx q[3];
rz(-1.5112359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9639637) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(0.56619823) q[2];
rz(-2.0171793) q[3];
sx q[3];
rz(-0.12523139) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2531256) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(0.26168564) q[1];
sx q[1];
rz(-1.2005922) q[1];
sx q[1];
rz(-2.190879) q[1];
rz(-0.28957257) q[2];
sx q[2];
rz(-3.0560383) q[2];
sx q[2];
rz(-0.4734931) q[2];
rz(0.64601267) q[3];
sx q[3];
rz(-2.1727242) q[3];
sx q[3];
rz(0.77632191) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
