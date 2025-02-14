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
rz(1.2122756) q[0];
sx q[0];
rz(-2.0097998) q[0];
sx q[0];
rz(1.3728859) q[0];
rz(-1.1276487) q[1];
sx q[1];
rz(-1.375066) q[1];
sx q[1];
rz(0.78627237) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72981156) q[0];
sx q[0];
rz(-1.3656817) q[0];
sx q[0];
rz(-1.1954444) q[0];
rz(0.55297466) q[2];
sx q[2];
rz(-2.1510501) q[2];
sx q[2];
rz(1.0667104) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7934301) q[1];
sx q[1];
rz(-1.600666) q[1];
sx q[1];
rz(-1.8601599) q[1];
rz(-pi) q[2];
rz(2.1350098) q[3];
sx q[3];
rz(-2.461883) q[3];
sx q[3];
rz(0.93384472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0218411) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(1.8915668) q[2];
rz(2.405808) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(-1.638394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4982872) q[0];
sx q[0];
rz(-1.9685638) q[0];
sx q[0];
rz(0.071320891) q[0];
rz(-1.1530863) q[1];
sx q[1];
rz(-2.3801897) q[1];
sx q[1];
rz(-0.52871314) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3708955) q[0];
sx q[0];
rz(-0.57790725) q[0];
sx q[0];
rz(-2.175287) q[0];
rz(0.18888338) q[2];
sx q[2];
rz(-1.3391277) q[2];
sx q[2];
rz(-0.29555368) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9829927) q[1];
sx q[1];
rz(-1.0550749) q[1];
sx q[1];
rz(2.2186568) q[1];
rz(1.9857446) q[3];
sx q[3];
rz(-1.7768716) q[3];
sx q[3];
rz(-1.2690488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4802287) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(-3.1129692) q[2];
rz(0.35401595) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(-0.55907512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.70327586) q[0];
sx q[0];
rz(-1.0030614) q[0];
sx q[0];
rz(2.2502374) q[0];
rz(-2.640653) q[1];
sx q[1];
rz(-1.5653862) q[1];
sx q[1];
rz(3.0827789) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5767305) q[0];
sx q[0];
rz(-1.4848827) q[0];
sx q[0];
rz(-1.6187173) q[0];
rz(-pi) q[1];
rz(0.28532736) q[2];
sx q[2];
rz(-2.0801736) q[2];
sx q[2];
rz(-2.2338108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0173024) q[1];
sx q[1];
rz(-2.2180386) q[1];
sx q[1];
rz(-2.1976297) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5427049) q[3];
sx q[3];
rz(-1.3904316) q[3];
sx q[3];
rz(1.2171868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(2.9050262) q[2];
rz(-2.2208354) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(-1.7304272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22223602) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(-2.2663569) q[0];
rz(1.596176) q[1];
sx q[1];
rz(-0.86219209) q[1];
sx q[1];
rz(3.0962871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047642853) q[0];
sx q[0];
rz(-1.6730621) q[0];
sx q[0];
rz(-0.065667466) q[0];
x q[1];
rz(0.50778054) q[2];
sx q[2];
rz(-2.1740401) q[2];
sx q[2];
rz(0.3321307) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2341174) q[1];
sx q[1];
rz(-1.6576515) q[1];
sx q[1];
rz(0.90076561) q[1];
x q[2];
rz(1.281146) q[3];
sx q[3];
rz(-0.60263205) q[3];
sx q[3];
rz(1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28216013) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(-1.947594) q[2];
rz(-2.2512839) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(-2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10876656) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(-2.4591675) q[0];
rz(-1.8531063) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(-1.8103745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2829943) q[0];
sx q[0];
rz(-2.4147075) q[0];
sx q[0];
rz(-1.5461393) q[0];
rz(-pi) q[1];
rz(2.8814828) q[2];
sx q[2];
rz(-1.1887728) q[2];
sx q[2];
rz(-2.0224689) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4466616) q[1];
sx q[1];
rz(-0.38685054) q[1];
sx q[1];
rz(0.37987654) q[1];
rz(-2.0985475) q[3];
sx q[3];
rz(-1.1748136) q[3];
sx q[3];
rz(2.2137303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0266626) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(0.66693532) q[2];
rz(0.55495787) q[3];
sx q[3];
rz(-2.4542377) q[3];
sx q[3];
rz(-0.2989029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894593) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(2.4378648) q[0];
rz(-0.16920371) q[1];
sx q[1];
rz(-1.5237944) q[1];
sx q[1];
rz(0.15883787) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7089091) q[0];
sx q[0];
rz(-1.7129717) q[0];
sx q[0];
rz(0.12169038) q[0];
rz(-pi) q[1];
rz(-0.67262291) q[2];
sx q[2];
rz(-2.9458698) q[2];
sx q[2];
rz(1.7253523) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4602214) q[1];
sx q[1];
rz(-1.4386144) q[1];
sx q[1];
rz(-0.48355196) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30384003) q[3];
sx q[3];
rz(-2.4812077) q[3];
sx q[3];
rz(0.25750289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.77330294) q[2];
sx q[2];
rz(-2.1998019) q[2];
sx q[2];
rz(-2.2840195) q[2];
rz(1.9722021) q[3];
sx q[3];
rz(-1.2811067) q[3];
sx q[3];
rz(2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17539772) q[0];
sx q[0];
rz(-2.3747787) q[0];
sx q[0];
rz(-2.4229557) q[0];
rz(2.8596558) q[1];
sx q[1];
rz(-1.3749296) q[1];
sx q[1];
rz(0.66863543) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2235093) q[0];
sx q[0];
rz(-1.232172) q[0];
sx q[0];
rz(-2.2302367) q[0];
rz(-2.0480506) q[2];
sx q[2];
rz(-0.68533939) q[2];
sx q[2];
rz(2.6743741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9033176) q[1];
sx q[1];
rz(-1.2438477) q[1];
sx q[1];
rz(-1.1823149) q[1];
rz(-1.5626146) q[3];
sx q[3];
rz(-2.2020209) q[3];
sx q[3];
rz(2.5142728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64688993) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(2.1051712) q[2];
rz(2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(-0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(2.4155937) q[0];
rz(1.9793824) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(-1.9451709) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1370173) q[0];
sx q[0];
rz(-1.5423623) q[0];
sx q[0];
rz(-1.397055) q[0];
x q[1];
rz(-0.047646626) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(0.36587151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5550674) q[1];
sx q[1];
rz(-1.0690926) q[1];
sx q[1];
rz(-2.1074159) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67025028) q[3];
sx q[3];
rz(-1.847451) q[3];
sx q[3];
rz(1.2069699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90765816) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(-2.2080803) q[2];
rz(2.524611) q[3];
sx q[3];
rz(-1.0022997) q[3];
sx q[3];
rz(-2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(0.053255178) q[0];
rz(2.8540197) q[1];
sx q[1];
rz(-2.2122999) q[1];
sx q[1];
rz(0.25921777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142562) q[0];
sx q[0];
rz(-1.5467318) q[0];
sx q[0];
rz(-1.5177478) q[0];
rz(-pi) q[1];
rz(-1.0561814) q[2];
sx q[2];
rz(-2.8764203) q[2];
sx q[2];
rz(-2.3725703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1500435) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(1.4808473) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9837228) q[3];
sx q[3];
rz(-2.3693858) q[3];
sx q[3];
rz(-1.7534353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5549434) q[2];
sx q[2];
rz(-1.8879075) q[2];
sx q[2];
rz(-0.1813691) q[2];
rz(-2.7465316) q[3];
sx q[3];
rz(-0.22160465) q[3];
sx q[3];
rz(-2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.3996537) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(-2.7594866) q[0];
rz(-0.29640472) q[1];
sx q[1];
rz(-1.0045241) q[1];
sx q[1];
rz(-1.872725) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4110376) q[0];
sx q[0];
rz(-2.0417622) q[0];
sx q[0];
rz(0.48226003) q[0];
rz(-pi) q[1];
rz(1.4051132) q[2];
sx q[2];
rz(-0.68000092) q[2];
sx q[2];
rz(0.14762893) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8812534) q[1];
sx q[1];
rz(-1.2311544) q[1];
sx q[1];
rz(-0.022625523) q[1];
rz(0.45370729) q[3];
sx q[3];
rz(-1.6804916) q[3];
sx q[3];
rz(-3.0806993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2962013) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(-0.38983795) q[2];
rz(-1.8440394) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(-1.3564431) q[1];
sx q[1];
rz(-1.0379797) q[1];
sx q[1];
rz(-1.5300068) q[1];
rz(-0.46089725) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(2.8036694) q[3];
sx q[3];
rz(-0.81339627) q[3];
sx q[3];
rz(0.42311121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
