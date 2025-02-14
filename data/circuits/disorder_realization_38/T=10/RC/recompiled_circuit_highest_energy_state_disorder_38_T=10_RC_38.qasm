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
rz(-0.28873697) q[0];
sx q[0];
rz(-2.4462235) q[0];
sx q[0];
rz(-0.26779548) q[0];
rz(-2.7195622) q[1];
sx q[1];
rz(-0.92075092) q[1];
sx q[1];
rz(-1.2738127) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.081163) q[0];
sx q[0];
rz(-1.3295242) q[0];
sx q[0];
rz(0.046242899) q[0];
rz(-pi) q[1];
rz(0.8795514) q[2];
sx q[2];
rz(-2.3143907) q[2];
sx q[2];
rz(-2.4776221) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63034025) q[1];
sx q[1];
rz(-0.66424886) q[1];
sx q[1];
rz(-3.0306007) q[1];
rz(-pi) q[2];
rz(1.7142322) q[3];
sx q[3];
rz(-1.6063476) q[3];
sx q[3];
rz(1.2811023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6196809) q[2];
sx q[2];
rz(-0.64138594) q[2];
sx q[2];
rz(-0.042595159) q[2];
rz(2.8604782) q[3];
sx q[3];
rz(-1.576985) q[3];
sx q[3];
rz(0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5113145) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(-1.0761155) q[0];
rz(2.1108744) q[1];
sx q[1];
rz(-1.1117671) q[1];
sx q[1];
rz(1.3105185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.020582) q[0];
sx q[0];
rz(-0.91269976) q[0];
sx q[0];
rz(3.115244) q[0];
rz(2.7380472) q[2];
sx q[2];
rz(-2.2650044) q[2];
sx q[2];
rz(-0.36118868) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4736997) q[1];
sx q[1];
rz(-1.4034499) q[1];
sx q[1];
rz(-1.4506063) q[1];
rz(-pi) q[2];
rz(2.5111273) q[3];
sx q[3];
rz(-2.5098636) q[3];
sx q[3];
rz(0.28030685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.217546) q[2];
sx q[2];
rz(-2.3895538) q[2];
sx q[2];
rz(0.95620608) q[2];
rz(1.8337967) q[3];
sx q[3];
rz(-0.81495133) q[3];
sx q[3];
rz(0.1213049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2995375) q[0];
sx q[0];
rz(-0.50899035) q[0];
sx q[0];
rz(-0.036238413) q[0];
rz(2.3402975) q[1];
sx q[1];
rz(-1.5943297) q[1];
sx q[1];
rz(1.3005728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1195923) q[0];
sx q[0];
rz(-1.4109857) q[0];
sx q[0];
rz(3.0335866) q[0];
rz(-pi) q[1];
rz(0.043134281) q[2];
sx q[2];
rz(-2.1675674) q[2];
sx q[2];
rz(1.5297001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3276931) q[1];
sx q[1];
rz(-2.3997953) q[1];
sx q[1];
rz(-0.90403907) q[1];
x q[2];
rz(-1.2040505) q[3];
sx q[3];
rz(-2.2270508) q[3];
sx q[3];
rz(-1.3964069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7654045) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(-0.21793951) q[2];
rz(-2.229522) q[3];
sx q[3];
rz(-1.4927161) q[3];
sx q[3];
rz(0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243778) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(0.81556129) q[0];
rz(-2.0888445) q[1];
sx q[1];
rz(-2.2595854) q[1];
sx q[1];
rz(-1.3515333) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4824351) q[0];
sx q[0];
rz(-1.1023942) q[0];
sx q[0];
rz(-1.4841561) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4250373) q[2];
sx q[2];
rz(-1.5045696) q[2];
sx q[2];
rz(-0.19342455) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4554418) q[1];
sx q[1];
rz(-1.0956005) q[1];
sx q[1];
rz(0.032217044) q[1];
rz(-pi) q[2];
rz(2.0729321) q[3];
sx q[3];
rz(-2.2229513) q[3];
sx q[3];
rz(-0.38092914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1286596) q[2];
sx q[2];
rz(-1.9482875) q[2];
sx q[2];
rz(2.7093757) q[2];
rz(0.84248078) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(-2.1459818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1603482) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(2.5904742) q[0];
rz(0.61800686) q[1];
sx q[1];
rz(-1.2851241) q[1];
sx q[1];
rz(1.0689703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38902125) q[0];
sx q[0];
rz(-0.71439904) q[0];
sx q[0];
rz(0.32984372) q[0];
x q[1];
rz(1.1177914) q[2];
sx q[2];
rz(-2.0108622) q[2];
sx q[2];
rz(-0.58811448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3772419) q[1];
sx q[1];
rz(-2.3016986) q[1];
sx q[1];
rz(-2.9725171) q[1];
x q[2];
rz(2.8884747) q[3];
sx q[3];
rz(-0.52682823) q[3];
sx q[3];
rz(-1.7071068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7908287) q[2];
sx q[2];
rz(-0.5841693) q[2];
sx q[2];
rz(2.186415) q[2];
rz(-2.5480934) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(-1.745863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3747568) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(-1.7497077) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(-1.3124189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29762938) q[0];
sx q[0];
rz(-2.2052551) q[0];
sx q[0];
rz(-2.5353801) q[0];
rz(-0.8719123) q[2];
sx q[2];
rz(-2.12247) q[2];
sx q[2];
rz(1.0629423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3692664) q[1];
sx q[1];
rz(-1.516894) q[1];
sx q[1];
rz(-0.51477706) q[1];
rz(1.7465318) q[3];
sx q[3];
rz(-0.32009691) q[3];
sx q[3];
rz(1.2552346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7890847) q[2];
sx q[2];
rz(-0.75560537) q[2];
sx q[2];
rz(1.7279203) q[2];
rz(2.6740354) q[3];
sx q[3];
rz(-0.65842015) q[3];
sx q[3];
rz(0.55268923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6308052) q[0];
sx q[0];
rz(-3.0785705) q[0];
sx q[0];
rz(0.68156534) q[0];
rz(-1.956578) q[1];
sx q[1];
rz(-0.76528913) q[1];
sx q[1];
rz(-0.78651816) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1543321) q[0];
sx q[0];
rz(-1.904878) q[0];
sx q[0];
rz(-0.41436899) q[0];
rz(-pi) q[1];
rz(0.94633905) q[2];
sx q[2];
rz(-1.985744) q[2];
sx q[2];
rz(-3.0111661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6326846) q[1];
sx q[1];
rz(-0.68636319) q[1];
sx q[1];
rz(-1.4014161) q[1];
rz(-pi) q[2];
rz(1.5249355) q[3];
sx q[3];
rz(-0.73165441) q[3];
sx q[3];
rz(-1.589244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6098392) q[2];
sx q[2];
rz(-0.24007758) q[2];
sx q[2];
rz(-1.5698203) q[2];
rz(-1.2872559) q[3];
sx q[3];
rz(-1.6183034) q[3];
sx q[3];
rz(2.1985445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.542273) q[0];
sx q[0];
rz(-2.4241408) q[0];
sx q[0];
rz(-1.9532816) q[0];
rz(-0.020847281) q[1];
sx q[1];
rz(-2.3884845) q[1];
sx q[1];
rz(-0.59897024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1513903) q[0];
sx q[0];
rz(-1.6693475) q[0];
sx q[0];
rz(-2.0771785) q[0];
rz(2.0161736) q[2];
sx q[2];
rz(-0.76596224) q[2];
sx q[2];
rz(-0.19419032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7341344) q[1];
sx q[1];
rz(-0.86454138) q[1];
sx q[1];
rz(0.64532537) q[1];
rz(-pi) q[2];
rz(-1.0858367) q[3];
sx q[3];
rz(-2.3502878) q[3];
sx q[3];
rz(-2.6216061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8591566) q[2];
sx q[2];
rz(-1.2500637) q[2];
sx q[2];
rz(0.51521987) q[2];
rz(2.6089) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(2.2044619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.027997967) q[0];
sx q[0];
rz(-2.4585215) q[0];
sx q[0];
rz(-2.7711476) q[0];
rz(2.0619552) q[1];
sx q[1];
rz(-0.41079435) q[1];
sx q[1];
rz(3.0830141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3315225) q[0];
sx q[0];
rz(-2.5096623) q[0];
sx q[0];
rz(-2.8721316) q[0];
rz(-pi) q[1];
rz(0.04444261) q[2];
sx q[2];
rz(-1.7979597) q[2];
sx q[2];
rz(-2.4054766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6296719) q[1];
sx q[1];
rz(-1.4751993) q[1];
sx q[1];
rz(2.3562576) q[1];
rz(-0.58038099) q[3];
sx q[3];
rz(-2.2182052) q[3];
sx q[3];
rz(-3.1064432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8687245) q[2];
sx q[2];
rz(-0.803002) q[2];
sx q[2];
rz(2.8105984) q[2];
rz(0.013414772) q[3];
sx q[3];
rz(-2.1881073) q[3];
sx q[3];
rz(1.0564055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57824221) q[0];
sx q[0];
rz(-0.42911068) q[0];
sx q[0];
rz(-0.44813928) q[0];
rz(-3.0478802) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(1.2154481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96998065) q[0];
sx q[0];
rz(-1.5413862) q[0];
sx q[0];
rz(1.1377115) q[0];
rz(-pi) q[1];
rz(2.3289069) q[2];
sx q[2];
rz(-1.865109) q[2];
sx q[2];
rz(-3.091623) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61291873) q[1];
sx q[1];
rz(-1.8380751) q[1];
sx q[1];
rz(-2.1211995) q[1];
x q[2];
rz(1.8849202) q[3];
sx q[3];
rz(-2.7249955) q[3];
sx q[3];
rz(-2.9521717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.965968) q[2];
sx q[2];
rz(-1.5529239) q[2];
sx q[2];
rz(-2.442181) q[2];
rz(1.2753963) q[3];
sx q[3];
rz(-1.1558775) q[3];
sx q[3];
rz(-1.2709966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4324343) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(-1.746183) q[1];
sx q[1];
rz(-1.9231053) q[1];
sx q[1];
rz(-1.3833192) q[1];
rz(2.4757842) q[2];
sx q[2];
rz(-1.5077017) q[2];
sx q[2];
rz(-2.2617819) q[2];
rz(2.6752409) q[3];
sx q[3];
rz(-1.6771355) q[3];
sx q[3];
rz(1.9775978) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
