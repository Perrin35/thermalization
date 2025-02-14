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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(-1.8921312) q[0];
rz(-0.14021048) q[1];
sx q[1];
rz(-1.6970716) q[1];
sx q[1];
rz(0.061847774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0573579) q[0];
sx q[0];
rz(-2.0664729) q[0];
sx q[0];
rz(-2.2593014) q[0];
x q[1];
rz(-3.05117) q[2];
sx q[2];
rz(-1.5704463) q[2];
sx q[2];
rz(1.4594913) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5995591) q[1];
sx q[1];
rz(-2.2890511) q[1];
sx q[1];
rz(0.12964779) q[1];
rz(-pi) q[2];
rz(-2.0560045) q[3];
sx q[3];
rz(-2.3923529) q[3];
sx q[3];
rz(-0.10252122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7141815) q[2];
sx q[2];
rz(-0.84250557) q[2];
sx q[2];
rz(0.25644914) q[2];
rz(2.9668258) q[3];
sx q[3];
rz(-1.9758965) q[3];
sx q[3];
rz(2.3037361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3278219) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-3.0178965) q[0];
rz(-0.23873121) q[1];
sx q[1];
rz(-2.3614466) q[1];
sx q[1];
rz(0.10496584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83451) q[0];
sx q[0];
rz(-1.7210467) q[0];
sx q[0];
rz(-2.2730877) q[0];
rz(0.41609515) q[2];
sx q[2];
rz(-1.7836092) q[2];
sx q[2];
rz(-2.0426322) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1233842) q[1];
sx q[1];
rz(-1.3088641) q[1];
sx q[1];
rz(-2.0831152) q[1];
rz(-2.4363748) q[3];
sx q[3];
rz(-1.2949484) q[3];
sx q[3];
rz(-2.5666756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.77439848) q[2];
sx q[2];
rz(-1.4906733) q[2];
sx q[2];
rz(1.4614089) q[2];
rz(-1.2857619) q[3];
sx q[3];
rz(-1.4444084) q[3];
sx q[3];
rz(0.27209601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937781) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(-1.0149957) q[0];
rz(-2.2442832) q[1];
sx q[1];
rz(-1.6474479) q[1];
sx q[1];
rz(-2.4651333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2059442) q[0];
sx q[0];
rz(-2.1876011) q[0];
sx q[0];
rz(0.031868462) q[0];
rz(-pi) q[1];
rz(-0.93134201) q[2];
sx q[2];
rz(-2.4081875) q[2];
sx q[2];
rz(2.1426107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30005042) q[1];
sx q[1];
rz(-1.9498697) q[1];
sx q[1];
rz(-2.9247523) q[1];
rz(-0.018142975) q[3];
sx q[3];
rz(-0.81644316) q[3];
sx q[3];
rz(1.3870365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2066388) q[2];
sx q[2];
rz(-0.97524869) q[2];
sx q[2];
rz(2.3812531) q[2];
rz(1.9350516) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(-2.2981203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7059785) q[0];
sx q[0];
rz(-0.62788457) q[0];
sx q[0];
rz(1.5465558) q[0];
rz(-0.71890038) q[1];
sx q[1];
rz(-1.8214106) q[1];
sx q[1];
rz(2.9100606) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0800832) q[0];
sx q[0];
rz(-0.69893796) q[0];
sx q[0];
rz(0.7028347) q[0];
rz(-pi) q[1];
rz(3.1377085) q[2];
sx q[2];
rz(-2.3489315) q[2];
sx q[2];
rz(2.5156227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0147103) q[1];
sx q[1];
rz(-1.0877123) q[1];
sx q[1];
rz(1.4654069) q[1];
rz(2.5367686) q[3];
sx q[3];
rz(-1.8993653) q[3];
sx q[3];
rz(0.25024763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4270758) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(0.60849774) q[2];
rz(-0.19615873) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(-2.1430446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8683559) q[0];
sx q[0];
rz(-0.65160692) q[0];
sx q[0];
rz(0.77907816) q[0];
rz(-2.211606) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(-1.9680061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55803199) q[0];
sx q[0];
rz(-2.8035218) q[0];
sx q[0];
rz(1.5589236) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7601682) q[2];
sx q[2];
rz(-2.5177285) q[2];
sx q[2];
rz(-0.80846918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0231578) q[1];
sx q[1];
rz(-1.506532) q[1];
sx q[1];
rz(1.2530909) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7125732) q[3];
sx q[3];
rz(-1.7439505) q[3];
sx q[3];
rz(0.28399434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.32213) q[2];
sx q[2];
rz(-0.16699114) q[2];
sx q[2];
rz(-0.67506153) q[2];
rz(-2.2760462) q[3];
sx q[3];
rz(-1.6277438) q[3];
sx q[3];
rz(1.9338098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753767) q[0];
sx q[0];
rz(-0.017711552) q[0];
sx q[0];
rz(2.2702763) q[0];
rz(-0.21698347) q[1];
sx q[1];
rz(-1.5419518) q[1];
sx q[1];
rz(-1.3380231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9674613) q[0];
sx q[0];
rz(-0.92087692) q[0];
sx q[0];
rz(1.4667257) q[0];
rz(-pi) q[1];
rz(-0.14071847) q[2];
sx q[2];
rz(-0.85007668) q[2];
sx q[2];
rz(-2.7776133) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29472953) q[1];
sx q[1];
rz(-1.4275075) q[1];
sx q[1];
rz(-3.0927977) q[1];
rz(-pi) q[2];
rz(2.7910978) q[3];
sx q[3];
rz(-0.51096254) q[3];
sx q[3];
rz(3.0929589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1285105) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(0.80005542) q[2];
rz(0.99177805) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(-0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8868788) q[0];
sx q[0];
rz(-2.7006221) q[0];
sx q[0];
rz(0.15175858) q[0];
rz(-1.4166547) q[1];
sx q[1];
rz(-2.2817426) q[1];
sx q[1];
rz(0.83622611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0043512303) q[0];
sx q[0];
rz(-2.3286331) q[0];
sx q[0];
rz(-1.0197958) q[0];
rz(-pi) q[1];
rz(-0.57668893) q[2];
sx q[2];
rz(-1.5637914) q[2];
sx q[2];
rz(-0.079886111) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9735896) q[1];
sx q[1];
rz(-0.26339809) q[1];
sx q[1];
rz(1.9067287) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8854675) q[3];
sx q[3];
rz(-0.47893347) q[3];
sx q[3];
rz(-1.1161982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6553216) q[2];
sx q[2];
rz(-0.063491193) q[2];
sx q[2];
rz(3.0625694) q[2];
rz(2.4231353) q[3];
sx q[3];
rz(-1.5114762) q[3];
sx q[3];
rz(0.92459905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.786161) q[0];
sx q[0];
rz(-2.5434255) q[0];
sx q[0];
rz(2.2736736) q[0];
rz(-0.68880853) q[1];
sx q[1];
rz(-0.30615607) q[1];
sx q[1];
rz(2.205663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.717632) q[0];
sx q[0];
rz(-2.0905295) q[0];
sx q[0];
rz(1.3613767) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9804968) q[2];
sx q[2];
rz(-1.2563224) q[2];
sx q[2];
rz(-2.1358228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3080146) q[1];
sx q[1];
rz(-1.8955827) q[1];
sx q[1];
rz(-0.056450162) q[1];
rz(-pi) q[2];
rz(0.74347382) q[3];
sx q[3];
rz(-1.5697362) q[3];
sx q[3];
rz(2.8618232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1615289) q[2];
sx q[2];
rz(-2.4825725) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(-1.0181381) q[3];
sx q[3];
rz(-1.1192106) q[3];
sx q[3];
rz(2.170678) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7629906) q[0];
sx q[0];
rz(-0.55736962) q[0];
sx q[0];
rz(-2.2879404) q[0];
rz(1.254982) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(-2.8010211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7434692) q[0];
sx q[0];
rz(-2.8883683) q[0];
sx q[0];
rz(0.88078334) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3236295) q[2];
sx q[2];
rz(-0.91193141) q[2];
sx q[2];
rz(0.62550046) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7813959) q[1];
sx q[1];
rz(-2.2193647) q[1];
sx q[1];
rz(-0.99332033) q[1];
x q[2];
rz(1.1084308) q[3];
sx q[3];
rz(-1.8183072) q[3];
sx q[3];
rz(2.4570297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7178932) q[2];
sx q[2];
rz(-1.7607949) q[2];
sx q[2];
rz(0.37810668) q[2];
rz(-2.9041491) q[3];
sx q[3];
rz(-2.0707371) q[3];
sx q[3];
rz(2.8606991) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0963999) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(-0.64724809) q[0];
rz(-1.1680394) q[1];
sx q[1];
rz(-2.2868575) q[1];
sx q[1];
rz(-0.38356575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76825324) q[0];
sx q[0];
rz(-3.1041234) q[0];
sx q[0];
rz(-0.61385198) q[0];
rz(-2.6099989) q[2];
sx q[2];
rz(-1.2832912) q[2];
sx q[2];
rz(0.93898857) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3473222) q[1];
sx q[1];
rz(-1.8350661) q[1];
sx q[1];
rz(1.3572973) q[1];
rz(-2.0862574) q[3];
sx q[3];
rz(-1.4773507) q[3];
sx q[3];
rz(-0.44778131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26455227) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(1.5569347) q[2];
rz(-1.5133739) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(-0.42678601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601111) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(1.1813286) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(-2.8988373) q[2];
sx q[2];
rz(-1.1052255) q[2];
sx q[2];
rz(0.33585264) q[2];
rz(-0.060130974) q[3];
sx q[3];
rz(-2.6562666) q[3];
sx q[3];
rz(-0.2851214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
