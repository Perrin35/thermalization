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
rz(1.3069557) q[0];
sx q[0];
rz(4.0210273) q[0];
sx q[0];
rz(9.9183912) q[0];
rz(1.8618795) q[1];
sx q[1];
rz(-0.6266098) q[1];
sx q[1];
rz(1.4785179) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8956229) q[0];
sx q[0];
rz(-2.1969271) q[0];
sx q[0];
rz(-2.5139721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25716146) q[2];
sx q[2];
rz(-0.92338054) q[2];
sx q[2];
rz(1.1360628) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50861102) q[1];
sx q[1];
rz(-0.87978432) q[1];
sx q[1];
rz(1.8553084) q[1];
rz(-pi) q[2];
rz(-0.36840393) q[3];
sx q[3];
rz(-0.96782902) q[3];
sx q[3];
rz(1.4908229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1349858) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(-1.1084278) q[2];
rz(-1.8290352) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(2.826622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38158622) q[0];
sx q[0];
rz(-2.3158323) q[0];
sx q[0];
rz(0.015856892) q[0];
rz(-0.99041692) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(2.2784065) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1970438) q[0];
sx q[0];
rz(-1.6672575) q[0];
sx q[0];
rz(2.405557) q[0];
rz(-pi) q[1];
rz(-1.241774) q[2];
sx q[2];
rz(-1.8890427) q[2];
sx q[2];
rz(1.1877738) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6241315) q[1];
sx q[1];
rz(-1.7449813) q[1];
sx q[1];
rz(2.6427173) q[1];
x q[2];
rz(0.14948577) q[3];
sx q[3];
rz(-0.85452291) q[3];
sx q[3];
rz(-3.1325185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4414759) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(0.46879834) q[2];
rz(-0.68764728) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0013292) q[0];
sx q[0];
rz(-2.2183473) q[0];
sx q[0];
rz(2.4053251) q[0];
rz(-2.7983792) q[1];
sx q[1];
rz(-0.88756573) q[1];
sx q[1];
rz(2.9578178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7385173) q[0];
sx q[0];
rz(-0.74664298) q[0];
sx q[0];
rz(-0.097642032) q[0];
rz(-pi) q[1];
rz(-0.60004514) q[2];
sx q[2];
rz(-2.349424) q[2];
sx q[2];
rz(-0.60952696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6334875) q[1];
sx q[1];
rz(-2.268666) q[1];
sx q[1];
rz(1.9644323) q[1];
rz(-1.0093498) q[3];
sx q[3];
rz(-0.61547503) q[3];
sx q[3];
rz(-1.4765679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78202128) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(-2.6256631) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(-1.1512604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16294031) q[0];
sx q[0];
rz(-1.2115703) q[0];
sx q[0];
rz(0.51280713) q[0];
rz(-2.6817952) q[1];
sx q[1];
rz(-1.5922981) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490246) q[0];
sx q[0];
rz(-1.6852753) q[0];
sx q[0];
rz(3.0106198) q[0];
rz(2.725179) q[2];
sx q[2];
rz(-1.8685307) q[2];
sx q[2];
rz(-2.1418051) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1000378) q[1];
sx q[1];
rz(-2.6943992) q[1];
sx q[1];
rz(-0.87765043) q[1];
x q[2];
rz(-0.38171347) q[3];
sx q[3];
rz(-0.36892051) q[3];
sx q[3];
rz(2.5526508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0189455) q[2];
sx q[2];
rz(-2.985869) q[2];
sx q[2];
rz(2.8141008) q[2];
rz(-0.28199768) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.17664385) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(0.97736812) q[0];
rz(0.36879677) q[1];
sx q[1];
rz(-0.86124033) q[1];
sx q[1];
rz(-1.4126973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644625) q[0];
sx q[0];
rz(-2.6106195) q[0];
sx q[0];
rz(1.9524379) q[0];
rz(-pi) q[1];
rz(1.0540038) q[2];
sx q[2];
rz(-1.5828504) q[2];
sx q[2];
rz(0.64753676) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.818644) q[1];
sx q[1];
rz(-1.6020157) q[1];
sx q[1];
rz(-1.853456) q[1];
x q[2];
rz(-1.2431954) q[3];
sx q[3];
rz(-2.1246111) q[3];
sx q[3];
rz(2.6462951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8480924) q[2];
sx q[2];
rz(-0.81587452) q[2];
sx q[2];
rz(0.14981848) q[2];
rz(-0.39854974) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77001101) q[0];
sx q[0];
rz(-0.12938975) q[0];
sx q[0];
rz(-2.5883801) q[0];
rz(-0.54168701) q[1];
sx q[1];
rz(-1.0463511) q[1];
sx q[1];
rz(-2.6936626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2509237) q[0];
sx q[0];
rz(-2.726993) q[0];
sx q[0];
rz(2.7756734) q[0];
x q[1];
rz(-0.42945736) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(-2.7384659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3975911) q[1];
sx q[1];
rz(-0.41644704) q[1];
sx q[1];
rz(-1.9201502) q[1];
rz(-pi) q[2];
rz(-2.8087696) q[3];
sx q[3];
rz(-1.5170005) q[3];
sx q[3];
rz(-1.4383565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64122671) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(2.3424303) q[2];
rz(-2.7568119) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(-2.1414115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36050972) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(-0.59180301) q[0];
rz(2.7599755) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(2.9891678) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29983175) q[0];
sx q[0];
rz(-0.61474568) q[0];
sx q[0];
rz(-2.1257945) q[0];
x q[1];
rz(0.32976361) q[2];
sx q[2];
rz(-0.50055365) q[2];
sx q[2];
rz(-0.19937521) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.287474) q[1];
sx q[1];
rz(-1.9175314) q[1];
sx q[1];
rz(-1.5902014) q[1];
rz(-pi) q[2];
rz(2.2493208) q[3];
sx q[3];
rz(-2.0990058) q[3];
sx q[3];
rz(1.5367374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3064208) q[2];
sx q[2];
rz(-0.99205899) q[2];
sx q[2];
rz(-1.6888118) q[2];
rz(-2.6460904) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(-2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.751916) q[0];
sx q[0];
rz(-3.1041807) q[0];
sx q[0];
rz(-2.9801242) q[0];
rz(0.031919315) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(1.8643103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91457462) q[0];
sx q[0];
rz(-2.380207) q[0];
sx q[0];
rz(-0.37933357) q[0];
x q[1];
rz(-2.0525682) q[2];
sx q[2];
rz(-2.6725997) q[2];
sx q[2];
rz(-0.8587786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1254221) q[1];
sx q[1];
rz(-1.6787663) q[1];
sx q[1];
rz(2.8362464) q[1];
rz(-0.038441258) q[3];
sx q[3];
rz(-1.2769967) q[3];
sx q[3];
rz(0.97803365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45234597) q[2];
sx q[2];
rz(-0.9114868) q[2];
sx q[2];
rz(-0.39608836) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(-0.8766492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006627) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(2.990429) q[0];
rz(-0.97686544) q[1];
sx q[1];
rz(-2.7604389) q[1];
sx q[1];
rz(-0.74473286) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1158651) q[0];
sx q[0];
rz(-2.7926499) q[0];
sx q[0];
rz(0.92744382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1579723) q[2];
sx q[2];
rz(-1.4399488) q[2];
sx q[2];
rz(-1.3589588) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3754282) q[1];
sx q[1];
rz(-0.64830983) q[1];
sx q[1];
rz(-3.1392155) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58988692) q[3];
sx q[3];
rz(-2.7253236) q[3];
sx q[3];
rz(1.2225162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5186844) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(0.94804478) q[2];
rz(-2.0840123) q[3];
sx q[3];
rz(-1.4761997) q[3];
sx q[3];
rz(2.0675366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093216151) q[0];
sx q[0];
rz(-2.8792448) q[0];
sx q[0];
rz(-2.9076305) q[0];
rz(1.9562862) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(1.2497466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3435182) q[0];
sx q[0];
rz(-1.8600704) q[0];
sx q[0];
rz(-0.57855655) q[0];
rz(-pi) q[1];
rz(0.62057497) q[2];
sx q[2];
rz(-1.7399551) q[2];
sx q[2];
rz(0.38693869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1162303) q[1];
sx q[1];
rz(-0.90060189) q[1];
sx q[1];
rz(2.4306562) q[1];
x q[2];
rz(-0.21548157) q[3];
sx q[3];
rz(-0.69284791) q[3];
sx q[3];
rz(0.97164916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.027111) q[2];
sx q[2];
rz(-1.9805084) q[2];
sx q[2];
rz(1.4410045) q[2];
rz(-0.13127413) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-2.89768) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048653614) q[0];
sx q[0];
rz(-2.1760512) q[0];
sx q[0];
rz(-2.0078134) q[0];
rz(1.0406021) q[1];
sx q[1];
rz(-0.84747172) q[1];
sx q[1];
rz(-0.52451959) q[1];
rz(0.94338633) q[2];
sx q[2];
rz(-1.8205943) q[2];
sx q[2];
rz(2.7216507) q[2];
rz(-1.7215988) q[3];
sx q[3];
rz(-2.1986113) q[3];
sx q[3];
rz(-2.1049706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
