OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7140952) q[0];
sx q[0];
rz(-2.5795955) q[0];
sx q[0];
rz(-0.23101097) q[0];
rz(-2.8958939) q[1];
sx q[1];
rz(-2.6872771) q[1];
sx q[1];
rz(1.8543724) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79345771) q[0];
sx q[0];
rz(-2.0525816) q[0];
sx q[0];
rz(-2.8346377) q[0];
rz(0.24536774) q[2];
sx q[2];
rz(-1.3431708) q[2];
sx q[2];
rz(-0.99492225) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8227777) q[1];
sx q[1];
rz(-2.9479369) q[1];
sx q[1];
rz(0.53656399) q[1];
rz(-pi) q[2];
rz(1.3955529) q[3];
sx q[3];
rz(-1.4966432) q[3];
sx q[3];
rz(1.7046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40453688) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(0.36049584) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(1.2688961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26038134) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(2.4847109) q[0];
rz(-2.8086713) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(-2.2064256) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2774076) q[0];
sx q[0];
rz(-1.5404697) q[0];
sx q[0];
rz(-3.0369989) q[0];
rz(-pi) q[1];
rz(-2.8043973) q[2];
sx q[2];
rz(-2.8697578) q[2];
sx q[2];
rz(1.7295009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.024927372) q[1];
sx q[1];
rz(-1.5313193) q[1];
sx q[1];
rz(2.9945606) q[1];
rz(-0.35547361) q[3];
sx q[3];
rz(-1.7661816) q[3];
sx q[3];
rz(0.030578407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8344581) q[2];
sx q[2];
rz(-1.5156526) q[2];
sx q[2];
rz(-1.9251941) q[2];
rz(0.52792102) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.5288178) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(-0.58498996) q[0];
rz(0.12475573) q[1];
sx q[1];
rz(-0.59586066) q[1];
sx q[1];
rz(-1.4488719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.751469) q[0];
sx q[0];
rz(-1.5737783) q[0];
sx q[0];
rz(-0.0026767038) q[0];
x q[1];
rz(2.9563006) q[2];
sx q[2];
rz(-3.0681562) q[2];
sx q[2];
rz(2.9543608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6022012) q[1];
sx q[1];
rz(-2.3219206) q[1];
sx q[1];
rz(2.17893) q[1];
rz(-pi) q[2];
rz(1.307906) q[3];
sx q[3];
rz(-1.2440727) q[3];
sx q[3];
rz(2.9430091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8831732) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(-1.4460571) q[2];
rz(2.3840733) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3676753) q[0];
sx q[0];
rz(-0.92905074) q[0];
sx q[0];
rz(1.0876592) q[0];
rz(-0.37519535) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(0.96022022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4488569) q[0];
sx q[0];
rz(-1.5310643) q[0];
sx q[0];
rz(-1.3055152) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2795936) q[2];
sx q[2];
rz(-3.0377977) q[2];
sx q[2];
rz(1.7218931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.474388) q[1];
sx q[1];
rz(-1.9992113) q[1];
sx q[1];
rz(-2.0640255) q[1];
rz(-2.8071515) q[3];
sx q[3];
rz(-0.28452793) q[3];
sx q[3];
rz(-0.65438017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.20597657) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(-1.6440294) q[2];
rz(0.86132541) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(-0.45708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1645666) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(2.5373051) q[0];
rz(-1.1445649) q[1];
sx q[1];
rz(-1.4860169) q[1];
sx q[1];
rz(0.42246517) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6910962) q[0];
sx q[0];
rz(-1.7393149) q[0];
sx q[0];
rz(0.082534747) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79549148) q[2];
sx q[2];
rz(-2.4078712) q[2];
sx q[2];
rz(0.41662859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3389272) q[1];
sx q[1];
rz(-2.4484854) q[1];
sx q[1];
rz(-0.1130123) q[1];
x q[2];
rz(0.38099184) q[3];
sx q[3];
rz(-1.89011) q[3];
sx q[3];
rz(1.8007985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3182688) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(0.65625119) q[2];
rz(1.5308135) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(-2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4055279) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(2.5166125) q[0];
rz(-2.3896353) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(-1.5350852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84941712) q[0];
sx q[0];
rz(-1.8596974) q[0];
sx q[0];
rz(-1.5194511) q[0];
rz(0.70930945) q[2];
sx q[2];
rz(-1.0934208) q[2];
sx q[2];
rz(-2.2717182) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2481164) q[1];
sx q[1];
rz(-1.3424338) q[1];
sx q[1];
rz(2.5121157) q[1];
rz(0.97134892) q[3];
sx q[3];
rz(-1.1068212) q[3];
sx q[3];
rz(-0.52677192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.357888) q[2];
sx q[2];
rz(-2.2124186) q[2];
rz(2.3164228) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850942) q[0];
sx q[0];
rz(-0.80699054) q[0];
sx q[0];
rz(1.7671385) q[0];
rz(-2.6311686) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(-0.62320954) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7316298) q[0];
sx q[0];
rz(-1.4193168) q[0];
sx q[0];
rz(1.9069457) q[0];
rz(-0.34131949) q[2];
sx q[2];
rz(-2.6067197) q[2];
sx q[2];
rz(0.90082263) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.12909266) q[1];
sx q[1];
rz(-2.9575037) q[1];
sx q[1];
rz(-2.2232008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21824746) q[3];
sx q[3];
rz(-0.75957662) q[3];
sx q[3];
rz(-1.4324566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0906543) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(-1.1673002) q[2];
rz(-0.022631571) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(1.1289271) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24564329) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(2.1242712) q[0];
rz(1.7474878) q[1];
sx q[1];
rz(-0.21109763) q[1];
sx q[1];
rz(0.62754935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.548259) q[0];
sx q[0];
rz(-2.152266) q[0];
sx q[0];
rz(2.6618746) q[0];
rz(1.4862138) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(1.5280452) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99929727) q[1];
sx q[1];
rz(-2.8311756) q[1];
sx q[1];
rz(0.69119549) q[1];
x q[2];
rz(-0.60365898) q[3];
sx q[3];
rz(-2.5839879) q[3];
sx q[3];
rz(-1.738036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62577406) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(-2.7308357) q[2];
rz(0.050203236) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(0.075411782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.582616) q[0];
sx q[0];
rz(-2.2447383) q[0];
sx q[0];
rz(0.26891747) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(3.0115829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1300505) q[0];
sx q[0];
rz(-0.39549144) q[0];
sx q[0];
rz(-0.96590913) q[0];
rz(-pi) q[1];
x q[1];
rz(1.214005) q[2];
sx q[2];
rz(-1.901682) q[2];
sx q[2];
rz(0.39168229) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8925608) q[1];
sx q[1];
rz(-1.7317803) q[1];
sx q[1];
rz(-2.575483) q[1];
rz(0.49874108) q[3];
sx q[3];
rz(-2.1864656) q[3];
sx q[3];
rz(1.9174674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2076063) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(-0.096573528) q[2];
rz(-1.6377595) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(-0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0132975) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(2.6085594) q[0];
rz(0.036272613) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(-1.4520377) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743458) q[0];
sx q[0];
rz(-1.1935992) q[0];
sx q[0];
rz(2.1826571) q[0];
rz(-1.6797941) q[2];
sx q[2];
rz(-1.2533292) q[2];
sx q[2];
rz(-0.017680971) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8126077) q[1];
sx q[1];
rz(-1.2115062) q[1];
sx q[1];
rz(1.6027662) q[1];
rz(2.96569) q[3];
sx q[3];
rz(-2.6697192) q[3];
sx q[3];
rz(-2.7681818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71296802) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(3.1089605) q[2];
rz(-0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(1.0802065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7964771) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(2.3445917) q[1];
sx q[1];
rz(-2.2140257) q[1];
sx q[1];
rz(0.066233403) q[1];
rz(2.6673139) q[2];
sx q[2];
rz(-2.5534292) q[2];
sx q[2];
rz(-3.0234887) q[2];
rz(0.12049051) q[3];
sx q[3];
rz(-2.3269666) q[3];
sx q[3];
rz(2.5757488) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
