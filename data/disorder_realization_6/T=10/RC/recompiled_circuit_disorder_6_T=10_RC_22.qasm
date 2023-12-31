OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(4.2545118) q[1];
sx q[1];
rz(1.7634044) q[1];
sx q[1];
rz(7.4982285) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5117447) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(1.3826136) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7772929) q[2];
sx q[2];
rz(-2.4476353) q[2];
sx q[2];
rz(-2.7147164) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3455968) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(0.3791581) q[1];
rz(-pi) q[2];
rz(-2.0753161) q[3];
sx q[3];
rz(-2.8404232) q[3];
sx q[3];
rz(2.1198213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(2.9585178) q[2];
rz(-2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(-2.8027957) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2167643) q[0];
sx q[0];
rz(-0.98470682) q[0];
sx q[0];
rz(-3.0490962) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23761959) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(1.2653637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51255608) q[1];
sx q[1];
rz(-1.7450383) q[1];
sx q[1];
rz(1.0838572) q[1];
rz(2.6529796) q[3];
sx q[3];
rz(-0.47227898) q[3];
sx q[3];
rz(-2.3924535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028458683) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-2.5862397) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1028324) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(1.0981512) q[2];
sx q[2];
rz(-2.522509) q[2];
sx q[2];
rz(1.667779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1173646) q[1];
sx q[1];
rz(-1.4091361) q[1];
sx q[1];
rz(-2.2092186) q[1];
rz(-pi) q[2];
rz(1.8954574) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(-0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(2.326791) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(-2.8864158) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6107451) q[0];
sx q[0];
rz(-1.8662211) q[0];
sx q[0];
rz(0.3770963) q[0];
rz(-pi) q[1];
rz(-2.7093676) q[2];
sx q[2];
rz(-2.4556293) q[2];
sx q[2];
rz(-1.9908817) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.081101) q[1];
sx q[1];
rz(-1.8685409) q[1];
sx q[1];
rz(-0.22025073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4820443) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(0.50597092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(-0.056093562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010904) q[0];
sx q[0];
rz(-0.70106693) q[0];
sx q[0];
rz(0.58347337) q[0];
x q[1];
rz(-0.34727879) q[2];
sx q[2];
rz(-2.0060853) q[2];
sx q[2];
rz(1.4594644) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1069378) q[1];
sx q[1];
rz(-1.5853303) q[1];
sx q[1];
rz(2.8388139) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0104996) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(1.2922985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(-0.8941935) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(0.9544968) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4650824) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(-0.17980534) q[0];
rz(-2.2013821) q[2];
sx q[2];
rz(-1.5427089) q[2];
sx q[2];
rz(2.0036151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94297385) q[1];
sx q[1];
rz(-1.0620411) q[1];
sx q[1];
rz(-1.0191304) q[1];
x q[2];
rz(-1.0682085) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(0.1474895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.59297562) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(0.61002237) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726495) q[0];
sx q[0];
rz(-2.5805051) q[0];
sx q[0];
rz(0.48604301) q[0];
rz(1.4201944) q[2];
sx q[2];
rz(-0.65956958) q[2];
sx q[2];
rz(1.2130376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2513189) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(3.1097079) q[1];
rz(-2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.6280744) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37305957) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(-1.2241227) q[0];
rz(-pi) q[1];
rz(1.2484776) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(-1.581574) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3031993) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(-1.4569015) q[1];
x q[2];
rz(1.1235808) q[3];
sx q[3];
rz(-1.2774602) q[3];
sx q[3];
rz(0.48188996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1071876) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.3945403) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(2.396778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9841524) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-2.0122583) q[0];
rz(-1.1698193) q[2];
sx q[2];
rz(-2.3224761) q[2];
sx q[2];
rz(-2.4405406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71036584) q[1];
sx q[1];
rz(-1.2214298) q[1];
sx q[1];
rz(2.9037895) q[1];
rz(2.3585988) q[3];
sx q[3];
rz(-1.8564312) q[3];
sx q[3];
rz(-2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.194681) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(1.9706479) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5491345) q[0];
sx q[0];
rz(-2.7146656) q[0];
x q[1];
rz(0.2756341) q[2];
sx q[2];
rz(-0.79489691) q[2];
sx q[2];
rz(-0.3030215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3552637) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(-1.5481871) q[1];
rz(-pi) q[2];
rz(0.42697866) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(-0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(-0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.6177572) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-0.70384937) q[2];
sx q[2];
rz(-2.2675632) q[2];
sx q[2];
rz(2.8013196) q[2];
rz(0.016146544) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
