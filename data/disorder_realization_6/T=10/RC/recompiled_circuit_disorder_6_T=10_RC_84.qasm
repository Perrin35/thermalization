OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6157827) q[0];
sx q[0];
rz(-1.4178185) q[0];
sx q[0];
rz(0.56086993) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(1.9265494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6998866) q[0];
sx q[0];
rz(-2.9268648) q[0];
sx q[0];
rz(1.0617274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(1.7125318) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4814264) q[1];
sx q[1];
rz(-2.5170442) q[1];
sx q[1];
rz(2.156483) q[1];
rz(1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(-0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.5391301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9248283) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(3.0490962) q[0];
rz(-1.7878754) q[2];
sx q[2];
rz(-0.84393822) q[2];
sx q[2];
rz(-2.196687) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7669249) q[1];
sx q[1];
rz(-0.51480773) q[1];
sx q[1];
rz(1.2109846) q[1];
x q[2];
rz(2.7178571) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(-2.7620897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-2.7084896) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12268513) q[0];
sx q[0];
rz(-2.2203831) q[0];
sx q[0];
rz(0.62223776) q[0];
rz(-1.0054587) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(-2.8440059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66545031) q[1];
sx q[1];
rz(-0.94201554) q[1];
sx q[1];
rz(-0.20035845) q[1];
rz(-pi) q[2];
x q[2];
rz(1.025612) q[3];
sx q[3];
rz(-1.7435929) q[3];
sx q[3];
rz(-1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4042523) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.2505442) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(-2.8864158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4677306) q[0];
sx q[0];
rz(-2.6669589) q[0];
sx q[0];
rz(0.69068308) q[0];
x q[1];
rz(-2.7093676) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(-1.150711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.081101) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(-2.9213419) q[1];
x q[2];
rz(1.049794) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(-1.9998159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-0.1023275) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(-3.0854991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027095196) q[0];
sx q[0];
rz(-1.2074911) q[0];
sx q[0];
rz(-0.61371213) q[0];
x q[1];
rz(-0.34727879) q[2];
sx q[2];
rz(-2.0060853) q[2];
sx q[2];
rz(-1.6821282) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1069378) q[1];
sx q[1];
rz(-1.5562623) q[1];
sx q[1];
rz(-0.30277877) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0104996) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(-1.2922985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-2.1870959) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6765103) q[0];
sx q[0];
rz(-2.0970779) q[0];
sx q[0];
rz(0.17980534) q[0];
rz(-2.2013821) q[2];
sx q[2];
rz(-1.5988837) q[2];
sx q[2];
rz(1.1379776) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.33659014) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(0.57979433) q[1];
rz(-pi) q[2];
rz(2.9783863) q[3];
sx q[3];
rz(-1.2833793) q[3];
sx q[3];
rz(-0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(-2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(0.61002237) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9305206) q[0];
sx q[0];
rz(-1.0809582) q[0];
sx q[0];
rz(1.8563489) q[0];
rz(-pi) q[1];
rz(-1.7213983) q[2];
sx q[2];
rz(-2.4820231) q[2];
sx q[2];
rz(1.9285551) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2513189) q[1];
sx q[1];
rz(-2.2379025) q[1];
sx q[1];
rz(3.1097079) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026222762) q[0];
sx q[0];
rz(-1.5407469) q[0];
sx q[0];
rz(3.1307334) q[0];
rz(-1.2543711) q[2];
sx q[2];
rz(-0.3379312) q[2];
sx q[2];
rz(0.28883176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8229586) q[1];
sx q[1];
rz(-1.6726603) q[1];
sx q[1];
rz(-0.46551367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(-0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(0.7448147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574402) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-1.1293344) q[0];
rz(-1.1698193) q[2];
sx q[2];
rz(-2.3224761) q[2];
sx q[2];
rz(0.70105201) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71036584) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(-0.23780312) q[1];
rz(-pi) q[2];
rz(0.39448491) q[3];
sx q[3];
rz(-0.82291616) q[3];
sx q[3];
rz(-0.2713954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.9706479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1146961) q[0];
sx q[0];
rz(-2.7141502) q[0];
sx q[0];
rz(-3.0893185) q[0];
x q[1];
rz(2.3659336) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(1.0722216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78466258) q[1];
sx q[1];
rz(-1.5934048) q[1];
sx q[1];
rz(3.1329586) q[1];
rz(2.5578299) q[3];
sx q[3];
rz(-2.6423892) q[3];
sx q[3];
rz(-1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(-0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(0.70384937) q[2];
sx q[2];
rz(-0.8740295) q[2];
sx q[2];
rz(-0.340273) q[2];
rz(-1.9042653) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
