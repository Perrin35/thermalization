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
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5117447) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(1.758979) q[0];
rz(0.66081337) q[2];
sx q[2];
rz(-1.3408957) q[2];
sx q[2];
rz(-1.4290609) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3455968) q[1];
sx q[1];
rz(-1.0618292) q[1];
sx q[1];
rz(2.7624346) q[1];
rz(2.0753161) q[3];
sx q[3];
rz(-2.8404232) q[3];
sx q[3];
rz(1.0217713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-0.077117292) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.6024626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2167643) q[0];
sx q[0];
rz(-0.98470682) q[0];
sx q[0];
rz(0.092496471) q[0];
x q[1];
rz(1.3537172) q[2];
sx q[2];
rz(-2.2976544) q[2];
sx q[2];
rz(-0.94490563) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(-2.944988) q[1];
rz(-pi) q[2];
rz(2.7178571) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(0.37950294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(0.65845931) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-0.55535299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12268513) q[0];
sx q[0];
rz(-0.92120954) q[0];
sx q[0];
rz(2.5193549) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31366445) q[2];
sx q[2];
rz(-1.0278388) q[2];
sx q[2];
rz(-2.0344337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1173646) q[1];
sx q[1];
rz(-1.7324565) q[1];
sx q[1];
rz(2.2092186) q[1];
x q[2];
rz(-1.025612) q[3];
sx q[3];
rz(-1.3979997) q[3];
sx q[3];
rz(-1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73734036) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8811532) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-0.81480169) q[0];
rz(-1.762215) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(2.8864158) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53084757) q[0];
sx q[0];
rz(-1.2753715) q[0];
sx q[0];
rz(0.3770963) q[0];
rz(-pi) q[1];
rz(-2.7093676) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(1.9908817) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0604917) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(0.22025073) q[1];
rz(1.6595483) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(-2.6356217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(-2.9684084) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(0.056093562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984021) q[0];
sx q[0];
rz(-1.0023596) q[0];
sx q[0];
rz(-1.1355023) q[0];
x q[1];
rz(0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(0.75013559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5406815) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(1.5555698) q[1];
x q[2];
rz(2.131093) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(1.8492941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-0.43219217) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.4661219) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(0.64754852) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(-2.1870959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1178841) q[0];
sx q[0];
rz(-2.5881898) q[0];
sx q[0];
rz(-1.8694359) q[0];
rz(-pi) q[1];
rz(-0.034770413) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(-2.6882753) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8442321) q[1];
sx q[1];
rz(-2.4095222) q[1];
sx q[1];
rz(-2.3872603) q[1];
x q[2];
rz(1.8618705) q[3];
sx q[3];
rz(-1.7272514) q[3];
sx q[3];
rz(2.1978956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(0.61002237) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9191223) q[0];
sx q[0];
rz(-1.3195992) q[0];
sx q[0];
rz(2.6343976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91673135) q[2];
sx q[2];
rz(-1.6628633) q[2];
sx q[2];
rz(2.9031861) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2513189) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(0.031884738) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(-0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.6280744) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37305957) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(1.91747) q[0];
x q[1];
rz(1.2484776) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(-1.581574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8229586) q[1];
sx q[1];
rz(-1.6726603) q[1];
sx q[1];
rz(-2.676079) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96111091) q[3];
sx q[3];
rz(-0.52934066) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(0.7448147) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72070044) q[0];
sx q[0];
rz(-1.8928796) q[0];
sx q[0];
rz(-0.78675227) q[0];
x q[1];
rz(-0.3955598) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.8849444) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.77764952) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(-1.9294444) q[1];
x q[2];
rz(-2.7471077) q[3];
sx q[3];
rz(-0.82291616) q[3];
sx q[3];
rz(-0.2713954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.1709447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6452713) q[0];
sx q[0];
rz(-1.5924581) q[0];
sx q[0];
rz(-2.7146656) q[0];
x q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(-2.4547581) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3552637) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(1.5934056) q[1];
rz(-pi) q[2];
rz(0.42697866) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(-0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(-0.5919624) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82407172) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(-3.042165) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-2.2291017) q[2];
sx q[2];
rz(-0.9463263) q[2];
sx q[2];
rz(1.8778388) q[2];
rz(-1.6173784) q[3];
sx q[3];
rz(-2.8077879) q[3];
sx q[3];
rz(2.2624364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
