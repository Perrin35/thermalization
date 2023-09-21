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
rz(-2.1808704) q[0];
sx q[0];
rz(-1.3836432) q[0];
sx q[0];
rz(-3.0357009) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.8006969) q[2];
sx q[2];
rz(-1.7125318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3455968) q[1];
sx q[1];
rz(-1.0618292) q[1];
sx q[1];
rz(-0.3791581) q[1];
x q[2];
rz(1.8362942) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(-3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8437682) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(-1.5391301) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30277006) q[0];
sx q[0];
rz(-1.4937703) q[0];
sx q[0];
rz(-0.98273116) q[0];
rz(-pi) q[1];
rz(-2.402926) q[2];
sx q[2];
rz(-1.4091485) q[2];
sx q[2];
rz(2.6612298) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1748845) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(-2.944988) q[1];
x q[2];
rz(2.7178571) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(0.37950294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12268513) q[0];
sx q[0];
rz(-0.92120954) q[0];
sx q[0];
rz(-0.62223776) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.136134) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(-0.29758673) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4761423) q[1];
sx q[1];
rz(-2.1995771) q[1];
sx q[1];
rz(-2.9412342) q[1];
rz(1.2461353) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811532) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2164453) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(-1.8871334) q[0];
rz(-1.9011263) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(1.4532879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42923388) q[1];
sx q[1];
rz(-2.7731967) q[1];
sx q[1];
rz(2.1894987) q[1];
rz(-pi) q[2];
rz(-1.049794) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(1.9998159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010904) q[0];
sx q[0];
rz(-0.70106693) q[0];
sx q[0];
rz(-2.5581193) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34727879) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(1.4594644) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6009112) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(1.5860228) q[1];
rz(-0.089245307) q[3];
sx q[3];
rz(-2.1292994) q[3];
sx q[3];
rz(-0.23111471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(-0.8941935) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.6754707) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(0.64754852) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.944825) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(-2.1041656) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5231832) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(-0.39436755) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1986188) q[1];
sx q[1];
rz(-2.0795515) q[1];
sx q[1];
rz(2.1224623) q[1];
x q[2];
rz(-1.2797221) q[3];
sx q[3];
rz(-1.7272514) q[3];
sx q[3];
rz(2.1978956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2110721) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.8563489) q[0];
rz(-0.1158175) q[2];
sx q[2];
rz(-2.2216185) q[2];
sx q[2];
rz(1.4027558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1998132) q[1];
sx q[1];
rz(-0.66775125) q[1];
sx q[1];
rz(-1.6112531) q[1];
rz(1.1915019) q[3];
sx q[3];
rz(-1.3790352) q[3];
sx q[3];
rz(0.8134884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(2.7526855) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7685331) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(1.91747) q[0];
rz(0.10891624) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(0.045217302) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.089162) q[1];
sx q[1];
rz(-0.4757291) q[1];
sx q[1];
rz(2.9176941) q[1];
x q[2];
rz(2.1804817) q[3];
sx q[3];
rz(-0.52934066) q[3];
sx q[3];
rz(1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0344051) q[2];
sx q[2];
rz(-1.9341058) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(-1.9125787) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(1.3945403) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(2.9558682) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(2.396778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72070044) q[0];
sx q[0];
rz(-1.2487131) q[0];
sx q[0];
rz(-0.78675227) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3486175) q[2];
sx q[2];
rz(-1.2816396) q[2];
sx q[2];
rz(2.5536429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8150755) q[1];
sx q[1];
rz(-0.41985598) q[1];
sx q[1];
rz(0.99680568) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9635779) q[3];
sx q[3];
rz(-2.3142356) q[3];
sx q[3];
rz(0.82069293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(2.7375896) q[0];
rz(3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.9706479) q[1];
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
rz(1.3002214) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(-2.4547581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78632894) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(-1.5481871) q[1];
x q[2];
rz(0.42697866) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.70384937) q[2];
sx q[2];
rz(-2.2675632) q[2];
sx q[2];
rz(2.8013196) q[2];
rz(1.2373274) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
