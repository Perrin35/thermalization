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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1808704) q[0];
sx q[0];
rz(-1.7579494) q[0];
sx q[0];
rz(-3.0357009) q[0];
rz(-pi) q[1];
rz(1.8589852) q[2];
sx q[2];
rz(-2.211314) q[2];
sx q[2];
rz(-0.033601947) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.558555) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(-1.0298883) q[1];
rz(-pi) q[2];
rz(-1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(2.9585178) q[2];
rz(-2.7637774) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-0.29418501) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(1.6024626) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7586655) q[0];
sx q[0];
rz(-2.5490952) q[0];
sx q[0];
rz(1.7090319) q[0];
rz(-pi) q[1];
rz(1.7878754) q[2];
sx q[2];
rz(-0.84393822) q[2];
sx q[2];
rz(2.196687) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(-0.19660463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6529796) q[3];
sx q[3];
rz(-2.6693137) q[3];
sx q[3];
rz(-0.74913914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(-1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(2.5862397) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(2.32248) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0434415) q[2];
sx q[2];
rz(-0.61908365) q[2];
sx q[2];
rz(-1.4738136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.33298102) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(1.3036742) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9402296) q[3];
sx q[3];
rz(-2.106973) q[3];
sx q[3];
rz(-0.45504967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92514738) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(1.2544592) q[0];
rz(0.43222506) q[2];
sx q[2];
rz(-0.6859633) q[2];
sx q[2];
rz(-1.150711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55528044) q[1];
sx q[1];
rz(-1.781207) q[1];
sx q[1];
rz(-1.8754688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6595483) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(-2.6356217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(-0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(3.0392652) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027095196) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(-2.5278805) q[0];
rz(-pi) q[1];
rz(-0.9390097) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(2.3914571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6519421) q[1];
sx q[1];
rz(-2.838476) q[1];
sx q[1];
rz(-3.0928844) q[1];
rz(-pi) q[2];
rz(-3.0523473) q[3];
sx q[3];
rz(-1.0122932) q[3];
sx q[3];
rz(-0.23111471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(-0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(2.1870959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6765103) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(-2.9617873) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5231832) q[2];
sx q[2];
rz(-2.5104668) q[2];
sx q[2];
rz(2.7472251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8050025) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(0.57979433) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8618705) q[3];
sx q[3];
rz(-1.7272514) q[3];
sx q[3];
rz(0.94369704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.59297562) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(1.0423638) q[2];
rz(0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(0.74321157) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(-2.5315703) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9305206) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(-1.2852438) q[0];
rz(-pi) q[1];
rz(-0.91673135) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(-2.9031861) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1998132) q[1];
sx q[1];
rz(-2.4738414) q[1];
sx q[1];
rz(-1.5303395) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20603541) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(2.2231893) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(1.6280744) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(2.4050074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442473) q[0];
sx q[0];
rz(-1.5816507) q[0];
sx q[0];
rz(1.5407451) q[0];
x q[1];
rz(3.0326764) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(3.0963754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8229586) q[1];
sx q[1];
rz(-1.4689323) q[1];
sx q[1];
rz(-2.676079) q[1];
rz(-pi) q[2];
rz(-0.96111091) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(-1.6314268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(1.9125787) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
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
rz(-1.5988388) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(2.396778) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9841524) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(2.0122583) q[0];
x q[1];
rz(2.7460329) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.8849444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4312268) q[1];
sx q[1];
rz(-1.2214298) q[1];
sx q[1];
rz(0.23780312) q[1];
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
rz(pi/2) q[1];
sx q[1];
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
x q[2];
sx q[3];
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
rz(-2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(2.7375896) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.9706479) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02689657) q[0];
sx q[0];
rz(-0.42744246) q[0];
sx q[0];
rz(-3.0893185) q[0];
x q[1];
rz(-1.3002214) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(2.4547581) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78466258) q[1];
sx q[1];
rz(-1.5934048) q[1];
sx q[1];
rz(0.008634062) q[1];
x q[2];
rz(1.2788494) q[3];
sx q[3];
rz(-1.9817096) q[3];
sx q[3];
rz(-1.3565855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(0.56636089) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(-1.6177572) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82407172) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(-0.099427632) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(0.70384937) q[2];
sx q[2];
rz(-0.8740295) q[2];
sx q[2];
rz(-0.340273) q[2];
rz(1.6173784) q[3];
sx q[3];
rz(-0.33380476) q[3];
sx q[3];
rz(-0.87915626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
