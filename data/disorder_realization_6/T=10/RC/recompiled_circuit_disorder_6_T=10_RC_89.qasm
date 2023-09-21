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
x q[1];
rz(1.8589852) q[2];
sx q[2];
rz(-2.211314) q[2];
sx q[2];
rz(-0.033601947) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.558555) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(-2.1117044) q[1];
x q[2];
rz(2.0753161) q[3];
sx q[3];
rz(-0.30116943) q[3];
sx q[3];
rz(-1.0217713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(-2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-3.0644754) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(1.6024626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30277006) q[0];
sx q[0];
rz(-1.4937703) q[0];
sx q[0];
rz(2.1588615) q[0];
rz(-pi) q[1];
rz(-2.402926) q[2];
sx q[2];
rz(-1.7324442) q[2];
sx q[2];
rz(-2.6612298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6290366) q[1];
sx q[1];
rz(-1.3965544) q[1];
sx q[1];
rz(-2.0577355) q[1];
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
x q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(2.7084896) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(0.55535299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0189075) q[0];
sx q[0];
rz(-0.92120954) q[0];
sx q[0];
rz(-2.5193549) q[0];
rz(-pi) q[1];
rz(-0.31366445) q[2];
sx q[2];
rz(-2.1137538) q[2];
sx q[2];
rz(-1.107159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33298102) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(1.3036742) q[1];
x q[2];
rz(-1.2461353) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(2.897443) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8811532) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92514738) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(-1.2544592) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9011263) q[2];
sx q[2];
rz(-2.1836046) q[2];
sx q[2];
rz(1.6883048) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0604917) q[1];
sx q[1];
rz(-1.8685409) q[1];
sx q[1];
rz(0.22025073) q[1];
rz(1.6595483) q[3];
sx q[3];
rz(-2.6188861) q[3];
sx q[3];
rz(2.6356217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(2.611768) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(0.056093562) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.010904) q[0];
sx q[0];
rz(-0.70106693) q[0];
sx q[0];
rz(2.5581193) q[0];
x q[1];
rz(-2.7943139) q[2];
sx q[2];
rz(-2.0060853) q[2];
sx q[2];
rz(1.6821282) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1069378) q[1];
sx q[1];
rz(-1.5853303) q[1];
sx q[1];
rz(2.8388139) q[1];
rz(1.7124743) q[3];
sx q[3];
rz(-2.5767527) q[3];
sx q[3];
rz(2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1405979) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.6754707) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(2.1870959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6765103) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(2.9617873) q[0];
rz(-pi) q[1];
rz(3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(0.45331732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33659014) q[1];
sx q[1];
rz(-2.046236) q[1];
sx q[1];
rz(2.5617983) q[1];
rz(0.16320634) q[3];
sx q[3];
rz(-1.2833793) q[3];
sx q[3];
rz(-2.4678469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59297562) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(1.0423638) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(-0.61002237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9191223) q[0];
sx q[0];
rz(-1.3195992) q[0];
sx q[0];
rz(-2.6343976) q[0];
x q[1];
rz(-0.91673135) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(0.23840657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1998132) q[1];
sx q[1];
rz(-0.66775125) q[1];
sx q[1];
rz(-1.6112531) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7685331) q[0];
sx q[0];
rz(-0.031950843) q[0];
sx q[0];
rz(1.91747) q[0];
x q[1];
rz(-1.2484776) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(1.5600187) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.089162) q[1];
sx q[1];
rz(-0.4757291) q[1];
sx q[1];
rz(-0.22389852) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1235808) q[3];
sx q[3];
rz(-1.8641324) q[3];
sx q[3];
rz(0.48188996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0344051) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-2.396778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1574402) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(2.0122583) q[0];
rz(0.3955598) q[2];
sx q[2];
rz(-0.8330847) q[2];
sx q[2];
rz(-1.8849444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4312268) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(2.9037895) q[1];
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
rz(-2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(2.1452346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.9706479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4963213) q[0];
sx q[0];
rz(-1.5491345) q[0];
sx q[0];
rz(0.42692703) q[0];
rz(-1.8413713) q[2];
sx q[2];
rz(-0.8136533) q[2];
sx q[2];
rz(-0.68683456) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7217692) q[1];
sx q[1];
rz(-0.024200736) q[1];
sx q[1];
rz(-1.9355378) q[1];
rz(2.5578299) q[3];
sx q[3];
rz(-2.6423892) q[3];
sx q[3];
rz(-1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(-3.042165) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-2.4026985) q[2];
sx q[2];
rz(-1.0514435) q[2];
sx q[2];
rz(-2.4098868) q[2];
rz(-1.2373274) q[3];
sx q[3];
rz(-1.5555391) q[3];
sx q[3];
rz(-2.4939668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];