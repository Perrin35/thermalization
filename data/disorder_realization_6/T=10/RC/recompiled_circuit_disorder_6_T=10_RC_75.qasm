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
rz(-2.5807227) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(-1.2150432) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96072223) q[0];
sx q[0];
rz(-1.3836432) q[0];
sx q[0];
rz(-0.10589177) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7772929) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(2.7147164) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3455968) q[1];
sx q[1];
rz(-1.0618292) q[1];
sx q[1];
rz(-0.3791581) q[1];
rz(-pi) q[2];
rz(-1.8362942) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(-0.063751566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9248283) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(-0.092496471) q[0];
rz(2.9039731) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(-1.876229) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51255608) q[1];
sx q[1];
rz(-1.3965544) q[1];
sx q[1];
rz(1.0838572) q[1];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(1.8542765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(0.55535299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1028324) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(-2.32248) q[0];
x q[1];
rz(2.136134) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(0.29758673) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4761423) q[1];
sx q[1];
rz(-2.1995771) q[1];
sx q[1];
rz(2.9412342) q[1];
rz(2.9402296) q[3];
sx q[3];
rz(-2.106973) q[3];
sx q[3];
rz(0.45504967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.2505442) q[2];
rz(2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811532) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(2.8864158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92514738) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(1.8871334) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63919477) q[2];
sx q[2];
rz(-1.3022458) q[2];
sx q[2];
rz(-3.064379) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.081101) q[1];
sx q[1];
rz(-1.8685409) q[1];
sx q[1];
rz(0.22025073) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0917986) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(-1.9998159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(2.611768) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027095196) q[0];
sx q[0];
rz(-1.2074911) q[0];
sx q[0];
rz(-0.61371213) q[0];
rz(-0.34727879) q[2];
sx q[2];
rz(-2.0060853) q[2];
sx q[2];
rz(-1.6821282) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.5860228) q[1];
rz(-1.0104996) q[3];
sx q[3];
rz(-1.6464525) q[3];
sx q[3];
rz(1.2922985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-0.43219217) q[2];
rz(-0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11319259) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(2.1870959) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0237085) q[0];
sx q[0];
rz(-0.55340289) q[0];
sx q[0];
rz(1.8694359) q[0];
rz(-pi) q[1];
rz(0.94021057) q[2];
sx q[2];
rz(-1.5427089) q[2];
sx q[2];
rz(-1.1379776) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.33659014) q[1];
sx q[1];
rz(-2.046236) q[1];
sx q[1];
rz(-0.57979433) q[1];
rz(-pi) q[2];
rz(1.0682085) q[3];
sx q[3];
rz(-0.32940255) q[3];
sx q[3];
rz(0.1474895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-0.61002237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3726495) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(-2.6555496) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1158175) q[2];
sx q[2];
rz(-2.2216185) q[2];
sx q[2];
rz(-1.4027558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.94177946) q[1];
sx q[1];
rz(-2.4738414) q[1];
sx q[1];
rz(1.6112531) q[1];
rz(-pi) q[2];
rz(-2.0537297) q[3];
sx q[3];
rz(-0.42290877) q[3];
sx q[3];
rz(1.9382697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-2.4050074) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5442473) q[0];
sx q[0];
rz(-1.559942) q[0];
sx q[0];
rz(1.6008475) q[0];
x q[1];
rz(1.2484776) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(1.581574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(-0.22389852) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32324507) q[3];
sx q[3];
rz(-1.1439699) q[3];
sx q[3];
rz(0.95110287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(2.4576808) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(0.7448147) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72070044) q[0];
sx q[0];
rz(-1.8928796) q[0];
sx q[0];
rz(0.78675227) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3955598) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(1.2566483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8150755) q[1];
sx q[1];
rz(-2.7217367) q[1];
sx q[1];
rz(-2.144787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9635779) q[3];
sx q[3];
rz(-0.82735705) q[3];
sx q[3];
rz(-0.82069293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1054489) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-2.1452346) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.9706479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0572646) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(-1.5469993) q[0];
rz(-2.8659586) q[2];
sx q[2];
rz(-0.79489691) q[2];
sx q[2];
rz(-0.3030215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3569301) q[1];
sx q[1];
rz(-1.5481879) q[1];
sx q[1];
rz(-3.1329586) q[1];
x q[2];
rz(-2.714614) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(-3.0468575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
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
