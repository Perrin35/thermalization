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
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62984798) q[0];
sx q[0];
rz(-1.6748322) q[0];
sx q[0];
rz(1.3826136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66081337) q[2];
sx q[2];
rz(-1.8006969) q[2];
sx q[2];
rz(1.7125318) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.79599586) q[1];
sx q[1];
rz(-1.0618292) q[1];
sx q[1];
rz(-0.3791581) q[1];
rz(0.14903544) q[3];
sx q[3];
rz(-1.3081074) q[3];
sx q[3];
rz(-1.5955773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-2.9585178) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(3.0644754) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(-1.5391301) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9248283) q[0];
sx q[0];
rz(-2.1568858) q[0];
sx q[0];
rz(-0.092496471) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23761959) q[2];
sx q[2];
rz(-0.75287205) q[2];
sx q[2];
rz(-1.876229) q[2];
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
rz(1.8061403) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(-1.8542765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(-2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(-1.9494879) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-2.5862397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(-0.81911266) q[0];
rz(-pi) q[1];
rz(1.0054587) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(-0.29758673) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33298102) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(-1.3036742) q[1];
rz(-pi) q[2];
rz(2.9402296) q[3];
sx q[3];
rz(-2.106973) q[3];
sx q[3];
rz(0.45504967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(-0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-2.6840211) q[0];
sx q[0];
rz(0.81480169) q[0];
rz(-1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(0.25517685) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53084757) q[0];
sx q[0];
rz(-1.8662211) q[0];
sx q[0];
rz(-0.3770963) q[0];
rz(-pi) q[1];
rz(1.2404664) q[2];
sx q[2];
rz(-2.1836046) q[2];
sx q[2];
rz(-1.4532879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.55528044) q[1];
sx q[1];
rz(-1.781207) q[1];
sx q[1];
rz(-1.2661238) q[1];
rz(-pi) q[2];
rz(1.6595483) q[3];
sx q[3];
rz(-0.52270652) q[3];
sx q[3];
rz(0.50597092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(2.9684084) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(-0.056093562) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1306886) q[0];
sx q[0];
rz(-2.4405257) q[0];
sx q[0];
rz(0.58347337) q[0];
x q[1];
rz(0.9390097) q[2];
sx q[2];
rz(-2.5917705) q[2];
sx q[2];
rz(-0.75013559) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.034654831) q[1];
sx q[1];
rz(-1.5853303) q[1];
sx q[1];
rz(-2.8388139) q[1];
rz(0.089245307) q[3];
sx q[3];
rz(-2.1292994) q[3];
sx q[3];
rz(-2.9104779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(0.43219217) q[2];
rz(0.8941935) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.4661219) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(-0.64754852) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19676767) q[0];
sx q[0];
rz(-1.7260572) q[0];
sx q[0];
rz(1.0374271) q[0];
x q[1];
rz(3.1068222) q[2];
sx q[2];
rz(-0.94049847) q[2];
sx q[2];
rz(2.6882753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1986188) q[1];
sx q[1];
rz(-2.0795515) q[1];
sx q[1];
rz(-1.0191304) q[1];
rz(-pi) q[2];
rz(-1.0682085) q[3];
sx q[3];
rz(-2.8121901) q[3];
sx q[3];
rz(-2.9941032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59297562) q[2];
sx q[2];
rz(-1.9081215) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-0.61002237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689432) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(-2.6555496) q[0];
x q[1];
rz(-0.1158175) q[2];
sx q[2];
rz(-0.91997416) q[2];
sx q[2];
rz(-1.4027558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89027379) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(0.031884738) q[1];
x q[2];
rz(-1.9500908) q[3];
sx q[3];
rz(-1.3790352) q[3];
sx q[3];
rz(0.8134884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8053749) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-0.91840333) q[2];
rz(1.5911128) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-2.7526855) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36088762) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(-2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-0.73658529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5973454) q[0];
sx q[0];
rz(-1.559942) q[0];
sx q[0];
rz(1.5407451) q[0];
rz(-pi) q[1];
rz(1.8872216) q[2];
sx q[2];
rz(-2.8036615) q[2];
sx q[2];
rz(2.8527609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0524307) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(-0.22389852) q[1];
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
rz(1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-2.396778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841524) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(1.1293344) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79297519) q[2];
sx q[2];
rz(-1.8599531) q[2];
sx q[2];
rz(2.5536429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3265171) q[1];
sx q[1];
rz(-2.7217367) q[1];
sx q[1];
rz(0.99680568) q[1];
rz(2.3585988) q[3];
sx q[3];
rz(-1.8564312) q[3];
sx q[3];
rz(-2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.194681) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-0.031127302) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.9706479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084328018) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(-1.5469993) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-0.8136533) q[2];
sx q[2];
rz(-0.68683456) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7217692) q[1];
sx q[1];
rz(-3.1173919) q[1];
sx q[1];
rz(-1.2060549) q[1];
rz(-1.2788494) q[3];
sx q[3];
rz(-1.1598831) q[3];
sx q[3];
rz(1.7850072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4460454) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(0.56636089) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82407172) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(2.4026985) q[2];
sx q[2];
rz(-2.0901491) q[2];
sx q[2];
rz(0.73170589) q[2];
rz(-3.1254461) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
