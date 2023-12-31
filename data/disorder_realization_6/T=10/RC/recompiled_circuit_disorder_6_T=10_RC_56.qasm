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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6998866) q[0];
sx q[0];
rz(-2.9268648) q[0];
sx q[0];
rz(-1.0617274) q[0];
rz(-2.7772929) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(0.42687624) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4814264) q[1];
sx q[1];
rz(-2.5170442) q[1];
sx q[1];
rz(0.98510965) q[1];
rz(-pi) q[2];
rz(0.14903544) q[3];
sx q[3];
rz(-1.3081074) q[3];
sx q[3];
rz(1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(2.9585178) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(-2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(1.6024626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7586655) q[0];
sx q[0];
rz(-2.5490952) q[0];
sx q[0];
rz(1.4325607) q[0];
rz(-pi) q[1];
x q[1];
rz(2.402926) q[2];
sx q[2];
rz(-1.4091485) q[2];
sx q[2];
rz(-2.6612298) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-1.0918573) q[1];
sx q[1];
rz(2.944988) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8061403) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(1.8542765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2960647) q[2];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(-0.55535299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99291891) q[0];
sx q[0];
rz(-0.86704555) q[0];
sx q[0];
rz(-0.91627319) q[0];
x q[1];
rz(-0.31366445) q[2];
sx q[2];
rz(-1.0278388) q[2];
sx q[2];
rz(-2.0344337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(-1.3036742) q[1];
rz(-pi) q[2];
rz(2.1159806) q[3];
sx q[3];
rz(-1.7435929) q[3];
sx q[3];
rz(1.9219414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-1.4677306) q[0];
sx q[0];
rz(-0.47463372) q[0];
sx q[0];
rz(-2.4509096) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5023979) q[2];
sx q[2];
rz(-1.3022458) q[2];
sx q[2];
rz(0.07721363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42923388) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(-2.1894987) q[1];
rz(-1.049794) q[3];
sx q[3];
rz(-1.61506) q[3];
sx q[3];
rz(1.1417768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8884376) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-2.9684084) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-2.9960222) q[3];
sx q[3];
rz(3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8431906) q[0];
sx q[0];
rz(-1.0023596) q[0];
sx q[0];
rz(1.1355023) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7943139) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(-1.6821282) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4896506) q[1];
sx q[1];
rz(-0.30311668) q[1];
sx q[1];
rz(-3.0928844) q[1];
rz(-3.0523473) q[3];
sx q[3];
rz(-2.1292994) q[3];
sx q[3];
rz(0.23111471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1405979) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(0.64754852) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(0.9544968) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19676767) q[0];
sx q[0];
rz(-1.7260572) q[0];
sx q[0];
rz(1.0374271) q[0];
x q[1];
rz(-3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(-0.45331732) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2973605) q[1];
sx q[1];
rz(-0.73207049) q[1];
sx q[1];
rz(-0.75433235) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0733842) q[3];
sx q[3];
rz(-0.32940255) q[3];
sx q[3];
rz(-0.1474895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(1.0423638) q[2];
rz(0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(-1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-2.3983811) q[0];
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
rz(1.9305206) q[0];
sx q[0];
rz(-2.0606344) q[0];
sx q[0];
rz(1.2852438) q[0];
rz(-0.91673135) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(0.23840657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4808018) q[1];
sx q[1];
rz(-1.5958438) q[1];
sx q[1];
rz(-2.2381496) q[1];
x q[2];
rz(-2.0537297) q[3];
sx q[3];
rz(-2.7186839) q[3];
sx q[3];
rz(1.203323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(1.5504799) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(1.5135182) q[0];
rz(0.52945119) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-0.73658529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7685331) q[0];
sx q[0];
rz(-0.031950843) q[0];
sx q[0];
rz(-1.91747) q[0];
rz(-pi) q[1];
rz(1.2543711) q[2];
sx q[2];
rz(-0.3379312) q[2];
sx q[2];
rz(-0.28883176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8383933) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(1.6846912) q[1];
x q[2];
rz(1.1235808) q[3];
sx q[3];
rz(-1.2774602) q[3];
sx q[3];
rz(0.48188996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0344051) q[2];
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
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443611) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-2.396778) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841524) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-1.1293344) q[0];
rz(1.9717734) q[2];
sx q[2];
rz(-0.81911659) q[2];
sx q[2];
rz(-0.70105201) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4312268) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(0.23780312) q[1];
rz(-pi) q[2];
rz(-1.9635779) q[3];
sx q[3];
rz(-2.3142356) q[3];
sx q[3];
rz(-2.3208997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(2.1452346) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74334082) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(3.1104654) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.9706479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084328018) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(1.5945934) q[0];
rz(-pi) q[1];
rz(0.2756341) q[2];
sx q[2];
rz(-2.3466957) q[2];
sx q[2];
rz(-2.8385712) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3569301) q[1];
sx q[1];
rz(-1.5481879) q[1];
sx q[1];
rz(-3.1329586) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5578299) q[3];
sx q[3];
rz(-0.49920344) q[3];
sx q[3];
rz(1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(0.56636089) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.099427632) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(2.4377433) q[2];
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
