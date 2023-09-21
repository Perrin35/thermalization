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
rz(1.9265494) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6998866) q[0];
sx q[0];
rz(-0.21472782) q[0];
sx q[0];
rz(2.0798652) q[0];
x q[1];
rz(-1.8589852) q[2];
sx q[2];
rz(-0.93027861) q[2];
sx q[2];
rz(3.1079907) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4814264) q[1];
sx q[1];
rz(-2.5170442) q[1];
sx q[1];
rz(2.156483) q[1];
x q[2];
rz(-0.14903544) q[3];
sx q[3];
rz(-1.3081074) q[3];
sx q[3];
rz(-1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(-2.7637774) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29782444) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(-0.077117292) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(-1.5391301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7586655) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(1.4325607) q[0];
x q[1];
rz(1.7878754) q[2];
sx q[2];
rz(-0.84393822) q[2];
sx q[2];
rz(-0.94490563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(-2.944988) q[1];
x q[2];
rz(2.6529796) q[3];
sx q[3];
rz(-0.47227898) q[3];
sx q[3];
rz(-2.3924535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.845528) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.5862397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(-0.31366445) q[2];
sx q[2];
rz(-1.0278388) q[2];
sx q[2];
rz(1.107159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-2.4858027) q[1];
sx q[1];
rz(1.8379184) q[1];
rz(2.9402296) q[3];
sx q[3];
rz(-1.0346197) q[3];
sx q[3];
rz(2.686543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(1.8910485) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.282225) q[3];
sx q[3];
rz(1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164453) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(-1.8871334) q[0];
rz(1.9011263) q[2];
sx q[2];
rz(-0.9579881) q[2];
sx q[2];
rz(1.6883048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0604917) q[1];
sx q[1];
rz(-1.8685409) q[1];
sx q[1];
rz(2.9213419) q[1];
rz(-pi) q[2];
rz(0.051024036) q[3];
sx q[3];
rz(-1.050356) q[3];
sx q[3];
rz(2.737962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2531551) q[2];
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
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859905) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(-1.7657071) q[0];
rz(-1.8638523) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(3.0854991) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2984021) q[0];
sx q[0];
rz(-1.0023596) q[0];
sx q[0];
rz(-2.0060904) q[0];
rz(-0.34727879) q[2];
sx q[2];
rz(-2.0060853) q[2];
sx q[2];
rz(1.4594644) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1069378) q[1];
sx q[1];
rz(-1.5562623) q[1];
sx q[1];
rz(-0.30277877) q[1];
rz(-pi) q[2];
rz(-1.0104996) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(1.8492941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.00099480199) q[2];
sx q[2];
rz(-0.91789118) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(0.8941935) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(-0.64754852) q[0];
rz(-1.2619069) q[1];
sx q[1];
rz(-1.4636661) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6765103) q[0];
sx q[0];
rz(-1.0445147) q[0];
sx q[0];
rz(2.9617873) q[0];
x q[1];
rz(0.94021057) q[2];
sx q[2];
rz(-1.5427089) q[2];
sx q[2];
rz(2.0036151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94297385) q[1];
sx q[1];
rz(-1.0620411) q[1];
sx q[1];
rz(2.1224623) q[1];
rz(-pi) q[2];
rz(2.9783863) q[3];
sx q[3];
rz(-1.2833793) q[3];
sx q[3];
rz(-0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(2.0992289) q[2];
rz(-2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(-1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(-0.74321157) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-2.5315703) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22247032) q[0];
sx q[0];
rz(-1.8219935) q[0];
sx q[0];
rz(0.50719502) q[0];
x q[1];
rz(3.0257752) q[2];
sx q[2];
rz(-0.91997416) q[2];
sx q[2];
rz(-1.4027558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2513189) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(-0.031884738) q[1];
rz(-pi) q[2];
rz(-0.20603541) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(-2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-0.91840333) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(0.38890719) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(-0.73658529) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026222762) q[0];
sx q[0];
rz(-1.6008458) q[0];
sx q[0];
rz(-0.010859246) q[0];
rz(1.893115) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(-1.581574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8383933) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(1.4569015) q[1];
x q[2];
rz(2.0180118) q[3];
sx q[3];
rz(-1.2774602) q[3];
sx q[3];
rz(2.6597027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(1.3945403) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-2.396778) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4208922) q[0];
sx q[0];
rz(-1.8928796) q[0];
sx q[0];
rz(-2.3548404) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79297519) q[2];
sx q[2];
rz(-1.2816396) q[2];
sx q[2];
rz(0.58794978) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3639431) q[1];
sx q[1];
rz(-1.793982) q[1];
sx q[1];
rz(-1.2121483) q[1];
rz(-1.1780147) q[3];
sx q[3];
rz(-2.3142356) q[3];
sx q[3];
rz(2.3208997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(1.9469117) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(3.1104654) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(-1.1709447) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6452713) q[0];
sx q[0];
rz(-1.5924581) q[0];
sx q[0];
rz(-2.7146656) q[0];
rz(-pi) q[1];
rz(1.8413713) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(2.4547581) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78632894) q[1];
sx q[1];
rz(-1.5794282) q[1];
sx q[1];
rz(1.5481871) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5578299) q[3];
sx q[3];
rz(-2.6423892) q[3];
sx q[3];
rz(-2.0016058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(2.5496303) q[2];
rz(-2.5752318) q[3];
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
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
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
rz(1.5242143) q[3];
sx q[3];
rz(-2.8077879) q[3];
sx q[3];
rz(2.2624364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
