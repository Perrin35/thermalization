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
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44170609) q[0];
sx q[0];
rz(-2.9268648) q[0];
sx q[0];
rz(-1.0617274) q[0];
rz(-pi) q[1];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.8006969) q[2];
sx q[2];
rz(-1.7125318) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3455968) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(0.3791581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14903544) q[3];
sx q[3];
rz(-1.8334853) q[3];
sx q[3];
rz(-1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(0.29418501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(0.077117292) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3829271) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(-1.7090319) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23761959) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(-1.2653637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6290366) q[1];
sx q[1];
rz(-1.7450383) q[1];
sx q[1];
rz(2.0577355) q[1];
rz(2.6529796) q[3];
sx q[3];
rz(-0.47227898) q[3];
sx q[3];
rz(0.74913914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.8639996) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-1.0226117) q[3];
sx q[3];
rz(0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.2116218) q[1];
sx q[1];
rz(2.5862397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(2.32248) q[0];
rz(-pi) q[1];
rz(-2.136134) q[2];
sx q[2];
rz(-1.8381422) q[2];
sx q[2];
rz(2.8440059) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33298102) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.3036742) q[1];
x q[2];
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
rz(-pi) q[1];
rz(-2.4042523) q[2];
sx q[2];
rz(-2.3534687) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.6916493) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(2.326791) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164453) q[0];
sx q[0];
rz(-1.2108004) q[0];
sx q[0];
rz(1.2544592) q[0];
x q[1];
rz(-2.7093676) q[2];
sx q[2];
rz(-2.4556293) q[2];
sx q[2];
rz(-1.9908817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42923388) q[1];
sx q[1];
rz(-2.7731967) q[1];
sx q[1];
rz(2.1894987) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.049794) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(-1.1417768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(3.0392652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.859905) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.7657071) q[0];
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
rz(0.027095196) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(0.61371213) q[0];
rz(-pi) q[1];
rz(-0.34727879) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(-1.4594644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4896506) q[1];
sx q[1];
rz(-0.30311668) q[1];
sx q[1];
rz(-3.0928844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.089245307) q[3];
sx q[3];
rz(-1.0122932) q[3];
sx q[3];
rz(-2.9104779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(-2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(1.6754707) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11319259) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(2.4940441) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(2.1870959) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19676767) q[0];
sx q[0];
rz(-1.7260572) q[0];
sx q[0];
rz(1.0374271) q[0];
rz(-1.5231832) q[2];
sx q[2];
rz(-0.63112586) q[2];
sx q[2];
rz(2.7472251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8050025) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(-0.57979433) q[1];
rz(2.9783863) q[3];
sx q[3];
rz(-1.8582134) q[3];
sx q[3];
rz(0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.548617) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-2.0992289) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-2.091566) q[3];
sx q[3];
rz(1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2838659) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(0.74321157) q[0];
rz(-1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(0.61002237) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689432) q[0];
sx q[0];
rz(-0.56108755) q[0];
sx q[0];
rz(0.48604301) q[0];
rz(-pi) q[1];
rz(1.7213983) q[2];
sx q[2];
rz(-0.65956958) q[2];
sx q[2];
rz(-1.2130376) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89027379) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(-3.1097079) q[1];
rz(-pi) q[2];
rz(0.20603541) q[3];
sx q[3];
rz(-1.9427951) q[3];
sx q[3];
rz(-2.4601065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(0.91840333) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(0.73658529) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026222762) q[0];
sx q[0];
rz(-1.5407469) q[0];
sx q[0];
rz(0.010859246) q[0];
rz(-pi) q[1];
rz(1.893115) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(-1.581574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(-0.22389852) q[1];
rz(-pi) q[2];
rz(-2.8183476) q[3];
sx q[3];
rz(-1.1439699) q[3];
sx q[3];
rz(-0.95110287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(-2.9558682) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-2.396778) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574402) q[0];
sx q[0];
rz(-2.3072349) q[0];
sx q[0];
rz(1.1293344) q[0];
rz(-pi) q[1];
rz(1.9717734) q[2];
sx q[2];
rz(-2.3224761) q[2];
sx q[2];
rz(0.70105201) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3265171) q[1];
sx q[1];
rz(-0.41985598) q[1];
sx q[1];
rz(-2.144787) q[1];
x q[2];
rz(2.3585988) q[3];
sx q[3];
rz(-1.2851614) q[3];
sx q[3];
rz(2.1180958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.9469117) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-0.99635807) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(2.7375896) q[0];
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
rz(-1.5491345) q[0];
sx q[0];
rz(-0.42692703) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3002214) q[2];
sx q[2];
rz(-0.8136533) q[2];
sx q[2];
rz(-0.68683456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78632894) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(-1.5934056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58376273) q[3];
sx q[3];
rz(-2.6423892) q[3];
sx q[3];
rz(1.1399869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-0.5919624) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-0.16470328) q[3];
sx q[3];
rz(1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-2.4377433) q[2];
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
