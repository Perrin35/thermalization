OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.69831508) q[0];
sx q[0];
rz(-0.3787711) q[0];
sx q[0];
rz(-1.9842499) q[0];
rz(-2.3565489) q[1];
sx q[1];
rz(-2.4554689) q[1];
sx q[1];
rz(-2.6148028) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136381) q[0];
sx q[0];
rz(-1.0109954) q[0];
sx q[0];
rz(2.7877121) q[0];
rz(-2.4073462) q[2];
sx q[2];
rz(-0.99319211) q[2];
sx q[2];
rz(-2.2478888) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4545645) q[1];
sx q[1];
rz(-2.6978081) q[1];
sx q[1];
rz(-0.084974809) q[1];
x q[2];
rz(1.0105114) q[3];
sx q[3];
rz(-1.5454195) q[3];
sx q[3];
rz(2.199774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36246768) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(-2.989952) q[2];
rz(-0.14196299) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4769984) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(2.3240996) q[0];
rz(3.0290161) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(2.2419825) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47682677) q[0];
sx q[0];
rz(-1.810084) q[0];
sx q[0];
rz(0.28061687) q[0];
x q[1];
rz(0.67018761) q[2];
sx q[2];
rz(-1.0600344) q[2];
sx q[2];
rz(2.0005884) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5499076) q[1];
sx q[1];
rz(-2.3909759) q[1];
sx q[1];
rz(-0.082521497) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51606222) q[3];
sx q[3];
rz(-2.1806697) q[3];
sx q[3];
rz(0.42760951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67109913) q[2];
sx q[2];
rz(-1.495139) q[2];
sx q[2];
rz(1.0279083) q[2];
rz(-0.73728621) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(-3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0055493) q[0];
sx q[0];
rz(-0.5558973) q[0];
sx q[0];
rz(1.1935724) q[0];
rz(-1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(1.4847635) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5201621) q[0];
sx q[0];
rz(-1.5957766) q[0];
sx q[0];
rz(-1.6170184) q[0];
rz(-pi) q[1];
rz(1.3593) q[2];
sx q[2];
rz(-1.6944318) q[2];
sx q[2];
rz(-2.5619626) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.130467) q[1];
sx q[1];
rz(-2.7724206) q[1];
sx q[1];
rz(-2.7586731) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.082529222) q[3];
sx q[3];
rz(-1.9325053) q[3];
sx q[3];
rz(-1.495273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9803311) q[2];
sx q[2];
rz(-1.223246) q[2];
sx q[2];
rz(-0.70183357) q[2];
rz(3.0724604) q[3];
sx q[3];
rz(-0.88874236) q[3];
sx q[3];
rz(-0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0858916) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(-0.62537801) q[0];
rz(-2.0471795) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(1.8249003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863281) q[0];
sx q[0];
rz(-1.0287971) q[0];
sx q[0];
rz(-2.9744451) q[0];
rz(2.8595692) q[2];
sx q[2];
rz(-2.8057753) q[2];
sx q[2];
rz(1.5329602) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3299371) q[1];
sx q[1];
rz(-1.7509169) q[1];
sx q[1];
rz(0.7094769) q[1];
rz(-pi) q[2];
rz(-0.98424498) q[3];
sx q[3];
rz(-0.86185019) q[3];
sx q[3];
rz(1.1701442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-2.4989765) q[2];
sx q[2];
rz(0.030108062) q[2];
rz(0.41664577) q[3];
sx q[3];
rz(-1.4285587) q[3];
sx q[3];
rz(3.0863975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3629214) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(1.2177421) q[0];
rz(-0.22449224) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(0.40103689) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3844135) q[0];
sx q[0];
rz(-1.9862439) q[0];
sx q[0];
rz(0.40796221) q[0];
rz(-0.65481477) q[2];
sx q[2];
rz(-1.8698955) q[2];
sx q[2];
rz(1.4465743) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68919888) q[1];
sx q[1];
rz(-1.3340063) q[1];
sx q[1];
rz(-1.522032) q[1];
x q[2];
rz(2.5878536) q[3];
sx q[3];
rz(-2.9512292) q[3];
sx q[3];
rz(-2.7680264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84247056) q[2];
sx q[2];
rz(-1.2229908) q[2];
sx q[2];
rz(-2.6341338) q[2];
rz(2.2211645) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(2.9612655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4949263) q[0];
sx q[0];
rz(-0.05519069) q[0];
sx q[0];
rz(-2.1221509) q[0];
rz(1.1722209) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(-1.2219465) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9393765) q[0];
sx q[0];
rz(-0.76986137) q[0];
sx q[0];
rz(-1.059993) q[0];
rz(0.035313531) q[2];
sx q[2];
rz(-2.0853634) q[2];
sx q[2];
rz(-0.071431486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5355365) q[1];
sx q[1];
rz(-2.5185985) q[1];
sx q[1];
rz(2.6583903) q[1];
rz(-2.0133408) q[3];
sx q[3];
rz(-2.5199515) q[3];
sx q[3];
rz(-2.3592202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6601861) q[2];
sx q[2];
rz(-0.6654827) q[2];
sx q[2];
rz(-2.0484203) q[2];
rz(2.6045351) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(0.66948906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2375803) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(1.681666) q[0];
rz(-1.6319252) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(-0.07930886) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8312165) q[0];
sx q[0];
rz(-2.8908911) q[0];
sx q[0];
rz(-1.901907) q[0];
rz(-pi) q[1];
rz(0.20360501) q[2];
sx q[2];
rz(-1.4362659) q[2];
sx q[2];
rz(0.30483887) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.053959) q[1];
sx q[1];
rz(-1.988639) q[1];
sx q[1];
rz(-1.2705951) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9856624) q[3];
sx q[3];
rz(-2.4178388) q[3];
sx q[3];
rz(-2.2035905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1377533) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(-2.7286781) q[2];
rz(-1.8462935) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(0.41954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9846648) q[0];
sx q[0];
rz(-2.4913737) q[0];
sx q[0];
rz(-0.28537634) q[0];
rz(-3.0889619) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(2.9579128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12626064) q[0];
sx q[0];
rz(-1.6063326) q[0];
sx q[0];
rz(-1.4592341) q[0];
x q[1];
rz(0.48641522) q[2];
sx q[2];
rz(-2.8210495) q[2];
sx q[2];
rz(-0.17580168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0252776) q[1];
sx q[1];
rz(-2.5417788) q[1];
sx q[1];
rz(-0.13529899) q[1];
rz(-0.64818188) q[3];
sx q[3];
rz(-1.8055918) q[3];
sx q[3];
rz(-2.8051493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24667428) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(2.0764009) q[2];
rz(1.4554321) q[3];
sx q[3];
rz(-1.8543517) q[3];
sx q[3];
rz(-2.5559032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0378549) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(-1.6023585) q[0];
rz(-2.1922951) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.7112214) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95887881) q[0];
sx q[0];
rz(-0.74107594) q[0];
sx q[0];
rz(2.4070021) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6550421) q[2];
sx q[2];
rz(-0.73053321) q[2];
sx q[2];
rz(2.4075395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0956968) q[1];
sx q[1];
rz(-1.2751082) q[1];
sx q[1];
rz(1.5398094) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1422906) q[3];
sx q[3];
rz(-1.3125889) q[3];
sx q[3];
rz(0.22120295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0091693) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(0.12651786) q[2];
rz(-2.8070519) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(-0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3009406) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(2.6244923) q[0];
rz(-3.0866947) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(0.089769207) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9835998) q[0];
sx q[0];
rz(-1.1674034) q[0];
sx q[0];
rz(2.1703815) q[0];
rz(-pi) q[1];
rz(-1.4710861) q[2];
sx q[2];
rz(-0.56312984) q[2];
sx q[2];
rz(-0.3143087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43347142) q[1];
sx q[1];
rz(-0.29954391) q[1];
sx q[1];
rz(-0.56510651) q[1];
x q[2];
rz(2.91342) q[3];
sx q[3];
rz(-1.2163278) q[3];
sx q[3];
rz(-1.6443524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6101997) q[2];
sx q[2];
rz(-2.0072082) q[2];
sx q[2];
rz(-0.72719491) q[2];
rz(0.7555035) q[3];
sx q[3];
rz(-2.754039) q[3];
sx q[3];
rz(0.21272794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1995734) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(1.749281) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(-0.10689312) q[2];
sx q[2];
rz(-0.86915599) q[2];
sx q[2];
rz(-0.55021777) q[2];
rz(-2.200442) q[3];
sx q[3];
rz(-1.0761906) q[3];
sx q[3];
rz(-3.0039345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
