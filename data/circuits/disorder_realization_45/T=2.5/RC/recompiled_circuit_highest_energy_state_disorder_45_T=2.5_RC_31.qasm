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
rz(-0.51195872) q[0];
sx q[0];
rz(-2.7380138) q[0];
sx q[0];
rz(3.1374748) q[0];
rz(-2.4601958) q[1];
sx q[1];
rz(-1.8713142) q[1];
sx q[1];
rz(-1.4407925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82872921) q[0];
sx q[0];
rz(-1.848319) q[0];
sx q[0];
rz(-2.9636431) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75136331) q[2];
sx q[2];
rz(-1.8247424) q[2];
sx q[2];
rz(2.6264409) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.550168) q[1];
sx q[1];
rz(-1.9240849) q[1];
sx q[1];
rz(-2.4387906) q[1];
x q[2];
rz(0.012083455) q[3];
sx q[3];
rz(-1.3971431) q[3];
sx q[3];
rz(-2.1326667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9422841) q[2];
sx q[2];
rz(-1.0465304) q[2];
sx q[2];
rz(-0.47041565) q[2];
rz(-2.2453902) q[3];
sx q[3];
rz(-1.6737409) q[3];
sx q[3];
rz(1.5907653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6386221) q[0];
sx q[0];
rz(-2.6574385) q[0];
sx q[0];
rz(2.7726987) q[0];
rz(2.6781354) q[1];
sx q[1];
rz(-1.2838485) q[1];
sx q[1];
rz(-2.8772112) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437816) q[0];
sx q[0];
rz(-2.4759295) q[0];
sx q[0];
rz(-1.3206318) q[0];
x q[1];
rz(2.5892841) q[2];
sx q[2];
rz(-2.8495516) q[2];
sx q[2];
rz(2.8949182) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7355613) q[1];
sx q[1];
rz(-1.8623061) q[1];
sx q[1];
rz(-0.31707615) q[1];
x q[2];
rz(1.6912429) q[3];
sx q[3];
rz(-0.79211881) q[3];
sx q[3];
rz(-1.219238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5276864) q[2];
sx q[2];
rz(-0.1656342) q[2];
sx q[2];
rz(1.3045093) q[2];
rz(2.8218609) q[3];
sx q[3];
rz(-0.78800646) q[3];
sx q[3];
rz(0.50362292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17427915) q[0];
sx q[0];
rz(-0.19198424) q[0];
sx q[0];
rz(-2.018003) q[0];
rz(-2.7536821) q[1];
sx q[1];
rz(-1.265641) q[1];
sx q[1];
rz(0.33043114) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215805) q[0];
sx q[0];
rz(-2.8224484) q[0];
sx q[0];
rz(-2.1681251) q[0];
x q[1];
rz(-3.0082012) q[2];
sx q[2];
rz(-0.61163227) q[2];
sx q[2];
rz(2.4470714) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36281113) q[1];
sx q[1];
rz(-1.8900202) q[1];
sx q[1];
rz(0.37324639) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64066842) q[3];
sx q[3];
rz(-2.0694149) q[3];
sx q[3];
rz(-2.9151268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.419751) q[2];
sx q[2];
rz(-1.9998877) q[2];
sx q[2];
rz(-1.9640131) q[2];
rz(1.8925331) q[3];
sx q[3];
rz(-1.306059) q[3];
sx q[3];
rz(-0.038399847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053442001) q[0];
sx q[0];
rz(-1.7120687) q[0];
sx q[0];
rz(0.088726774) q[0];
rz(-0.42452043) q[1];
sx q[1];
rz(-2.4219234) q[1];
sx q[1];
rz(-1.4094062) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8363894) q[0];
sx q[0];
rz(-0.38581784) q[0];
sx q[0];
rz(-1.057339) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1698977) q[2];
sx q[2];
rz(-1.9847775) q[2];
sx q[2];
rz(-1.5091015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3864435) q[1];
sx q[1];
rz(-2.1514822) q[1];
sx q[1];
rz(-2.0232852) q[1];
x q[2];
rz(-0.32281422) q[3];
sx q[3];
rz(-2.0695311) q[3];
sx q[3];
rz(0.43969595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38922629) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(-2.7092095) q[2];
rz(2.6835175) q[3];
sx q[3];
rz(-1.4832486) q[3];
sx q[3];
rz(-2.1336011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78419375) q[0];
sx q[0];
rz(-2.9965239) q[0];
sx q[0];
rz(-0.31132895) q[0];
rz(-1.1773102) q[1];
sx q[1];
rz(-1.9639587) q[1];
sx q[1];
rz(2.6572773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70894079) q[0];
sx q[0];
rz(-1.6659032) q[0];
sx q[0];
rz(1.4115788) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2782477) q[2];
sx q[2];
rz(-1.496409) q[2];
sx q[2];
rz(1.1400346) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6020585) q[1];
sx q[1];
rz(-1.9440398) q[1];
sx q[1];
rz(0.15284027) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9919162) q[3];
sx q[3];
rz(-1.6972741) q[3];
sx q[3];
rz(1.5095131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8208661) q[2];
sx q[2];
rz(-1.7118688) q[2];
sx q[2];
rz(0.42496625) q[2];
rz(3.037437) q[3];
sx q[3];
rz(-2.7724373) q[3];
sx q[3];
rz(-2.4764376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.7729618) q[0];
sx q[0];
rz(-0.63987982) q[0];
sx q[0];
rz(0.1567008) q[0];
rz(-2.8796097) q[1];
sx q[1];
rz(-1.9888839) q[1];
sx q[1];
rz(1.461747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91297594) q[0];
sx q[0];
rz(-1.3973087) q[0];
sx q[0];
rz(-0.96488379) q[0];
x q[1];
rz(-0.22905519) q[2];
sx q[2];
rz(-2.2543) q[2];
sx q[2];
rz(2.7102269) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9380381) q[1];
sx q[1];
rz(-2.9239836) q[1];
sx q[1];
rz(-0.59534351) q[1];
rz(-pi) q[2];
rz(-0.42243345) q[3];
sx q[3];
rz(-2.3452873) q[3];
sx q[3];
rz(-0.56047201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4713952) q[2];
sx q[2];
rz(-1.6452226) q[2];
sx q[2];
rz(0.71693286) q[2];
rz(2.9251621) q[3];
sx q[3];
rz(-1.2856893) q[3];
sx q[3];
rz(-1.9560248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7316498) q[0];
sx q[0];
rz(-0.2922295) q[0];
sx q[0];
rz(1.2662079) q[0];
rz(0.75824291) q[1];
sx q[1];
rz(-1.1781324) q[1];
sx q[1];
rz(-2.0936802) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3558725) q[0];
sx q[0];
rz(-1.2730755) q[0];
sx q[0];
rz(-0.96426086) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4518544) q[2];
sx q[2];
rz(-2.0819217) q[2];
sx q[2];
rz(3.1060807) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0430773) q[1];
sx q[1];
rz(-1.9334353) q[1];
sx q[1];
rz(-3.0671544) q[1];
x q[2];
rz(-2.6664958) q[3];
sx q[3];
rz(-0.49834278) q[3];
sx q[3];
rz(-1.4628177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8083814) q[2];
sx q[2];
rz(-1.0071249) q[2];
sx q[2];
rz(-3.132931) q[2];
rz(2.1374785) q[3];
sx q[3];
rz(-2.6486371) q[3];
sx q[3];
rz(-2.139411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36326161) q[0];
sx q[0];
rz(-2.9982428) q[0];
sx q[0];
rz(-2.1891201) q[0];
rz(1.6078423) q[1];
sx q[1];
rz(-1.2865103) q[1];
sx q[1];
rz(1.0362157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0506315) q[0];
sx q[0];
rz(-2.1078175) q[0];
sx q[0];
rz(2.0799252) q[0];
rz(-pi) q[1];
rz(0.47889809) q[2];
sx q[2];
rz(-0.69456813) q[2];
sx q[2];
rz(-0.50146363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.10401512) q[1];
sx q[1];
rz(-2.0682145) q[1];
sx q[1];
rz(0.90105547) q[1];
x q[2];
rz(0.50004543) q[3];
sx q[3];
rz(-2.3498021) q[3];
sx q[3];
rz(-2.4337976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2967534) q[2];
sx q[2];
rz(-0.71594683) q[2];
sx q[2];
rz(-2.1412264) q[2];
rz(1.2760466) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(1.067151) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412096) q[0];
sx q[0];
rz(-2.2333133) q[0];
sx q[0];
rz(2.8874183) q[0];
rz(1.8036448) q[1];
sx q[1];
rz(-0.86054069) q[1];
sx q[1];
rz(-0.77802229) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36150071) q[0];
sx q[0];
rz(-0.92467148) q[0];
sx q[0];
rz(-1.9801937) q[0];
x q[1];
rz(-2.4619589) q[2];
sx q[2];
rz(-0.25654116) q[2];
sx q[2];
rz(-2.7367965) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5601156) q[1];
sx q[1];
rz(-0.72598493) q[1];
sx q[1];
rz(-2.5037161) q[1];
x q[2];
rz(-1.9970787) q[3];
sx q[3];
rz(-0.23565764) q[3];
sx q[3];
rz(2.896275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0622056) q[2];
sx q[2];
rz(-1.5176682) q[2];
sx q[2];
rz(3.0089231) q[2];
rz(2.0343871) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(1.4440822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4001813) q[0];
sx q[0];
rz(-1.1962698) q[0];
sx q[0];
rz(-0.74914002) q[0];
rz(0.49508849) q[1];
sx q[1];
rz(-1.9063213) q[1];
sx q[1];
rz(1.3912158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86319727) q[0];
sx q[0];
rz(-0.94587284) q[0];
sx q[0];
rz(1.004659) q[0];
x q[1];
rz(-2.7374908) q[2];
sx q[2];
rz(-1.1617336) q[2];
sx q[2];
rz(2.7210195) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0616845) q[1];
sx q[1];
rz(-2.7322953) q[1];
sx q[1];
rz(3.0319935) q[1];
x q[2];
rz(2.6364348) q[3];
sx q[3];
rz(-1.0978175) q[3];
sx q[3];
rz(1.5468171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5007925) q[2];
sx q[2];
rz(-1.4747138) q[2];
sx q[2];
rz(2.2643209) q[2];
rz(-2.6194465) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(1.6845901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7397639) q[0];
sx q[0];
rz(-0.54720989) q[0];
sx q[0];
rz(-1.3450958) q[0];
rz(1.2783891) q[1];
sx q[1];
rz(-2.8291193) q[1];
sx q[1];
rz(3.1183174) q[1];
rz(0.042365778) q[2];
sx q[2];
rz(-0.63127098) q[2];
sx q[2];
rz(1.9628738) q[2];
rz(1.9587026) q[3];
sx q[3];
rz(-1.1454875) q[3];
sx q[3];
rz(-2.3844908) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
