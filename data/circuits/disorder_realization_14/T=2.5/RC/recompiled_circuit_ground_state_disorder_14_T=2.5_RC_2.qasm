OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.1190342) q[0];
sx q[0];
rz(-2.0329539) q[0];
sx q[0];
rz(-0.57060567) q[0];
rz(1.624149) q[1];
sx q[1];
rz(-0.94944209) q[1];
sx q[1];
rz(1.5784664) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9304912) q[0];
sx q[0];
rz(-1.6864713) q[0];
sx q[0];
rz(1.132908) q[0];
x q[1];
rz(-2.3307595) q[2];
sx q[2];
rz(-1.8301415) q[2];
sx q[2];
rz(-1.3400067) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7540463) q[1];
sx q[1];
rz(-1.7472194) q[1];
sx q[1];
rz(-2.7642718) q[1];
rz(-pi) q[2];
rz(2.6622979) q[3];
sx q[3];
rz(-0.66761469) q[3];
sx q[3];
rz(-1.3541753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40330517) q[2];
sx q[2];
rz(-0.83633509) q[2];
sx q[2];
rz(-2.9264012) q[2];
rz(0.0067986851) q[3];
sx q[3];
rz(-2.2577622) q[3];
sx q[3];
rz(-0.96338898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(2.6007518) q[0];
sx q[0];
rz(-2.231926) q[0];
sx q[0];
rz(-1.2700861) q[0];
rz(-1.0441682) q[1];
sx q[1];
rz(-2.7296598) q[1];
sx q[1];
rz(2.5792436) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017496271) q[0];
sx q[0];
rz(-1.3096333) q[0];
sx q[0];
rz(-3.0246252) q[0];
rz(-1.1435881) q[2];
sx q[2];
rz(-1.0173804) q[2];
sx q[2];
rz(1.1737905) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96404379) q[1];
sx q[1];
rz(-2.4052087) q[1];
sx q[1];
rz(-0.089929637) q[1];
rz(2.8329912) q[3];
sx q[3];
rz(-1.2720593) q[3];
sx q[3];
rz(-1.5200391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1336141) q[2];
sx q[2];
rz(-2.9972711) q[2];
sx q[2];
rz(-0.14145429) q[2];
rz(-2.134363) q[3];
sx q[3];
rz(-1.4814582) q[3];
sx q[3];
rz(-0.13834113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7815642) q[0];
sx q[0];
rz(-1.6369632) q[0];
sx q[0];
rz(0.052852782) q[0];
rz(-1.7618893) q[1];
sx q[1];
rz(-2.3638937) q[1];
sx q[1];
rz(-2.8676829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6099458) q[0];
sx q[0];
rz(-1.393712) q[0];
sx q[0];
rz(-3.1041932) q[0];
rz(0.81676151) q[2];
sx q[2];
rz(-1.5026983) q[2];
sx q[2];
rz(0.20055972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63800784) q[1];
sx q[1];
rz(-1.5440739) q[1];
sx q[1];
rz(0.77049945) q[1];
rz(-pi) q[2];
rz(-0.11577167) q[3];
sx q[3];
rz(-1.8325509) q[3];
sx q[3];
rz(-3.0512471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86458796) q[2];
sx q[2];
rz(-1.4952345) q[2];
sx q[2];
rz(0.60025364) q[2];
rz(-2.4872335) q[3];
sx q[3];
rz(-2.8221966) q[3];
sx q[3];
rz(-1.4734242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.3510503) q[0];
sx q[0];
rz(-1.496614) q[0];
sx q[0];
rz(-0.45781621) q[0];
rz(-0.59902016) q[1];
sx q[1];
rz(-0.51282239) q[1];
sx q[1];
rz(-2.8694428) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3783589) q[0];
sx q[0];
rz(-2.1714179) q[0];
sx q[0];
rz(-1.2637704) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6105) q[2];
sx q[2];
rz(-2.0390687) q[2];
sx q[2];
rz(0.35781579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49644955) q[1];
sx q[1];
rz(-1.6890059) q[1];
sx q[1];
rz(-0.18564863) q[1];
rz(0.64025819) q[3];
sx q[3];
rz(-1.7096412) q[3];
sx q[3];
rz(-1.1306896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.51364556) q[2];
sx q[2];
rz(-0.38868913) q[2];
sx q[2];
rz(1.2316616) q[2];
rz(-1.7197459) q[3];
sx q[3];
rz(-1.809241) q[3];
sx q[3];
rz(-0.96760145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.642639) q[0];
sx q[0];
rz(-0.068878219) q[0];
sx q[0];
rz(-2.4054085) q[0];
rz(0.63756293) q[1];
sx q[1];
rz(-1.2022737) q[1];
sx q[1];
rz(0.41086248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7124568) q[0];
sx q[0];
rz(-0.68198689) q[0];
sx q[0];
rz(0.76759591) q[0];
x q[1];
rz(-1.246894) q[2];
sx q[2];
rz(-2.5482087) q[2];
sx q[2];
rz(1.3932799) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.73130979) q[1];
sx q[1];
rz(-1.5640794) q[1];
sx q[1];
rz(2.5120458) q[1];
rz(-0.75848364) q[3];
sx q[3];
rz(-2.3173769) q[3];
sx q[3];
rz(0.099140204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8132482) q[2];
sx q[2];
rz(-1.5466377) q[2];
sx q[2];
rz(2.6206139) q[2];
rz(-3.0327435) q[3];
sx q[3];
rz(-2.1284926) q[3];
sx q[3];
rz(2.525135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2206889) q[0];
sx q[0];
rz(-1.9158798) q[0];
sx q[0];
rz(-1.2589681) q[0];
rz(2.068223) q[1];
sx q[1];
rz(-1.6683038) q[1];
sx q[1];
rz(-2.4980256) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94144017) q[0];
sx q[0];
rz(-2.6905031) q[0];
sx q[0];
rz(3.1300225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3544092) q[2];
sx q[2];
rz(-1.2020276) q[2];
sx q[2];
rz(2.5390196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7593708) q[1];
sx q[1];
rz(-0.77414162) q[1];
sx q[1];
rz(1.8060807) q[1];
rz(2.4197641) q[3];
sx q[3];
rz(-2.4284644) q[3];
sx q[3];
rz(-1.6588746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7966938) q[2];
sx q[2];
rz(-2.0774136) q[2];
sx q[2];
rz(-1.3382443) q[2];
rz(-1.4517339) q[3];
sx q[3];
rz(-1.2029878) q[3];
sx q[3];
rz(-0.10045997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.79784262) q[0];
sx q[0];
rz(-1.8666973) q[0];
sx q[0];
rz(-0.95329681) q[0];
rz(-0.15996179) q[1];
sx q[1];
rz(-2.5123031) q[1];
sx q[1];
rz(0.86872855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4788268) q[0];
sx q[0];
rz(-1.7870169) q[0];
sx q[0];
rz(-0.89333306) q[0];
rz(0.73448278) q[2];
sx q[2];
rz(-0.96436635) q[2];
sx q[2];
rz(-2.7766428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27426961) q[1];
sx q[1];
rz(-1.549701) q[1];
sx q[1];
rz(2.8460851) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89959156) q[3];
sx q[3];
rz(-1.7760522) q[3];
sx q[3];
rz(1.156607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4538883) q[2];
sx q[2];
rz(-1.1331465) q[2];
sx q[2];
rz(-2.647661) q[2];
rz(1.4880873) q[3];
sx q[3];
rz(-2.2736214) q[3];
sx q[3];
rz(0.021763703) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22773753) q[0];
sx q[0];
rz(-1.4233002) q[0];
sx q[0];
rz(-0.29105759) q[0];
rz(0.13045467) q[1];
sx q[1];
rz(-2.7619402) q[1];
sx q[1];
rz(-1.5694654) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30877844) q[0];
sx q[0];
rz(-0.052292112) q[0];
sx q[0];
rz(-1.3175658) q[0];
rz(-pi) q[1];
x q[1];
rz(1.356691) q[2];
sx q[2];
rz(-2.4567025) q[2];
sx q[2];
rz(-2.8321526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2184184) q[1];
sx q[1];
rz(-1.6863135) q[1];
sx q[1];
rz(0.61489132) q[1];
rz(0.30867851) q[3];
sx q[3];
rz(-1.4490845) q[3];
sx q[3];
rz(0.50060035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36720413) q[2];
sx q[2];
rz(-2.9640894) q[2];
sx q[2];
rz(0.11167488) q[2];
rz(-0.81237927) q[3];
sx q[3];
rz(-1.8492536) q[3];
sx q[3];
rz(-1.6151379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3647301) q[0];
sx q[0];
rz(-2.9193945) q[0];
sx q[0];
rz(0.36866933) q[0];
rz(-0.61019507) q[1];
sx q[1];
rz(-1.5851603) q[1];
sx q[1];
rz(2.1048224) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42175779) q[0];
sx q[0];
rz(-1.0204691) q[0];
sx q[0];
rz(1.3070413) q[0];
x q[1];
rz(3.1366903) q[2];
sx q[2];
rz(-2.0837726) q[2];
sx q[2];
rz(0.42286985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6994991) q[1];
sx q[1];
rz(-1.7101062) q[1];
sx q[1];
rz(-2.3087835) q[1];
rz(-pi) q[2];
rz(1.0843363) q[3];
sx q[3];
rz(-1.4245302) q[3];
sx q[3];
rz(-1.4649966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6250352) q[2];
sx q[2];
rz(-2.3380307) q[2];
sx q[2];
rz(2.6885923) q[2];
rz(0.9520483) q[3];
sx q[3];
rz(-2.0249849) q[3];
sx q[3];
rz(3.0208352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9348738) q[0];
sx q[0];
rz(-3.0937338) q[0];
sx q[0];
rz(-2.5913443) q[0];
rz(-1.0847367) q[1];
sx q[1];
rz(-2.5526498) q[1];
sx q[1];
rz(-1.750607) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15875951) q[0];
sx q[0];
rz(-1.5313494) q[0];
sx q[0];
rz(-1.0692014) q[0];
x q[1];
rz(2.9736217) q[2];
sx q[2];
rz(-1.0412626) q[2];
sx q[2];
rz(-2.2574097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.7470737) q[1];
sx q[1];
rz(-2.2929774) q[1];
sx q[1];
rz(2.4492521) q[1];
rz(-pi) q[2];
x q[2];
rz(3.080322) q[3];
sx q[3];
rz(-1.6271504) q[3];
sx q[3];
rz(2.0298093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37798205) q[2];
sx q[2];
rz(-0.72145975) q[2];
sx q[2];
rz(-1.1627496) q[2];
rz(-2.5871596) q[3];
sx q[3];
rz(-0.88397637) q[3];
sx q[3];
rz(2.673431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32600317) q[0];
sx q[0];
rz(-1.8392039) q[0];
sx q[0];
rz(-1.4890672) q[0];
rz(2.9266657) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(-2.8561572) q[2];
sx q[2];
rz(-0.68842059) q[2];
sx q[2];
rz(-0.60981087) q[2];
rz(-1.636408) q[3];
sx q[3];
rz(-0.8498522) q[3];
sx q[3];
rz(-0.1838799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
