OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0225585) q[0];
sx q[0];
rz(2.0329539) q[0];
sx q[0];
rz(11.995765) q[0];
rz(-1.5174436) q[1];
sx q[1];
rz(-2.1921506) q[1];
sx q[1];
rz(1.5631262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9304912) q[0];
sx q[0];
rz(-1.6864713) q[0];
sx q[0];
rz(2.0086847) q[0];
rz(1.9384266) q[2];
sx q[2];
rz(-2.3470633) q[2];
sx q[2];
rz(-0.49434915) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38754639) q[1];
sx q[1];
rz(-1.7472194) q[1];
sx q[1];
rz(-0.37732084) q[1];
rz(-pi) q[2];
rz(1.9195032) q[3];
sx q[3];
rz(-0.98920663) q[3];
sx q[3];
rz(-0.76954704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7382875) q[2];
sx q[2];
rz(-2.3052576) q[2];
sx q[2];
rz(-2.9264012) q[2];
rz(-3.134794) q[3];
sx q[3];
rz(-0.88383043) q[3];
sx q[3];
rz(0.96338898) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5408408) q[0];
sx q[0];
rz(-0.90966666) q[0];
sx q[0];
rz(-1.2700861) q[0];
rz(-2.0974244) q[1];
sx q[1];
rz(-2.7296598) q[1];
sx q[1];
rz(-2.5792436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017496271) q[0];
sx q[0];
rz(-1.3096333) q[0];
sx q[0];
rz(3.0246252) q[0];
rz(-pi) q[1];
rz(-2.5508443) q[2];
sx q[2];
rz(-2.4563654) q[2];
sx q[2];
rz(0.45999377) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6015502) q[1];
sx q[1];
rz(-1.6311495) q[1];
sx q[1];
rz(-2.4072231) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8834194) q[3];
sx q[3];
rz(-1.2762831) q[3];
sx q[3];
rz(-0.14430639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0079785138) q[2];
sx q[2];
rz(-0.14432159) q[2];
sx q[2];
rz(-3.0001384) q[2];
rz(2.134363) q[3];
sx q[3];
rz(-1.4814582) q[3];
sx q[3];
rz(-3.0032515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7815642) q[0];
sx q[0];
rz(-1.5046295) q[0];
sx q[0];
rz(-0.052852782) q[0];
rz(1.3797034) q[1];
sx q[1];
rz(-2.3638937) q[1];
sx q[1];
rz(-2.8676829) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045740728) q[0];
sx q[0];
rz(-1.533982) q[0];
sx q[0];
rz(1.748002) q[0];
x q[1];
rz(-1.4714965) q[2];
sx q[2];
rz(-0.75649951) q[2];
sx q[2];
rz(1.4425636) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5035848) q[1];
sx q[1];
rz(-1.5975188) q[1];
sx q[1];
rz(-0.77049945) q[1];
rz(-3.025821) q[3];
sx q[3];
rz(-1.3090418) q[3];
sx q[3];
rz(-3.0512471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86458796) q[2];
sx q[2];
rz(-1.6463582) q[2];
sx q[2];
rz(0.60025364) q[2];
rz(-2.4872335) q[3];
sx q[3];
rz(-0.31939605) q[3];
sx q[3];
rz(-1.6681685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3510503) q[0];
sx q[0];
rz(-1.6449787) q[0];
sx q[0];
rz(-2.6837764) q[0];
rz(-0.59902016) q[1];
sx q[1];
rz(-0.51282239) q[1];
sx q[1];
rz(0.27214989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7632338) q[0];
sx q[0];
rz(-0.97017479) q[0];
sx q[0];
rz(1.8778223) q[0];
x q[1];
rz(-2.3568677) q[2];
sx q[2];
rz(-0.69284495) q[2];
sx q[2];
rz(1.8681519) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49644955) q[1];
sx q[1];
rz(-1.6890059) q[1];
sx q[1];
rz(-2.955944) q[1];
rz(1.3982716) q[3];
sx q[3];
rz(-2.2039045) q[3];
sx q[3];
rz(-2.5987491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6279471) q[2];
sx q[2];
rz(-2.7529035) q[2];
sx q[2];
rz(-1.909931) q[2];
rz(-1.7197459) q[3];
sx q[3];
rz(-1.809241) q[3];
sx q[3];
rz(-0.96760145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.642639) q[0];
sx q[0];
rz(-0.068878219) q[0];
sx q[0];
rz(0.73618412) q[0];
rz(0.63756293) q[1];
sx q[1];
rz(-1.2022737) q[1];
sx q[1];
rz(0.41086248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4291358) q[0];
sx q[0];
rz(-0.68198689) q[0];
sx q[0];
rz(0.76759591) q[0];
rz(-2.9301398) q[2];
sx q[2];
rz(-2.1295068) q[2];
sx q[2];
rz(1.3635456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3113271) q[1];
sx q[1];
rz(-2.5120148) q[1];
sx q[1];
rz(0.011407995) q[1];
x q[2];
rz(2.383109) q[3];
sx q[3];
rz(-2.3173769) q[3];
sx q[3];
rz(0.099140204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32834443) q[2];
sx q[2];
rz(-1.594955) q[2];
sx q[2];
rz(-2.6206139) q[2];
rz(-3.0327435) q[3];
sx q[3];
rz(-2.1284926) q[3];
sx q[3];
rz(-0.61645761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2206889) q[0];
sx q[0];
rz(-1.2257129) q[0];
sx q[0];
rz(-1.8826245) q[0];
rz(2.068223) q[1];
sx q[1];
rz(-1.4732889) q[1];
sx q[1];
rz(2.4980256) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2001525) q[0];
sx q[0];
rz(-2.6905031) q[0];
sx q[0];
rz(3.1300225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7648206) q[2];
sx q[2];
rz(-1.7724282) q[2];
sx q[2];
rz(-2.0942935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3583329) q[1];
sx q[1];
rz(-1.734501) q[1];
sx q[1];
rz(2.3309784) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72182853) q[3];
sx q[3];
rz(-2.4284644) q[3];
sx q[3];
rz(-1.4827181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7966938) q[2];
sx q[2];
rz(-1.0641791) q[2];
sx q[2];
rz(-1.8033484) q[2];
rz(-1.6898588) q[3];
sx q[3];
rz(-1.9386049) q[3];
sx q[3];
rz(3.0411327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79784262) q[0];
sx q[0];
rz(-1.2748953) q[0];
sx q[0];
rz(0.95329681) q[0];
rz(-0.15996179) q[1];
sx q[1];
rz(-0.62928951) q[1];
sx q[1];
rz(2.2728641) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9729444) q[0];
sx q[0];
rz(-2.4356844) q[0];
sx q[0];
rz(1.9078518) q[0];
x q[1];
rz(-2.4071099) q[2];
sx q[2];
rz(-2.1772263) q[2];
sx q[2];
rz(-0.36494985) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.27426961) q[1];
sx q[1];
rz(-1.549701) q[1];
sx q[1];
rz(0.29550754) q[1];
rz(1.8938167) q[3];
sx q[3];
rz(-2.4443808) q[3];
sx q[3];
rz(2.4761971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4538883) q[2];
sx q[2];
rz(-2.0084461) q[2];
sx q[2];
rz(-0.49393168) q[2];
rz(-1.6535053) q[3];
sx q[3];
rz(-0.86797124) q[3];
sx q[3];
rz(3.119829) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9138551) q[0];
sx q[0];
rz(-1.4233002) q[0];
sx q[0];
rz(-0.29105759) q[0];
rz(0.13045467) q[1];
sx q[1];
rz(-0.37965241) q[1];
sx q[1];
rz(1.5694654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5792517) q[0];
sx q[0];
rz(-1.6214193) q[0];
sx q[0];
rz(-0.013112092) q[0];
rz(0.89719551) q[2];
sx q[2];
rz(-1.7056124) q[2];
sx q[2];
rz(1.7133985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92317428) q[1];
sx q[1];
rz(-1.6863135) q[1];
sx q[1];
rz(0.61489132) q[1];
rz(-1.4431105) q[3];
sx q[3];
rz(-1.2644759) q[3];
sx q[3];
rz(1.1088913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36720413) q[2];
sx q[2];
rz(-2.9640894) q[2];
sx q[2];
rz(3.0299178) q[2];
rz(-2.3292134) q[3];
sx q[3];
rz(-1.8492536) q[3];
sx q[3];
rz(-1.5264548) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3647301) q[0];
sx q[0];
rz(-2.9193945) q[0];
sx q[0];
rz(0.36866933) q[0];
rz(-2.5313976) q[1];
sx q[1];
rz(-1.5851603) q[1];
sx q[1];
rz(-2.1048224) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198349) q[0];
sx q[0];
rz(-1.0204691) q[0];
sx q[0];
rz(1.8345513) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0837777) q[2];
sx q[2];
rz(-1.566525) q[2];
sx q[2];
rz(1.1503324) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2543176) q[1];
sx q[1];
rz(-0.84158449) q[1];
sx q[1];
rz(-0.18730728) q[1];
rz(-pi) q[2];
rz(0.16513326) q[3];
sx q[3];
rz(-2.0516178) q[3];
sx q[3];
rz(2.9588678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5165575) q[2];
sx q[2];
rz(-0.80356193) q[2];
sx q[2];
rz(0.45300031) q[2];
rz(-0.9520483) q[3];
sx q[3];
rz(-2.0249849) q[3];
sx q[3];
rz(0.12075748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20671885) q[0];
sx q[0];
rz(-3.0937338) q[0];
sx q[0];
rz(-2.5913443) q[0];
rz(2.0568559) q[1];
sx q[1];
rz(-0.58894283) q[1];
sx q[1];
rz(-1.3909856) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3904141) q[0];
sx q[0];
rz(-1.0696279) q[0];
sx q[0];
rz(-3.0966109) q[0];
rz(-1.2925658) q[2];
sx q[2];
rz(-2.5884853) q[2];
sx q[2];
rz(1.9335374) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3251614) q[1];
sx q[1];
rz(-2.0702989) q[1];
sx q[1];
rz(-2.4234467) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3971294) q[3];
sx q[3];
rz(-3.0583707) q[3];
sx q[3];
rz(-1.2017488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37798205) q[2];
sx q[2];
rz(-0.72145975) q[2];
sx q[2];
rz(1.9788431) q[2];
rz(-0.55443305) q[3];
sx q[3];
rz(-0.88397637) q[3];
sx q[3];
rz(0.46816167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.32600317) q[0];
sx q[0];
rz(-1.8392039) q[0];
sx q[0];
rz(-1.4890672) q[0];
rz(-0.21492699) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(1.3431637) q[2];
sx q[2];
rz(-0.91522436) q[2];
sx q[2];
rz(2.1686423) q[2];
rz(3.0671185) q[3];
sx q[3];
rz(-2.4182033) q[3];
sx q[3];
rz(-0.084666336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
