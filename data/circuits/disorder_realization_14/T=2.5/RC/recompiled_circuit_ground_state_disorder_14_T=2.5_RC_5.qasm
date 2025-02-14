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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9304912) q[0];
sx q[0];
rz(-1.6864713) q[0];
sx q[0];
rz(-1.132908) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9384266) q[2];
sx q[2];
rz(-0.7945294) q[2];
sx q[2];
rz(0.49434915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5414428) q[1];
sx q[1];
rz(-2.7268616) q[1];
sx q[1];
rz(2.6909237) q[1];
rz(-pi) q[2];
rz(1.2220895) q[3];
sx q[3];
rz(-0.98920663) q[3];
sx q[3];
rz(0.76954704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40330517) q[2];
sx q[2];
rz(-2.3052576) q[2];
sx q[2];
rz(0.21519145) q[2];
rz(3.134794) q[3];
sx q[3];
rz(-2.2577622) q[3];
sx q[3];
rz(-2.1782037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6007518) q[0];
sx q[0];
rz(-2.231926) q[0];
sx q[0];
rz(-1.2700861) q[0];
rz(2.0974244) q[1];
sx q[1];
rz(-2.7296598) q[1];
sx q[1];
rz(-0.56234908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320181) q[0];
sx q[0];
rz(-0.28561297) q[0];
sx q[0];
rz(1.1591039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9980045) q[2];
sx q[2];
rz(-2.1242122) q[2];
sx q[2];
rz(1.1737905) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1775489) q[1];
sx q[1];
rz(-2.4052087) q[1];
sx q[1];
rz(0.089929637) q[1];
rz(-1.8834194) q[3];
sx q[3];
rz(-1.2762831) q[3];
sx q[3];
rz(2.9972863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.0079785138) q[2];
sx q[2];
rz(-0.14432159) q[2];
sx q[2];
rz(-0.14145429) q[2];
rz(-1.0072297) q[3];
sx q[3];
rz(-1.6601345) q[3];
sx q[3];
rz(3.0032515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36002845) q[0];
sx q[0];
rz(-1.5046295) q[0];
sx q[0];
rz(-3.0887399) q[0];
rz(-1.3797034) q[1];
sx q[1];
rz(-2.3638937) q[1];
sx q[1];
rz(-0.27390972) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045740728) q[0];
sx q[0];
rz(-1.533982) q[0];
sx q[0];
rz(-1.3935907) q[0];
x q[1];
rz(-1.6700961) q[2];
sx q[2];
rz(-0.75649951) q[2];
sx q[2];
rz(-1.4425636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63800784) q[1];
sx q[1];
rz(-1.5440739) q[1];
sx q[1];
rz(0.77049945) q[1];
rz(3.025821) q[3];
sx q[3];
rz(-1.8325509) q[3];
sx q[3];
rz(-3.0512471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2770047) q[2];
sx q[2];
rz(-1.4952345) q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.6837764) q[0];
rz(0.59902016) q[1];
sx q[1];
rz(-2.6287703) q[1];
sx q[1];
rz(0.27214989) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7632338) q[0];
sx q[0];
rz(-2.1714179) q[0];
sx q[0];
rz(1.8778223) q[0];
x q[1];
rz(2.3568677) q[2];
sx q[2];
rz(-2.4487477) q[2];
sx q[2];
rz(1.8681519) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49644955) q[1];
sx q[1];
rz(-1.4525868) q[1];
sx q[1];
rz(2.955944) q[1];
x q[2];
rz(1.7433211) q[3];
sx q[3];
rz(-0.93768812) q[3];
sx q[3];
rz(-2.5987491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.51364556) q[2];
sx q[2];
rz(-2.7529035) q[2];
sx q[2];
rz(-1.2316616) q[2];
rz(1.4218467) q[3];
sx q[3];
rz(-1.809241) q[3];
sx q[3];
rz(-0.96760145) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49895367) q[0];
sx q[0];
rz(-0.068878219) q[0];
sx q[0];
rz(2.4054085) q[0];
rz(0.63756293) q[1];
sx q[1];
rz(-1.2022737) q[1];
sx q[1];
rz(-2.7307302) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7124568) q[0];
sx q[0];
rz(-0.68198689) q[0];
sx q[0];
rz(-0.76759591) q[0];
x q[1];
rz(1.246894) q[2];
sx q[2];
rz(-0.59338394) q[2];
sx q[2];
rz(1.3932799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83026559) q[1];
sx q[1];
rz(-2.5120148) q[1];
sx q[1];
rz(0.011407995) q[1];
rz(-0.66524283) q[3];
sx q[3];
rz(-1.0415631) q[3];
sx q[3];
rz(-1.0981263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32834443) q[2];
sx q[2];
rz(-1.594955) q[2];
sx q[2];
rz(-0.52097875) q[2];
rz(-0.10884918) q[3];
sx q[3];
rz(-2.1284926) q[3];
sx q[3];
rz(0.61645761) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2206889) q[0];
sx q[0];
rz(-1.9158798) q[0];
sx q[0];
rz(1.2589681) q[0];
rz(-1.0733696) q[1];
sx q[1];
rz(-1.4732889) q[1];
sx q[1];
rz(-0.64356709) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130085) q[0];
sx q[0];
rz(-1.1197392) q[0];
sx q[0];
rz(-1.5764007) q[0];
rz(-pi) q[1];
rz(-2.6344755) q[2];
sx q[2];
rz(-2.7165453) q[2];
sx q[2];
rz(0.054946446) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7832598) q[1];
sx q[1];
rz(-1.734501) q[1];
sx q[1];
rz(-0.81061426) q[1];
rz(-pi) q[2];
rz(-2.0900299) q[3];
sx q[3];
rz(-2.084084) q[3];
sx q[3];
rz(2.5198873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3448989) q[2];
sx q[2];
rz(-2.0774136) q[2];
sx q[2];
rz(1.8033484) q[2];
rz(1.4517339) q[3];
sx q[3];
rz(-1.2029878) q[3];
sx q[3];
rz(0.10045997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.79784262) q[0];
sx q[0];
rz(-1.8666973) q[0];
sx q[0];
rz(0.95329681) q[0];
rz(-2.9816309) q[1];
sx q[1];
rz(-2.5123031) q[1];
sx q[1];
rz(-0.86872855) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1686483) q[0];
sx q[0];
rz(-2.4356844) q[0];
sx q[0];
rz(1.9078518) q[0];
x q[1];
rz(2.4071099) q[2];
sx q[2];
rz(-0.96436635) q[2];
sx q[2];
rz(-0.36494985) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27426961) q[1];
sx q[1];
rz(-1.5918917) q[1];
sx q[1];
rz(2.8460851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2420011) q[3];
sx q[3];
rz(-1.7760522) q[3];
sx q[3];
rz(-1.9849856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6877044) q[2];
sx q[2];
rz(-1.1331465) q[2];
sx q[2];
rz(-0.49393168) q[2];
rz(-1.6535053) q[3];
sx q[3];
rz(-2.2736214) q[3];
sx q[3];
rz(-3.119829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22773753) q[0];
sx q[0];
rz(-1.4233002) q[0];
sx q[0];
rz(-2.8505351) q[0];
rz(0.13045467) q[1];
sx q[1];
rz(-0.37965241) q[1];
sx q[1];
rz(1.5694654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8328142) q[0];
sx q[0];
rz(-3.0893005) q[0];
sx q[0];
rz(-1.8240269) q[0];
rz(1.7849017) q[2];
sx q[2];
rz(-2.4567025) q[2];
sx q[2];
rz(-0.30944007) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56640139) q[1];
sx q[1];
rz(-2.1809887) q[1];
sx q[1];
rz(-1.4296878) q[1];
rz(1.6984822) q[3];
sx q[3];
rz(-1.8771168) q[3];
sx q[3];
rz(-1.1088913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7743885) q[2];
sx q[2];
rz(-0.17750326) q[2];
sx q[2];
rz(-3.0299178) q[2];
rz(2.3292134) q[3];
sx q[3];
rz(-1.2923391) q[3];
sx q[3];
rz(-1.5264548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3647301) q[0];
sx q[0];
rz(-0.22219816) q[0];
sx q[0];
rz(-2.7729233) q[0];
rz(-2.5313976) q[1];
sx q[1];
rz(-1.5851603) q[1];
sx q[1];
rz(1.0367702) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0087413) q[0];
sx q[0];
rz(-1.3467107) q[0];
sx q[0];
rz(-0.56613825) q[0];
x q[1];
rz(-2.0837777) q[2];
sx q[2];
rz(-1.566525) q[2];
sx q[2];
rz(1.9912602) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1644124) q[1];
sx q[1];
rz(-2.3930139) q[1];
sx q[1];
rz(1.776265) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9764594) q[3];
sx q[3];
rz(-1.0899749) q[3];
sx q[3];
rz(-0.1827249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.6250352) q[2];
sx q[2];
rz(-2.3380307) q[2];
sx q[2];
rz(-0.45300031) q[2];
rz(-0.9520483) q[3];
sx q[3];
rz(-2.0249849) q[3];
sx q[3];
rz(0.12075748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20671885) q[0];
sx q[0];
rz(-0.047858866) q[0];
sx q[0];
rz(2.5913443) q[0];
rz(-1.0847367) q[1];
sx q[1];
rz(-0.58894283) q[1];
sx q[1];
rz(1.750607) q[1];
rz(-pi) q[2];
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
rz(-2.10033) q[2];
sx q[2];
rz(-0.88418294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9908438) q[1];
sx q[1];
rz(-2.1865784) q[1];
sx q[1];
rz(-0.94373871) q[1];
x q[2];
rz(0.061270631) q[3];
sx q[3];
rz(-1.6271504) q[3];
sx q[3];
rz(-2.0298093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.37798205) q[2];
sx q[2];
rz(-0.72145975) q[2];
sx q[2];
rz(-1.1627496) q[2];
rz(2.5871596) q[3];
sx q[3];
rz(-0.88397637) q[3];
sx q[3];
rz(0.46816167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32600317) q[0];
sx q[0];
rz(-1.3023888) q[0];
sx q[0];
rz(1.6525255) q[0];
rz(2.9266657) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(0.28543546) q[2];
sx q[2];
rz(-0.68842059) q[2];
sx q[2];
rz(-0.60981087) q[2];
rz(-1.5051846) q[3];
sx q[3];
rz(-2.2917404) q[3];
sx q[3];
rz(2.9577128) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
