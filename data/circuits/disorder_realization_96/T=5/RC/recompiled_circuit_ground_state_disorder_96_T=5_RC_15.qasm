OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.191303) q[0];
sx q[0];
rz(-2.8714955) q[0];
sx q[0];
rz(-0.88859963) q[0];
rz(1.8285881) q[1];
sx q[1];
rz(-1.5421901) q[1];
sx q[1];
rz(1.3808274) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.439405) q[0];
sx q[0];
rz(-2.8719465) q[0];
sx q[0];
rz(-2.3764971) q[0];
rz(-1.0384473) q[2];
sx q[2];
rz(-0.84538904) q[2];
sx q[2];
rz(-0.40276819) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7270679) q[1];
sx q[1];
rz(-1.5611851) q[1];
sx q[1];
rz(0.33551402) q[1];
x q[2];
rz(2.1112135) q[3];
sx q[3];
rz(-0.4183397) q[3];
sx q[3];
rz(2.6650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31972739) q[2];
sx q[2];
rz(-1.8165908) q[2];
sx q[2];
rz(-2.3925609) q[2];
rz(2.9876515) q[3];
sx q[3];
rz(-1.0356244) q[3];
sx q[3];
rz(0.099893959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.446949) q[0];
sx q[0];
rz(-1.4780937) q[0];
sx q[0];
rz(1.2167759) q[0];
rz(1.0617537) q[1];
sx q[1];
rz(-2.3371425) q[1];
sx q[1];
rz(0.42713508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7593133) q[0];
sx q[0];
rz(-1.0045035) q[0];
sx q[0];
rz(-1.1575538) q[0];
x q[1];
rz(-1.2522302) q[2];
sx q[2];
rz(-0.67020352) q[2];
sx q[2];
rz(-2.2549652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6500351) q[1];
sx q[1];
rz(-1.2833529) q[1];
sx q[1];
rz(0.14658714) q[1];
x q[2];
rz(-1.3668381) q[3];
sx q[3];
rz(-1.6220304) q[3];
sx q[3];
rz(-2.8683087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0624258) q[2];
sx q[2];
rz(-1.7288952) q[2];
sx q[2];
rz(-1.742935) q[2];
rz(2.9178197) q[3];
sx q[3];
rz(-2.2327435) q[3];
sx q[3];
rz(1.5435425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021521213) q[0];
sx q[0];
rz(-1.3503617) q[0];
sx q[0];
rz(3.0960826) q[0];
rz(-1.9728164) q[1];
sx q[1];
rz(-1.2380506) q[1];
sx q[1];
rz(0.57560903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94852266) q[0];
sx q[0];
rz(-2.2098477) q[0];
sx q[0];
rz(1.5568514) q[0];
rz(0.017150684) q[2];
sx q[2];
rz(-0.91020012) q[2];
sx q[2];
rz(1.2397546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82390416) q[1];
sx q[1];
rz(-2.1145193) q[1];
sx q[1];
rz(2.8871817) q[1];
rz(-1.5381673) q[3];
sx q[3];
rz(-2.4885119) q[3];
sx q[3];
rz(0.075270502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0098972926) q[2];
sx q[2];
rz(-0.74247777) q[2];
sx q[2];
rz(2.1843145) q[2];
rz(-0.57473985) q[3];
sx q[3];
rz(-1.6114019) q[3];
sx q[3];
rz(1.8360229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9855373) q[0];
sx q[0];
rz(-0.95273459) q[0];
sx q[0];
rz(-1.8804469) q[0];
rz(3.0592697) q[1];
sx q[1];
rz(-1.0876834) q[1];
sx q[1];
rz(-1.3538768) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9808423) q[0];
sx q[0];
rz(-2.7880362) q[0];
sx q[0];
rz(0.73676957) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48575966) q[2];
sx q[2];
rz(-1.4076715) q[2];
sx q[2];
rz(0.32147929) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9685843) q[1];
sx q[1];
rz(-2.6436596) q[1];
sx q[1];
rz(-0.11346101) q[1];
rz(0.4889826) q[3];
sx q[3];
rz(-0.95150286) q[3];
sx q[3];
rz(1.3448546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7202619) q[2];
sx q[2];
rz(-2.223184) q[2];
sx q[2];
rz(1.8850231) q[2];
rz(-0.71150696) q[3];
sx q[3];
rz(-1.3308176) q[3];
sx q[3];
rz(2.8594657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8307777) q[0];
sx q[0];
rz(-0.84596914) q[0];
sx q[0];
rz(0.34314439) q[0];
rz(0.061773069) q[1];
sx q[1];
rz(-2.1715178) q[1];
sx q[1];
rz(-1.4168581) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34300229) q[0];
sx q[0];
rz(-1.8483642) q[0];
sx q[0];
rz(-2.9647602) q[0];
rz(-pi) q[1];
rz(2.4409962) q[2];
sx q[2];
rz(-1.3178409) q[2];
sx q[2];
rz(2.7299936) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38450634) q[1];
sx q[1];
rz(-1.0913186) q[1];
sx q[1];
rz(-0.15285413) q[1];
x q[2];
rz(2.0569818) q[3];
sx q[3];
rz(-2.668368) q[3];
sx q[3];
rz(1.5718232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5938277) q[2];
sx q[2];
rz(-1.6739028) q[2];
sx q[2];
rz(2.055638) q[2];
rz(0.91935277) q[3];
sx q[3];
rz(-1.7707526) q[3];
sx q[3];
rz(2.5277188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3980961) q[0];
sx q[0];
rz(-1.2092104) q[0];
sx q[0];
rz(-2.3486775) q[0];
rz(-2.4010557) q[1];
sx q[1];
rz(-0.99895993) q[1];
sx q[1];
rz(0.83121306) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7949267) q[0];
sx q[0];
rz(-1.7788609) q[0];
sx q[0];
rz(-1.1870334) q[0];
rz(-pi) q[1];
rz(-3.1341427) q[2];
sx q[2];
rz(-1.2937355) q[2];
sx q[2];
rz(-1.7817093) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38922849) q[1];
sx q[1];
rz(-0.78654754) q[1];
sx q[1];
rz(2.0081372) q[1];
rz(-pi) q[2];
rz(1.100086) q[3];
sx q[3];
rz(-1.5522955) q[3];
sx q[3];
rz(-2.5262566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0703766) q[2];
sx q[2];
rz(-1.8698317) q[2];
sx q[2];
rz(-1.1191204) q[2];
rz(1.8708771) q[3];
sx q[3];
rz(-1.1071353) q[3];
sx q[3];
rz(1.7601815) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374461) q[0];
sx q[0];
rz(-2.9747712) q[0];
sx q[0];
rz(1.5995837) q[0];
rz(2.1265325) q[1];
sx q[1];
rz(-1.5559745) q[1];
sx q[1];
rz(-2.9687845) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8625582) q[0];
sx q[0];
rz(-2.1365215) q[0];
sx q[0];
rz(2.5507257) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2020993) q[2];
sx q[2];
rz(-2.2208344) q[2];
sx q[2];
rz(1.4024613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1041278) q[1];
sx q[1];
rz(-1.8191254) q[1];
sx q[1];
rz(2.009997) q[1];
x q[2];
rz(-0.76877131) q[3];
sx q[3];
rz(-1.1535346) q[3];
sx q[3];
rz(-0.64266838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.479594) q[2];
sx q[2];
rz(-0.84344232) q[2];
sx q[2];
rz(2.1638828) q[2];
rz(0.57502037) q[3];
sx q[3];
rz(-1.2049371) q[3];
sx q[3];
rz(-1.2303111) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8136895) q[0];
sx q[0];
rz(-0.29569018) q[0];
sx q[0];
rz(3.1258702) q[0];
rz(1.9533336) q[1];
sx q[1];
rz(-0.63242042) q[1];
sx q[1];
rz(-1.9421008) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2000113) q[0];
sx q[0];
rz(-2.0111548) q[0];
sx q[0];
rz(2.8811127) q[0];
rz(-pi) q[1];
rz(-1.472961) q[2];
sx q[2];
rz(-2.6487158) q[2];
sx q[2];
rz(-2.5081483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75592985) q[1];
sx q[1];
rz(-2.0634868) q[1];
sx q[1];
rz(1.0316959) q[1];
rz(-1.1114798) q[3];
sx q[3];
rz(-0.62255854) q[3];
sx q[3];
rz(2.7186269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99379313) q[2];
sx q[2];
rz(-1.9027998) q[2];
sx q[2];
rz(1.7564868) q[2];
rz(2.5177453) q[3];
sx q[3];
rz(-1.2522937) q[3];
sx q[3];
rz(-1.0998211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95595908) q[0];
sx q[0];
rz(-0.82056844) q[0];
sx q[0];
rz(2.4483335) q[0];
rz(-2.8765826) q[1];
sx q[1];
rz(-2.3131504) q[1];
sx q[1];
rz(-1.6185435) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5582568) q[0];
sx q[0];
rz(-0.77363217) q[0];
sx q[0];
rz(-3.1257854) q[0];
rz(-0.62317836) q[2];
sx q[2];
rz(-2.0668536) q[2];
sx q[2];
rz(0.11962275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48623006) q[1];
sx q[1];
rz(-1.910348) q[1];
sx q[1];
rz(1.5961673) q[1];
x q[2];
rz(1.6783488) q[3];
sx q[3];
rz(-2.5578024) q[3];
sx q[3];
rz(1.8369758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30274621) q[2];
sx q[2];
rz(-1.6281717) q[2];
sx q[2];
rz(1.4576853) q[2];
rz(0.95528209) q[3];
sx q[3];
rz(-2.6421319) q[3];
sx q[3];
rz(2.3063851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38462287) q[0];
sx q[0];
rz(-0.56461016) q[0];
sx q[0];
rz(-2.6589174) q[0];
rz(-2.2118498) q[1];
sx q[1];
rz(-1.0044121) q[1];
sx q[1];
rz(-2.7453056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38448725) q[0];
sx q[0];
rz(-1.9845909) q[0];
sx q[0];
rz(1.1767469) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17896374) q[2];
sx q[2];
rz(-1.71993) q[2];
sx q[2];
rz(-2.3863132) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0925818) q[1];
sx q[1];
rz(-1.3900458) q[1];
sx q[1];
rz(0.84047079) q[1];
x q[2];
rz(-1.2108874) q[3];
sx q[3];
rz(-1.9511838) q[3];
sx q[3];
rz(1.4402657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3489909) q[2];
sx q[2];
rz(-1.4395809) q[2];
sx q[2];
rz(-2.8144042) q[2];
rz(2.3606825) q[3];
sx q[3];
rz(-1.1612929) q[3];
sx q[3];
rz(-1.6646632) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1031716) q[0];
sx q[0];
rz(-2.1506943) q[0];
sx q[0];
rz(-2.8839169) q[0];
rz(-2.7720263) q[1];
sx q[1];
rz(-2.259544) q[1];
sx q[1];
rz(2.4824711) q[1];
rz(-0.69330458) q[2];
sx q[2];
rz(-1.5360166) q[2];
sx q[2];
rz(1.7770384) q[2];
rz(0.9624858) q[3];
sx q[3];
rz(-1.2506093) q[3];
sx q[3];
rz(0.067230732) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
