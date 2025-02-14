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
rz(-1.3402101) q[0];
sx q[0];
rz(-2.8647381) q[0];
sx q[0];
rz(1.0387596) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(-2.3050397) q[1];
sx q[1];
rz(0.41681448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3562389) q[0];
sx q[0];
rz(-0.84945852) q[0];
sx q[0];
rz(1.4291886) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7567467) q[2];
sx q[2];
rz(-2.550808) q[2];
sx q[2];
rz(0.42921517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3089247) q[1];
sx q[1];
rz(-1.0005669) q[1];
sx q[1];
rz(2.4491688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5384212) q[3];
sx q[3];
rz(-1.4093295) q[3];
sx q[3];
rz(-2.3700455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(1.2539585) q[2];
rz(1.5597255) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.9480243) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(-0.76876202) q[0];
rz(-1.7747152) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(1.0381402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65377849) q[0];
sx q[0];
rz(-0.013049203) q[0];
sx q[0];
rz(2.292146) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72121303) q[2];
sx q[2];
rz(-1.291847) q[2];
sx q[2];
rz(-2.9245289) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1662054) q[1];
sx q[1];
rz(-2.176099) q[1];
sx q[1];
rz(2.0209842) q[1];
x q[2];
rz(-1.6065381) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(-1.3066178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64608964) q[2];
sx q[2];
rz(-1.8547408) q[2];
sx q[2];
rz(2.3213279) q[2];
rz(2.8881554) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(-0.36111116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609817) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(2.2542727) q[0];
rz(-2.6112828) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(-1.1393772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4654008) q[0];
sx q[0];
rz(-1.4290819) q[0];
sx q[0];
rz(2.7928674) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0372889) q[2];
sx q[2];
rz(-0.59743687) q[2];
sx q[2];
rz(3.1410599) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4344221) q[1];
sx q[1];
rz(-1.5853197) q[1];
sx q[1];
rz(2.8959031) q[1];
x q[2];
rz(0.091413012) q[3];
sx q[3];
rz(-0.52052906) q[3];
sx q[3];
rz(2.6176069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.016971074) q[2];
sx q[2];
rz(-1.4554224) q[2];
sx q[2];
rz(2.7023081) q[2];
rz(0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4029694) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(-0.11548197) q[0];
rz(0.11784095) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(1.8353362) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7457434) q[0];
sx q[0];
rz(-1.1538528) q[0];
sx q[0];
rz(-2.8265619) q[0];
x q[1];
rz(-2.6682786) q[2];
sx q[2];
rz(-0.71719786) q[2];
sx q[2];
rz(-1.0953449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0688248) q[1];
sx q[1];
rz(-2.1842274) q[1];
sx q[1];
rz(2.3062506) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6870858) q[3];
sx q[3];
rz(-2.5657885) q[3];
sx q[3];
rz(0.35279122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8370886) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(-2.9869249) q[2];
rz(1.539544) q[3];
sx q[3];
rz(-1.3402901) q[3];
sx q[3];
rz(0.60607564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57347572) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(0.83531761) q[0];
rz(0.71594816) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(0.11944019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855921) q[0];
sx q[0];
rz(-1.5876666) q[0];
sx q[0];
rz(-1.5149679) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30952366) q[2];
sx q[2];
rz(-2.672827) q[2];
sx q[2];
rz(-3.0075207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.844567) q[1];
sx q[1];
rz(-0.78999619) q[1];
sx q[1];
rz(0.67761919) q[1];
rz(-0.53262696) q[3];
sx q[3];
rz(-1.9386734) q[3];
sx q[3];
rz(-1.0581349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2401838) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(-0.060997941) q[2];
rz(0.048246233) q[3];
sx q[3];
rz(-1.9082853) q[3];
sx q[3];
rz(2.4129996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.401684) q[0];
sx q[0];
rz(-0.70867276) q[0];
sx q[0];
rz(-2.0106864) q[0];
rz(2.6629958) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(-1.7344249) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78188183) q[0];
sx q[0];
rz(-3.0880083) q[0];
sx q[0];
rz(-3.1126541) q[0];
rz(-pi) q[1];
rz(0.84848225) q[2];
sx q[2];
rz(-1.3568078) q[2];
sx q[2];
rz(-3.1297562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9764912) q[1];
sx q[1];
rz(-2.8860984) q[1];
sx q[1];
rz(1.1063904) q[1];
x q[2];
rz(-2.4396517) q[3];
sx q[3];
rz(-0.27595584) q[3];
sx q[3];
rz(-0.70358932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56167928) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(-1.9631867) q[2];
rz(2.0922349) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(-1.8572846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9310164) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(-1.806102) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(0.30805045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0548693) q[0];
sx q[0];
rz(-2.6941774) q[0];
sx q[0];
rz(-0.50680508) q[0];
rz(2.5045583) q[2];
sx q[2];
rz(-0.61372988) q[2];
sx q[2];
rz(1.5368652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96770331) q[1];
sx q[1];
rz(-0.99882616) q[1];
sx q[1];
rz(0.58764761) q[1];
rz(1.9568029) q[3];
sx q[3];
rz(-1.6148477) q[3];
sx q[3];
rz(0.06202997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7293952) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(0.13128734) q[2];
rz(-0.73733759) q[3];
sx q[3];
rz(-1.3056583) q[3];
sx q[3];
rz(2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8312296) q[0];
sx q[0];
rz(-1.0897626) q[0];
sx q[0];
rz(-2.6112153) q[0];
rz(-1.7165548) q[1];
sx q[1];
rz(-2.0030463) q[1];
sx q[1];
rz(2.4200965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24151267) q[0];
sx q[0];
rz(-2.3519313) q[0];
sx q[0];
rz(-0.65171839) q[0];
rz(2.8558735) q[2];
sx q[2];
rz(-2.6556394) q[2];
sx q[2];
rz(3.1348117) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65039448) q[1];
sx q[1];
rz(-1.7868687) q[1];
sx q[1];
rz(-2.9465527) q[1];
rz(0.11105342) q[3];
sx q[3];
rz(-2.3559582) q[3];
sx q[3];
rz(1.3291886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(-1.4777769) q[2];
rz(-1.9780698) q[3];
sx q[3];
rz(-1.082837) q[3];
sx q[3];
rz(-1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32996938) q[0];
sx q[0];
rz(-1.4412619) q[0];
sx q[0];
rz(-2.5323618) q[0];
rz(-1.5513264) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.4564266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1066172) q[0];
sx q[0];
rz(-1.9272695) q[0];
sx q[0];
rz(0.26654213) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11498835) q[2];
sx q[2];
rz(-0.48795944) q[2];
sx q[2];
rz(-1.8166627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8943772) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(-1.2213329) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84734579) q[3];
sx q[3];
rz(-2.06978) q[3];
sx q[3];
rz(2.8049198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41813254) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(-0.88031236) q[2];
rz(-0.20600016) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77077615) q[0];
sx q[0];
rz(-2.4801065) q[0];
sx q[0];
rz(-3.1296375) q[0];
rz(1.5097584) q[1];
sx q[1];
rz(-0.44523528) q[1];
sx q[1];
rz(0.59648046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76145455) q[0];
sx q[0];
rz(-1.589847) q[0];
sx q[0];
rz(0.46940243) q[0];
rz(2.7786461) q[2];
sx q[2];
rz(-1.3273113) q[2];
sx q[2];
rz(-1.7044662) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0833019) q[1];
sx q[1];
rz(-1.4920248) q[1];
sx q[1];
rz(-1.0942671) q[1];
rz(-2.6634253) q[3];
sx q[3];
rz(-2.8088125) q[3];
sx q[3];
rz(1.5733583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(-0.19700024) q[2];
rz(-1.152285) q[3];
sx q[3];
rz(-2.9983493) q[3];
sx q[3];
rz(-0.44767374) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780509) q[0];
sx q[0];
rz(-1.0095689) q[0];
sx q[0];
rz(2.9472245) q[0];
rz(-2.9327783) q[1];
sx q[1];
rz(-1.6046235) q[1];
sx q[1];
rz(-1.0135289) q[1];
rz(-1.6848736) q[2];
sx q[2];
rz(-2.4488505) q[2];
sx q[2];
rz(1.746576) q[2];
rz(-0.98255929) q[3];
sx q[3];
rz(-2.0053902) q[3];
sx q[3];
rz(-2.7872661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
