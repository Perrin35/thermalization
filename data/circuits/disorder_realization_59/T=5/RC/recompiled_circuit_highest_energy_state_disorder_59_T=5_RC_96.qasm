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
rz(0.94035971) q[0];
sx q[0];
rz(4.263989) q[0];
sx q[0];
rz(12.372237) q[0];
rz(-2.0464719) q[1];
sx q[1];
rz(-1.4580589) q[1];
sx q[1];
rz(-2.3966052) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8391957) q[0];
sx q[0];
rz(-1.6212362) q[0];
sx q[0];
rz(-3.0390059) q[0];
x q[1];
rz(-2.216624) q[2];
sx q[2];
rz(-2.2708974) q[2];
sx q[2];
rz(-2.2440804) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.68406635) q[1];
sx q[1];
rz(-1.4722595) q[1];
sx q[1];
rz(-0.24344488) q[1];
rz(-pi) q[2];
rz(-2.2017415) q[3];
sx q[3];
rz(-2.5999024) q[3];
sx q[3];
rz(1.9909137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3377043) q[2];
sx q[2];
rz(-1.3362198) q[2];
sx q[2];
rz(-1.3275006) q[2];
rz(1.7727857) q[3];
sx q[3];
rz(-2.8076706) q[3];
sx q[3];
rz(2.9201065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.1526445) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(0.055135559) q[0];
rz(-2.4368743) q[1];
sx q[1];
rz(-1.5780508) q[1];
sx q[1];
rz(2.096874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25356141) q[0];
sx q[0];
rz(-0.46875254) q[0];
sx q[0];
rz(-0.6387438) q[0];
x q[1];
rz(0.80514812) q[2];
sx q[2];
rz(-0.74294801) q[2];
sx q[2];
rz(-2.2066903) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56480184) q[1];
sx q[1];
rz(-1.163244) q[1];
sx q[1];
rz(2.3840586) q[1];
x q[2];
rz(-2.5651188) q[3];
sx q[3];
rz(-0.11291355) q[3];
sx q[3];
rz(-0.1422595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.81950554) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(1.3219924) q[2];
rz(-2.3540438) q[3];
sx q[3];
rz(-1.7371663) q[3];
sx q[3];
rz(-1.0549841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2128485) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(-0.71841946) q[0];
rz(0.33572117) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(0.95300037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18291423) q[0];
sx q[0];
rz(-2.1116858) q[0];
sx q[0];
rz(0.042634115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0510873) q[2];
sx q[2];
rz(-1.698709) q[2];
sx q[2];
rz(-3.0232883) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.73317161) q[1];
sx q[1];
rz(-1.1984899) q[1];
sx q[1];
rz(2.2835963) q[1];
x q[2];
rz(2.0805667) q[3];
sx q[3];
rz(-1.022664) q[3];
sx q[3];
rz(-0.47993101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1373875) q[2];
sx q[2];
rz(-1.1342099) q[2];
sx q[2];
rz(-2.4887264) q[2];
rz(-0.79814923) q[3];
sx q[3];
rz(-1.8316725) q[3];
sx q[3];
rz(-2.4979512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.882944) q[0];
sx q[0];
rz(-0.84027165) q[0];
sx q[0];
rz(-0.81664455) q[0];
rz(2.4500997) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(0.74849558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2740606) q[0];
sx q[0];
rz(-0.36991773) q[0];
sx q[0];
rz(2.5648613) q[0];
rz(-pi) q[1];
rz(-1.9709936) q[2];
sx q[2];
rz(-2.9098401) q[2];
sx q[2];
rz(1.0585001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2556127) q[1];
sx q[1];
rz(-1.1062262) q[1];
sx q[1];
rz(2.2963312) q[1];
x q[2];
rz(-2.9660712) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(2.1057749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3889435) q[2];
sx q[2];
rz(-2.3708673) q[2];
sx q[2];
rz(0.79461092) q[2];
rz(-1.4258344) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(-2.7690601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36885095) q[0];
sx q[0];
rz(-1.0798825) q[0];
sx q[0];
rz(0.4976196) q[0];
rz(-2.3720062) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(-2.41113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4896547) q[0];
sx q[0];
rz(-1.3562855) q[0];
sx q[0];
rz(0.69846054) q[0];
x q[1];
rz(-2.1986352) q[2];
sx q[2];
rz(-1.7779609) q[2];
sx q[2];
rz(-0.10016537) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9047315) q[1];
sx q[1];
rz(-1.7805459) q[1];
sx q[1];
rz(0.81717234) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0570083) q[3];
sx q[3];
rz(-1.9745805) q[3];
sx q[3];
rz(-2.9069834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2490425) q[2];
sx q[2];
rz(-1.6232792) q[2];
sx q[2];
rz(2.6882233) q[2];
rz(1.8306277) q[3];
sx q[3];
rz(-0.1736719) q[3];
sx q[3];
rz(-3.0176676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675568) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(1.7932844) q[0];
rz(-1.096161) q[1];
sx q[1];
rz(-1.6588914) q[1];
sx q[1];
rz(0.62166628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28866119) q[0];
sx q[0];
rz(-0.71446361) q[0];
sx q[0];
rz(-0.024152569) q[0];
x q[1];
rz(-1.4716343) q[2];
sx q[2];
rz(-1.9648203) q[2];
sx q[2];
rz(-2.9539915) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7335947) q[1];
sx q[1];
rz(-0.80469614) q[1];
sx q[1];
rz(-0.3478653) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0249025) q[3];
sx q[3];
rz(-0.98034795) q[3];
sx q[3];
rz(1.6131372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2094476) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(-1.7032334) q[2];
rz(-2.1608593) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88406968) q[0];
sx q[0];
rz(-1.3781837) q[0];
sx q[0];
rz(-2.8164465) q[0];
rz(0.50666058) q[1];
sx q[1];
rz(-2.1279361) q[1];
sx q[1];
rz(-1.8258757) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7804805) q[0];
sx q[0];
rz(-2.4473815) q[0];
sx q[0];
rz(-1.3785465) q[0];
rz(-pi) q[1];
rz(2.1998134) q[2];
sx q[2];
rz(-0.67275362) q[2];
sx q[2];
rz(-2.7655809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3218477) q[1];
sx q[1];
rz(-0.79230601) q[1];
sx q[1];
rz(-2.5755432) q[1];
rz(-pi) q[2];
rz(2.2030284) q[3];
sx q[3];
rz(-1.5364416) q[3];
sx q[3];
rz(1.1704579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.331984) q[2];
sx q[2];
rz(-2.1029682) q[2];
sx q[2];
rz(1.6161551) q[2];
rz(1.7841313) q[3];
sx q[3];
rz(-1.4493891) q[3];
sx q[3];
rz(-1.6239032) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656533) q[0];
sx q[0];
rz(-2.0015367) q[0];
sx q[0];
rz(2.4923988) q[0];
rz(0.80750418) q[1];
sx q[1];
rz(-2.0048001) q[1];
sx q[1];
rz(1.3884707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20304414) q[0];
sx q[0];
rz(-1.5309835) q[0];
sx q[0];
rz(0.96588246) q[0];
rz(-pi) q[1];
rz(-1.2818579) q[2];
sx q[2];
rz(-1.7035489) q[2];
sx q[2];
rz(-1.8366739) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.786504) q[1];
sx q[1];
rz(-1.5032217) q[1];
sx q[1];
rz(1.3994292) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47775538) q[3];
sx q[3];
rz(-2.7057608) q[3];
sx q[3];
rz(-1.2660668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2041152) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(2.333763) q[2];
rz(2.4701123) q[3];
sx q[3];
rz(-1.5059794) q[3];
sx q[3];
rz(-1.5427264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9143739) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(0.45289034) q[0];
rz(0.29531404) q[1];
sx q[1];
rz(-1.9120049) q[1];
sx q[1];
rz(1.7426851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3810661) q[0];
sx q[0];
rz(-0.84080836) q[0];
sx q[0];
rz(0.11294079) q[0];
rz(-pi) q[1];
rz(1.8434502) q[2];
sx q[2];
rz(-2.6116707) q[2];
sx q[2];
rz(2.8255445) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4537482) q[1];
sx q[1];
rz(-1.1198938) q[1];
sx q[1];
rz(-1.1503673) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7679498) q[3];
sx q[3];
rz(-1.592336) q[3];
sx q[3];
rz(2.7113999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7704775) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(-0.57658833) q[2];
rz(2.9234486) q[3];
sx q[3];
rz(-0.79272565) q[3];
sx q[3];
rz(1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52467728) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(3.0531378) q[0];
rz(-2.6103919) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(1.6204576) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5514497) q[0];
sx q[0];
rz(-1.9634754) q[0];
sx q[0];
rz(-0.090242437) q[0];
x q[1];
rz(-1.0934866) q[2];
sx q[2];
rz(-1.6132406) q[2];
sx q[2];
rz(2.1316656) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.056235751) q[1];
sx q[1];
rz(-1.7952982) q[1];
sx q[1];
rz(-1.6457998) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5456387) q[3];
sx q[3];
rz(-1.9203382) q[3];
sx q[3];
rz(0.28448018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8497808) q[2];
sx q[2];
rz(-1.3498638) q[2];
sx q[2];
rz(-1.6144217) q[2];
rz(-1.269086) q[3];
sx q[3];
rz(-2.1554155) q[3];
sx q[3];
rz(-2.0932978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0328746) q[0];
sx q[0];
rz(-1.1058818) q[0];
sx q[0];
rz(-0.82431128) q[0];
rz(2.6678008) q[1];
sx q[1];
rz(-1.5722678) q[1];
sx q[1];
rz(-1.5617465) q[1];
rz(1.9971581) q[2];
sx q[2];
rz(-2.1254267) q[2];
sx q[2];
rz(2.8993901) q[2];
rz(0.99226034) q[3];
sx q[3];
rz(-1.0836061) q[3];
sx q[3];
rz(1.5988812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
