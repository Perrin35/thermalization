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
rz(-2.2012329) q[0];
sx q[0];
rz(-1.1223963) q[0];
sx q[0];
rz(0.19413343) q[0];
rz(1.0951207) q[1];
sx q[1];
rz(-1.6835338) q[1];
sx q[1];
rz(2.3966052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7238104) q[0];
sx q[0];
rz(-0.11427721) q[0];
sx q[0];
rz(2.6835915) q[0];
rz(2.216624) q[2];
sx q[2];
rz(-0.87069521) q[2];
sx q[2];
rz(-2.2440804) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2304334) q[1];
sx q[1];
rz(-1.8130366) q[1];
sx q[1];
rz(-1.6723067) q[1];
rz(2.0230832) q[3];
sx q[3];
rz(-1.8798401) q[3];
sx q[3];
rz(-3.0024101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8038883) q[2];
sx q[2];
rz(-1.3362198) q[2];
sx q[2];
rz(-1.814092) q[2];
rz(-1.7727857) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(2.9201065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.1526445) q[0];
sx q[0];
rz(-1.7247609) q[0];
sx q[0];
rz(0.055135559) q[0];
rz(-0.70471835) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(2.096874) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8880312) q[0];
sx q[0];
rz(-0.46875254) q[0];
sx q[0];
rz(2.5028489) q[0];
rz(-pi) q[1];
rz(-0.56684871) q[2];
sx q[2];
rz(-2.0802312) q[2];
sx q[2];
rz(1.8519362) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56480184) q[1];
sx q[1];
rz(-1.163244) q[1];
sx q[1];
rz(-2.3840586) q[1];
x q[2];
rz(2.5651188) q[3];
sx q[3];
rz(-0.11291355) q[3];
sx q[3];
rz(0.1422595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3220871) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(1.8196003) q[2];
rz(2.3540438) q[3];
sx q[3];
rz(-1.7371663) q[3];
sx q[3];
rz(-2.0866086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2128485) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(-0.71841946) q[0];
rz(-2.8058715) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(0.95300037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18291423) q[0];
sx q[0];
rz(-2.1116858) q[0];
sx q[0];
rz(0.042634115) q[0];
rz(-pi) q[1];
rz(-1.8242055) q[2];
sx q[2];
rz(-2.6077787) q[2];
sx q[2];
rz(1.6718503) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.408421) q[1];
sx q[1];
rz(-1.9431027) q[1];
sx q[1];
rz(-2.2835963) q[1];
x q[2];
rz(-2.4673053) q[3];
sx q[3];
rz(-2.4113048) q[3];
sx q[3];
rz(-1.3004608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.0042051729) q[2];
sx q[2];
rz(-1.1342099) q[2];
sx q[2];
rz(-0.65286621) q[2];
rz(-0.79814923) q[3];
sx q[3];
rz(-1.8316725) q[3];
sx q[3];
rz(-2.4979512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.2586486) q[0];
sx q[0];
rz(-2.301321) q[0];
sx q[0];
rz(-2.3249481) q[0];
rz(2.4500997) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(0.74849558) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25835055) q[0];
sx q[0];
rz(-1.2628947) q[0];
sx q[0];
rz(1.3624205) q[0];
rz(-pi) q[1];
rz(-1.7848133) q[2];
sx q[2];
rz(-1.4811918) q[2];
sx q[2];
rz(-2.2387308) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69341993) q[1];
sx q[1];
rz(-0.93575555) q[1];
sx q[1];
rz(-2.5513812) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3393551) q[3];
sx q[3];
rz(-1.7417296) q[3];
sx q[3];
rz(0.49491301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3889435) q[2];
sx q[2];
rz(-2.3708673) q[2];
sx q[2];
rz(2.3469817) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36885095) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(-0.4976196) q[0];
rz(0.76958641) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(-2.41113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896547) q[0];
sx q[0];
rz(-1.7853072) q[0];
sx q[0];
rz(-2.4431321) q[0];
rz(1.914417) q[2];
sx q[2];
rz(-0.65672749) q[2];
sx q[2];
rz(1.7467787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5268908) q[1];
sx q[1];
rz(-2.3040471) q[1];
sx q[1];
rz(-0.28403838) q[1];
x q[2];
rz(0.45017193) q[3];
sx q[3];
rz(-2.014959) q[3];
sx q[3];
rz(1.1314363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2490425) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(-0.45336938) q[2];
rz(1.310965) q[3];
sx q[3];
rz(-2.9679208) q[3];
sx q[3];
rz(0.12392509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675568) q[0];
sx q[0];
rz(-1.0834162) q[0];
sx q[0];
rz(-1.7932844) q[0];
rz(1.096161) q[1];
sx q[1];
rz(-1.4827012) q[1];
sx q[1];
rz(-2.5199264) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2638877) q[0];
sx q[0];
rz(-1.5549721) q[0];
sx q[0];
rz(0.71431922) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39577272) q[2];
sx q[2];
rz(-1.6623375) q[2];
sx q[2];
rz(1.7965732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7324896) q[1];
sx q[1];
rz(-1.81899) q[1];
sx q[1];
rz(0.77381882) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4824941) q[3];
sx q[3];
rz(-0.78135282) q[3];
sx q[3];
rz(2.3574061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2094476) q[2];
sx q[2];
rz(-0.79938447) q[2];
sx q[2];
rz(-1.7032334) q[2];
rz(-0.98073331) q[3];
sx q[3];
rz(-1.2659729) q[3];
sx q[3];
rz(-1.929662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88406968) q[0];
sx q[0];
rz(-1.3781837) q[0];
sx q[0];
rz(-2.8164465) q[0];
rz(-2.6349321) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(-1.315717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36111212) q[0];
sx q[0];
rz(-2.4473815) q[0];
sx q[0];
rz(-1.7630462) q[0];
x q[1];
rz(-0.99847128) q[2];
sx q[2];
rz(-1.1954167) q[2];
sx q[2];
rz(-1.7121512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.81974492) q[1];
sx q[1];
rz(-2.3492866) q[1];
sx q[1];
rz(2.5755432) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5127038) q[3];
sx q[3];
rz(-2.5085555) q[3];
sx q[3];
rz(2.7881088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.331984) q[2];
sx q[2];
rz(-2.1029682) q[2];
sx q[2];
rz(-1.6161551) q[2];
rz(1.7841313) q[3];
sx q[3];
rz(-1.4493891) q[3];
sx q[3];
rz(1.5176895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1656533) q[0];
sx q[0];
rz(-1.140056) q[0];
sx q[0];
rz(2.4923988) q[0];
rz(0.80750418) q[1];
sx q[1];
rz(-1.1367926) q[1];
sx q[1];
rz(1.753122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8313468) q[0];
sx q[0];
rz(-0.60605907) q[0];
sx q[0];
rz(1.640727) q[0];
rz(-pi) q[1];
rz(-1.2818579) q[2];
sx q[2];
rz(-1.4380437) q[2];
sx q[2];
rz(1.8366739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35508868) q[1];
sx q[1];
rz(-1.638371) q[1];
sx q[1];
rz(1.7421634) q[1];
rz(-pi) q[2];
rz(1.7817333) q[3];
sx q[3];
rz(-1.1865133) q[3];
sx q[3];
rz(2.3944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93747741) q[2];
sx q[2];
rz(-2.9948339) q[2];
sx q[2];
rz(0.80782962) q[2];
rz(-0.67148036) q[3];
sx q[3];
rz(-1.5059794) q[3];
sx q[3];
rz(1.5988662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143739) q[0];
sx q[0];
rz(-0.53633538) q[0];
sx q[0];
rz(0.45289034) q[0];
rz(0.29531404) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(1.3989075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59205306) q[0];
sx q[0];
rz(-0.73707923) q[0];
sx q[0];
rz(1.6960742) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0844551) q[2];
sx q[2];
rz(-1.4342564) q[2];
sx q[2];
rz(-1.650102) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6878444) q[1];
sx q[1];
rz(-2.0216989) q[1];
sx q[1];
rz(-1.9912254) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3736428) q[3];
sx q[3];
rz(-1.5492566) q[3];
sx q[3];
rz(-0.43019274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7704775) q[2];
sx q[2];
rz(-2.2687843) q[2];
sx q[2];
rz(0.57658833) q[2];
rz(2.9234486) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(-1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52467728) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(0.088454811) q[0];
rz(-2.6103919) q[1];
sx q[1];
rz(-0.4041268) q[1];
sx q[1];
rz(1.6204576) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94604055) q[0];
sx q[0];
rz(-1.4874391) q[0];
sx q[0];
rz(1.964919) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4786104) q[2];
sx q[2];
rz(-0.47904821) q[2];
sx q[2];
rz(2.6625815) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38177165) q[1];
sx q[1];
rz(-0.23649892) q[1];
sx q[1];
rz(0.31707724) q[1];
rz(-pi) q[2];
rz(-0.57595595) q[3];
sx q[3];
rz(-0.68000845) q[3];
sx q[3];
rz(0.81871225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8497808) q[2];
sx q[2];
rz(-1.7917289) q[2];
sx q[2];
rz(1.527171) q[2];
rz(-1.8725066) q[3];
sx q[3];
rz(-0.98617712) q[3];
sx q[3];
rz(-2.0932978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.108718) q[0];
sx q[0];
rz(-1.1058818) q[0];
sx q[0];
rz(-0.82431128) q[0];
rz(0.4737919) q[1];
sx q[1];
rz(-1.5693249) q[1];
sx q[1];
rz(1.5798461) q[1];
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
