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
rz(-0.33997384) q[0];
sx q[0];
rz(-0.61859328) q[0];
sx q[0];
rz(0.11237385) q[0];
rz(-2.1984341) q[1];
sx q[1];
rz(-1.4961493) q[1];
sx q[1];
rz(1.5785616) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0267082) q[0];
sx q[0];
rz(-2.0100532) q[0];
sx q[0];
rz(-0.54217302) q[0];
x q[1];
rz(-1.4077827) q[2];
sx q[2];
rz(-2.2634677) q[2];
sx q[2];
rz(2.7977347) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1670363) q[1];
sx q[1];
rz(-0.99460973) q[1];
sx q[1];
rz(-0.78414708) q[1];
rz(0.89157414) q[3];
sx q[3];
rz(-1.8603311) q[3];
sx q[3];
rz(1.139977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9872687) q[2];
sx q[2];
rz(-1.6538041) q[2];
sx q[2];
rz(2.7538015) q[2];
rz(-0.44529861) q[3];
sx q[3];
rz(-2.6058091) q[3];
sx q[3];
rz(2.9545412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33291373) q[0];
sx q[0];
rz(-0.50030047) q[0];
sx q[0];
rz(2.152541) q[0];
rz(1.5736846) q[1];
sx q[1];
rz(-0.75850073) q[1];
sx q[1];
rz(3.0737976) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2382375) q[0];
sx q[0];
rz(-0.35404917) q[0];
sx q[0];
rz(-2.7992008) q[0];
x q[1];
rz(-2.8591462) q[2];
sx q[2];
rz(-1.2722655) q[2];
sx q[2];
rz(1.1242497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47956672) q[1];
sx q[1];
rz(-2.5218039) q[1];
sx q[1];
rz(2.8940043) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12884295) q[3];
sx q[3];
rz(-2.2097603) q[3];
sx q[3];
rz(2.2966677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13061358) q[2];
sx q[2];
rz(-2.4108672) q[2];
sx q[2];
rz(2.0558426) q[2];
rz(-2.7340414) q[3];
sx q[3];
rz(-1.4980039) q[3];
sx q[3];
rz(-1.7727859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8959592) q[0];
sx q[0];
rz(-0.85773724) q[0];
sx q[0];
rz(-1.6813543) q[0];
rz(-1.1145837) q[1];
sx q[1];
rz(-2.3268301) q[1];
sx q[1];
rz(-2.1688555) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7290608) q[0];
sx q[0];
rz(-1.5648989) q[0];
sx q[0];
rz(-0.062882857) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8778828) q[2];
sx q[2];
rz(-0.59520856) q[2];
sx q[2];
rz(-0.021364621) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.430065) q[1];
sx q[1];
rz(-1.0257755) q[1];
sx q[1];
rz(-1.3793634) q[1];
rz(-pi) q[2];
rz(0.89935715) q[3];
sx q[3];
rz(-1.4345508) q[3];
sx q[3];
rz(-1.5929102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14791791) q[2];
sx q[2];
rz(-2.6213054) q[2];
sx q[2];
rz(-0.49368039) q[2];
rz(2.4364021) q[3];
sx q[3];
rz(-0.54809904) q[3];
sx q[3];
rz(1.3956453) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.482835) q[0];
sx q[0];
rz(-2.5107497) q[0];
sx q[0];
rz(0.054542907) q[0];
rz(1.6627809) q[1];
sx q[1];
rz(-0.80816591) q[1];
sx q[1];
rz(0.85173839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211713) q[0];
sx q[0];
rz(-1.2454438) q[0];
sx q[0];
rz(-2.6634548) q[0];
rz(1.4657945) q[2];
sx q[2];
rz(-1.9622905) q[2];
sx q[2];
rz(-0.016648559) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1143491) q[1];
sx q[1];
rz(-1.7173697) q[1];
sx q[1];
rz(-0.40302009) q[1];
rz(-pi) q[2];
rz(2.4486831) q[3];
sx q[3];
rz(-1.9468782) q[3];
sx q[3];
rz(-0.90510891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30107522) q[2];
sx q[2];
rz(-3.0121351) q[2];
sx q[2];
rz(1.2811309) q[2];
rz(1.8649795) q[3];
sx q[3];
rz(-1.3161148) q[3];
sx q[3];
rz(-2.2006456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304826) q[0];
sx q[0];
rz(-3.0597882) q[0];
sx q[0];
rz(1.4962037) q[0];
rz(1.243535) q[1];
sx q[1];
rz(-1.4356109) q[1];
sx q[1];
rz(0.03874716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4604038) q[0];
sx q[0];
rz(-0.28875586) q[0];
sx q[0];
rz(2.0577752) q[0];
x q[1];
rz(-2.4960802) q[2];
sx q[2];
rz(-0.61290462) q[2];
sx q[2];
rz(0.77790368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8643724) q[1];
sx q[1];
rz(-0.38063875) q[1];
sx q[1];
rz(0.71613272) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0729573) q[3];
sx q[3];
rz(-0.58311948) q[3];
sx q[3];
rz(2.0770962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1425928) q[2];
sx q[2];
rz(-0.6730538) q[2];
sx q[2];
rz(-2.1633945) q[2];
rz(0.081196872) q[3];
sx q[3];
rz(-1.2673763) q[3];
sx q[3];
rz(0.80011884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7112432) q[0];
sx q[0];
rz(-2.5177903) q[0];
sx q[0];
rz(-1.2687564) q[0];
rz(2.2637892) q[1];
sx q[1];
rz(-1.1962079) q[1];
sx q[1];
rz(0.21003221) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9531404) q[0];
sx q[0];
rz(-2.3658381) q[0];
sx q[0];
rz(2.1908568) q[0];
rz(-2.7619475) q[2];
sx q[2];
rz(-1.2894508) q[2];
sx q[2];
rz(-0.18923551) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.91466838) q[1];
sx q[1];
rz(-1.6579501) q[1];
sx q[1];
rz(2.8512849) q[1];
rz(-pi) q[2];
rz(0.29573567) q[3];
sx q[3];
rz(-1.1533168) q[3];
sx q[3];
rz(-1.118286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36748019) q[2];
sx q[2];
rz(-1.3051278) q[2];
sx q[2];
rz(-1.0779862) q[2];
rz(0.20400253) q[3];
sx q[3];
rz(-1.9670468) q[3];
sx q[3];
rz(-1.7417804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.792895) q[0];
sx q[0];
rz(-0.63454151) q[0];
sx q[0];
rz(0.10575159) q[0];
rz(1.3144685) q[1];
sx q[1];
rz(-2.1038838) q[1];
sx q[1];
rz(2.7076142) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26555222) q[0];
sx q[0];
rz(-1.5469451) q[0];
sx q[0];
rz(2.0562459) q[0];
rz(0.94900382) q[2];
sx q[2];
rz(-0.45897537) q[2];
sx q[2];
rz(-1.768809) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.909653) q[1];
sx q[1];
rz(-2.3756333) q[1];
sx q[1];
rz(-3.0161117) q[1];
rz(-pi) q[2];
rz(-1.4922797) q[3];
sx q[3];
rz(-1.8338876) q[3];
sx q[3];
rz(-0.81180868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5964261) q[2];
sx q[2];
rz(-0.20238987) q[2];
sx q[2];
rz(-1.8983967) q[2];
rz(-0.09178129) q[3];
sx q[3];
rz(-1.4877078) q[3];
sx q[3];
rz(1.6091326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6169601) q[0];
sx q[0];
rz(-1.4240823) q[0];
sx q[0];
rz(2.4272163) q[0];
rz(0.83348918) q[1];
sx q[1];
rz(-2.1099213) q[1];
sx q[1];
rz(2.3918236) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3692255) q[0];
sx q[0];
rz(-1.1208659) q[0];
sx q[0];
rz(0.97226592) q[0];
rz(-pi) q[1];
rz(0.60192666) q[2];
sx q[2];
rz(-1.9052124) q[2];
sx q[2];
rz(-2.189881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0088391) q[1];
sx q[1];
rz(-2.3633966) q[1];
sx q[1];
rz(1.9254033) q[1];
rz(1.8014471) q[3];
sx q[3];
rz(-2.8845571) q[3];
sx q[3];
rz(1.5046985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82911503) q[2];
sx q[2];
rz(-1.95936) q[2];
sx q[2];
rz(1.0324837) q[2];
rz(-1.1517582) q[3];
sx q[3];
rz(-1.0843029) q[3];
sx q[3];
rz(1.9657709) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6308003) q[0];
sx q[0];
rz(-1.3277338) q[0];
sx q[0];
rz(3.0603141) q[0];
rz(2.6649113) q[1];
sx q[1];
rz(-1.3414693) q[1];
sx q[1];
rz(-0.016990677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9067334) q[0];
sx q[0];
rz(-1.181266) q[0];
sx q[0];
rz(-0.46789668) q[0];
rz(-pi) q[1];
rz(-2.4912253) q[2];
sx q[2];
rz(-1.128049) q[2];
sx q[2];
rz(-2.3897417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.022804755) q[1];
sx q[1];
rz(-2.9937016) q[1];
sx q[1];
rz(2.1639362) q[1];
rz(-3.105427) q[3];
sx q[3];
rz(-1.6828575) q[3];
sx q[3];
rz(-1.6074635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2396669) q[2];
sx q[2];
rz(-1.6096175) q[2];
sx q[2];
rz(1.3942744) q[2];
rz(1.4782921) q[3];
sx q[3];
rz(-0.52634382) q[3];
sx q[3];
rz(-0.31005508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8095241) q[0];
sx q[0];
rz(-2.9991751) q[0];
sx q[0];
rz(-0.19752565) q[0];
rz(1.7212414) q[1];
sx q[1];
rz(-2.0887801) q[1];
sx q[1];
rz(-1.24409) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0164393) q[0];
sx q[0];
rz(-2.2028515) q[0];
sx q[0];
rz(1.7470361) q[0];
rz(0.19322531) q[2];
sx q[2];
rz(-2.3086289) q[2];
sx q[2];
rz(-2.0173313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1411775) q[1];
sx q[1];
rz(-1.1338935) q[1];
sx q[1];
rz(-1.2440597) q[1];
rz(-0.029989244) q[3];
sx q[3];
rz(-0.99195671) q[3];
sx q[3];
rz(-2.8473583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55599219) q[2];
sx q[2];
rz(-2.8828794) q[2];
sx q[2];
rz(0.82536215) q[2];
rz(-2.182492) q[3];
sx q[3];
rz(-0.93950713) q[3];
sx q[3];
rz(1.8382898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319204) q[0];
sx q[0];
rz(-0.29255022) q[0];
sx q[0];
rz(-1.1462611) q[0];
rz(1.7068901) q[1];
sx q[1];
rz(-2.2271894) q[1];
sx q[1];
rz(0.31417876) q[1];
rz(1.0028253) q[2];
sx q[2];
rz(-0.88063721) q[2];
sx q[2];
rz(0.70692393) q[2];
rz(-2.0758177) q[3];
sx q[3];
rz(-1.4151971) q[3];
sx q[3];
rz(-0.74756037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
