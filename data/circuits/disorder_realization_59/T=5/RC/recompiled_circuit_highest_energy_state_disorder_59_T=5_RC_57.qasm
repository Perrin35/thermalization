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
rz(-2.0464719) q[1];
sx q[1];
rz(-1.4580589) q[1];
sx q[1];
rz(-2.3966052) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4177822) q[0];
sx q[0];
rz(-0.11427721) q[0];
sx q[0];
rz(-2.6835915) q[0];
x q[1];
rz(2.5212633) q[2];
sx q[2];
rz(-0.91378586) q[2];
sx q[2];
rz(0.034053085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50959831) q[1];
sx q[1];
rz(-2.879329) q[1];
sx q[1];
rz(0.38919501) q[1];
rz(-pi) q[2];
rz(-0.3410913) q[3];
sx q[3];
rz(-2.00019) q[3];
sx q[3];
rz(1.856696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8038883) q[2];
sx q[2];
rz(-1.8053728) q[2];
sx q[2];
rz(-1.3275006) q[2];
rz(-1.7727857) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(2.9201065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889481) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(-3.0864571) q[0];
rz(2.4368743) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(2.096874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2392512) q[0];
sx q[0];
rz(-1.2980868) q[0];
sx q[0];
rz(-0.38614778) q[0];
x q[1];
rz(-2.3364445) q[2];
sx q[2];
rz(-2.3986446) q[2];
sx q[2];
rz(2.2066903) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4942532) q[1];
sx q[1];
rz(-2.2534465) q[1];
sx q[1];
rz(-1.0346336) q[1];
rz(-pi) q[2];
rz(0.09478507) q[3];
sx q[3];
rz(-1.5093424) q[3];
sx q[3];
rz(-2.2866142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.81950554) q[2];
sx q[2];
rz(-2.4884067) q[2];
sx q[2];
rz(1.8196003) q[2];
rz(0.78754887) q[3];
sx q[3];
rz(-1.7371663) q[3];
sx q[3];
rz(-1.0549841) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2128485) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(0.71841946) q[0];
rz(2.8058715) q[1];
sx q[1];
rz(-1.4664783) q[1];
sx q[1];
rz(0.95300037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3659201) q[0];
sx q[0];
rz(-1.5342511) q[0];
sx q[0];
rz(-2.1120872) q[0];
x q[1];
rz(-0.14710958) q[2];
sx q[2];
rz(-1.0557555) q[2];
sx q[2];
rz(-1.7619606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2359537) q[1];
sx q[1];
rz(-0.78879005) q[1];
sx q[1];
rz(1.0324638) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0805667) q[3];
sx q[3];
rz(-2.1189287) q[3];
sx q[3];
rz(0.47993101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1373875) q[2];
sx q[2];
rz(-2.0073828) q[2];
sx q[2];
rz(2.4887264) q[2];
rz(-2.3434434) q[3];
sx q[3];
rz(-1.3099202) q[3];
sx q[3];
rz(-2.4979512) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.882944) q[0];
sx q[0];
rz(-2.301321) q[0];
sx q[0];
rz(-2.3249481) q[0];
rz(-0.69149292) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(0.74849558) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2740606) q[0];
sx q[0];
rz(-0.36991773) q[0];
sx q[0];
rz(0.57673133) q[0];
x q[1];
rz(-0.091684999) q[2];
sx q[2];
rz(-1.7839415) q[2];
sx q[2];
rz(-0.64848778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.88598) q[1];
sx q[1];
rz(-2.0353664) q[1];
sx q[1];
rz(0.8452615) q[1];
rz(-2.9660712) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(-1.0358178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7526492) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(-0.79461092) q[2];
rz(-1.4258344) q[3];
sx q[3];
rz(-1.040193) q[3];
sx q[3];
rz(-0.37253255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36885095) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(2.6439731) q[0];
rz(-0.76958641) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(-0.73046267) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4896547) q[0];
sx q[0];
rz(-1.7853072) q[0];
sx q[0];
rz(2.4431321) q[0];
rz(-pi) q[1];
rz(0.25409296) q[2];
sx q[2];
rz(-2.1832018) q[2];
sx q[2];
rz(1.3224365) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0260057) q[1];
sx q[1];
rz(-0.77669826) q[1];
sx q[1];
rz(1.8724043) q[1];
rz(1.0845844) q[3];
sx q[3];
rz(-1.9745805) q[3];
sx q[3];
rz(-0.23460925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.89255014) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(-2.6882233) q[2];
rz(-1.8306277) q[3];
sx q[3];
rz(-2.9679208) q[3];
sx q[3];
rz(-3.0176676) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3675568) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(-1.7932844) q[0];
rz(1.096161) q[1];
sx q[1];
rz(-1.4827012) q[1];
sx q[1];
rz(0.62166628) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2638877) q[0];
sx q[0];
rz(-1.5866205) q[0];
sx q[0];
rz(-0.71431922) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4716343) q[2];
sx q[2];
rz(-1.1767724) q[2];
sx q[2];
rz(2.9539915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.409103) q[1];
sx q[1];
rz(-1.81899) q[1];
sx q[1];
rz(-0.77381882) q[1];
x q[2];
rz(-0.65909855) q[3];
sx q[3];
rz(-2.3602398) q[3];
sx q[3];
rz(2.3574061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2094476) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(-1.7032334) q[2];
rz(-0.98073331) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.257523) q[0];
sx q[0];
rz(-1.7634089) q[0];
sx q[0];
rz(0.32514611) q[0];
rz(-2.6349321) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(1.8258757) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5324253) q[0];
sx q[0];
rz(-2.2497674) q[0];
sx q[0];
rz(-2.9838613) q[0];
x q[1];
rz(-2.1431214) q[2];
sx q[2];
rz(-1.9461759) q[2];
sx q[2];
rz(-1.7121512) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3218477) q[1];
sx q[1];
rz(-2.3492866) q[1];
sx q[1];
rz(2.5755432) q[1];
rz(1.6288888) q[3];
sx q[3];
rz(-0.63303715) q[3];
sx q[3];
rz(2.7881088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.331984) q[2];
sx q[2];
rz(-2.1029682) q[2];
sx q[2];
rz(-1.5254376) q[2];
rz(1.7841313) q[3];
sx q[3];
rz(-1.4493891) q[3];
sx q[3];
rz(-1.6239032) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759393) q[0];
sx q[0];
rz(-2.0015367) q[0];
sx q[0];
rz(-0.64919382) q[0];
rz(2.3340885) q[1];
sx q[1];
rz(-2.0048001) q[1];
sx q[1];
rz(1.753122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20304414) q[0];
sx q[0];
rz(-1.5309835) q[0];
sx q[0];
rz(-0.96588246) q[0];
rz(-1.2818579) q[2];
sx q[2];
rz(-1.4380437) q[2];
sx q[2];
rz(-1.3049187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9375704) q[1];
sx q[1];
rz(-1.7417684) q[1];
sx q[1];
rz(-0.06857605) q[1];
x q[2];
rz(-2.7494568) q[3];
sx q[3];
rz(-1.7661422) q[3];
sx q[3];
rz(0.74356642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.93747741) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(0.80782962) q[2];
rz(2.4701123) q[3];
sx q[3];
rz(-1.5059794) q[3];
sx q[3];
rz(1.5988662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272187) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(-2.6887023) q[0];
rz(-2.8462786) q[1];
sx q[1];
rz(-1.9120049) q[1];
sx q[1];
rz(1.7426851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5495396) q[0];
sx q[0];
rz(-0.73707923) q[0];
sx q[0];
rz(1.4455185) q[0];
rz(-1.0571376) q[2];
sx q[2];
rz(-1.7073362) q[2];
sx q[2];
rz(-1.650102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4860515) q[1];
sx q[1];
rz(-0.60645804) q[1];
sx q[1];
rz(2.441179) q[1];
rz(-1.3736428) q[3];
sx q[3];
rz(-1.592336) q[3];
sx q[3];
rz(0.43019274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37111515) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(2.5650043) q[2];
rz(0.21814403) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(-1.2208337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(3.0531378) q[0];
rz(-0.5312008) q[1];
sx q[1];
rz(-2.7374659) q[1];
sx q[1];
rz(1.6204576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94604055) q[0];
sx q[0];
rz(-1.6541535) q[0];
sx q[0];
rz(-1.1766737) q[0];
rz(-pi) q[1];
rz(2.0481061) q[2];
sx q[2];
rz(-1.6132406) q[2];
sx q[2];
rz(2.1316656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6437592) q[1];
sx q[1];
rz(-1.4976785) q[1];
sx q[1];
rz(0.22511367) q[1];
x q[2];
rz(1.1559256) q[3];
sx q[3];
rz(-2.1262827) q[3];
sx q[3];
rz(-1.6270571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8497808) q[2];
sx q[2];
rz(-1.3498638) q[2];
sx q[2];
rz(1.6144217) q[2];
rz(-1.8725066) q[3];
sx q[3];
rz(-0.98617712) q[3];
sx q[3];
rz(1.0482949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.108718) q[0];
sx q[0];
rz(-1.1058818) q[0];
sx q[0];
rz(-0.82431128) q[0];
rz(2.6678008) q[1];
sx q[1];
rz(-1.5722678) q[1];
sx q[1];
rz(-1.5617465) q[1];
rz(0.59745477) q[2];
sx q[2];
rz(-1.2115527) q[2];
sx q[2];
rz(1.5633898) q[2];
rz(2.1493323) q[3];
sx q[3];
rz(-2.0579865) q[3];
sx q[3];
rz(-1.5427114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
