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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632089) q[0];
sx q[0];
rz(-1.4683405) q[0];
sx q[0];
rz(1.5200903) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92496867) q[2];
sx q[2];
rz(-2.2708974) q[2];
sx q[2];
rz(-2.2440804) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91115921) q[1];
sx q[1];
rz(-1.8130366) q[1];
sx q[1];
rz(-1.6723067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1185094) q[3];
sx q[3];
rz(-1.8798401) q[3];
sx q[3];
rz(-3.0024101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3377043) q[2];
sx q[2];
rz(-1.8053728) q[2];
sx q[2];
rz(-1.814092) q[2];
rz(1.368807) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(2.9201065) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9889481) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(3.0864571) q[0];
rz(2.4368743) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(-1.0447186) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8880312) q[0];
sx q[0];
rz(-0.46875254) q[0];
sx q[0];
rz(-2.5028489) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80514812) q[2];
sx q[2];
rz(-2.3986446) q[2];
sx q[2];
rz(2.2066903) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5767908) q[1];
sx q[1];
rz(-1.9783486) q[1];
sx q[1];
rz(2.3840586) q[1];
x q[2];
rz(1.6325266) q[3];
sx q[3];
rz(-1.6654019) q[3];
sx q[3];
rz(2.4199362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3220871) q[2];
sx q[2];
rz(-2.4884067) q[2];
sx q[2];
rz(-1.8196003) q[2];
rz(-2.3540438) q[3];
sx q[3];
rz(-1.7371663) q[3];
sx q[3];
rz(2.0866086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2128485) q[0];
sx q[0];
rz(-2.4725914) q[0];
sx q[0];
rz(0.71841946) q[0];
rz(2.8058715) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(-0.95300037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9586784) q[0];
sx q[0];
rz(-1.0299068) q[0];
sx q[0];
rz(-0.042634115) q[0];
rz(-pi) q[1];
rz(2.0905054) q[2];
sx q[2];
rz(-1.698709) q[2];
sx q[2];
rz(-3.0232883) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.53297199) q[1];
sx q[1];
rz(-0.91580456) q[1];
sx q[1];
rz(-0.47651948) q[1];
rz(-pi) q[2];
rz(0.6742874) q[3];
sx q[3];
rz(-2.4113048) q[3];
sx q[3];
rz(1.8411318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0042051729) q[2];
sx q[2];
rz(-1.1342099) q[2];
sx q[2];
rz(-0.65286621) q[2];
rz(-0.79814923) q[3];
sx q[3];
rz(-1.3099202) q[3];
sx q[3];
rz(2.4979512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2586486) q[0];
sx q[0];
rz(-2.301321) q[0];
sx q[0];
rz(-0.81664455) q[0];
rz(2.4500997) q[1];
sx q[1];
rz(-1.7934727) q[1];
sx q[1];
rz(-0.74849558) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2484528) q[0];
sx q[0];
rz(-1.7692385) q[0];
sx q[0];
rz(2.8273185) q[0];
rz(-1.9709936) q[2];
sx q[2];
rz(-0.23175254) q[2];
sx q[2];
rz(-1.0585001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4481727) q[1];
sx q[1];
rz(-2.2058371) q[1];
sx q[1];
rz(-2.5513812) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17552142) q[3];
sx q[3];
rz(-1.3427882) q[3];
sx q[3];
rz(2.1057749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7526492) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(-2.3469817) q[2];
rz(1.7157582) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(0.37253255) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36885095) q[0];
sx q[0];
rz(-2.0617101) q[0];
sx q[0];
rz(2.6439731) q[0];
rz(-2.3720062) q[1];
sx q[1];
rz(-1.1520569) q[1];
sx q[1];
rz(2.41113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.045864) q[0];
sx q[0];
rz(-2.2501643) q[0];
sx q[0];
rz(-1.2936399) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94295741) q[2];
sx q[2];
rz(-1.3636317) q[2];
sx q[2];
rz(3.0414273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9047315) q[1];
sx q[1];
rz(-1.3610468) q[1];
sx q[1];
rz(0.81717234) q[1];
rz(-pi) q[2];
rz(-2.6914207) q[3];
sx q[3];
rz(-2.014959) q[3];
sx q[3];
rz(1.1314363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89255014) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(-2.6882233) q[2];
rz(1.310965) q[3];
sx q[3];
rz(-2.9679208) q[3];
sx q[3];
rz(-3.0176676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675568) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(-1.3483082) q[0];
rz(1.096161) q[1];
sx q[1];
rz(-1.4827012) q[1];
sx q[1];
rz(-2.5199264) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.877705) q[0];
sx q[0];
rz(-1.5866205) q[0];
sx q[0];
rz(0.71431922) q[0];
rz(-1.4716343) q[2];
sx q[2];
rz(-1.1767724) q[2];
sx q[2];
rz(2.9539915) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.073879378) q[1];
sx q[1];
rz(-2.3150958) q[1];
sx q[1];
rz(-1.9113051) q[1];
x q[2];
rz(2.4824941) q[3];
sx q[3];
rz(-2.3602398) q[3];
sx q[3];
rz(2.3574061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93214503) q[2];
sx q[2];
rz(-0.79938447) q[2];
sx q[2];
rz(-1.4383593) q[2];
rz(2.1608593) q[3];
sx q[3];
rz(-1.2659729) q[3];
sx q[3];
rz(1.2119306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.257523) q[0];
sx q[0];
rz(-1.7634089) q[0];
sx q[0];
rz(-2.8164465) q[0];
rz(-2.6349321) q[1];
sx q[1];
rz(-2.1279361) q[1];
sx q[1];
rz(1.315717) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611826) q[0];
sx q[0];
rz(-1.4482486) q[0];
sx q[0];
rz(-2.2558801) q[0];
rz(-pi) q[1];
rz(-0.94177926) q[2];
sx q[2];
rz(-0.67275362) q[2];
sx q[2];
rz(0.37601177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1707304) q[1];
sx q[1];
rz(-1.9625753) q[1];
sx q[1];
rz(2.4337593) q[1];
x q[2];
rz(2.2030284) q[3];
sx q[3];
rz(-1.5364416) q[3];
sx q[3];
rz(-1.9711348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.331984) q[2];
sx q[2];
rz(-1.0386244) q[2];
sx q[2];
rz(1.6161551) q[2];
rz(-1.3574613) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-2.0048001) q[1];
sx q[1];
rz(1.3884707) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3102458) q[0];
sx q[0];
rz(-0.60605907) q[0];
sx q[0];
rz(1.5008657) q[0];
rz(-pi) q[1];
rz(-3.0031707) q[2];
sx q[2];
rz(-1.8571203) q[2];
sx q[2];
rz(0.22655205) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35508868) q[1];
sx q[1];
rz(-1.638371) q[1];
sx q[1];
rz(1.3994292) q[1];
rz(-1.7817333) q[3];
sx q[3];
rz(-1.9550793) q[3];
sx q[3];
rz(2.3944642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93747741) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(-0.80782962) q[2];
rz(2.4701123) q[3];
sx q[3];
rz(-1.6356133) q[3];
sx q[3];
rz(-1.5988662) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143739) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(2.6887023) q[0];
rz(-0.29531404) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(1.7426851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88576372) q[0];
sx q[0];
rz(-1.4867146) q[0];
sx q[0];
rz(-0.83763116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8434502) q[2];
sx q[2];
rz(-2.6116707) q[2];
sx q[2];
rz(2.8255445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6878444) q[1];
sx q[1];
rz(-2.0216989) q[1];
sx q[1];
rz(-1.1503673) q[1];
x q[2];
rz(-0.021965071) q[3];
sx q[3];
rz(-1.3736892) q[3];
sx q[3];
rz(-2.0052912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7704775) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(2.5650043) q[2];
rz(-0.21814403) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(-1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52467728) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(-3.0531378) q[0];
rz(0.5312008) q[1];
sx q[1];
rz(-2.7374659) q[1];
sx q[1];
rz(-1.6204576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3192465) q[0];
sx q[0];
rz(-2.7392031) q[0];
sx q[0];
rz(-1.3565543) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.047777339) q[2];
sx q[2];
rz(-2.0476404) q[2];
sx q[2];
rz(-2.5587815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4978334) q[1];
sx q[1];
rz(-1.6439142) q[1];
sx q[1];
rz(2.916479) q[1];
rz(-pi) q[2];
rz(-0.57595595) q[3];
sx q[3];
rz(-0.68000845) q[3];
sx q[3];
rz(0.81871225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8497808) q[2];
sx q[2];
rz(-1.3498638) q[2];
sx q[2];
rz(1.527171) q[2];
rz(1.269086) q[3];
sx q[3];
rz(-2.1554155) q[3];
sx q[3];
rz(-1.0482949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0328746) q[0];
sx q[0];
rz(-2.0357108) q[0];
sx q[0];
rz(2.3172814) q[0];
rz(2.6678008) q[1];
sx q[1];
rz(-1.5722678) q[1];
sx q[1];
rz(-1.5617465) q[1];
rz(-1.9971581) q[2];
sx q[2];
rz(-1.016166) q[2];
sx q[2];
rz(-0.24220256) q[2];
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
