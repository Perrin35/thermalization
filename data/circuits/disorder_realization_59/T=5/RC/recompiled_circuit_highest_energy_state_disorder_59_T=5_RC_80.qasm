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
rz(-2.0191963) q[0];
sx q[0];
rz(2.9474592) q[0];
rz(-2.0464719) q[1];
sx q[1];
rz(-1.4580589) q[1];
sx q[1];
rz(-2.3966052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632089) q[0];
sx q[0];
rz(-1.4683405) q[0];
sx q[0];
rz(-1.5200903) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.92496867) q[2];
sx q[2];
rz(-2.2708974) q[2];
sx q[2];
rz(-0.89751228) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4575263) q[1];
sx q[1];
rz(-1.4722595) q[1];
sx q[1];
rz(0.24344488) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3410913) q[3];
sx q[3];
rz(-1.1414027) q[3];
sx q[3];
rz(1.2848967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3377043) q[2];
sx q[2];
rz(-1.3362198) q[2];
sx q[2];
rz(1.814092) q[2];
rz(-1.368807) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(-2.9201065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1526445) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(-3.0864571) q[0];
rz(-0.70471835) q[1];
sx q[1];
rz(-1.5635419) q[1];
sx q[1];
rz(2.096874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7009788) q[0];
sx q[0];
rz(-1.9419646) q[0];
sx q[0];
rz(1.2775902) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98590322) q[2];
sx q[2];
rz(-2.0587181) q[2];
sx q[2];
rz(-0.019854644) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5767908) q[1];
sx q[1];
rz(-1.163244) q[1];
sx q[1];
rz(-0.75753404) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.509066) q[3];
sx q[3];
rz(-1.6654019) q[3];
sx q[3];
rz(-0.7216565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3220871) q[2];
sx q[2];
rz(-0.6531859) q[2];
sx q[2];
rz(1.3219924) q[2];
rz(0.78754887) q[3];
sx q[3];
rz(-1.7371663) q[3];
sx q[3];
rz(-1.0549841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9287441) q[0];
sx q[0];
rz(-2.4725914) q[0];
sx q[0];
rz(0.71841946) q[0];
rz(-0.33572117) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(2.1885923) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26557661) q[0];
sx q[0];
rz(-0.5424005) q[0];
sx q[0];
rz(1.6416373) q[0];
rz(1.0510873) q[2];
sx q[2];
rz(-1.4428836) q[2];
sx q[2];
rz(0.11830439) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6086207) q[1];
sx q[1];
rz(-2.2257881) q[1];
sx q[1];
rz(-2.6650732) q[1];
rz(-pi) q[2];
rz(-0.6742874) q[3];
sx q[3];
rz(-2.4113048) q[3];
sx q[3];
rz(-1.8411318) q[3];
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
rz(0.79814923) q[3];
sx q[3];
rz(-1.3099202) q[3];
sx q[3];
rz(-2.4979512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.7934727) q[1];
sx q[1];
rz(-0.74849558) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2740606) q[0];
sx q[0];
rz(-2.7716749) q[0];
sx q[0];
rz(-2.5648613) q[0];
x q[1];
rz(3.0499077) q[2];
sx q[2];
rz(-1.3576512) q[2];
sx q[2];
rz(-2.4931049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.15259057) q[1];
sx q[1];
rz(-2.3034597) q[1];
sx q[1];
rz(0.92392606) q[1];
x q[2];
rz(1.8022376) q[3];
sx q[3];
rz(-1.7417296) q[3];
sx q[3];
rz(-2.6466796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3889435) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(2.3469817) q[2];
rz(1.7157582) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(0.37253255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-0.36885095) q[0];
sx q[0];
rz(-1.0798825) q[0];
sx q[0];
rz(2.6439731) q[0];
rz(0.76958641) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(-2.41113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8121634) q[0];
sx q[0];
rz(-0.72532988) q[0];
sx q[0];
rz(-0.32666387) q[0];
rz(0.94295741) q[2];
sx q[2];
rz(-1.7779609) q[2];
sx q[2];
rz(3.0414273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0260057) q[1];
sx q[1];
rz(-2.3648944) q[1];
sx q[1];
rz(1.8724043) q[1];
rz(-pi) q[2];
rz(-1.0845844) q[3];
sx q[3];
rz(-1.1670122) q[3];
sx q[3];
rz(2.9069834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2490425) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(0.45336938) q[2];
rz(1.310965) q[3];
sx q[3];
rz(-2.9679208) q[3];
sx q[3];
rz(0.12392509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675568) q[0];
sx q[0];
rz(-1.0834162) q[0];
sx q[0];
rz(1.3483082) q[0];
rz(-2.0454316) q[1];
sx q[1];
rz(-1.6588914) q[1];
sx q[1];
rz(2.5199264) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2638877) q[0];
sx q[0];
rz(-1.5549721) q[0];
sx q[0];
rz(-0.71431922) q[0];
rz(-pi) q[1];
rz(2.7458199) q[2];
sx q[2];
rz(-1.6623375) q[2];
sx q[2];
rz(-1.3450194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7335947) q[1];
sx q[1];
rz(-0.80469614) q[1];
sx q[1];
rz(2.7937274) q[1];
rz(0.65909855) q[3];
sx q[3];
rz(-2.3602398) q[3];
sx q[3];
rz(-2.3574061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2094476) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(-1.4383593) q[2];
rz(-2.1608593) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.929662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88406968) q[0];
sx q[0];
rz(-1.7634089) q[0];
sx q[0];
rz(-0.32514611) q[0];
rz(0.50666058) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(1.8258757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60916736) q[0];
sx q[0];
rz(-0.89182527) q[0];
sx q[0];
rz(-2.9838613) q[0];
x q[1];
rz(0.43834941) q[2];
sx q[2];
rz(-2.0988771) q[2];
sx q[2];
rz(-0.37330369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.81974492) q[1];
sx q[1];
rz(-2.3492866) q[1];
sx q[1];
rz(-0.56604947) q[1];
x q[2];
rz(-1.6288888) q[3];
sx q[3];
rz(-0.63303715) q[3];
sx q[3];
rz(0.35348383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80960861) q[2];
sx q[2];
rz(-1.0386244) q[2];
sx q[2];
rz(1.5254376) q[2];
rz(-1.7841313) q[3];
sx q[3];
rz(-1.6922035) q[3];
sx q[3];
rz(1.5176895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9759393) q[0];
sx q[0];
rz(-2.0015367) q[0];
sx q[0];
rz(-0.64919382) q[0];
rz(-0.80750418) q[1];
sx q[1];
rz(-1.1367926) q[1];
sx q[1];
rz(1.3884707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3952636) q[0];
sx q[0];
rz(-2.1751624) q[0];
sx q[0];
rz(3.0932032) q[0];
rz(-pi) q[1];
rz(-1.8597348) q[2];
sx q[2];
rz(-1.4380437) q[2];
sx q[2];
rz(-1.8366739) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5538881) q[1];
sx q[1];
rz(-0.18408751) q[1];
sx q[1];
rz(-1.9485997) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47775538) q[3];
sx q[3];
rz(-0.43583187) q[3];
sx q[3];
rz(-1.8755258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93747741) q[2];
sx q[2];
rz(-2.9948339) q[2];
sx q[2];
rz(-0.80782962) q[2];
rz(0.67148036) q[3];
sx q[3];
rz(-1.5059794) q[3];
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
rz(-1.9143739) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(-2.6887023) q[0];
rz(-2.8462786) q[1];
sx q[1];
rz(-1.2295877) q[1];
sx q[1];
rz(-1.7426851) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2558289) q[0];
sx q[0];
rz(-1.4867146) q[0];
sx q[0];
rz(2.3039615) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0571376) q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-0.30946839) q[1];
sx q[1];
rz(-1.9469643) q[1];
sx q[1];
rz(-2.6539567) q[1];
rz(0.021965071) q[3];
sx q[3];
rz(-1.7679035) q[3];
sx q[3];
rz(-2.0052912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37111515) q[2];
sx q[2];
rz(-2.2687843) q[2];
sx q[2];
rz(-0.57658833) q[2];
rz(-0.21814403) q[3];
sx q[3];
rz(-2.348867) q[3];
sx q[3];
rz(-1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169154) q[0];
sx q[0];
rz(-2.8941532) q[0];
sx q[0];
rz(0.088454811) q[0];
rz(-0.5312008) q[1];
sx q[1];
rz(-2.7374659) q[1];
sx q[1];
rz(1.6204576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82234612) q[0];
sx q[0];
rz(-2.7392031) q[0];
sx q[0];
rz(1.3565543) q[0];
rz(-pi) q[1];
rz(-1.6629822) q[2];
sx q[2];
rz(-2.6625444) q[2];
sx q[2];
rz(-0.47901115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.056235751) q[1];
sx q[1];
rz(-1.3462945) q[1];
sx q[1];
rz(1.6457998) q[1];
rz(-pi) q[2];
rz(-2.5656367) q[3];
sx q[3];
rz(-0.68000845) q[3];
sx q[3];
rz(2.3228804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29181186) q[2];
sx q[2];
rz(-1.7917289) q[2];
sx q[2];
rz(-1.527171) q[2];
rz(1.8725066) q[3];
sx q[3];
rz(-2.1554155) q[3];
sx q[3];
rz(1.0482949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(2.5441379) q[2];
sx q[2];
rz(-1.93004) q[2];
sx q[2];
rz(-1.5782028) q[2];
rz(2.3403917) q[3];
sx q[3];
rz(-2.4036435) q[3];
sx q[3];
rz(2.5477464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
