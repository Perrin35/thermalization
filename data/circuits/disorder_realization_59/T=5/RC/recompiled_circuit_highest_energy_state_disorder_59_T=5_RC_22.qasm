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
rz(-2.9474592) q[0];
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
rz(2.8391957) q[0];
sx q[0];
rz(-1.6212362) q[0];
sx q[0];
rz(3.0390059) q[0];
rz(-pi) q[1];
rz(-0.62032932) q[2];
sx q[2];
rz(-2.2278068) q[2];
sx q[2];
rz(3.1075396) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4575263) q[1];
sx q[1];
rz(-1.6693331) q[1];
sx q[1];
rz(0.24344488) q[1];
rz(0.93985112) q[3];
sx q[3];
rz(-0.54169023) q[3];
sx q[3];
rz(-1.9909137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3377043) q[2];
sx q[2];
rz(-1.8053728) q[2];
sx q[2];
rz(-1.814092) q[2];
rz(1.368807) q[3];
sx q[3];
rz(-0.33392206) q[3];
sx q[3];
rz(-0.22148618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1526445) q[0];
sx q[0];
rz(-1.4168318) q[0];
sx q[0];
rz(-3.0864571) q[0];
rz(2.4368743) q[1];
sx q[1];
rz(-1.5780508) q[1];
sx q[1];
rz(-2.096874) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8880312) q[0];
sx q[0];
rz(-2.6728401) q[0];
sx q[0];
rz(-0.6387438) q[0];
rz(-0.98590322) q[2];
sx q[2];
rz(-2.0587181) q[2];
sx q[2];
rz(3.121738) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56480184) q[1];
sx q[1];
rz(-1.163244) q[1];
sx q[1];
rz(2.3840586) q[1];
x q[2];
rz(3.0468076) q[3];
sx q[3];
rz(-1.6322502) q[3];
sx q[3];
rz(-2.2866142) q[3];
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
rz(-1.4044263) q[3];
sx q[3];
rz(-2.0866086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2128485) q[0];
sx q[0];
rz(-0.66900122) q[0];
sx q[0];
rz(-0.71841946) q[0];
rz(2.8058715) q[1];
sx q[1];
rz(-1.6751143) q[1];
sx q[1];
rz(-0.95300037) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.876016) q[0];
sx q[0];
rz(-0.5424005) q[0];
sx q[0];
rz(-1.6416373) q[0];
rz(-pi) q[1];
rz(1.3173872) q[2];
sx q[2];
rz(-0.53381397) q[2];
sx q[2];
rz(-1.6718503) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.408421) q[1];
sx q[1];
rz(-1.1984899) q[1];
sx q[1];
rz(2.2835963) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5312214) q[3];
sx q[3];
rz(-1.1412176) q[3];
sx q[3];
rz(2.3342032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1373875) q[2];
sx q[2];
rz(-2.0073828) q[2];
sx q[2];
rz(0.65286621) q[2];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.882944) q[0];
sx q[0];
rz(-2.301321) q[0];
sx q[0];
rz(2.3249481) q[0];
rz(-0.69149292) q[1];
sx q[1];
rz(-1.34812) q[1];
sx q[1];
rz(0.74849558) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25835055) q[0];
sx q[0];
rz(-1.2628947) q[0];
sx q[0];
rz(-1.3624205) q[0];
rz(-1.7848133) q[2];
sx q[2];
rz(-1.6604009) q[2];
sx q[2];
rz(-0.90286189) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.88598) q[1];
sx q[1];
rz(-1.1062262) q[1];
sx q[1];
rz(-0.8452615) q[1];
x q[2];
rz(-2.2159202) q[3];
sx q[3];
rz(-2.8547849) q[3];
sx q[3];
rz(-1.4405026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7526492) q[2];
sx q[2];
rz(-0.77072531) q[2];
sx q[2];
rz(0.79461092) q[2];
rz(-1.7157582) q[3];
sx q[3];
rz(-2.1013997) q[3];
sx q[3];
rz(2.7690601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(0.76958641) q[1];
sx q[1];
rz(-1.9895357) q[1];
sx q[1];
rz(0.73046267) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896547) q[0];
sx q[0];
rz(-1.3562855) q[0];
sx q[0];
rz(-0.69846054) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.914417) q[2];
sx q[2];
rz(-0.65672749) q[2];
sx q[2];
rz(-1.7467787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2368612) q[1];
sx q[1];
rz(-1.3610468) q[1];
sx q[1];
rz(-2.3244203) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6914207) q[3];
sx q[3];
rz(-2.014959) q[3];
sx q[3];
rz(-2.0101563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2490425) q[2];
sx q[2];
rz(-1.5183134) q[2];
sx q[2];
rz(-2.6882233) q[2];
rz(-1.8306277) q[3];
sx q[3];
rz(-0.1736719) q[3];
sx q[3];
rz(-0.12392509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740358) q[0];
sx q[0];
rz(-2.0581764) q[0];
sx q[0];
rz(1.3483082) q[0];
rz(1.096161) q[1];
sx q[1];
rz(-1.6588914) q[1];
sx q[1];
rz(-0.62166628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28866119) q[0];
sx q[0];
rz(-0.71446361) q[0];
sx q[0];
rz(3.1174401) q[0];
rz(-pi) q[1];
rz(1.4716343) q[2];
sx q[2];
rz(-1.1767724) q[2];
sx q[2];
rz(-2.9539915) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0677133) q[1];
sx q[1];
rz(-2.3150958) q[1];
sx q[1];
rz(1.2302876) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(0.93214503) q[2];
sx q[2];
rz(-2.3422082) q[2];
sx q[2];
rz(-1.7032334) q[2];
rz(0.98073331) q[3];
sx q[3];
rz(-1.8756198) q[3];
sx q[3];
rz(-1.929662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88406968) q[0];
sx q[0];
rz(-1.3781837) q[0];
sx q[0];
rz(-0.32514611) q[0];
rz(-2.6349321) q[1];
sx q[1];
rz(-1.0136565) q[1];
sx q[1];
rz(-1.315717) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0611826) q[0];
sx q[0];
rz(-1.6933441) q[0];
sx q[0];
rz(2.2558801) q[0];
x q[1];
rz(-0.43834941) q[2];
sx q[2];
rz(-1.0427156) q[2];
sx q[2];
rz(2.768289) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81974492) q[1];
sx q[1];
rz(-2.3492866) q[1];
sx q[1];
rz(-0.56604947) q[1];
rz(0.93856424) q[3];
sx q[3];
rz(-1.5364416) q[3];
sx q[3];
rz(1.9711348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.331984) q[2];
sx q[2];
rz(-1.0386244) q[2];
sx q[2];
rz(1.6161551) q[2];
rz(-1.7841313) q[3];
sx q[3];
rz(-1.6922035) q[3];
sx q[3];
rz(-1.6239032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9759393) q[0];
sx q[0];
rz(-1.140056) q[0];
sx q[0];
rz(-0.64919382) q[0];
rz(0.80750418) q[1];
sx q[1];
rz(-2.0048001) q[1];
sx q[1];
rz(1.3884707) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3102458) q[0];
sx q[0];
rz(-2.5355336) q[0];
sx q[0];
rz(-1.5008657) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0090603) q[2];
sx q[2];
rz(-2.8243941) q[2];
sx q[2];
rz(2.4567921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5877046) q[1];
sx q[1];
rz(-0.18408751) q[1];
sx q[1];
rz(-1.192993) q[1];
x q[2];
rz(2.6638373) q[3];
sx q[3];
rz(-0.43583187) q[3];
sx q[3];
rz(-1.2660668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93747741) q[2];
sx q[2];
rz(-0.14675879) q[2];
sx q[2];
rz(2.333763) q[2];
rz(-0.67148036) q[3];
sx q[3];
rz(-1.6356133) q[3];
sx q[3];
rz(-1.5988662) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9143739) q[0];
sx q[0];
rz(-2.6052573) q[0];
sx q[0];
rz(-2.6887023) q[0];
rz(2.8462786) q[1];
sx q[1];
rz(-1.9120049) q[1];
sx q[1];
rz(1.3989075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2558289) q[0];
sx q[0];
rz(-1.4867146) q[0];
sx q[0];
rz(-0.83763116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0844551) q[2];
sx q[2];
rz(-1.4342564) q[2];
sx q[2];
rz(1.4914907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4860515) q[1];
sx q[1];
rz(-2.5351346) q[1];
sx q[1];
rz(0.70041366) q[1];
rz(-3.1196276) q[3];
sx q[3];
rz(-1.7679035) q[3];
sx q[3];
rz(1.1363014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37111515) q[2];
sx q[2];
rz(-0.8728084) q[2];
sx q[2];
rz(0.57658833) q[2];
rz(2.9234486) q[3];
sx q[3];
rz(-0.79272565) q[3];
sx q[3];
rz(1.920759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6169154) q[0];
sx q[0];
rz(-0.24743947) q[0];
sx q[0];
rz(-0.088454811) q[0];
rz(0.5312008) q[1];
sx q[1];
rz(-2.7374659) q[1];
sx q[1];
rz(1.521135) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82234612) q[0];
sx q[0];
rz(-2.7392031) q[0];
sx q[0];
rz(1.3565543) q[0];
x q[1];
rz(-3.0938153) q[2];
sx q[2];
rz(-1.0939523) q[2];
sx q[2];
rz(-2.5587815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6437592) q[1];
sx q[1];
rz(-1.4976785) q[1];
sx q[1];
rz(-2.916479) q[1];
rz(-pi) q[2];
rz(0.59595396) q[3];
sx q[3];
rz(-1.2212545) q[3];
sx q[3];
rz(-2.8571125) q[3];
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
rz(-1.269086) q[3];
sx q[3];
rz(-0.98617712) q[3];
sx q[3];
rz(2.0932978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
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
rz(1.0328746) q[0];
sx q[0];
rz(-2.0357108) q[0];
sx q[0];
rz(2.3172814) q[0];
rz(-0.4737919) q[1];
sx q[1];
rz(-1.5722678) q[1];
sx q[1];
rz(-1.5617465) q[1];
rz(2.5529591) q[2];
sx q[2];
rz(-2.4559173) q[2];
sx q[2];
rz(0.46951175) q[2];
rz(-0.80120091) q[3];
sx q[3];
rz(-2.4036435) q[3];
sx q[3];
rz(2.5477464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
