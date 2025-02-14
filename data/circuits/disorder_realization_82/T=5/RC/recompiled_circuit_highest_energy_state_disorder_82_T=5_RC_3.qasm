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
rz(-2.1190116) q[0];
sx q[0];
rz(-2.1562161) q[0];
sx q[0];
rz(2.0289542) q[0];
rz(0.39482173) q[1];
sx q[1];
rz(-1.1345175) q[1];
sx q[1];
rz(-1.1910103) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9562839) q[0];
sx q[0];
rz(-0.6260159) q[0];
sx q[0];
rz(0.52395384) q[0];
rz(-pi) q[1];
rz(-1.979492) q[2];
sx q[2];
rz(-1.6922277) q[2];
sx q[2];
rz(1.8885776) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8748693) q[1];
sx q[1];
rz(-1.6376347) q[1];
sx q[1];
rz(0.066796227) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84221091) q[3];
sx q[3];
rz(-2.5274383) q[3];
sx q[3];
rz(2.24005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4224008) q[2];
sx q[2];
rz(-2.9802608) q[2];
sx q[2];
rz(-0.42616978) q[2];
rz(-0.69914493) q[3];
sx q[3];
rz(-1.9482502) q[3];
sx q[3];
rz(-2.7307811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8461269) q[0];
sx q[0];
rz(-2.182425) q[0];
sx q[0];
rz(-2.3058983) q[0];
rz(2.6169547) q[1];
sx q[1];
rz(-1.6181889) q[1];
sx q[1];
rz(1.6297278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1439813) q[0];
sx q[0];
rz(-1.8845692) q[0];
sx q[0];
rz(2.1812516) q[0];
rz(0.76194069) q[2];
sx q[2];
rz(-2.2968869) q[2];
sx q[2];
rz(-1.8619137) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8162093) q[1];
sx q[1];
rz(-1.1126425) q[1];
sx q[1];
rz(1.3558595) q[1];
x q[2];
rz(1.2542975) q[3];
sx q[3];
rz(-2.0997556) q[3];
sx q[3];
rz(-1.0366755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52874804) q[2];
sx q[2];
rz(-1.560644) q[2];
sx q[2];
rz(0.35231248) q[2];
rz(2.6490372) q[3];
sx q[3];
rz(-2.0217321) q[3];
sx q[3];
rz(-2.742384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50017971) q[0];
sx q[0];
rz(-0.68920511) q[0];
sx q[0];
rz(2.8885762) q[0];
rz(0.5223271) q[1];
sx q[1];
rz(-0.42799196) q[1];
sx q[1];
rz(-2.7740251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0899058) q[0];
sx q[0];
rz(-1.5977706) q[0];
sx q[0];
rz(1.5083675) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39607145) q[2];
sx q[2];
rz(-2.1537915) q[2];
sx q[2];
rz(-2.1427659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8525703) q[1];
sx q[1];
rz(-1.3636359) q[1];
sx q[1];
rz(-3.0604355) q[1];
rz(-pi) q[2];
rz(3.0250164) q[3];
sx q[3];
rz(-1.7128782) q[3];
sx q[3];
rz(0.97654479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14105178) q[2];
sx q[2];
rz(-2.8512569) q[2];
sx q[2];
rz(2.7153437) q[2];
rz(-2.352377) q[3];
sx q[3];
rz(-1.1169249) q[3];
sx q[3];
rz(-1.5392019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(2.5696647) q[0];
sx q[0];
rz(-2.8043788) q[0];
sx q[0];
rz(-0.1871044) q[0];
rz(-1.3693753) q[1];
sx q[1];
rz(-1.2792842) q[1];
sx q[1];
rz(0.61417907) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3423796) q[0];
sx q[0];
rz(-2.6465017) q[0];
sx q[0];
rz(-1.6922047) q[0];
rz(-pi) q[1];
rz(-2.1000142) q[2];
sx q[2];
rz(-0.61162156) q[2];
sx q[2];
rz(-0.84866103) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.56127292) q[1];
sx q[1];
rz(-2.4497708) q[1];
sx q[1];
rz(-2.0136334) q[1];
rz(-pi) q[2];
x q[2];
rz(2.92431) q[3];
sx q[3];
rz(-2.7585895) q[3];
sx q[3];
rz(0.10620761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46089178) q[2];
sx q[2];
rz(-0.9047752) q[2];
sx q[2];
rz(1.2659849) q[2];
rz(-2.8547817) q[3];
sx q[3];
rz(-0.81227842) q[3];
sx q[3];
rz(0.12107818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1665961) q[0];
sx q[0];
rz(-0.91355211) q[0];
sx q[0];
rz(-1.1443369) q[0];
rz(0.76885778) q[1];
sx q[1];
rz(-0.91173333) q[1];
sx q[1];
rz(-1.2164046) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74690565) q[0];
sx q[0];
rz(-1.1228859) q[0];
sx q[0];
rz(-2.7120717) q[0];
rz(-pi) q[1];
rz(-1.7038962) q[2];
sx q[2];
rz(-1.6665227) q[2];
sx q[2];
rz(-1.7246703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8647841) q[1];
sx q[1];
rz(-1.6549161) q[1];
sx q[1];
rz(1.2963956) q[1];
rz(-pi) q[2];
rz(-2.7796916) q[3];
sx q[3];
rz(-1.6503157) q[3];
sx q[3];
rz(-1.9504296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6151109) q[2];
sx q[2];
rz(-2.5418126) q[2];
sx q[2];
rz(2.6748924) q[2];
rz(-0.99272054) q[3];
sx q[3];
rz(-1.7388742) q[3];
sx q[3];
rz(-1.7464975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4863131) q[0];
sx q[0];
rz(-1.7236973) q[0];
sx q[0];
rz(-0.41743761) q[0];
rz(0.4256658) q[1];
sx q[1];
rz(-2.3047431) q[1];
sx q[1];
rz(0.72798896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5363184) q[0];
sx q[0];
rz(-2.0173376) q[0];
sx q[0];
rz(2.1552318) q[0];
rz(-pi) q[1];
rz(1.7547771) q[2];
sx q[2];
rz(-2.3021491) q[2];
sx q[2];
rz(1.9860742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83367741) q[1];
sx q[1];
rz(-2.6628072) q[1];
sx q[1];
rz(0.91955863) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40188222) q[3];
sx q[3];
rz(-1.0392351) q[3];
sx q[3];
rz(-1.8762003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41842469) q[2];
sx q[2];
rz(-1.6614513) q[2];
sx q[2];
rz(-1.9815014) q[2];
rz(-0.55231071) q[3];
sx q[3];
rz(-0.6902802) q[3];
sx q[3];
rz(-2.6753329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1138678) q[0];
sx q[0];
rz(-2.213527) q[0];
sx q[0];
rz(3.0357251) q[0];
rz(-1.8766807) q[1];
sx q[1];
rz(-0.55095878) q[1];
sx q[1];
rz(-1.6884621) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4111076) q[0];
sx q[0];
rz(-0.8115226) q[0];
sx q[0];
rz(1.1280035) q[0];
rz(-pi) q[1];
rz(-1.7488453) q[2];
sx q[2];
rz(-1.4499003) q[2];
sx q[2];
rz(-2.4089782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8034722) q[1];
sx q[1];
rz(-2.8104467) q[1];
sx q[1];
rz(-1.1263333) q[1];
rz(-2.9960521) q[3];
sx q[3];
rz(-0.91314935) q[3];
sx q[3];
rz(-1.4484792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.037131) q[2];
sx q[2];
rz(-1.1300602) q[2];
sx q[2];
rz(-1.5360443) q[2];
rz(1.7724841) q[3];
sx q[3];
rz(-2.5684999) q[3];
sx q[3];
rz(1.8039186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8417514) q[0];
sx q[0];
rz(-1.890269) q[0];
sx q[0];
rz(1.1773671) q[0];
rz(-1.3031225) q[1];
sx q[1];
rz(-1.6603371) q[1];
sx q[1];
rz(-0.67828137) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2934351) q[0];
sx q[0];
rz(-0.57852902) q[0];
sx q[0];
rz(1.0941157) q[0];
rz(-1.6839333) q[2];
sx q[2];
rz(-1.8461421) q[2];
sx q[2];
rz(-2.6983698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5694414) q[1];
sx q[1];
rz(-1.2673693) q[1];
sx q[1];
rz(-1.7904758) q[1];
x q[2];
rz(2.5362064) q[3];
sx q[3];
rz(-1.6035214) q[3];
sx q[3];
rz(-2.721019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9699041) q[2];
sx q[2];
rz(-2.6320612) q[2];
sx q[2];
rz(2.2129538) q[2];
rz(-1.6396133) q[3];
sx q[3];
rz(-2.00311) q[3];
sx q[3];
rz(-1.5750711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.11012) q[0];
sx q[0];
rz(-2.0125772) q[0];
sx q[0];
rz(1.336115) q[0];
rz(-2.8307092) q[1];
sx q[1];
rz(-1.5300405) q[1];
sx q[1];
rz(-1.329604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3361137) q[0];
sx q[0];
rz(-0.74148899) q[0];
sx q[0];
rz(0.51261897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2822658) q[2];
sx q[2];
rz(-2.2999249) q[2];
sx q[2];
rz(-2.2657889) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3327766) q[1];
sx q[1];
rz(-0.13012281) q[1];
sx q[1];
rz(-0.4904186) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0692762) q[3];
sx q[3];
rz(-2.2304689) q[3];
sx q[3];
rz(-0.31859627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10799321) q[2];
sx q[2];
rz(-0.42069837) q[2];
sx q[2];
rz(-1.9795974) q[2];
rz(1.1437931) q[3];
sx q[3];
rz(-1.9547209) q[3];
sx q[3];
rz(2.7403045) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6859632) q[0];
sx q[0];
rz(-0.81448737) q[0];
sx q[0];
rz(-2.234835) q[0];
rz(-2.7470159) q[1];
sx q[1];
rz(-0.60293174) q[1];
sx q[1];
rz(-1.1366064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0011463116) q[0];
sx q[0];
rz(-1.3147722) q[0];
sx q[0];
rz(1.2885068) q[0];
x q[1];
rz(-1.8086484) q[2];
sx q[2];
rz(-2.8360711) q[2];
sx q[2];
rz(-0.50070465) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26173466) q[1];
sx q[1];
rz(-2.2764808) q[1];
sx q[1];
rz(2.9565587) q[1];
x q[2];
rz(-1.2796938) q[3];
sx q[3];
rz(-1.2754692) q[3];
sx q[3];
rz(0.81084033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0236728) q[2];
sx q[2];
rz(-1.6912141) q[2];
sx q[2];
rz(-1.0517906) q[2];
rz(-3.1045095) q[3];
sx q[3];
rz(-0.6784234) q[3];
sx q[3];
rz(1.3858494) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7294075) q[0];
sx q[0];
rz(-2.2426028) q[0];
sx q[0];
rz(-2.7329024) q[0];
rz(0.15569923) q[1];
sx q[1];
rz(-1.0380048) q[1];
sx q[1];
rz(0.17539594) q[1];
rz(2.7338244) q[2];
sx q[2];
rz(-2.2121439) q[2];
sx q[2];
rz(-0.045814093) q[2];
rz(-1.2721522) q[3];
sx q[3];
rz(-2.1437672) q[3];
sx q[3];
rz(-1.7199316) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
