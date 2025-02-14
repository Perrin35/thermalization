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
rz(0.8402549) q[0];
sx q[0];
rz(4.1794887) q[0];
sx q[0];
rz(10.269796) q[0];
rz(0.71212274) q[1];
sx q[1];
rz(-2.1399838) q[1];
sx q[1];
rz(1.4955624) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7094249) q[0];
sx q[0];
rz(-2.3948095) q[0];
sx q[0];
rz(-0.92632697) q[0];
x q[1];
rz(1.7208112) q[2];
sx q[2];
rz(-1.9505462) q[2];
sx q[2];
rz(-0.94204547) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7993303) q[1];
sx q[1];
rz(-1.0940219) q[1];
sx q[1];
rz(-2.083198) q[1];
x q[2];
rz(1.846662) q[3];
sx q[3];
rz(-1.952269) q[3];
sx q[3];
rz(1.4700397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.738395) q[2];
sx q[2];
rz(-1.9591816) q[2];
sx q[2];
rz(-0.5564059) q[2];
rz(-0.49499908) q[3];
sx q[3];
rz(-2.7904816) q[3];
sx q[3];
rz(-0.88465148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(0.86968386) q[0];
sx q[0];
rz(-3.0392635) q[0];
sx q[0];
rz(-2.5129357) q[0];
rz(2.2643845) q[1];
sx q[1];
rz(-2.6429206) q[1];
sx q[1];
rz(-0.25310755) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019397639) q[0];
sx q[0];
rz(-3.1240648) q[0];
sx q[0];
rz(-2.0445834) q[0];
rz(-0.57665373) q[2];
sx q[2];
rz(-2.6284784) q[2];
sx q[2];
rz(-0.54158825) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69324609) q[1];
sx q[1];
rz(-2.9102059) q[1];
sx q[1];
rz(0.53129249) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83929245) q[3];
sx q[3];
rz(-1.7613693) q[3];
sx q[3];
rz(-2.1803602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6796598) q[2];
sx q[2];
rz(-1.9340645) q[2];
sx q[2];
rz(-2.0750462) q[2];
rz(1.3072394) q[3];
sx q[3];
rz(-2.6756838) q[3];
sx q[3];
rz(0.90827847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2147373) q[0];
sx q[0];
rz(-2.1934788) q[0];
sx q[0];
rz(-0.98534775) q[0];
rz(-0.39380479) q[1];
sx q[1];
rz(-2.1486798) q[1];
sx q[1];
rz(-0.52678144) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39296374) q[0];
sx q[0];
rz(-2.606578) q[0];
sx q[0];
rz(-1.4070542) q[0];
rz(-0.71957876) q[2];
sx q[2];
rz(-2.4256896) q[2];
sx q[2];
rz(0.86052513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5722428) q[1];
sx q[1];
rz(-1.7399638) q[1];
sx q[1];
rz(-1.1430686) q[1];
rz(-pi) q[2];
rz(2.1680896) q[3];
sx q[3];
rz(-1.9466578) q[3];
sx q[3];
rz(-1.5099389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7963316) q[2];
sx q[2];
rz(-0.76097208) q[2];
sx q[2];
rz(-0.14656466) q[2];
rz(-2.2740299) q[3];
sx q[3];
rz(-0.87351322) q[3];
sx q[3];
rz(-0.7578907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72916156) q[0];
sx q[0];
rz(-0.49598345) q[0];
sx q[0];
rz(0.55369401) q[0];
rz(-1.6665392) q[1];
sx q[1];
rz(-0.29860425) q[1];
sx q[1];
rz(-2.8972304) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2840773) q[0];
sx q[0];
rz(-1.2173442) q[0];
sx q[0];
rz(-1.2156603) q[0];
x q[1];
rz(1.9072794) q[2];
sx q[2];
rz(-0.76533106) q[2];
sx q[2];
rz(-0.31631563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3267388) q[1];
sx q[1];
rz(-2.5757901) q[1];
sx q[1];
rz(2.4292355) q[1];
rz(-3.1126106) q[3];
sx q[3];
rz(-0.71235114) q[3];
sx q[3];
rz(-0.94635743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.060140572) q[2];
sx q[2];
rz(-2.0144561) q[2];
sx q[2];
rz(0.78090182) q[2];
rz(2.7447356) q[3];
sx q[3];
rz(-0.01440993) q[3];
sx q[3];
rz(1.074033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9020554) q[0];
sx q[0];
rz(-2.9750415) q[0];
sx q[0];
rz(0.2567513) q[0];
rz(2.9638839) q[1];
sx q[1];
rz(-0.54568988) q[1];
sx q[1];
rz(0.94388747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092800822) q[0];
sx q[0];
rz(-1.9328797) q[0];
sx q[0];
rz(-2.0440153) q[0];
x q[1];
rz(1.2805528) q[2];
sx q[2];
rz(-1.7566881) q[2];
sx q[2];
rz(-3.0772532) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0066932) q[1];
sx q[1];
rz(-1.0311145) q[1];
sx q[1];
rz(1.5657019) q[1];
x q[2];
rz(-1.0141171) q[3];
sx q[3];
rz(-1.1667487) q[3];
sx q[3];
rz(1.7168644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2029767) q[2];
sx q[2];
rz(-2.1111033) q[2];
sx q[2];
rz(-2.9317648) q[2];
rz(0.046791568) q[3];
sx q[3];
rz(-1.4593461) q[3];
sx q[3];
rz(-0.34745026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059134722) q[0];
sx q[0];
rz(-2.3248398) q[0];
sx q[0];
rz(1.9385852) q[0];
rz(-1.5164392) q[1];
sx q[1];
rz(-0.37767437) q[1];
sx q[1];
rz(-3.1086521) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0982978) q[0];
sx q[0];
rz(-1.8445208) q[0];
sx q[0];
rz(1.0874463) q[0];
rz(-pi) q[1];
rz(0.85196544) q[2];
sx q[2];
rz(-2.4174066) q[2];
sx q[2];
rz(0.23978309) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6194544) q[1];
sx q[1];
rz(-0.65548766) q[1];
sx q[1];
rz(-2.5966899) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5046963) q[3];
sx q[3];
rz(-0.99352443) q[3];
sx q[3];
rz(-0.12337109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5861627) q[2];
sx q[2];
rz(-1.0032434) q[2];
sx q[2];
rz(0.45247751) q[2];
rz(-2.9943976) q[3];
sx q[3];
rz(-2.9087524) q[3];
sx q[3];
rz(1.2991306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.074987) q[0];
sx q[0];
rz(-2.7050278) q[0];
sx q[0];
rz(0.35705745) q[0];
rz(2.1287411) q[1];
sx q[1];
rz(-1.169299) q[1];
sx q[1];
rz(0.036651932) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22003156) q[0];
sx q[0];
rz(-2.4274094) q[0];
sx q[0];
rz(0.69252391) q[0];
x q[1];
rz(-2.9190262) q[2];
sx q[2];
rz(-0.43912008) q[2];
sx q[2];
rz(2.665208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4067024) q[1];
sx q[1];
rz(-1.6631366) q[1];
sx q[1];
rz(-0.86812302) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5753079) q[3];
sx q[3];
rz(-1.6016212) q[3];
sx q[3];
rz(0.84336057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6645633) q[2];
sx q[2];
rz(-2.1881115) q[2];
sx q[2];
rz(-0.089740962) q[2];
rz(-1.8803895) q[3];
sx q[3];
rz(-1.829105) q[3];
sx q[3];
rz(-1.7173654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66202128) q[0];
sx q[0];
rz(-0.57858545) q[0];
sx q[0];
rz(1.4805967) q[0];
rz(-1.4501694) q[1];
sx q[1];
rz(-2.4822576) q[1];
sx q[1];
rz(-2.3816542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63045365) q[0];
sx q[0];
rz(-1.5604815) q[0];
sx q[0];
rz(1.6869839) q[0];
x q[1];
rz(2.5133695) q[2];
sx q[2];
rz(-1.4951653) q[2];
sx q[2];
rz(3.0043907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7260598) q[1];
sx q[1];
rz(-2.5431255) q[1];
sx q[1];
rz(2.3319671) q[1];
x q[2];
rz(-1.8283056) q[3];
sx q[3];
rz(-2.0638568) q[3];
sx q[3];
rz(0.3230394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30457589) q[2];
sx q[2];
rz(-1.872007) q[2];
sx q[2];
rz(3.0820091) q[2];
rz(-0.42851055) q[3];
sx q[3];
rz(-0.84463745) q[3];
sx q[3];
rz(-2.5394411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8932327) q[0];
sx q[0];
rz(-0.28286523) q[0];
sx q[0];
rz(2.7981113) q[0];
rz(0.19613014) q[1];
sx q[1];
rz(-2.2562512) q[1];
sx q[1];
rz(-0.79998618) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72092743) q[0];
sx q[0];
rz(-1.1625097) q[0];
sx q[0];
rz(0.83939206) q[0];
rz(-2.7550789) q[2];
sx q[2];
rz(-2.5441351) q[2];
sx q[2];
rz(-2.6253485) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5243036) q[1];
sx q[1];
rz(-2.5802748) q[1];
sx q[1];
rz(-0.08331068) q[1];
rz(2.2535747) q[3];
sx q[3];
rz(-1.5170485) q[3];
sx q[3];
rz(-2.2537504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.34667748) q[2];
sx q[2];
rz(-2.3778043) q[2];
sx q[2];
rz(-2.9128892) q[2];
rz(1.4165357) q[3];
sx q[3];
rz(-1.4226457) q[3];
sx q[3];
rz(-0.67030877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59920853) q[0];
sx q[0];
rz(-2.6238361) q[0];
sx q[0];
rz(1.9589348) q[0];
rz(2.5987437) q[1];
sx q[1];
rz(-0.12195568) q[1];
sx q[1];
rz(-0.45546946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96129721) q[0];
sx q[0];
rz(-1.1396798) q[0];
sx q[0];
rz(0.53133734) q[0];
rz(-2.9331839) q[2];
sx q[2];
rz(-1.128607) q[2];
sx q[2];
rz(-0.65763523) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7389981) q[1];
sx q[1];
rz(-1.4345503) q[1];
sx q[1];
rz(2.6722277) q[1];
x q[2];
rz(0.72584589) q[3];
sx q[3];
rz(-2.5739087) q[3];
sx q[3];
rz(2.3633702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21696422) q[2];
sx q[2];
rz(-2.2483726) q[2];
sx q[2];
rz(0.42579892) q[2];
rz(-0.18481542) q[3];
sx q[3];
rz(-2.5033689) q[3];
sx q[3];
rz(3.1137915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4375147) q[0];
sx q[0];
rz(-1.6117493) q[0];
sx q[0];
rz(3.1060863) q[0];
rz(-2.6620445) q[1];
sx q[1];
rz(-1.5794812) q[1];
sx q[1];
rz(-1.5522122) q[1];
rz(0.59500792) q[2];
sx q[2];
rz(-1.7562661) q[2];
sx q[2];
rz(-1.4944639) q[2];
rz(0.95994031) q[3];
sx q[3];
rz(-1.452594) q[3];
sx q[3];
rz(0.95127524) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
