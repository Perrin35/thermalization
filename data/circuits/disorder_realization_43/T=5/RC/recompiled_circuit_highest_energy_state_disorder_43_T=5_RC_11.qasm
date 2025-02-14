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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(0.55846941) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(1.37473) q[1];
sx q[1];
rz(9.3461499) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0538875) q[0];
sx q[0];
rz(-1.0358521) q[0];
sx q[0];
rz(2.22552) q[0];
x q[1];
rz(1.4457277) q[2];
sx q[2];
rz(-0.89587051) q[2];
sx q[2];
rz(-1.9927911) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.47657644) q[1];
sx q[1];
rz(-1.4890097) q[1];
sx q[1];
rz(2.8881914) q[1];
rz(-1.2492248) q[3];
sx q[3];
rz(-1.9059637) q[3];
sx q[3];
rz(-0.17091076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89995304) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(-2.107035) q[2];
rz(-0.11040802) q[3];
sx q[3];
rz(-1.4805099) q[3];
sx q[3];
rz(-0.30606562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0050874) q[0];
sx q[0];
rz(-1.445048) q[0];
sx q[0];
rz(-1.9553631) q[0];
rz(-1.8738481) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(-1.3892106) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88436172) q[0];
sx q[0];
rz(-0.62998929) q[0];
sx q[0];
rz(-2.1068186) q[0];
rz(2.9302338) q[2];
sx q[2];
rz(-0.49079259) q[2];
sx q[2];
rz(2.9826814) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8858365) q[1];
sx q[1];
rz(-2.6955018) q[1];
sx q[1];
rz(-2.2913833) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0458353) q[3];
sx q[3];
rz(-2.1928722) q[3];
sx q[3];
rz(-0.84703895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.17239751) q[2];
sx q[2];
rz(-2.0266504) q[2];
sx q[2];
rz(1.4364852) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-0.45537046) q[3];
sx q[3];
rz(0.026738515) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7773975) q[0];
sx q[0];
rz(-1.3854249) q[0];
sx q[0];
rz(-1.2170323) q[0];
rz(0.2150391) q[1];
sx q[1];
rz(-1.8937078) q[1];
sx q[1];
rz(1.4395813) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7584981) q[0];
sx q[0];
rz(-1.4466084) q[0];
sx q[0];
rz(0.66207768) q[0];
rz(-pi) q[1];
x q[1];
rz(0.093482253) q[2];
sx q[2];
rz(-1.0614402) q[2];
sx q[2];
rz(2.2619154) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90104687) q[1];
sx q[1];
rz(-0.26296294) q[1];
sx q[1];
rz(2.0790624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2974627) q[3];
sx q[3];
rz(-1.36051) q[3];
sx q[3];
rz(0.027829011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8000468) q[2];
sx q[2];
rz(-1.6365106) q[2];
sx q[2];
rz(-2.3254507) q[2];
rz(0.0084361313) q[3];
sx q[3];
rz(-0.71788994) q[3];
sx q[3];
rz(1.5754023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208045) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(2.8724443) q[0];
rz(-1.3828269) q[1];
sx q[1];
rz(-2.5803284) q[1];
sx q[1];
rz(0.47163481) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29971805) q[0];
sx q[0];
rz(-0.73900676) q[0];
sx q[0];
rz(1.3625381) q[0];
rz(0.083375562) q[2];
sx q[2];
rz(-0.40503392) q[2];
sx q[2];
rz(2.0138559) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2216404) q[1];
sx q[1];
rz(-2.5121009) q[1];
sx q[1];
rz(0.55947495) q[1];
x q[2];
rz(-0.35393039) q[3];
sx q[3];
rz(-2.2056863) q[3];
sx q[3];
rz(0.77087444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4857594) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(-1.887623) q[2];
rz(-0.29821011) q[3];
sx q[3];
rz(-1.8658274) q[3];
sx q[3];
rz(1.5378753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1042079) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(-2.3846159) q[0];
rz(-1.0589927) q[1];
sx q[1];
rz(-1.4451566) q[1];
sx q[1];
rz(-0.34245488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.379102) q[0];
sx q[0];
rz(-0.44371334) q[0];
sx q[0];
rz(1.1136454) q[0];
rz(2.8169698) q[2];
sx q[2];
rz(-0.89541905) q[2];
sx q[2];
rz(-1.7243232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7389512) q[1];
sx q[1];
rz(-0.61719751) q[1];
sx q[1];
rz(2.4715241) q[1];
rz(-pi) q[2];
rz(0.60635546) q[3];
sx q[3];
rz(-2.6487051) q[3];
sx q[3];
rz(2.3541401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.69514889) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(0.62015074) q[2];
rz(-2.5946963) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3243489) q[0];
sx q[0];
rz(-0.20340782) q[0];
sx q[0];
rz(1.2212344) q[0];
rz(-2.0381894) q[1];
sx q[1];
rz(-1.1627448) q[1];
sx q[1];
rz(-0.81072909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2317144) q[0];
sx q[0];
rz(-2.8389769) q[0];
sx q[0];
rz(0.55993863) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76580183) q[2];
sx q[2];
rz(-2.3277936) q[2];
sx q[2];
rz(1.9794996) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.4524881) q[1];
sx q[1];
rz(-2.9924722) q[1];
sx q[1];
rz(-1.945687) q[1];
rz(1.9493413) q[3];
sx q[3];
rz(-1.2929799) q[3];
sx q[3];
rz(-1.6343821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.80060426) q[2];
sx q[2];
rz(-1.1977414) q[2];
sx q[2];
rz(-1.9761337) q[2];
rz(-0.10945877) q[3];
sx q[3];
rz(-2.1200924) q[3];
sx q[3];
rz(-3.0808595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8282181) q[0];
sx q[0];
rz(-3.1309541) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(0.24318801) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(2.5004255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0571447) q[0];
sx q[0];
rz(-2.085745) q[0];
sx q[0];
rz(-0.65339974) q[0];
x q[1];
rz(-2.7566822) q[2];
sx q[2];
rz(-0.78854568) q[2];
sx q[2];
rz(0.53764082) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9037364) q[1];
sx q[1];
rz(-2.3163539) q[1];
sx q[1];
rz(-2.0261057) q[1];
rz(1.0585467) q[3];
sx q[3];
rz(-0.44971684) q[3];
sx q[3];
rz(-1.095677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0082561) q[2];
sx q[2];
rz(-1.1587605) q[2];
sx q[2];
rz(-2.1978417) q[2];
rz(-2.4933695) q[3];
sx q[3];
rz(-0.74130487) q[3];
sx q[3];
rz(-0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42565313) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-0.42914036) q[0];
rz(0.40402135) q[1];
sx q[1];
rz(-2.7262913) q[1];
sx q[1];
rz(-0.9187575) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63916535) q[0];
sx q[0];
rz(-0.60898655) q[0];
sx q[0];
rz(-2.7091685) q[0];
x q[1];
rz(-1.8317675) q[2];
sx q[2];
rz(-1.4497641) q[2];
sx q[2];
rz(2.3132035) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1707216) q[1];
sx q[1];
rz(-1.7890463) q[1];
sx q[1];
rz(-1.2531325) q[1];
rz(-pi) q[2];
x q[2];
rz(1.346896) q[3];
sx q[3];
rz(-2.8200275) q[3];
sx q[3];
rz(-1.079815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.018844133) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(1.3502236) q[2];
rz(-0.33251479) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(-1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74828446) q[0];
sx q[0];
rz(-0.64993334) q[0];
sx q[0];
rz(-2.6848324) q[0];
rz(-1.7204174) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(-3.0866887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47726705) q[0];
sx q[0];
rz(-2.3033963) q[0];
sx q[0];
rz(-3.102836) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3308067) q[2];
sx q[2];
rz(-0.56466757) q[2];
sx q[2];
rz(1.1155903) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0662288) q[1];
sx q[1];
rz(-1.7875331) q[1];
sx q[1];
rz(-2.0677057) q[1];
rz(0.53590448) q[3];
sx q[3];
rz(-1.2782989) q[3];
sx q[3];
rz(-1.3939987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88825893) q[2];
sx q[2];
rz(-0.71037018) q[2];
sx q[2];
rz(-0.19676512) q[2];
rz(1.0737859) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(2.4055068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9842904) q[0];
sx q[0];
rz(-2.2881916) q[0];
sx q[0];
rz(2.2421457) q[0];
rz(-2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(1.4003632) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.78681) q[0];
sx q[0];
rz(-1.7799095) q[0];
sx q[0];
rz(2.1378921) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1323523) q[2];
sx q[2];
rz(-2.6137527) q[2];
sx q[2];
rz(-0.49293016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33161641) q[1];
sx q[1];
rz(-0.73996937) q[1];
sx q[1];
rz(-0.59438594) q[1];
x q[2];
rz(-2.1651044) q[3];
sx q[3];
rz(-2.1955238) q[3];
sx q[3];
rz(-0.67052602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0805936) q[2];
sx q[2];
rz(-2.521535) q[2];
sx q[2];
rz(-1.9798123) q[2];
rz(-1.0787841) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(-2.259528) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921191) q[0];
sx q[0];
rz(-1.5220806) q[0];
sx q[0];
rz(1.4694389) q[0];
rz(3.0081765) q[1];
sx q[1];
rz(-1.3795556) q[1];
sx q[1];
rz(2.1605927) q[1];
rz(-1.6242003) q[2];
sx q[2];
rz(-1.5697865) q[2];
sx q[2];
rz(-3.1278232) q[2];
rz(2.6399163) q[3];
sx q[3];
rz(-1.6335084) q[3];
sx q[3];
rz(-0.15214534) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
