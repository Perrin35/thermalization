OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4053722) q[0];
sx q[0];
rz(-0.020981941) q[0];
sx q[0];
rz(-1.1809281) q[0];
rz(2.6612072) q[1];
sx q[1];
rz(-1.1552224) q[1];
sx q[1];
rz(0.46673271) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3213171) q[0];
sx q[0];
rz(-1.2877687) q[0];
sx q[0];
rz(0.77518605) q[0];
rz(-pi) q[1];
rz(-2.0601022) q[2];
sx q[2];
rz(-2.0976411) q[2];
sx q[2];
rz(-0.88050264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3148823) q[1];
sx q[1];
rz(-2.9032384) q[1];
sx q[1];
rz(1.5520067) q[1];
rz(-pi) q[2];
rz(0.67812293) q[3];
sx q[3];
rz(-2.0270542) q[3];
sx q[3];
rz(2.9564136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0537609) q[2];
sx q[2];
rz(-1.6873282) q[2];
sx q[2];
rz(1.9161179) q[2];
rz(0.37442225) q[3];
sx q[3];
rz(-2.008308) q[3];
sx q[3];
rz(-2.5882914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9541009) q[0];
sx q[0];
rz(-0.88045374) q[0];
sx q[0];
rz(-0.53400293) q[0];
rz(-2.7502637) q[1];
sx q[1];
rz(-0.44328872) q[1];
sx q[1];
rz(0.88723976) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0011883) q[0];
sx q[0];
rz(-2.6713704) q[0];
sx q[0];
rz(0.73064877) q[0];
rz(-pi) q[1];
rz(0.75669719) q[2];
sx q[2];
rz(-1.9644418) q[2];
sx q[2];
rz(1.4243038) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6897517) q[1];
sx q[1];
rz(-1.6997758) q[1];
sx q[1];
rz(-0.20492985) q[1];
rz(-1.0126566) q[3];
sx q[3];
rz(-2.5502656) q[3];
sx q[3];
rz(-1.8048353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5947764) q[2];
sx q[2];
rz(-1.3204601) q[2];
sx q[2];
rz(2.7776264) q[2];
rz(-1.724203) q[3];
sx q[3];
rz(-0.28984362) q[3];
sx q[3];
rz(2.3072306) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.904838) q[0];
sx q[0];
rz(-0.098243864) q[0];
sx q[0];
rz(-0.33016095) q[0];
rz(-2.5579021) q[1];
sx q[1];
rz(-0.81283641) q[1];
sx q[1];
rz(2.6469753) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1415633) q[0];
sx q[0];
rz(-0.38262472) q[0];
sx q[0];
rz(0.79866807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0614659) q[2];
sx q[2];
rz(-1.3441836) q[2];
sx q[2];
rz(-1.3148215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45060748) q[1];
sx q[1];
rz(-2.4276884) q[1];
sx q[1];
rz(2.3074478) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0298585) q[3];
sx q[3];
rz(-0.98971266) q[3];
sx q[3];
rz(1.5257208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3459449) q[2];
sx q[2];
rz(-0.66850418) q[2];
sx q[2];
rz(3.0188959) q[2];
rz(0.042938558) q[3];
sx q[3];
rz(-1.0791082) q[3];
sx q[3];
rz(-0.19449657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73404679) q[0];
sx q[0];
rz(-2.6341944) q[0];
sx q[0];
rz(-0.88974446) q[0];
rz(0.45097688) q[1];
sx q[1];
rz(-2.273592) q[1];
sx q[1];
rz(2.6873592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9376603) q[0];
sx q[0];
rz(-1.2302006) q[0];
sx q[0];
rz(1.457859) q[0];
x q[1];
rz(-1.5689911) q[2];
sx q[2];
rz(-1.1434525) q[2];
sx q[2];
rz(-1.3722562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52773038) q[1];
sx q[1];
rz(-2.3017028) q[1];
sx q[1];
rz(-1.1682603) q[1];
x q[2];
rz(2.9353981) q[3];
sx q[3];
rz(-2.7795994) q[3];
sx q[3];
rz(-2.9235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98895994) q[2];
sx q[2];
rz(-0.6344513) q[2];
sx q[2];
rz(0.75759849) q[2];
rz(0.18975137) q[3];
sx q[3];
rz(-1.8228143) q[3];
sx q[3];
rz(2.5943713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309288) q[0];
sx q[0];
rz(-0.78721109) q[0];
sx q[0];
rz(-1.2934562) q[0];
rz(2.5594607) q[1];
sx q[1];
rz(-1.4385834) q[1];
sx q[1];
rz(0.17939803) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6027733) q[0];
sx q[0];
rz(-1.7742549) q[0];
sx q[0];
rz(-2.957983) q[0];
x q[1];
rz(1.3264234) q[2];
sx q[2];
rz(-2.2269911) q[2];
sx q[2];
rz(2.6442621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.363626) q[1];
sx q[1];
rz(-0.95172626) q[1];
sx q[1];
rz(2.9113063) q[1];
x q[2];
rz(-0.87406335) q[3];
sx q[3];
rz(-1.1723601) q[3];
sx q[3];
rz(1.5521727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77357117) q[2];
sx q[2];
rz(-0.37006912) q[2];
sx q[2];
rz(-0.5698815) q[2];
rz(-0.44024769) q[3];
sx q[3];
rz(-1.2443685) q[3];
sx q[3];
rz(-2.7888035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.2920947) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(2.0000892) q[0];
rz(-1.3453206) q[1];
sx q[1];
rz(-0.99015403) q[1];
sx q[1];
rz(-2.9703531) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9123403) q[0];
sx q[0];
rz(-2.0496589) q[0];
sx q[0];
rz(1.4570057) q[0];
rz(2.0948579) q[2];
sx q[2];
rz(-1.9643133) q[2];
sx q[2];
rz(1.6167906) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0192193) q[1];
sx q[1];
rz(-1.4046486) q[1];
sx q[1];
rz(2.5801587) q[1];
rz(-pi) q[2];
rz(-2.1275131) q[3];
sx q[3];
rz(-1.5814335) q[3];
sx q[3];
rz(2.6680846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77219999) q[2];
sx q[2];
rz(-2.2056396) q[2];
sx q[2];
rz(3.0041223) q[2];
rz(1.8933659) q[3];
sx q[3];
rz(-0.56709254) q[3];
sx q[3];
rz(1.9198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0565979) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(-1.4884663) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.4005125) q[1];
sx q[1];
rz(-1.3135501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1689455) q[0];
sx q[0];
rz(-2.9914476) q[0];
sx q[0];
rz(-2.1779968) q[0];
x q[1];
rz(-2.9503472) q[2];
sx q[2];
rz(-1.4594263) q[2];
sx q[2];
rz(1.065606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.32282) q[1];
sx q[1];
rz(-1.5427151) q[1];
sx q[1];
rz(-0.089046003) q[1];
rz(-1.5093027) q[3];
sx q[3];
rz(-0.69685059) q[3];
sx q[3];
rz(-1.9887672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58275783) q[2];
sx q[2];
rz(-2.1119327) q[2];
sx q[2];
rz(0.47490698) q[2];
rz(0.20259914) q[3];
sx q[3];
rz(-2.30862) q[3];
sx q[3];
rz(-1.9995662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982034) q[0];
sx q[0];
rz(-0.13309637) q[0];
sx q[0];
rz(-0.30580795) q[0];
rz(0.43931475) q[1];
sx q[1];
rz(-2.3190934) q[1];
sx q[1];
rz(1.3154715) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1179349) q[0];
sx q[0];
rz(-2.2542037) q[0];
sx q[0];
rz(1.5853496) q[0];
x q[1];
rz(-1.015918) q[2];
sx q[2];
rz(-1.3279997) q[2];
sx q[2];
rz(1.0958874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0076239) q[1];
sx q[1];
rz(-2.3465912) q[1];
sx q[1];
rz(-2.9717173) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7301808) q[3];
sx q[3];
rz(-1.0209382) q[3];
sx q[3];
rz(-0.362277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7266015) q[2];
sx q[2];
rz(-1.5982268) q[2];
sx q[2];
rz(-0.42465633) q[2];
rz(-0.031938227) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(0.64938515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6040819) q[0];
sx q[0];
rz(-2.4896955) q[0];
sx q[0];
rz(-0.55484581) q[0];
rz(-2.8765053) q[1];
sx q[1];
rz(-1.6731429) q[1];
sx q[1];
rz(-1.9404985) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62702589) q[0];
sx q[0];
rz(-1.1794588) q[0];
sx q[0];
rz(-2.5364801) q[0];
rz(-pi) q[1];
rz(-1.3987819) q[2];
sx q[2];
rz(-1.2988678) q[2];
sx q[2];
rz(0.045814117) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7857901) q[1];
sx q[1];
rz(-1.5340163) q[1];
sx q[1];
rz(2.793502) q[1];
rz(0.74984925) q[3];
sx q[3];
rz(-1.5886652) q[3];
sx q[3];
rz(2.8211879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25192866) q[2];
sx q[2];
rz(-0.98968518) q[2];
sx q[2];
rz(-2.2970693) q[2];
rz(1.1318413) q[3];
sx q[3];
rz(-2.2364538) q[3];
sx q[3];
rz(1.3593146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.1691386) q[0];
sx q[0];
rz(-2.7511981) q[0];
sx q[0];
rz(1.1727232) q[0];
rz(-0.68924618) q[1];
sx q[1];
rz(-1.0968364) q[1];
sx q[1];
rz(0.26395878) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57676753) q[0];
sx q[0];
rz(-1.1297313) q[0];
sx q[0];
rz(0.70617981) q[0];
x q[1];
rz(0.9969292) q[2];
sx q[2];
rz(-0.80748122) q[2];
sx q[2];
rz(-1.5308183) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2745299) q[1];
sx q[1];
rz(-1.2559051) q[1];
sx q[1];
rz(-0.19537433) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92688074) q[3];
sx q[3];
rz(-1.6686693) q[3];
sx q[3];
rz(-2.904908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3802203) q[2];
sx q[2];
rz(-1.994588) q[2];
sx q[2];
rz(1.1205565) q[2];
rz(1.5348966) q[3];
sx q[3];
rz(-0.94529072) q[3];
sx q[3];
rz(1.5103316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071038889) q[0];
sx q[0];
rz(-2.4865535) q[0];
sx q[0];
rz(3.0336663) q[0];
rz(-2.4212266) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(0.82577827) q[2];
sx q[2];
rz(-2.0407531) q[2];
sx q[2];
rz(-2.4498418) q[2];
rz(-1.646044) q[3];
sx q[3];
rz(-0.9552707) q[3];
sx q[3];
rz(0.09494119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
