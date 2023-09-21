OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7802785) q[0];
sx q[0];
rz(-1.8457992) q[0];
sx q[0];
rz(-1.5048024) q[0];
rz(-1.7668031) q[2];
sx q[2];
rz(-0.73050806) q[2];
sx q[2];
rz(1.0175878) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3589188) q[1];
sx q[1];
rz(-2.3896857) q[1];
sx q[1];
rz(-2.9208675) q[1];
x q[2];
rz(-1.3189089) q[3];
sx q[3];
rz(-0.10676521) q[3];
sx q[3];
rz(1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(-2.1271558) q[2];
rz(2.9075918) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(0.28449374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(-2.1388334) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013872) q[0];
sx q[0];
rz(-1.1367072) q[0];
sx q[0];
rz(0.65625221) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1728671) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(-2.6181521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1042852) q[1];
sx q[1];
rz(-0.98854317) q[1];
sx q[1];
rz(-2.9557455) q[1];
rz(-pi) q[2];
rz(-2.6456804) q[3];
sx q[3];
rz(-2.4818588) q[3];
sx q[3];
rz(-1.0752614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(2.9906452) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(-3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2484444) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-0.77899581) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(0.88358203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8944091) q[0];
sx q[0];
rz(-2.7054225) q[0];
sx q[0];
rz(-2.7437074) q[0];
rz(-0.94647879) q[2];
sx q[2];
rz(-2.6636071) q[2];
sx q[2];
rz(-1.9566655) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5501432) q[1];
sx q[1];
rz(-1.5933697) q[1];
sx q[1];
rz(-1.8557465) q[1];
x q[2];
rz(-1.6885353) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8399923) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-2.278573) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(2.5298932) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94565369) q[0];
sx q[0];
rz(-1.6473856) q[0];
sx q[0];
rz(2.4038195) q[0];
rz(-pi) q[1];
rz(2.8407211) q[2];
sx q[2];
rz(-2.8340325) q[2];
sx q[2];
rz(0.077066271) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78754567) q[1];
sx q[1];
rz(-2.5552632) q[1];
sx q[1];
rz(1.3613308) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9228966) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(0.72363879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(2.8660529) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-0.87669796) q[0];
rz(-2.450401) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-0.91526389) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4782151) q[0];
sx q[0];
rz(-1.6179248) q[0];
sx q[0];
rz(2.1992654) q[0];
rz(-pi) q[1];
x q[1];
rz(0.077276262) q[2];
sx q[2];
rz(-0.52646598) q[2];
sx q[2];
rz(-2.4174945) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1716869) q[1];
sx q[1];
rz(-1.415442) q[1];
sx q[1];
rz(0.032226493) q[1];
rz(-2.4004585) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22121945) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(-2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(1.0790496) q[0];
rz(0.59533978) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(-2.7450096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3298033) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(1.1008218) q[0];
rz(0.47607143) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(1.6518041) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.035316) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(-1.4186526) q[1];
x q[2];
rz(1.8144238) q[3];
sx q[3];
rz(-1.0001839) q[3];
sx q[3];
rz(-0.43587303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(-0.78480762) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9046017) q[0];
sx q[0];
rz(-0.31123529) q[0];
sx q[0];
rz(-2.8471332) q[0];
x q[1];
rz(2.9992505) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(3.1359283) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3583402) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.7794442) q[1];
x q[2];
rz(1.3133615) q[3];
sx q[3];
rz(-2.7830887) q[3];
sx q[3];
rz(0.076971019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37041935) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(-2.1408634) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(-3.1088366) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(-0.89950048) q[0];
rz(2.2135127) q[2];
sx q[2];
rz(-2.4734481) q[2];
sx q[2];
rz(1.1754787) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9084839) q[1];
sx q[1];
rz(-1.9875803) q[1];
sx q[1];
rz(-0.064518708) q[1];
x q[2];
rz(-0.0891536) q[3];
sx q[3];
rz(-1.9490064) q[3];
sx q[3];
rz(-2.0034727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(3.0416097) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.6149678) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582726) q[0];
sx q[0];
rz(-2.9968046) q[0];
sx q[0];
rz(-0.38082122) q[0];
rz(1.7301241) q[2];
sx q[2];
rz(-2.1630641) q[2];
sx q[2];
rz(-1.9221905) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.607) q[1];
sx q[1];
rz(-1.8796762) q[1];
sx q[1];
rz(-0.81307462) q[1];
x q[2];
rz(-1.1205348) q[3];
sx q[3];
rz(-2.1053208) q[3];
sx q[3];
rz(3.0578606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0316524) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(-3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.7074701) q[0];
rz(-1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(-0.5982582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0043250672) q[0];
sx q[0];
rz(-0.30170545) q[0];
sx q[0];
rz(-2.3803677) q[0];
rz(-0.21226378) q[2];
sx q[2];
rz(-1.4258254) q[2];
sx q[2];
rz(-0.14130172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88565247) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(3.0401858) q[1];
rz(-pi) q[2];
rz(-1.7353021) q[3];
sx q[3];
rz(-1.7269772) q[3];
sx q[3];
rz(2.3567049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(0.59371289) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(-2.0098856) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3020637) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.2364173) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-0.61959735) q[2];
sx q[2];
rz(-1.4705114) q[2];
sx q[2];
rz(0.32497139) q[2];
rz(0.89647722) q[3];
sx q[3];
rz(-1.8835246) q[3];
sx q[3];
rz(-2.6261331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
