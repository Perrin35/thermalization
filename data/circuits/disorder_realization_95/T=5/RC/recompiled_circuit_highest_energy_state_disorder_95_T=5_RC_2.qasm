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
rz(1.3746102) q[0];
sx q[0];
rz(4.1144028) q[0];
sx q[0];
rz(10.65539) q[0];
rz(1.5988916) q[1];
sx q[1];
rz(3.4945421) q[1];
sx q[1];
rz(8.9897692) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1327788) q[0];
sx q[0];
rz(-2.7933196) q[0];
sx q[0];
rz(-1.2754557) q[0];
x q[1];
rz(2.6224766) q[2];
sx q[2];
rz(-0.26072219) q[2];
sx q[2];
rz(2.1671425) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5953345) q[1];
sx q[1];
rz(-1.2915242) q[1];
sx q[1];
rz(2.3550849) q[1];
x q[2];
rz(-0.013641274) q[3];
sx q[3];
rz(-1.6023926) q[3];
sx q[3];
rz(2.7073366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9445442) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(-0.91828263) q[2];
rz(2.7783172) q[3];
sx q[3];
rz(-1.2493635) q[3];
sx q[3];
rz(1.3242807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0488224) q[0];
sx q[0];
rz(-0.73960441) q[0];
sx q[0];
rz(-1.1356461) q[0];
rz(-2.8593235) q[1];
sx q[1];
rz(-1.5156563) q[1];
sx q[1];
rz(-2.938882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7736063) q[0];
sx q[0];
rz(-1.3895036) q[0];
sx q[0];
rz(1.3884991) q[0];
rz(1.5957256) q[2];
sx q[2];
rz(-2.6597616) q[2];
sx q[2];
rz(-0.78311759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2308826) q[1];
sx q[1];
rz(-2.6697914) q[1];
sx q[1];
rz(2.8698026) q[1];
x q[2];
rz(-1.1003157) q[3];
sx q[3];
rz(-0.011547877) q[3];
sx q[3];
rz(-1.5401767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2477766) q[2];
sx q[2];
rz(-0.42422366) q[2];
sx q[2];
rz(0.068537863) q[2];
rz(-2.2872772) q[3];
sx q[3];
rz(-1.2645384) q[3];
sx q[3];
rz(-3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0900367) q[0];
sx q[0];
rz(-0.55110252) q[0];
sx q[0];
rz(2.0447482) q[0];
rz(0.54496533) q[1];
sx q[1];
rz(-2.4577591) q[1];
sx q[1];
rz(-2.9473238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.504963) q[0];
sx q[0];
rz(-1.2218814) q[0];
sx q[0];
rz(-1.1303085) q[0];
x q[1];
rz(-1.0610007) q[2];
sx q[2];
rz(-0.61697996) q[2];
sx q[2];
rz(1.4088907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.0021469963) q[1];
sx q[1];
rz(-1.0771828) q[1];
sx q[1];
rz(1.1612372) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2195593) q[3];
sx q[3];
rz(-0.3742758) q[3];
sx q[3];
rz(0.55340278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7668309) q[2];
sx q[2];
rz(-1.3730048) q[2];
sx q[2];
rz(2.4135446) q[2];
rz(-2.9211365) q[3];
sx q[3];
rz(-0.15153344) q[3];
sx q[3];
rz(0.63173405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24882889) q[0];
sx q[0];
rz(-1.4735824) q[0];
sx q[0];
rz(-1.0695176) q[0];
rz(0.88691521) q[1];
sx q[1];
rz(-0.506217) q[1];
sx q[1];
rz(1.8198397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3794083) q[0];
sx q[0];
rz(-0.43566049) q[0];
sx q[0];
rz(0.092677516) q[0];
x q[1];
rz(1.555205) q[2];
sx q[2];
rz(-1.331845) q[2];
sx q[2];
rz(-0.76214253) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48044887) q[1];
sx q[1];
rz(-0.44434198) q[1];
sx q[1];
rz(1.4037057) q[1];
x q[2];
rz(-2.6245313) q[3];
sx q[3];
rz(-1.0365465) q[3];
sx q[3];
rz(-1.2599961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.47361031) q[2];
sx q[2];
rz(-2.0986291) q[2];
sx q[2];
rz(1.854151) q[2];
rz(2.1227664) q[3];
sx q[3];
rz(-0.5580709) q[3];
sx q[3];
rz(-2.4655925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.92893112) q[0];
sx q[0];
rz(-0.23290817) q[0];
sx q[0];
rz(-0.88053298) q[0];
rz(-3.0048043) q[1];
sx q[1];
rz(-2.9158178) q[1];
sx q[1];
rz(-2.5313964) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5870415) q[0];
sx q[0];
rz(-0.64233883) q[0];
sx q[0];
rz(-1.5498398) q[0];
x q[1];
rz(2.7099508) q[2];
sx q[2];
rz(-1.7100309) q[2];
sx q[2];
rz(2.8406258) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5194441) q[1];
sx q[1];
rz(-0.3560027) q[1];
sx q[1];
rz(-2.555833) q[1];
x q[2];
rz(-0.43936748) q[3];
sx q[3];
rz(-1.6939079) q[3];
sx q[3];
rz(1.1151552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5337164) q[2];
sx q[2];
rz(-2.4166962) q[2];
sx q[2];
rz(1.2302715) q[2];
rz(0.88165927) q[3];
sx q[3];
rz(-1.449838) q[3];
sx q[3];
rz(1.8264044) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45516685) q[0];
sx q[0];
rz(-1.5436341) q[0];
sx q[0];
rz(-2.0879188) q[0];
rz(-1.3764489) q[1];
sx q[1];
rz(-2.0339298) q[1];
sx q[1];
rz(2.460316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4453955) q[0];
sx q[0];
rz(-2.4413709) q[0];
sx q[0];
rz(-0.025431319) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37364423) q[2];
sx q[2];
rz(-2.3609997) q[2];
sx q[2];
rz(-0.045750387) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5743915) q[1];
sx q[1];
rz(-1.8387477) q[1];
sx q[1];
rz(1.059157) q[1];
x q[2];
rz(-1.525649) q[3];
sx q[3];
rz(-2.4490926) q[3];
sx q[3];
rz(-0.11603234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.072558746) q[2];
sx q[2];
rz(-0.37798887) q[2];
sx q[2];
rz(0.47522137) q[2];
rz(-1.2230988) q[3];
sx q[3];
rz(-2.1971072) q[3];
sx q[3];
rz(-0.16330115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.908919) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(-0.21937823) q[0];
rz(2.1351922) q[1];
sx q[1];
rz(-1.4234797) q[1];
sx q[1];
rz(-1.2054319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74586263) q[0];
sx q[0];
rz(-0.59915483) q[0];
sx q[0];
rz(1.5320734) q[0];
rz(-pi) q[1];
rz(2.5616165) q[2];
sx q[2];
rz(-2.6784228) q[2];
sx q[2];
rz(-1.620174) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0675065) q[1];
sx q[1];
rz(-1.5879022) q[1];
sx q[1];
rz(2.9422804) q[1];
rz(-pi) q[2];
rz(-3.0053776) q[3];
sx q[3];
rz(-2.0117617) q[3];
sx q[3];
rz(-1.7441526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.22126234) q[2];
sx q[2];
rz(-0.39322501) q[2];
sx q[2];
rz(0.66366759) q[2];
rz(-0.71781939) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(0.17620152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32888907) q[0];
sx q[0];
rz(-2.1825574) q[0];
sx q[0];
rz(1.9004199) q[0];
rz(-2.7136956) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(1.7180299) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9374838) q[0];
sx q[0];
rz(-0.23074575) q[0];
sx q[0];
rz(-2.8088914) q[0];
rz(-pi) q[1];
rz(2.2888118) q[2];
sx q[2];
rz(-1.6907534) q[2];
sx q[2];
rz(2.6991685) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9291722) q[1];
sx q[1];
rz(-0.66832405) q[1];
sx q[1];
rz(1.6051488) q[1];
rz(0.33496952) q[3];
sx q[3];
rz(-1.6050361) q[3];
sx q[3];
rz(-2.2881495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9972035) q[2];
sx q[2];
rz(-0.68789566) q[2];
sx q[2];
rz(-0.29590657) q[2];
rz(1.0812673) q[3];
sx q[3];
rz(-1.6176977) q[3];
sx q[3];
rz(-0.14779873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57701552) q[0];
sx q[0];
rz(-0.58203375) q[0];
sx q[0];
rz(-0.70593315) q[0];
rz(-2.9544746) q[1];
sx q[1];
rz(-2.3046604) q[1];
sx q[1];
rz(2.6954938) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69639403) q[0];
sx q[0];
rz(-0.15226752) q[0];
sx q[0];
rz(-3.1294786) q[0];
rz(0.31766625) q[2];
sx q[2];
rz(-1.1199754) q[2];
sx q[2];
rz(2.7433155) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6793788) q[1];
sx q[1];
rz(-2.129038) q[1];
sx q[1];
rz(2.7427196) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8817552) q[3];
sx q[3];
rz(-0.59229718) q[3];
sx q[3];
rz(-0.01334503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17911653) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(1.4013438) q[2];
rz(0.60351795) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(0.13874273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1097581) q[0];
sx q[0];
rz(-1.7367481) q[0];
sx q[0];
rz(-3.0314714) q[0];
rz(1.8203863) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(0.22470156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53407828) q[0];
sx q[0];
rz(-0.67994962) q[0];
sx q[0];
rz(2.2282766) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0761833) q[2];
sx q[2];
rz(-0.32534625) q[2];
sx q[2];
rz(-2.4869652) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4555295) q[1];
sx q[1];
rz(-1.7357329) q[1];
sx q[1];
rz(-2.4071724) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7431926) q[3];
sx q[3];
rz(-1.5972923) q[3];
sx q[3];
rz(-2.4717028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8297742) q[2];
sx q[2];
rz(-2.8303787) q[2];
sx q[2];
rz(-0.52660006) q[2];
rz(2.7569568) q[3];
sx q[3];
rz(-1.8943818) q[3];
sx q[3];
rz(-2.9068936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(0.93152355) q[0];
sx q[0];
rz(-0.43545224) q[0];
sx q[0];
rz(-0.42238105) q[0];
rz(0.79413636) q[1];
sx q[1];
rz(-1.4310478) q[1];
sx q[1];
rz(-1.4337883) q[1];
rz(0.88243816) q[2];
sx q[2];
rz(-1.0461764) q[2];
sx q[2];
rz(-1.7027693) q[2];
rz(-1.9392813) q[3];
sx q[3];
rz(-1.9574584) q[3];
sx q[3];
rz(-3.0493469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
