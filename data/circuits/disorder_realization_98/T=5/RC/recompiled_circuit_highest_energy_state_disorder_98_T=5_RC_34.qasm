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
rz(-0.67796081) q[0];
sx q[0];
rz(6.0089587) q[0];
sx q[0];
rz(11.275509) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8052657) q[0];
sx q[0];
rz(-1.8417497) q[0];
sx q[0];
rz(1.8150369) q[0];
rz(-pi) q[1];
rz(1.8617282) q[2];
sx q[2];
rz(-0.81255823) q[2];
sx q[2];
rz(-1.1464553) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95924458) q[1];
sx q[1];
rz(-2.3492491) q[1];
sx q[1];
rz(0.46087973) q[1];
x q[2];
rz(2.6628482) q[3];
sx q[3];
rz(-2.1447721) q[3];
sx q[3];
rz(0.51306242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(-0.901326) q[2];
rz(0.46402913) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(-3.071781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0483911) q[0];
sx q[0];
rz(-3.1049325) q[0];
sx q[0];
rz(-1.8638336) q[0];
rz(3.0473862) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(-1.6069848) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8746895) q[0];
sx q[0];
rz(-2.1703566) q[0];
sx q[0];
rz(-0.039681704) q[0];
x q[1];
rz(-0.056939967) q[2];
sx q[2];
rz(-1.3888265) q[2];
sx q[2];
rz(1.6873311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6834641) q[1];
sx q[1];
rz(-0.56741112) q[1];
sx q[1];
rz(2.6415884) q[1];
rz(-0.40462287) q[3];
sx q[3];
rz(-1.957486) q[3];
sx q[3];
rz(2.2819732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4216807) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-0.16033944) q[2];
rz(-2.4028589) q[3];
sx q[3];
rz(-2.2952047) q[3];
sx q[3];
rz(1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
rz(-0.62683231) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(2.6318188) q[0];
rz(1.7104507) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(0.74657718) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7187481) q[0];
sx q[0];
rz(-1.5908518) q[0];
sx q[0];
rz(1.0182747) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66303894) q[2];
sx q[2];
rz(-2.2929077) q[2];
sx q[2];
rz(2.7977365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.67184) q[1];
sx q[1];
rz(-0.79974175) q[1];
sx q[1];
rz(2.0861162) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88991092) q[3];
sx q[3];
rz(-1.8783675) q[3];
sx q[3];
rz(-0.5336844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9868682) q[2];
sx q[2];
rz(-2.8747323) q[2];
sx q[2];
rz(1.9179087) q[2];
rz(-1.0736505) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735485) q[0];
sx q[0];
rz(-2.2548964) q[0];
sx q[0];
rz(2.7622188) q[0];
rz(0.1768449) q[1];
sx q[1];
rz(-1.481448) q[1];
sx q[1];
rz(0.79536974) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75126266) q[0];
sx q[0];
rz(-0.38026938) q[0];
sx q[0];
rz(-1.8230536) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8537117) q[2];
sx q[2];
rz(-2.3672315) q[2];
sx q[2];
rz(1.1224358) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25177342) q[1];
sx q[1];
rz(-1.3135513) q[1];
sx q[1];
rz(-2.1104269) q[1];
rz(-pi) q[2];
rz(0.42511149) q[3];
sx q[3];
rz(-3.0038341) q[3];
sx q[3];
rz(2.3099096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7103601) q[2];
sx q[2];
rz(-1.3453588) q[2];
sx q[2];
rz(1.7809407) q[2];
rz(-2.1121173) q[3];
sx q[3];
rz(-1.4808713) q[3];
sx q[3];
rz(1.3220538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(1.5486451) q[0];
rz(-0.40052888) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.4097479) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9841138) q[0];
sx q[0];
rz(-1.3712008) q[0];
sx q[0];
rz(-0.73961135) q[0];
rz(-0.39938853) q[2];
sx q[2];
rz(-0.3872954) q[2];
sx q[2];
rz(-3.0687817) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6741989) q[1];
sx q[1];
rz(-1.129303) q[1];
sx q[1];
rz(0.91798325) q[1];
x q[2];
rz(1.665776) q[3];
sx q[3];
rz(-0.62415571) q[3];
sx q[3];
rz(-2.399866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3132402) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(0.43133119) q[2];
rz(3.0865772) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(-0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677652) q[0];
sx q[0];
rz(-0.25280935) q[0];
sx q[0];
rz(1.025169) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-2.5621474) q[1];
sx q[1];
rz(2.2377009) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9867803) q[0];
sx q[0];
rz(-2.203699) q[0];
sx q[0];
rz(0.38582071) q[0];
rz(-pi) q[1];
rz(-2.2270062) q[2];
sx q[2];
rz(-1.0433955) q[2];
sx q[2];
rz(1.630115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55906193) q[1];
sx q[1];
rz(-0.85147038) q[1];
sx q[1];
rz(1.933631) q[1];
rz(-pi) q[2];
rz(-0.76670209) q[3];
sx q[3];
rz(-0.23582102) q[3];
sx q[3];
rz(-2.9028149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8338833) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-0.21635381) q[2];
rz(-0.89573914) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(-1.5007796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(-0.69881451) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(0.85174495) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9601933) q[0];
sx q[0];
rz(-1.6087247) q[0];
sx q[0];
rz(-1.0344124) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5965271) q[2];
sx q[2];
rz(-0.99359578) q[2];
sx q[2];
rz(2.4943697) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56543437) q[1];
sx q[1];
rz(-1.2093739) q[1];
sx q[1];
rz(0.71908497) q[1];
rz(-pi) q[2];
rz(2.0533356) q[3];
sx q[3];
rz(-0.98615328) q[3];
sx q[3];
rz(-2.6934153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.8375887) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(2.6336929) q[2];
rz(-2.1128283) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(0.60104162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69865882) q[0];
sx q[0];
rz(-1.826094) q[0];
sx q[0];
rz(-3.1319295) q[0];
rz(1.6161605) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(-2.756871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54668173) q[0];
sx q[0];
rz(-1.718504) q[0];
sx q[0];
rz(1.8099996) q[0];
rz(-pi) q[1];
rz(-2.9392249) q[2];
sx q[2];
rz(-1.9500537) q[2];
sx q[2];
rz(-2.8059562) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.87054756) q[1];
sx q[1];
rz(-1.3636075) q[1];
sx q[1];
rz(-1.8533587) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99920239) q[3];
sx q[3];
rz(-0.88547844) q[3];
sx q[3];
rz(-0.1911605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.559451) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(0.01586308) q[2];
rz(-1.2265685) q[3];
sx q[3];
rz(-2.3358986) q[3];
sx q[3];
rz(-0.62172833) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7568307) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(2.050198) q[0];
rz(0.81820828) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(-2.4726726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9364939) q[0];
sx q[0];
rz(-0.97081796) q[0];
sx q[0];
rz(2.4121802) q[0];
rz(-1.0878272) q[2];
sx q[2];
rz(-2.3685799) q[2];
sx q[2];
rz(0.20394606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.131625) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(1.4373006) q[1];
rz(-pi) q[2];
rz(-2.8112765) q[3];
sx q[3];
rz(-1.013375) q[3];
sx q[3];
rz(-1.5239609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5074671) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(-2.1779306) q[2];
rz(1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(-1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1249579) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(-1.6868663) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.5364372) q[1];
sx q[1];
rz(0.32807168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0963225) q[0];
sx q[0];
rz(-3.0131786) q[0];
sx q[0];
rz(1.1044109) q[0];
rz(-pi) q[1];
rz(-2.5727082) q[2];
sx q[2];
rz(-2.9279104) q[2];
sx q[2];
rz(-1.3485707) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0218483) q[1];
sx q[1];
rz(-1.1895864) q[1];
sx q[1];
rz(0.66701835) q[1];
x q[2];
rz(-3.0343247) q[3];
sx q[3];
rz(-1.7431431) q[3];
sx q[3];
rz(1.8037667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6803117) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(0.4168365) q[2];
rz(2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(2.7509287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9417435) q[0];
sx q[0];
rz(-1.2171634) q[0];
sx q[0];
rz(-1.4456277) q[0];
rz(2.4327714) q[1];
sx q[1];
rz(-2.2292021) q[1];
sx q[1];
rz(3.0880047) q[1];
rz(0.20119737) q[2];
sx q[2];
rz(-2.1977949) q[2];
sx q[2];
rz(2.0044873) q[2];
rz(0.96434595) q[3];
sx q[3];
rz(-2.3026005) q[3];
sx q[3];
rz(2.8540924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
