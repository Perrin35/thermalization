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
rz(0.024700392) q[0];
sx q[0];
rz(-1.2588809) q[0];
sx q[0];
rz(1.2727241) q[0];
rz(1.3920353) q[1];
sx q[1];
rz(-1.4104383) q[1];
sx q[1];
rz(2.2478204) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0964684) q[0];
sx q[0];
rz(-0.91334263) q[0];
sx q[0];
rz(-0.7492926) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29869334) q[2];
sx q[2];
rz(-1.325404) q[2];
sx q[2];
rz(-1.1358062) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53613816) q[1];
sx q[1];
rz(-1.0404603) q[1];
sx q[1];
rz(-2.082389) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36811604) q[3];
sx q[3];
rz(-1.5541346) q[3];
sx q[3];
rz(-2.1429495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.36665234) q[2];
sx q[2];
rz(-1.6795936) q[2];
sx q[2];
rz(2.6095663) q[2];
rz(-2.4074647) q[3];
sx q[3];
rz(-2.8668154) q[3];
sx q[3];
rz(-1.1067357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4183913) q[0];
sx q[0];
rz(-0.98576236) q[0];
sx q[0];
rz(3.0611839) q[0];
rz(2.6037727) q[1];
sx q[1];
rz(-2.0729013) q[1];
sx q[1];
rz(2.417876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88392576) q[0];
sx q[0];
rz(-0.91086713) q[0];
sx q[0];
rz(-2.1377449) q[0];
x q[1];
rz(-2.2494456) q[2];
sx q[2];
rz(-0.79024678) q[2];
sx q[2];
rz(2.9228766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.557505) q[1];
sx q[1];
rz(-2.1592916) q[1];
sx q[1];
rz(3.1397538) q[1];
rz(-pi) q[2];
rz(-2.6234954) q[3];
sx q[3];
rz(-2.2541917) q[3];
sx q[3];
rz(-2.4404614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.85372743) q[2];
sx q[2];
rz(-2.008805) q[2];
sx q[2];
rz(0.90616027) q[2];
rz(2.7791038) q[3];
sx q[3];
rz(-1.5972842) q[3];
sx q[3];
rz(3.0947065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2694117) q[0];
sx q[0];
rz(-1.325664) q[0];
sx q[0];
rz(1.0915225) q[0];
rz(0.12256924) q[1];
sx q[1];
rz(-1.7054319) q[1];
sx q[1];
rz(-1.7950119) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5751981) q[0];
sx q[0];
rz(-1.9267531) q[0];
sx q[0];
rz(-0.031513799) q[0];
rz(-pi) q[1];
rz(-0.41287072) q[2];
sx q[2];
rz(-0.70253583) q[2];
sx q[2];
rz(-1.2324049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.31981644) q[1];
sx q[1];
rz(-0.21142658) q[1];
sx q[1];
rz(-1.8199304) q[1];
x q[2];
rz(2.6123206) q[3];
sx q[3];
rz(-1.4733757) q[3];
sx q[3];
rz(-0.71089632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71318212) q[2];
sx q[2];
rz(-1.1298263) q[2];
sx q[2];
rz(-0.23207363) q[2];
rz(-1.170916) q[3];
sx q[3];
rz(-2.5210095) q[3];
sx q[3];
rz(-0.28461972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93037) q[0];
sx q[0];
rz(-2.3486597) q[0];
sx q[0];
rz(-1.5027745) q[0];
rz(0.59208313) q[1];
sx q[1];
rz(-1.1555669) q[1];
sx q[1];
rz(2.2519462) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8439595) q[0];
sx q[0];
rz(-1.0464764) q[0];
sx q[0];
rz(2.6056933) q[0];
x q[1];
rz(2.0060894) q[2];
sx q[2];
rz(-0.5331299) q[2];
sx q[2];
rz(-2.1338303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1210416) q[1];
sx q[1];
rz(-2.7747323) q[1];
sx q[1];
rz(-1.0479142) q[1];
x q[2];
rz(-0.38460807) q[3];
sx q[3];
rz(-1.1879564) q[3];
sx q[3];
rz(-2.2205381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4952116) q[2];
sx q[2];
rz(-1.4633598) q[2];
sx q[2];
rz(0.67592534) q[2];
rz(-1.4365139) q[3];
sx q[3];
rz(-2.8016475) q[3];
sx q[3];
rz(0.16726141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2109569) q[0];
sx q[0];
rz(-2.8093331) q[0];
sx q[0];
rz(0.19530547) q[0];
rz(3.0209814) q[1];
sx q[1];
rz(-2.9258525) q[1];
sx q[1];
rz(-2.9511071) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99283907) q[0];
sx q[0];
rz(-1.5061814) q[0];
sx q[0];
rz(-3.0535327) q[0];
rz(-pi) q[1];
rz(-2.4649925) q[2];
sx q[2];
rz(-1.5034505) q[2];
sx q[2];
rz(-1.2852033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45791679) q[1];
sx q[1];
rz(-1.4701587) q[1];
sx q[1];
rz(1.5186054) q[1];
rz(-pi) q[2];
rz(1.9917914) q[3];
sx q[3];
rz(-0.99830571) q[3];
sx q[3];
rz(-1.5508625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6327989) q[2];
sx q[2];
rz(-1.6570647) q[2];
sx q[2];
rz(-1.9602027) q[2];
rz(1.7975636) q[3];
sx q[3];
rz(-0.7414147) q[3];
sx q[3];
rz(0.77176362) q[3];
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
rz(pi/2) q[3];
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
rz(0.33264273) q[0];
sx q[0];
rz(-1.7985666) q[0];
sx q[0];
rz(1.1395662) q[0];
rz(0.66967669) q[1];
sx q[1];
rz(-0.91595903) q[1];
sx q[1];
rz(1.12961) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78167262) q[0];
sx q[0];
rz(-1.8501213) q[0];
sx q[0];
rz(-1.2646227) q[0];
x q[1];
rz(1.8367581) q[2];
sx q[2];
rz(-0.87087357) q[2];
sx q[2];
rz(1.0476607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5796389) q[1];
sx q[1];
rz(-0.41746751) q[1];
sx q[1];
rz(-2.6446656) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4127168) q[3];
sx q[3];
rz(-1.7493093) q[3];
sx q[3];
rz(-1.3016537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7149522) q[2];
sx q[2];
rz(-1.4098097) q[2];
sx q[2];
rz(2.7362774) q[2];
rz(-2.3146368) q[3];
sx q[3];
rz(-0.6260286) q[3];
sx q[3];
rz(2.2820182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.7760794) q[0];
sx q[0];
rz(-0.022553355) q[0];
sx q[0];
rz(-2.205701) q[0];
rz(2.2731958) q[1];
sx q[1];
rz(-2.6135542) q[1];
sx q[1];
rz(1.7220928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0763797) q[0];
sx q[0];
rz(-1.2166436) q[0];
sx q[0];
rz(-2.5962418) q[0];
rz(-2.9090857) q[2];
sx q[2];
rz(-2.0662796) q[2];
sx q[2];
rz(1.8956309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.81379997) q[1];
sx q[1];
rz(-0.84974242) q[1];
sx q[1];
rz(-0.93443296) q[1];
rz(-0.21162947) q[3];
sx q[3];
rz(-0.0018001477) q[3];
sx q[3];
rz(1.3789919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3203656) q[2];
sx q[2];
rz(-1.3996404) q[2];
sx q[2];
rz(2.121675) q[2];
rz(-0.50926456) q[3];
sx q[3];
rz(-1.441322) q[3];
sx q[3];
rz(-0.078484623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34198636) q[0];
sx q[0];
rz(-1.5155563) q[0];
sx q[0];
rz(-1.9827783) q[0];
rz(-0.47053567) q[1];
sx q[1];
rz(-1.7262986) q[1];
sx q[1];
rz(-1.6758957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19150951) q[0];
sx q[0];
rz(-1.5198764) q[0];
sx q[0];
rz(1.5945192) q[0];
rz(2.7297425) q[2];
sx q[2];
rz(-0.73397103) q[2];
sx q[2];
rz(-3.0628772) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.012124) q[1];
sx q[1];
rz(-2.2290342) q[1];
sx q[1];
rz(1.4736498) q[1];
x q[2];
rz(-0.37181446) q[3];
sx q[3];
rz(-0.68228693) q[3];
sx q[3];
rz(-0.90250795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1450242) q[2];
sx q[2];
rz(-1.8137167) q[2];
sx q[2];
rz(0.81929755) q[2];
rz(0.79646349) q[3];
sx q[3];
rz(-0.4070681) q[3];
sx q[3];
rz(-1.4888633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040937) q[0];
sx q[0];
rz(-0.096385328) q[0];
sx q[0];
rz(0.82164422) q[0];
rz(0.67283982) q[1];
sx q[1];
rz(-2.5655589) q[1];
sx q[1];
rz(-1.0427262) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0053528) q[0];
sx q[0];
rz(-1.7757105) q[0];
sx q[0];
rz(-1.4968275) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94093948) q[2];
sx q[2];
rz(-0.90412882) q[2];
sx q[2];
rz(-0.43064865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71127001) q[1];
sx q[1];
rz(-1.9414328) q[1];
sx q[1];
rz(2.1881585) q[1];
rz(-pi) q[2];
rz(1.4185227) q[3];
sx q[3];
rz(-1.8209753) q[3];
sx q[3];
rz(-2.6891249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78160703) q[2];
sx q[2];
rz(-1.2722509) q[2];
sx q[2];
rz(-2.2992415) q[2];
rz(2.602747) q[3];
sx q[3];
rz(-1.6539961) q[3];
sx q[3];
rz(0.52411383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.46776029) q[0];
sx q[0];
rz(-0.14586511) q[0];
sx q[0];
rz(-1.2354596) q[0];
rz(-0.94888672) q[1];
sx q[1];
rz(-2.3088539) q[1];
sx q[1];
rz(1.6611151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2713) q[0];
sx q[0];
rz(-0.90134951) q[0];
sx q[0];
rz(2.7842194) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83809488) q[2];
sx q[2];
rz(-0.75586719) q[2];
sx q[2];
rz(-0.91516337) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9292641) q[1];
sx q[1];
rz(-0.18106279) q[1];
sx q[1];
rz(0.4275112) q[1];
rz(1.3289024) q[3];
sx q[3];
rz(-1.8902361) q[3];
sx q[3];
rz(-2.0177359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.71198717) q[2];
sx q[2];
rz(-1.7179855) q[2];
sx q[2];
rz(0.31022662) q[2];
rz(-0.45983908) q[3];
sx q[3];
rz(-2.583677) q[3];
sx q[3];
rz(1.3742113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9643758) q[0];
sx q[0];
rz(-2.0073267) q[0];
sx q[0];
rz(1.0649756) q[0];
rz(-0.89827697) q[1];
sx q[1];
rz(-0.54294642) q[1];
sx q[1];
rz(-0.44566659) q[1];
rz(-2.0736135) q[2];
sx q[2];
rz(-0.72952727) q[2];
sx q[2];
rz(-2.9222957) q[2];
rz(0.38705714) q[3];
sx q[3];
rz(-2.195812) q[3];
sx q[3];
rz(0.14433911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
