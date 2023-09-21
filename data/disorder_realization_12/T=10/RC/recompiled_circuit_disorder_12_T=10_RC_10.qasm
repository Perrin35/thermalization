OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(-2.4175329) q[0];
sx q[0];
rz(1.6568503) q[0];
rz(-1.9384664) q[1];
sx q[1];
rz(-2.6180747) q[1];
sx q[1];
rz(0.88820052) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0523383) q[0];
sx q[0];
rz(-1.5262145) q[0];
sx q[0];
rz(2.9705439) q[0];
rz(-pi) q[1];
rz(0.16562478) q[2];
sx q[2];
rz(-2.1016444) q[2];
sx q[2];
rz(2.3373375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8656617) q[1];
sx q[1];
rz(-1.7377503) q[1];
sx q[1];
rz(1.5713515) q[1];
x q[2];
rz(-0.039460823) q[3];
sx q[3];
rz(-2.4426115) q[3];
sx q[3];
rz(-0.82074245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.858294) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(-2.0236012) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(2.4867687) q[0];
rz(-1.9251992) q[1];
sx q[1];
rz(-1.1647859) q[1];
sx q[1];
rz(-0.12589802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9585796) q[0];
sx q[0];
rz(-2.5218681) q[0];
sx q[0];
rz(-2.241914) q[0];
x q[1];
rz(-0.6217896) q[2];
sx q[2];
rz(-0.99025531) q[2];
sx q[2];
rz(2.7948126) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44015917) q[1];
sx q[1];
rz(-1.5408181) q[1];
sx q[1];
rz(-1.6975801) q[1];
rz(-pi) q[2];
rz(-0.84729362) q[3];
sx q[3];
rz(-1.5788659) q[3];
sx q[3];
rz(1.1524259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(2.753567) q[2];
rz(1.7175425) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(0.084005984) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.1598587) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.738049) q[0];
sx q[0];
rz(-1.4883211) q[0];
sx q[0];
rz(-2.0820046) q[0];
rz(-0.94295393) q[2];
sx q[2];
rz(-1.9938333) q[2];
sx q[2];
rz(-1.1689651) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68418903) q[1];
sx q[1];
rz(-2.0261129) q[1];
sx q[1];
rz(-1.9940358) q[1];
rz(-pi) q[2];
rz(2.1163164) q[3];
sx q[3];
rz(-1.6189515) q[3];
sx q[3];
rz(-2.2930876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40538654) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(-2.5260177) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(1.6595586) q[0];
rz(-0.7011134) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(2.8569417) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2157001) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(1.0924933) q[0];
x q[1];
rz(1.2830164) q[2];
sx q[2];
rz(-1.7680401) q[2];
sx q[2];
rz(1.4892727) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22420158) q[1];
sx q[1];
rz(-2.9926139) q[1];
sx q[1];
rz(2.1267183) q[1];
rz(-pi) q[2];
rz(0.986226) q[3];
sx q[3];
rz(-2.3826736) q[3];
sx q[3];
rz(-0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(0.56751928) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(2.0295985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48859566) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(-1.1524221) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(0.93793905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6803166) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(1.3374431) q[0];
x q[1];
rz(-2.1300975) q[2];
sx q[2];
rz(-1.1202295) q[2];
sx q[2];
rz(0.089103854) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8590947) q[1];
sx q[1];
rz(-1.5364093) q[1];
sx q[1];
rz(-3.1153468) q[1];
rz(-pi) q[2];
rz(0.33850833) q[3];
sx q[3];
rz(-2.7673628) q[3];
sx q[3];
rz(-0.81774536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(-1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43907169) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51150409) q[0];
sx q[0];
rz(-0.43404365) q[0];
sx q[0];
rz(-0.50891288) q[0];
rz(-pi) q[1];
rz(1.7700023) q[2];
sx q[2];
rz(-0.7253941) q[2];
sx q[2];
rz(0.2142011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6186657) q[1];
sx q[1];
rz(-1.3892801) q[1];
sx q[1];
rz(1.8261441) q[1];
rz(-pi) q[2];
rz(3.0895124) q[3];
sx q[3];
rz(-0.77643231) q[3];
sx q[3];
rz(-1.6805502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6094728) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(0.042908948) q[0];
rz(2.2242916) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(0.50061217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4681049) q[0];
sx q[0];
rz(-1.3930495) q[0];
sx q[0];
rz(-0.034981473) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31502864) q[2];
sx q[2];
rz(-1.050596) q[2];
sx q[2];
rz(-0.5859642) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.73270479) q[1];
sx q[1];
rz(-1.9703456) q[1];
sx q[1];
rz(1.7016181) q[1];
rz(-2.72243) q[3];
sx q[3];
rz(-0.8014285) q[3];
sx q[3];
rz(-0.49405801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6033972) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(-2.4659757) q[2];
rz(0.44089857) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(-0.52136695) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.065141) q[0];
sx q[0];
rz(-0.32442176) q[0];
sx q[0];
rz(1.0674397) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.4656461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7536403) q[0];
sx q[0];
rz(-2.6590829) q[0];
sx q[0];
rz(-1.8637191) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6797811) q[2];
sx q[2];
rz(-1.7562859) q[2];
sx q[2];
rz(-2.2456004) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9379942) q[1];
sx q[1];
rz(-1.9059056) q[1];
sx q[1];
rz(2.0957698) q[1];
x q[2];
rz(2.397714) q[3];
sx q[3];
rz(-2.3489967) q[3];
sx q[3];
rz(2.511123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33891588) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(-0.87351292) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2575689) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(0.73721686) q[0];
rz(3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-2.2163056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5697445) q[0];
sx q[0];
rz(-2.0664555) q[0];
sx q[0];
rz(-0.75832383) q[0];
rz(0.16889062) q[2];
sx q[2];
rz(-2.4714111) q[2];
sx q[2];
rz(2.9549213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9463897) q[1];
sx q[1];
rz(-0.22559799) q[1];
sx q[1];
rz(2.5220847) q[1];
x q[2];
rz(2.7221189) q[3];
sx q[3];
rz(-2.2443716) q[3];
sx q[3];
rz(-2.2485965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(0.090099661) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(-1.1130921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(2.8758077) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(0.75751799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.950338) q[0];
sx q[0];
rz(-1.3779103) q[0];
sx q[0];
rz(0.14001503) q[0];
rz(-pi) q[1];
rz(2.1610356) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(-0.57934258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0150891) q[1];
sx q[1];
rz(-1.9941829) q[1];
sx q[1];
rz(-0.24366118) q[1];
rz(-pi) q[2];
rz(3.1086139) q[3];
sx q[3];
rz(-2.4313201) q[3];
sx q[3];
rz(-0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1511128) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(2.347836) q[2];
rz(0.63888597) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6486075) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-2.7040504) q[2];
sx q[2];
rz(-1.2699288) q[2];
sx q[2];
rz(-1.5581836) q[2];
rz(-0.99132514) q[3];
sx q[3];
rz(-1.8425059) q[3];
sx q[3];
rz(-0.89210638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];