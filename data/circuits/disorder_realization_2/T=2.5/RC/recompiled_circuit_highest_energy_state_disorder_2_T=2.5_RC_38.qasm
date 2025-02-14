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
rz(-0.19718117) q[0];
sx q[0];
rz(-2.0694759) q[0];
sx q[0];
rz(1.5972692) q[0];
rz(-2.7838669) q[1];
sx q[1];
rz(-1.8341213) q[1];
sx q[1];
rz(2.7657685) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8944112) q[0];
sx q[0];
rz(-2.3355749) q[0];
sx q[0];
rz(-1.6707808) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2959458) q[2];
sx q[2];
rz(-1.0604727) q[2];
sx q[2];
rz(-0.21462552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.42106146) q[1];
sx q[1];
rz(-2.1481247) q[1];
sx q[1];
rz(2.4666822) q[1];
x q[2];
rz(1.1268257) q[3];
sx q[3];
rz(-2.4587963) q[3];
sx q[3];
rz(1.7449774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69349217) q[2];
sx q[2];
rz(-0.38303146) q[2];
sx q[2];
rz(2.7310272) q[2];
rz(-0.52452123) q[3];
sx q[3];
rz(-2.9295242) q[3];
sx q[3];
rz(-1.9922403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18325026) q[0];
sx q[0];
rz(-0.14160937) q[0];
sx q[0];
rz(-0.70567065) q[0];
rz(1.2880098) q[1];
sx q[1];
rz(-2.5649004) q[1];
sx q[1];
rz(-2.0075331) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3243038) q[0];
sx q[0];
rz(-1.1686443) q[0];
sx q[0];
rz(2.3020846) q[0];
rz(-pi) q[1];
rz(-1.7524598) q[2];
sx q[2];
rz(-1.2328892) q[2];
sx q[2];
rz(-2.365571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6780939) q[1];
sx q[1];
rz(-1.948852) q[1];
sx q[1];
rz(-0.75611703) q[1];
x q[2];
rz(0.56523486) q[3];
sx q[3];
rz(-1.4197403) q[3];
sx q[3];
rz(-2.8961669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80295339) q[2];
sx q[2];
rz(-2.6889375) q[2];
sx q[2];
rz(2.6111641) q[2];
rz(-0.28702921) q[3];
sx q[3];
rz(-2.65286) q[3];
sx q[3];
rz(0.73634017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257783) q[0];
sx q[0];
rz(-1.1017355) q[0];
sx q[0];
rz(2.0871479) q[0];
rz(-1.4796481) q[1];
sx q[1];
rz(-2.2190861) q[1];
sx q[1];
rz(1.066801) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22196348) q[0];
sx q[0];
rz(-2.9502075) q[0];
sx q[0];
rz(-1.3758977) q[0];
x q[1];
rz(0.11458204) q[2];
sx q[2];
rz(-2.4620585) q[2];
sx q[2];
rz(0.058341786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.40854) q[1];
sx q[1];
rz(-2.618737) q[1];
sx q[1];
rz(1.1001292) q[1];
rz(-pi) q[2];
rz(-1.9122804) q[3];
sx q[3];
rz(-0.74401281) q[3];
sx q[3];
rz(0.0014425576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58831424) q[2];
sx q[2];
rz(-1.30043) q[2];
sx q[2];
rz(1.1021357) q[2];
rz(-2.6831902) q[3];
sx q[3];
rz(-1.5285243) q[3];
sx q[3];
rz(-0.63157356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.691953) q[0];
sx q[0];
rz(-2.3821558) q[0];
sx q[0];
rz(0.48680437) q[0];
rz(-1.5454769) q[1];
sx q[1];
rz(-2.7524452) q[1];
sx q[1];
rz(-2.5130689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8299943) q[0];
sx q[0];
rz(-1.7033368) q[0];
sx q[0];
rz(1.9674942) q[0];
rz(-pi) q[1];
rz(-1.7968463) q[2];
sx q[2];
rz(-0.58494431) q[2];
sx q[2];
rz(1.9247885) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13884096) q[1];
sx q[1];
rz(-1.3831487) q[1];
sx q[1];
rz(-1.1739338) q[1];
rz(0.071672385) q[3];
sx q[3];
rz(-2.1183156) q[3];
sx q[3];
rz(-1.6184834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5635166) q[2];
sx q[2];
rz(-3.0264644) q[2];
sx q[2];
rz(-0.76187491) q[2];
rz(2.9231807) q[3];
sx q[3];
rz(-0.29774791) q[3];
sx q[3];
rz(1.0485605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31576425) q[0];
sx q[0];
rz(-0.35218069) q[0];
sx q[0];
rz(2.5370362) q[0];
rz(-1.8479337) q[1];
sx q[1];
rz(-2.3217521) q[1];
sx q[1];
rz(0.34436071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0492951) q[0];
sx q[0];
rz(-1.5504748) q[0];
sx q[0];
rz(-1.4882213) q[0];
rz(-1.9042964) q[2];
sx q[2];
rz(-2.1715559) q[2];
sx q[2];
rz(1.5983901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1923848) q[1];
sx q[1];
rz(-2.4864962) q[1];
sx q[1];
rz(0.75836436) q[1];
x q[2];
rz(-2.133338) q[3];
sx q[3];
rz(-0.75762094) q[3];
sx q[3];
rz(-1.7239095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1206426) q[2];
sx q[2];
rz(-0.63434333) q[2];
sx q[2];
rz(2.9316588) q[2];
rz(2.0338992) q[3];
sx q[3];
rz(-1.451713) q[3];
sx q[3];
rz(-0.2618739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27860677) q[0];
sx q[0];
rz(-2.5513273) q[0];
sx q[0];
rz(-3.0721831) q[0];
rz(-3.1076) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(-0.45190826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44563189) q[0];
sx q[0];
rz(-1.1813786) q[0];
sx q[0];
rz(1.7444064) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68695171) q[2];
sx q[2];
rz(-1.0712306) q[2];
sx q[2];
rz(0.14909185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.375651) q[1];
sx q[1];
rz(-3.1187291) q[1];
sx q[1];
rz(1.2506865) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8256559) q[3];
sx q[3];
rz(-1.3075124) q[3];
sx q[3];
rz(1.7039434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41428337) q[2];
sx q[2];
rz(-2.3392129) q[2];
sx q[2];
rz(1.9963973) q[2];
rz(-2.1833873) q[3];
sx q[3];
rz(-1.9327791) q[3];
sx q[3];
rz(-2.8632979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77022839) q[0];
sx q[0];
rz(-2.6537277) q[0];
sx q[0];
rz(-2.2801953) q[0];
rz(2.1791012) q[1];
sx q[1];
rz(-0.95314127) q[1];
sx q[1];
rz(-0.15693396) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.033534) q[0];
sx q[0];
rz(-2.3367051) q[0];
sx q[0];
rz(-1.0548841) q[0];
rz(-pi) q[1];
rz(-0.10617039) q[2];
sx q[2];
rz(-1.9261126) q[2];
sx q[2];
rz(-2.4715142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2409901) q[1];
sx q[1];
rz(-2.1442736) q[1];
sx q[1];
rz(2.1585224) q[1];
rz(-pi) q[2];
rz(-3.0928218) q[3];
sx q[3];
rz(-1.3579821) q[3];
sx q[3];
rz(-1.770894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.018205) q[2];
sx q[2];
rz(-0.79424477) q[2];
sx q[2];
rz(-2.0920853) q[2];
rz(-2.6617995) q[3];
sx q[3];
rz(-0.19194651) q[3];
sx q[3];
rz(0.16201365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4341693) q[0];
sx q[0];
rz(-2.8987085) q[0];
sx q[0];
rz(-2.6450787) q[0];
rz(-1.4855344) q[1];
sx q[1];
rz(-0.86692202) q[1];
sx q[1];
rz(0.022580126) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2537947) q[0];
sx q[0];
rz(-0.27844515) q[0];
sx q[0];
rz(0.87775567) q[0];
rz(-pi) q[1];
rz(-2.0099925) q[2];
sx q[2];
rz(-0.15838693) q[2];
sx q[2];
rz(1.8902376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3564313) q[1];
sx q[1];
rz(-1.3884434) q[1];
sx q[1];
rz(-2.3556832) q[1];
rz(-pi) q[2];
rz(-0.19666861) q[3];
sx q[3];
rz(-1.0447352) q[3];
sx q[3];
rz(-2.9064309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61985832) q[2];
sx q[2];
rz(-1.1405742) q[2];
sx q[2];
rz(-0.015259585) q[2];
rz(-2.7742625) q[3];
sx q[3];
rz(-2.6945249) q[3];
sx q[3];
rz(0.36146155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5527282) q[0];
sx q[0];
rz(-2.9155154) q[0];
sx q[0];
rz(0.077762522) q[0];
rz(-3.0231158) q[1];
sx q[1];
rz(-2.8006554) q[1];
sx q[1];
rz(-2.4854787) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828147) q[0];
sx q[0];
rz(-1.4442632) q[0];
sx q[0];
rz(-1.2650498) q[0];
rz(-pi) q[1];
rz(-1.3030197) q[2];
sx q[2];
rz(-1.736044) q[2];
sx q[2];
rz(2.6337773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1826577) q[1];
sx q[1];
rz(-1.9709341) q[1];
sx q[1];
rz(-1.047664) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3967152) q[3];
sx q[3];
rz(-0.69760579) q[3];
sx q[3];
rz(0.11817237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6180827) q[2];
sx q[2];
rz(-1.9205576) q[2];
sx q[2];
rz(-2.6005884) q[2];
rz(-0.4506166) q[3];
sx q[3];
rz(-1.6559947) q[3];
sx q[3];
rz(-0.97644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4377874) q[0];
sx q[0];
rz(-1.3492275) q[0];
sx q[0];
rz(-3.034814) q[0];
rz(-1.59185) q[1];
sx q[1];
rz(-0.50250643) q[1];
sx q[1];
rz(1.7639311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6984682) q[0];
sx q[0];
rz(-1.5266935) q[0];
sx q[0];
rz(-2.0206142) q[0];
x q[1];
rz(-0.3129199) q[2];
sx q[2];
rz(-1.7831137) q[2];
sx q[2];
rz(-2.2231399) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71860524) q[1];
sx q[1];
rz(-1.3235301) q[1];
sx q[1];
rz(0.24691564) q[1];
x q[2];
rz(2.315134) q[3];
sx q[3];
rz(-0.71751201) q[3];
sx q[3];
rz(-0.96271587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0926853) q[2];
sx q[2];
rz(-1.7203628) q[2];
sx q[2];
rz(0.29472026) q[2];
rz(-0.60485351) q[3];
sx q[3];
rz(-1.2302783) q[3];
sx q[3];
rz(-0.10943432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43815628) q[0];
sx q[0];
rz(-1.3180757) q[0];
sx q[0];
rz(-0.99676589) q[0];
rz(3.1238212) q[1];
sx q[1];
rz(-1.2635063) q[1];
sx q[1];
rz(1.8020188) q[1];
rz(-0.80452917) q[2];
sx q[2];
rz(-1.1297516) q[2];
sx q[2];
rz(-0.17937112) q[2];
rz(0.83497077) q[3];
sx q[3];
rz(-2.2899262) q[3];
sx q[3];
rz(-0.038577608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
