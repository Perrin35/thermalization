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
rz(0.83956194) q[0];
sx q[0];
rz(3.8008939) q[0];
sx q[0];
rz(10.492926) q[0];
rz(0.90574342) q[1];
sx q[1];
rz(-1.4248166) q[1];
sx q[1];
rz(-0.6583156) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5323344) q[0];
sx q[0];
rz(-2.242385) q[0];
sx q[0];
rz(0.39686578) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3320902) q[2];
sx q[2];
rz(-0.90478071) q[2];
sx q[2];
rz(1.314338) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5707455) q[1];
sx q[1];
rz(-2.1553504) q[1];
sx q[1];
rz(-2.5258614) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2194276) q[3];
sx q[3];
rz(-2.0982705) q[3];
sx q[3];
rz(-1.4722401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1540404) q[2];
sx q[2];
rz(-1.1215569) q[2];
sx q[2];
rz(-1.7169607) q[2];
rz(2.3512261) q[3];
sx q[3];
rz(-2.3873886) q[3];
sx q[3];
rz(-2.8823901) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454999) q[0];
sx q[0];
rz(-0.51829618) q[0];
sx q[0];
rz(-1.8081283) q[0];
rz(-1.2893527) q[1];
sx q[1];
rz(-2.5080296) q[1];
sx q[1];
rz(-0.11817558) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2424392) q[0];
sx q[0];
rz(-2.1964873) q[0];
sx q[0];
rz(-3.1345128) q[0];
rz(-pi) q[1];
rz(-2.1344654) q[2];
sx q[2];
rz(-1.7810797) q[2];
sx q[2];
rz(-1.0894094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3505248) q[1];
sx q[1];
rz(-1.9842897) q[1];
sx q[1];
rz(0.86422635) q[1];
rz(-pi) q[2];
rz(2.4746861) q[3];
sx q[3];
rz(-0.42508891) q[3];
sx q[3];
rz(0.57890427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3893343) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(0.83750677) q[2];
rz(2.5181455) q[3];
sx q[3];
rz(-2.0313171) q[3];
sx q[3];
rz(0.71201223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74152827) q[0];
sx q[0];
rz(-0.15151227) q[0];
sx q[0];
rz(2.9140299) q[0];
rz(-0.96861068) q[1];
sx q[1];
rz(-2.249735) q[1];
sx q[1];
rz(-1.6516986) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9084599) q[0];
sx q[0];
rz(-1.8709749) q[0];
sx q[0];
rz(1.079193) q[0];
rz(-pi) q[1];
x q[1];
rz(0.030185791) q[2];
sx q[2];
rz(-1.7329114) q[2];
sx q[2];
rz(2.8679304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35830733) q[1];
sx q[1];
rz(-2.6926125) q[1];
sx q[1];
rz(2.994356) q[1];
rz(-pi) q[2];
rz(1.8996542) q[3];
sx q[3];
rz(-0.76995459) q[3];
sx q[3];
rz(-2.9640523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0861133) q[2];
sx q[2];
rz(-1.5429292) q[2];
sx q[2];
rz(2.1231269) q[2];
rz(-0.57824072) q[3];
sx q[3];
rz(-2.6522804) q[3];
sx q[3];
rz(-1.9448634) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8483491) q[0];
sx q[0];
rz(-1.542955) q[0];
sx q[0];
rz(1.6795213) q[0];
rz(-2.2734185) q[1];
sx q[1];
rz(-1.9055007) q[1];
sx q[1];
rz(-0.47913924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7443607) q[0];
sx q[0];
rz(-1.8313837) q[0];
sx q[0];
rz(0.67407383) q[0];
rz(-pi) q[1];
rz(-2.5136679) q[2];
sx q[2];
rz(-0.52573813) q[2];
sx q[2];
rz(3.1319654) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9769765) q[1];
sx q[1];
rz(-2.6880126) q[1];
sx q[1];
rz(-0.11670392) q[1];
rz(1.3105743) q[3];
sx q[3];
rz(-1.3129915) q[3];
sx q[3];
rz(1.5797404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9118871) q[2];
sx q[2];
rz(-2.2356326) q[2];
sx q[2];
rz(0.68944302) q[2];
rz(0.42129579) q[3];
sx q[3];
rz(-0.7553941) q[3];
sx q[3];
rz(-1.9740483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(0.01920779) q[0];
sx q[0];
rz(-0.27007073) q[0];
sx q[0];
rz(-0.58575678) q[0];
rz(-2.6981804) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(2.401039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6218005) q[0];
sx q[0];
rz(-0.92042887) q[0];
sx q[0];
rz(-0.67861329) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29467543) q[2];
sx q[2];
rz(-2.1496219) q[2];
sx q[2];
rz(-0.097172849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.003518) q[1];
sx q[1];
rz(-1.8205531) q[1];
sx q[1];
rz(-1.300092) q[1];
rz(-pi) q[2];
rz(-1.8216474) q[3];
sx q[3];
rz(-1.2870868) q[3];
sx q[3];
rz(-1.914806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1195141) q[2];
sx q[2];
rz(-1.7064648) q[2];
sx q[2];
rz(-0.73149991) q[2];
rz(-0.59705192) q[3];
sx q[3];
rz(-2.4686333) q[3];
sx q[3];
rz(1.2520242) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6325833) q[0];
sx q[0];
rz(-2.3409797) q[0];
sx q[0];
rz(-1.5923694) q[0];
rz(2.1005311) q[1];
sx q[1];
rz(-0.58716232) q[1];
sx q[1];
rz(-0.36137533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072424732) q[0];
sx q[0];
rz(-2.4033045) q[0];
sx q[0];
rz(-0.17743547) q[0];
rz(1.560426) q[2];
sx q[2];
rz(-0.25744312) q[2];
sx q[2];
rz(2.4598222) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6238072) q[1];
sx q[1];
rz(-1.1248765) q[1];
sx q[1];
rz(1.2257027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18852461) q[3];
sx q[3];
rz(-2.162121) q[3];
sx q[3];
rz(0.44805474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8063712) q[2];
sx q[2];
rz(-2.1101895) q[2];
sx q[2];
rz(-1.6356989) q[2];
rz(-2.7527572) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(-0.65765643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036670551) q[0];
sx q[0];
rz(-2.4668283) q[0];
sx q[0];
rz(-2.2210333) q[0];
rz(2.8432644) q[1];
sx q[1];
rz(-1.0519271) q[1];
sx q[1];
rz(-1.9783609) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0759493) q[0];
sx q[0];
rz(-1.8952184) q[0];
sx q[0];
rz(1.5604273) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3506868) q[2];
sx q[2];
rz(-1.8703572) q[2];
sx q[2];
rz(-1.1096481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.071347941) q[1];
sx q[1];
rz(-1.5694968) q[1];
sx q[1];
rz(-2.595129) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1441394) q[3];
sx q[3];
rz(-2.8152044) q[3];
sx q[3];
rz(-0.24392621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7982911) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(-0.11246559) q[2];
rz(0.28137842) q[3];
sx q[3];
rz(-2.3795542) q[3];
sx q[3];
rz(-1.5828097) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52044582) q[0];
sx q[0];
rz(-2.1805094) q[0];
sx q[0];
rz(3.0302826) q[0];
rz(0.81980199) q[1];
sx q[1];
rz(-0.83378053) q[1];
sx q[1];
rz(0.05096635) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726377) q[0];
sx q[0];
rz(-1.8678471) q[0];
sx q[0];
rz(-2.9651276) q[0];
rz(-1.3459008) q[2];
sx q[2];
rz(-0.95976613) q[2];
sx q[2];
rz(-2.261206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.513235) q[1];
sx q[1];
rz(-1.7019148) q[1];
sx q[1];
rz(0.031946957) q[1];
rz(-pi) q[2];
rz(2.1858864) q[3];
sx q[3];
rz(-0.73266163) q[3];
sx q[3];
rz(0.77377711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.082197949) q[2];
sx q[2];
rz(-1.1092564) q[2];
sx q[2];
rz(1.294403) q[2];
rz(-2.1829677) q[3];
sx q[3];
rz(-0.38161033) q[3];
sx q[3];
rz(-0.81058782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628919) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(0.94619757) q[0];
rz(1.9323438) q[1];
sx q[1];
rz(-0.47018662) q[1];
sx q[1];
rz(1.6016003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9935051) q[0];
sx q[0];
rz(-2.3984841) q[0];
sx q[0];
rz(1.1412727) q[0];
rz(0.49175434) q[2];
sx q[2];
rz(-1.7778313) q[2];
sx q[2];
rz(1.4018261) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18064776) q[1];
sx q[1];
rz(-1.0156609) q[1];
sx q[1];
rz(1.8349232) q[1];
rz(-pi) q[2];
rz(-0.33312914) q[3];
sx q[3];
rz(-0.58099834) q[3];
sx q[3];
rz(-2.6521386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9340747) q[2];
sx q[2];
rz(-1.8002276) q[2];
sx q[2];
rz(1.5156457) q[2];
rz(-2.3385284) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(-1.1481185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.76252) q[0];
sx q[0];
rz(-1.2305434) q[0];
sx q[0];
rz(1.3364963) q[0];
rz(-1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(1.9660827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7911691) q[0];
sx q[0];
rz(-0.8514815) q[0];
sx q[0];
rz(-1.8940303) q[0];
rz(-0.78091518) q[2];
sx q[2];
rz(-1.7405207) q[2];
sx q[2];
rz(0.11955027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6503295) q[1];
sx q[1];
rz(-1.6230738) q[1];
sx q[1];
rz(0.13904528) q[1];
x q[2];
rz(-1.8776348) q[3];
sx q[3];
rz(-2.1458907) q[3];
sx q[3];
rz(1.8317312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.056957873) q[2];
sx q[2];
rz(-1.6677758) q[2];
sx q[2];
rz(-1.3062306) q[2];
rz(-2.7133387) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(0.38124198) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.320095) q[0];
sx q[0];
rz(-1.4023254) q[0];
sx q[0];
rz(-1.0925972) q[0];
rz(-0.26662695) q[1];
sx q[1];
rz(-1.0933924) q[1];
sx q[1];
rz(-0.46270121) q[1];
rz(-0.54556432) q[2];
sx q[2];
rz(-0.64020276) q[2];
sx q[2];
rz(1.1545867) q[2];
rz(-1.7206031) q[3];
sx q[3];
rz(-1.5021742) q[3];
sx q[3];
rz(-1.4175137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
