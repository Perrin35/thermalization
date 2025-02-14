OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(-2.3999441) q[1];
sx q[1];
rz(2.9902966) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79729743) q[0];
sx q[0];
rz(-1.5452641) q[0];
sx q[0];
rz(-1.6961369) q[0];
rz(1.0213576) q[2];
sx q[2];
rz(-1.4714699) q[2];
sx q[2];
rz(-1.316357) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94231168) q[1];
sx q[1];
rz(-2.1122871) q[1];
sx q[1];
rz(-0.6813867) q[1];
rz(-pi) q[2];
rz(2.862418) q[3];
sx q[3];
rz(-0.82133365) q[3];
sx q[3];
rz(0.73842919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0155045) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(-0.31910953) q[2];
rz(-1.1110405) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(-1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023271712) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(-0.58854377) q[0];
rz(0.59666807) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4174735) q[0];
sx q[0];
rz(-0.80097526) q[0];
sx q[0];
rz(2.1885314) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0701635) q[2];
sx q[2];
rz(-2.0703346) q[2];
sx q[2];
rz(-0.75354924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10728285) q[1];
sx q[1];
rz(-1.4713381) q[1];
sx q[1];
rz(-2.0080272) q[1];
rz(0.94542687) q[3];
sx q[3];
rz(-1.6540048) q[3];
sx q[3];
rz(1.041846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1179463) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(-0.092197593) q[2];
rz(-0.89208952) q[3];
sx q[3];
rz(-0.8725608) q[3];
sx q[3];
rz(1.0649072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43111619) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(-0.098467501) q[0];
rz(-2.5473728) q[1];
sx q[1];
rz(-1.0585982) q[1];
sx q[1];
rz(2.1509511) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38483175) q[0];
sx q[0];
rz(-2.7568026) q[0];
sx q[0];
rz(1.8033474) q[0];
rz(2.3805101) q[2];
sx q[2];
rz(-1.154261) q[2];
sx q[2];
rz(2.0628945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9746863) q[1];
sx q[1];
rz(-0.94407996) q[1];
sx q[1];
rz(2.8381366) q[1];
x q[2];
rz(-2.0001171) q[3];
sx q[3];
rz(-2.6344382) q[3];
sx q[3];
rz(2.0961268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6777665) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(2.3393935) q[2];
rz(-1.5944611) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(0.86435634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39740729) q[0];
sx q[0];
rz(-2.7825401) q[0];
sx q[0];
rz(1.3209976) q[0];
rz(1.6162704) q[1];
sx q[1];
rz(-0.95322144) q[1];
sx q[1];
rz(-2.9811409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4686543) q[0];
sx q[0];
rz(-0.78862353) q[0];
sx q[0];
rz(-1.9140052) q[0];
rz(-0.58893369) q[2];
sx q[2];
rz(-2.9351165) q[2];
sx q[2];
rz(-1.4331417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7509979) q[1];
sx q[1];
rz(-0.92727755) q[1];
sx q[1];
rz(2.949763) q[1];
x q[2];
rz(-2.5009551) q[3];
sx q[3];
rz(-2.5036219) q[3];
sx q[3];
rz(-2.9393794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3204331) q[2];
sx q[2];
rz(-1.4904138) q[2];
sx q[2];
rz(0.56905812) q[2];
rz(1.2218366) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4959167) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(2.3107279) q[0];
rz(-1.9310541) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(2.1499706) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9308335) q[0];
sx q[0];
rz(-0.88085876) q[0];
sx q[0];
rz(2.7372975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8610271) q[2];
sx q[2];
rz(-2.2380851) q[2];
sx q[2];
rz(1.9633121) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6779532) q[1];
sx q[1];
rz(-1.7262184) q[1];
sx q[1];
rz(0.83828204) q[1];
rz(-0.95007054) q[3];
sx q[3];
rz(-0.80214989) q[3];
sx q[3];
rz(-0.42250326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7603989) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(-0.35923108) q[2];
rz(0.71581101) q[3];
sx q[3];
rz(-0.78407136) q[3];
sx q[3];
rz(-1.951096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0642218) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(-2.6203058) q[0];
rz(-2.298666) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(0.11046031) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4058286) q[0];
sx q[0];
rz(-0.89019201) q[0];
sx q[0];
rz(-1.8981147) q[0];
x q[1];
rz(-3.1119124) q[2];
sx q[2];
rz(-0.92845193) q[2];
sx q[2];
rz(-0.50291598) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46347324) q[1];
sx q[1];
rz(-1.3729401) q[1];
sx q[1];
rz(-1.0956647) q[1];
rz(-pi) q[2];
rz(-2.8824948) q[3];
sx q[3];
rz(-2.0784335) q[3];
sx q[3];
rz(-0.51578427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9342186) q[2];
sx q[2];
rz(-1.5418345) q[2];
sx q[2];
rz(-0.99384394) q[2];
rz(-2.5908616) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(-1.5229092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9024502) q[0];
sx q[0];
rz(-0.37343326) q[0];
sx q[0];
rz(2.9504839) q[0];
rz(-0.36901078) q[1];
sx q[1];
rz(-1.399682) q[1];
sx q[1];
rz(0.63327995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12145081) q[0];
sx q[0];
rz(-0.25853863) q[0];
sx q[0];
rz(2.7516737) q[0];
rz(0.72004135) q[2];
sx q[2];
rz(-1.9332464) q[2];
sx q[2];
rz(-1.3316278) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5577199) q[1];
sx q[1];
rz(-0.28350949) q[1];
sx q[1];
rz(-2.3260172) q[1];
rz(2.4426887) q[3];
sx q[3];
rz(-1.4443384) q[3];
sx q[3];
rz(1.5497643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17948991) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(-2.7607259) q[2];
rz(-1.3880091) q[3];
sx q[3];
rz(-1.0264779) q[3];
sx q[3];
rz(1.4306205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4717344) q[0];
sx q[0];
rz(-2.7405881) q[0];
sx q[0];
rz(-1.5995455) q[0];
rz(-1.3767287) q[1];
sx q[1];
rz(-1.3841261) q[1];
sx q[1];
rz(-0.9309887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5409398) q[0];
sx q[0];
rz(-0.64221817) q[0];
sx q[0];
rz(-0.14677958) q[0];
rz(2.6561894) q[2];
sx q[2];
rz(-0.50386643) q[2];
sx q[2];
rz(-0.47270838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37759763) q[1];
sx q[1];
rz(-0.74255172) q[1];
sx q[1];
rz(1.5686036) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92967195) q[3];
sx q[3];
rz(-2.3768209) q[3];
sx q[3];
rz(1.1906033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.57637438) q[2];
sx q[2];
rz(-2.5810676) q[2];
sx q[2];
rz(3.0726688) q[2];
rz(-1.8779514) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(0.69839683) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7428335) q[0];
sx q[0];
rz(-0.95785207) q[0];
sx q[0];
rz(0.051890705) q[0];
rz(2.9715111) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(-0.35194078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1387685) q[0];
sx q[0];
rz(-1.616298) q[0];
sx q[0];
rz(-0.36624927) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1743618) q[2];
sx q[2];
rz(-1.033342) q[2];
sx q[2];
rz(-3.1307901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7406977) q[1];
sx q[1];
rz(-2.4603421) q[1];
sx q[1];
rz(-0.91886144) q[1];
rz(-1.3424344) q[3];
sx q[3];
rz(-0.31410892) q[3];
sx q[3];
rz(2.4921162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6841782) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(-2.7049098) q[2];
rz(-0.39143482) q[3];
sx q[3];
rz(-1.4623564) q[3];
sx q[3];
rz(-2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6201685) q[0];
sx q[0];
rz(-2.4346209) q[0];
sx q[0];
rz(0.11908764) q[0];
rz(-1.8424312) q[1];
sx q[1];
rz(-1.3928587) q[1];
sx q[1];
rz(-1.3577168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1022259) q[0];
sx q[0];
rz(-2.358317) q[0];
sx q[0];
rz(-2.164515) q[0];
rz(-2.3215649) q[2];
sx q[2];
rz(-1.8402647) q[2];
sx q[2];
rz(-2.7633689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57798856) q[1];
sx q[1];
rz(-1.23659) q[1];
sx q[1];
rz(2.5086705) q[1];
rz(2.6068654) q[3];
sx q[3];
rz(-0.89294723) q[3];
sx q[3];
rz(1.1191561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.208821) q[2];
sx q[2];
rz(-1.5522771) q[2];
sx q[2];
rz(0.61441747) q[2];
rz(1.5116073) q[3];
sx q[3];
rz(-0.67174086) q[3];
sx q[3];
rz(-3.0577799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83196249) q[0];
sx q[0];
rz(-0.93292581) q[0];
sx q[0];
rz(0.25136872) q[0];
rz(2.108719) q[1];
sx q[1];
rz(-1.947247) q[1];
sx q[1];
rz(-1.4364545) q[1];
rz(3.0713785) q[2];
sx q[2];
rz(-1.4660942) q[2];
sx q[2];
rz(1.6906665) q[2];
rz(-2.4182416) q[3];
sx q[3];
rz(-2.1925329) q[3];
sx q[3];
rz(-0.88919269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
