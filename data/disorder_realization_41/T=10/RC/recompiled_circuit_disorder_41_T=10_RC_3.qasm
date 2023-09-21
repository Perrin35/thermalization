OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37786814) q[0];
sx q[0];
rz(-1.2963429) q[0];
sx q[0];
rz(-0.26835693) q[0];
rz(-0.29834892) q[2];
sx q[2];
rz(-0.48765182) q[2];
sx q[2];
rz(-1.2717441) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5948407) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(-0.91083281) q[1];
x q[2];
rz(-1.0825726) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(-2.0522576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39820406) q[0];
sx q[0];
rz(-2.2687015) q[0];
sx q[0];
rz(-2.8845805) q[0];
rz(-pi) q[1];
rz(0.40541655) q[2];
sx q[2];
rz(-2.9581796) q[2];
sx q[2];
rz(-2.237052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0577382) q[1];
sx q[1];
rz(-1.6946304) q[1];
sx q[1];
rz(2.3480575) q[1];
x q[2];
rz(1.1284626) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(-0.80004063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-2.7868328) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86706485) q[0];
sx q[0];
rz(-1.9221677) q[0];
sx q[0];
rz(-1.685313) q[0];
rz(1.7180195) q[2];
sx q[2];
rz(-1.8185116) q[2];
sx q[2];
rz(-1.9927646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55089009) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(0.18443702) q[1];
x q[2];
rz(-2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(-2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53326398) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(0.91745013) q[0];
rz(-2.4687185) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-0.26487574) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0565856) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(2.4823275) q[0];
rz(-0.27819602) q[2];
sx q[2];
rz(-1.2831266) q[2];
sx q[2];
rz(-0.96565914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59135624) q[1];
sx q[1];
rz(-1.1657506) q[1];
sx q[1];
rz(0.77781271) q[1];
rz(-pi) q[2];
rz(2.3817252) q[3];
sx q[3];
rz(-2.2025975) q[3];
sx q[3];
rz(-2.2894273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-2.662861) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(2.1889401) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53341502) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(-1.7081225) q[0];
x q[1];
rz(2.7258337) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(1.8005467) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9597185) q[1];
sx q[1];
rz(-1.8078783) q[1];
sx q[1];
rz(-1.167776) q[1];
rz(-1.6126552) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(-2.4822513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7197363) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(2.9690202) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1738152) q[0];
sx q[0];
rz(-0.41668188) q[0];
sx q[0];
rz(1.5370876) q[0];
rz(-1.0748323) q[2];
sx q[2];
rz(-2.6325588) q[2];
sx q[2];
rz(-2.3376655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8166568) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(-1.1549969) q[1];
rz(1.3401838) q[3];
sx q[3];
rz(-0.95313822) q[3];
sx q[3];
rz(-2.2889083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(2.3849934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7998357) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(-2.1205649) q[0];
rz(-pi) q[1];
rz(-0.41140326) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(1.411737) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3488256) q[1];
sx q[1];
rz(-1.4139237) q[1];
sx q[1];
rz(1.3385593) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6301304) q[3];
sx q[3];
rz(-1.7620798) q[3];
sx q[3];
rz(0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(2.1527122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7635927) q[0];
sx q[0];
rz(-0.53056301) q[0];
sx q[0];
rz(-0.87688045) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093401508) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(-1.8708558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7372072) q[1];
sx q[1];
rz(-0.42006856) q[1];
sx q[1];
rz(-1.7350446) q[1];
x q[2];
rz(0.51560651) q[3];
sx q[3];
rz(-2.0022087) q[3];
sx q[3];
rz(-0.069375667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(0.60920238) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0926889) q[0];
sx q[0];
rz(-0.19177076) q[0];
sx q[0];
rz(-1.2094686) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74553211) q[2];
sx q[2];
rz(-2.9420256) q[2];
sx q[2];
rz(1.0996639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1210291) q[1];
sx q[1];
rz(-2.4596446) q[1];
sx q[1];
rz(0.49973439) q[1];
x q[2];
rz(-2.14823) q[3];
sx q[3];
rz(-1.1953029) q[3];
sx q[3];
rz(-1.4162228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(-0.8264181) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.7451161) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012750082) q[0];
sx q[0];
rz(-1.2484974) q[0];
sx q[0];
rz(-1.661257) q[0];
x q[1];
rz(2.4834677) q[2];
sx q[2];
rz(-1.3326367) q[2];
sx q[2];
rz(-0.70891526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9979981) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(2.8423611) q[1];
x q[2];
rz(2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(-2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(0.940154) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(0.084074323) q[3];
sx q[3];
rz(-2.6064059) q[3];
sx q[3];
rz(2.4168766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];