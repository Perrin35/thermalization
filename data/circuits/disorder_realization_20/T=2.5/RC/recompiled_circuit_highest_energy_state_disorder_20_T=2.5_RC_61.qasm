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
rz(-0.72467726) q[0];
sx q[0];
rz(-1.8301355) q[0];
sx q[0];
rz(2.3367982) q[0];
rz(-2.5465487) q[1];
sx q[1];
rz(-0.19785985) q[1];
sx q[1];
rz(-1.2208389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9210584) q[0];
sx q[0];
rz(-3.0166114) q[0];
sx q[0];
rz(0.28470953) q[0];
rz(-pi) q[1];
rz(2.6788459) q[2];
sx q[2];
rz(-1.695206) q[2];
sx q[2];
rz(-0.98906987) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2455622) q[1];
sx q[1];
rz(-0.82791019) q[1];
sx q[1];
rz(1.6309392) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5086256) q[3];
sx q[3];
rz(-1.4598914) q[3];
sx q[3];
rz(-1.9680915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7120984) q[2];
sx q[2];
rz(-2.6338989) q[2];
sx q[2];
rz(-1.3362159) q[2];
rz(-1.3974238) q[3];
sx q[3];
rz(-1.6349399) q[3];
sx q[3];
rz(1.5857182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812688) q[0];
sx q[0];
rz(-1.6345432) q[0];
sx q[0];
rz(0.67833483) q[0];
rz(-2.5256269) q[1];
sx q[1];
rz(-2.1149642) q[1];
sx q[1];
rz(2.9982627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6852138) q[0];
sx q[0];
rz(-0.36093484) q[0];
sx q[0];
rz(1.1120615) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4392339) q[2];
sx q[2];
rz(-2.4189197) q[2];
sx q[2];
rz(-0.45718788) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2611609) q[1];
sx q[1];
rz(-1.6182233) q[1];
sx q[1];
rz(-2.6359184) q[1];
x q[2];
rz(2.7160104) q[3];
sx q[3];
rz(-1.2493361) q[3];
sx q[3];
rz(-1.5150573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6501288) q[2];
sx q[2];
rz(-2.3626707) q[2];
sx q[2];
rz(1.2196563) q[2];
rz(-0.31798142) q[3];
sx q[3];
rz(-0.82428437) q[3];
sx q[3];
rz(-0.73941755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7431444) q[0];
sx q[0];
rz(-2.2205181) q[0];
sx q[0];
rz(0.63968023) q[0];
rz(-1.3321053) q[1];
sx q[1];
rz(-2.2001241) q[1];
sx q[1];
rz(-2.926362) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71986854) q[0];
sx q[0];
rz(-0.54066896) q[0];
sx q[0];
rz(2.7175277) q[0];
rz(-pi) q[1];
rz(0.0052506487) q[2];
sx q[2];
rz(-1.8551146) q[2];
sx q[2];
rz(-2.9470987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0063246) q[1];
sx q[1];
rz(-2.746114) q[1];
sx q[1];
rz(1.5979824) q[1];
x q[2];
rz(-1.1138737) q[3];
sx q[3];
rz(-1.6789241) q[3];
sx q[3];
rz(0.34764744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8362063) q[2];
sx q[2];
rz(-0.47577327) q[2];
sx q[2];
rz(2.9212908) q[2];
rz(2.3588755) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(0.14557423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43133217) q[0];
sx q[0];
rz(-0.98589698) q[0];
sx q[0];
rz(-1.8248935) q[0];
rz(0.45007625) q[1];
sx q[1];
rz(-1.4349667) q[1];
sx q[1];
rz(-0.11134527) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24303007) q[0];
sx q[0];
rz(-0.70625171) q[0];
sx q[0];
rz(-2.7468922) q[0];
x q[1];
rz(-1.8315122) q[2];
sx q[2];
rz(-2.8088154) q[2];
sx q[2];
rz(-2.8944601) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7584472) q[1];
sx q[1];
rz(-1.0186334) q[1];
sx q[1];
rz(-2.8124185) q[1];
rz(-2.0599635) q[3];
sx q[3];
rz(-2.0738154) q[3];
sx q[3];
rz(-0.67219998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7340118) q[2];
sx q[2];
rz(-1.13052) q[2];
sx q[2];
rz(0.13047516) q[2];
rz(-0.16034165) q[3];
sx q[3];
rz(-2.4800143) q[3];
sx q[3];
rz(0.17016889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57123667) q[0];
sx q[0];
rz(-2.9354876) q[0];
sx q[0];
rz(-3.0233132) q[0];
rz(-2.2453399) q[1];
sx q[1];
rz(-1.4786485) q[1];
sx q[1];
rz(2.5884657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61643415) q[0];
sx q[0];
rz(-1.8983952) q[0];
sx q[0];
rz(-2.6382951) q[0];
x q[1];
rz(0.48001473) q[2];
sx q[2];
rz(-0.79593796) q[2];
sx q[2];
rz(0.98668879) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.432888) q[1];
sx q[1];
rz(-0.9917534) q[1];
sx q[1];
rz(-0.20198532) q[1];
x q[2];
rz(2.7378169) q[3];
sx q[3];
rz(-2.194768) q[3];
sx q[3];
rz(1.7221019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3502189) q[2];
sx q[2];
rz(-2.4198136) q[2];
sx q[2];
rz(-2.1283894) q[2];
rz(-0.31747097) q[3];
sx q[3];
rz(-1.443913) q[3];
sx q[3];
rz(2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1182227) q[0];
sx q[0];
rz(-0.88339266) q[0];
sx q[0];
rz(-0.50514847) q[0];
rz(0.39816868) q[1];
sx q[1];
rz(-1.265637) q[1];
sx q[1];
rz(1.8123951) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4666207) q[0];
sx q[0];
rz(-1.734991) q[0];
sx q[0];
rz(1.0999099) q[0];
rz(-pi) q[1];
rz(-0.5661613) q[2];
sx q[2];
rz(-2.5306411) q[2];
sx q[2];
rz(-0.47537714) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8201645) q[1];
sx q[1];
rz(-1.0099995) q[1];
sx q[1];
rz(0.1295156) q[1];
x q[2];
rz(-1.2990713) q[3];
sx q[3];
rz(-1.8193805) q[3];
sx q[3];
rz(-2.9523073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.087611) q[2];
sx q[2];
rz(-1.3302646) q[2];
sx q[2];
rz(2.3114253) q[2];
rz(-1.6603445) q[3];
sx q[3];
rz(-1.0426499) q[3];
sx q[3];
rz(1.164485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8692577) q[0];
sx q[0];
rz(-2.2255958) q[0];
sx q[0];
rz(-1.2087615) q[0];
rz(0.2441949) q[1];
sx q[1];
rz(-1.1714275) q[1];
sx q[1];
rz(0.29921439) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0335145) q[0];
sx q[0];
rz(-0.21768269) q[0];
sx q[0];
rz(1.449145) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0396378) q[2];
sx q[2];
rz(-2.1130145) q[2];
sx q[2];
rz(1.416666) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4296885) q[1];
sx q[1];
rz(-1.619408) q[1];
sx q[1];
rz(-0.91829388) q[1];
rz(-1.9817623) q[3];
sx q[3];
rz(-0.91051447) q[3];
sx q[3];
rz(3.0294047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19305912) q[2];
sx q[2];
rz(-1.7391354) q[2];
sx q[2];
rz(-2.194727) q[2];
rz(1.7638505) q[3];
sx q[3];
rz(-0.54072127) q[3];
sx q[3];
rz(0.64259678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4544025) q[0];
sx q[0];
rz(-0.43455046) q[0];
sx q[0];
rz(-0.53892556) q[0];
rz(-2.9680805) q[1];
sx q[1];
rz(-0.75650802) q[1];
sx q[1];
rz(2.7573746) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.656216) q[0];
sx q[0];
rz(-1.1278825) q[0];
sx q[0];
rz(-1.3963267) q[0];
rz(2.8649533) q[2];
sx q[2];
rz(-0.62046548) q[2];
sx q[2];
rz(-2.7850604) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7795664) q[1];
sx q[1];
rz(-3.0770731) q[1];
sx q[1];
rz(1.1532591) q[1];
x q[2];
rz(-0.027111091) q[3];
sx q[3];
rz(-2.1246006) q[3];
sx q[3];
rz(-1.1983271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3197799) q[2];
sx q[2];
rz(-0.76924789) q[2];
sx q[2];
rz(0.5963076) q[2];
rz(-2.1042306) q[3];
sx q[3];
rz(-2.3682902) q[3];
sx q[3];
rz(-1.1843225) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768815) q[0];
sx q[0];
rz(-2.3887971) q[0];
sx q[0];
rz(0.14213128) q[0];
rz(2.0243952) q[1];
sx q[1];
rz(-1.486472) q[1];
sx q[1];
rz(-1.2275009) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8637674) q[0];
sx q[0];
rz(-1.8466878) q[0];
sx q[0];
rz(0.49580144) q[0];
rz(-pi) q[1];
rz(-1.6864482) q[2];
sx q[2];
rz(-3.0031548) q[2];
sx q[2];
rz(1.9659245) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.043316226) q[1];
sx q[1];
rz(-2.195561) q[1];
sx q[1];
rz(-3.0618145) q[1];
rz(-pi) q[2];
rz(-1.0957509) q[3];
sx q[3];
rz(-0.79986279) q[3];
sx q[3];
rz(2.0020655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0907937) q[2];
sx q[2];
rz(-0.24238452) q[2];
sx q[2];
rz(0.79831115) q[2];
rz(1.5618886) q[3];
sx q[3];
rz(-1.7589933) q[3];
sx q[3];
rz(-0.17670512) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.056219) q[0];
sx q[0];
rz(-2.4985785) q[0];
sx q[0];
rz(3.1363078) q[0];
rz(-0.082848631) q[1];
sx q[1];
rz(-1.1033892) q[1];
sx q[1];
rz(2.8740035) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0326313) q[0];
sx q[0];
rz(-1.259045) q[0];
sx q[0];
rz(-0.70272081) q[0];
rz(-2.204699) q[2];
sx q[2];
rz(-0.46910252) q[2];
sx q[2];
rz(-1.3718951) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5983683) q[1];
sx q[1];
rz(-2.7562047) q[1];
sx q[1];
rz(2.1699127) q[1];
x q[2];
rz(-2.0844719) q[3];
sx q[3];
rz(-1.1506211) q[3];
sx q[3];
rz(-0.021180245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.573632) q[2];
sx q[2];
rz(-2.6601514) q[2];
sx q[2];
rz(-1.9285704) q[2];
rz(1.2447478) q[3];
sx q[3];
rz(-0.48143482) q[3];
sx q[3];
rz(1.1904233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5382814) q[0];
sx q[0];
rz(-2.17157) q[0];
sx q[0];
rz(-3.0441913) q[0];
rz(-3.1311323) q[1];
sx q[1];
rz(-2.2757826) q[1];
sx q[1];
rz(-0.59715685) q[1];
rz(-1.4821178) q[2];
sx q[2];
rz(-2.1046998) q[2];
sx q[2];
rz(1.7490341) q[2];
rz(-0.17531323) q[3];
sx q[3];
rz(-0.98876035) q[3];
sx q[3];
rz(0.50553346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
