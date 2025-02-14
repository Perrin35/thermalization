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
rz(0.24231237) q[0];
sx q[0];
rz(4.0153541) q[0];
sx q[0];
rz(10.913475) q[0];
rz(1.0701264) q[1];
sx q[1];
rz(-2.5843599) q[1];
sx q[1];
rz(-0.4692404) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0996286) q[0];
sx q[0];
rz(-1.6828186) q[0];
sx q[0];
rz(-0.69667453) q[0];
rz(-pi) q[1];
rz(2.8618545) q[2];
sx q[2];
rz(-1.7953897) q[2];
sx q[2];
rz(-0.43817929) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9201913) q[1];
sx q[1];
rz(-1.2556352) q[1];
sx q[1];
rz(1.0454476) q[1];
rz(-pi) q[2];
rz(-2.4116573) q[3];
sx q[3];
rz(-1.1570108) q[3];
sx q[3];
rz(-1.034193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0327518) q[2];
sx q[2];
rz(-2.8665172) q[2];
sx q[2];
rz(-1.3230422) q[2];
rz(-0.9465341) q[3];
sx q[3];
rz(-2.5338379) q[3];
sx q[3];
rz(0.058145903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7846947) q[0];
sx q[0];
rz(-0.31814831) q[0];
sx q[0];
rz(1.963105) q[0];
rz(-2.5937041) q[1];
sx q[1];
rz(-1.9375487) q[1];
sx q[1];
rz(1.9080124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89900333) q[0];
sx q[0];
rz(-1.543643) q[0];
sx q[0];
rz(-1.8544556) q[0];
rz(0.080670653) q[2];
sx q[2];
rz(-1.4051132) q[2];
sx q[2];
rz(2.5625429) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4405304) q[1];
sx q[1];
rz(-0.077210285) q[1];
sx q[1];
rz(-2.3838897) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61290755) q[3];
sx q[3];
rz(-0.41348413) q[3];
sx q[3];
rz(0.72847937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2445688) q[2];
sx q[2];
rz(-2.343101) q[2];
sx q[2];
rz(1.7986253) q[2];
rz(-3.1161984) q[3];
sx q[3];
rz(-1.3586905) q[3];
sx q[3];
rz(0.079518147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19567604) q[0];
sx q[0];
rz(-1.8657277) q[0];
sx q[0];
rz(0.50862408) q[0];
rz(1.7018082) q[1];
sx q[1];
rz(-1.516781) q[1];
sx q[1];
rz(1.5135099) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1726097) q[0];
sx q[0];
rz(-1.727761) q[0];
sx q[0];
rz(2.5353372) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2824487) q[2];
sx q[2];
rz(-2.9133688) q[2];
sx q[2];
rz(1.3068401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73838961) q[1];
sx q[1];
rz(-1.566267) q[1];
sx q[1];
rz(-0.036935115) q[1];
rz(-pi) q[2];
rz(-2.3044488) q[3];
sx q[3];
rz(-0.52609936) q[3];
sx q[3];
rz(0.21080454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.14579183) q[2];
sx q[2];
rz(-2.3172816) q[2];
sx q[2];
rz(-0.15131797) q[2];
rz(0.96210903) q[3];
sx q[3];
rz(-1.5113219) q[3];
sx q[3];
rz(-2.7675653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91604084) q[0];
sx q[0];
rz(-2.2489838) q[0];
sx q[0];
rz(1.4501866) q[0];
rz(-2.080503) q[1];
sx q[1];
rz(-0.42409813) q[1];
sx q[1];
rz(1.1294956) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2573525) q[0];
sx q[0];
rz(-1.7654618) q[0];
sx q[0];
rz(-2.3759222) q[0];
rz(-pi) q[1];
rz(0.53217066) q[2];
sx q[2];
rz(-1.9111562) q[2];
sx q[2];
rz(0.81102358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68431932) q[1];
sx q[1];
rz(-2.4114747) q[1];
sx q[1];
rz(1.3616379) q[1];
rz(1.2673504) q[3];
sx q[3];
rz(-0.11622322) q[3];
sx q[3];
rz(0.35467142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9811161) q[2];
sx q[2];
rz(-0.99157292) q[2];
sx q[2];
rz(1.1264832) q[2];
rz(2.9142761) q[3];
sx q[3];
rz(-1.7132297) q[3];
sx q[3];
rz(-3.0912002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8166167) q[0];
sx q[0];
rz(-0.80642527) q[0];
sx q[0];
rz(0.6977914) q[0];
rz(-1.6119488) q[1];
sx q[1];
rz(-0.35747129) q[1];
sx q[1];
rz(2.7059817) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2078081) q[0];
sx q[0];
rz(-1.6949886) q[0];
sx q[0];
rz(-1.3072312) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0514652) q[2];
sx q[2];
rz(-0.70106912) q[2];
sx q[2];
rz(2.4034479) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10098329) q[1];
sx q[1];
rz(-0.48982778) q[1];
sx q[1];
rz(0.71160729) q[1];
rz(-pi) q[2];
rz(0.27542476) q[3];
sx q[3];
rz(-1.6425624) q[3];
sx q[3];
rz(0.83930086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7464253) q[2];
sx q[2];
rz(-0.94234157) q[2];
sx q[2];
rz(-2.49474) q[2];
rz(2.4344889) q[3];
sx q[3];
rz(-3.0349858) q[3];
sx q[3];
rz(-2.148237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3498722) q[0];
sx q[0];
rz(-2.7339022) q[0];
sx q[0];
rz(-0.83734018) q[0];
rz(0.95036858) q[1];
sx q[1];
rz(-1.252754) q[1];
sx q[1];
rz(0.18384917) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1342574) q[0];
sx q[0];
rz(-2.7401972) q[0];
sx q[0];
rz(2.4069277) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9230336) q[2];
sx q[2];
rz(-1.248385) q[2];
sx q[2];
rz(-0.36568322) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.819471) q[1];
sx q[1];
rz(-2.2985704) q[1];
sx q[1];
rz(2.6202109) q[1];
rz(-2.9086439) q[3];
sx q[3];
rz(-1.4714575) q[3];
sx q[3];
rz(2.8976822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6763372) q[2];
sx q[2];
rz(-2.4691171) q[2];
sx q[2];
rz(-0.23784168) q[2];
rz(3.0658718) q[3];
sx q[3];
rz(-3.0140311) q[3];
sx q[3];
rz(0.41847509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7521562) q[0];
sx q[0];
rz(-0.51920033) q[0];
sx q[0];
rz(-0.12553781) q[0];
rz(0.43352661) q[1];
sx q[1];
rz(-2.5018689) q[1];
sx q[1];
rz(1.7009521) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5880116) q[0];
sx q[0];
rz(-1.5097678) q[0];
sx q[0];
rz(-0.048179772) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7309181) q[2];
sx q[2];
rz(-2.0537964) q[2];
sx q[2];
rz(-1.4087189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.10769612) q[1];
sx q[1];
rz(-1.0748861) q[1];
sx q[1];
rz(-1.1065844) q[1];
rz(-pi) q[2];
rz(-0.85236211) q[3];
sx q[3];
rz(-0.96489513) q[3];
sx q[3];
rz(-0.62921274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0354249) q[2];
sx q[2];
rz(-1.1664392) q[2];
sx q[2];
rz(-2.9336145) q[2];
rz(2.4584127) q[3];
sx q[3];
rz(-2.4222789) q[3];
sx q[3];
rz(-1.4815559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55027562) q[0];
sx q[0];
rz(-2.2324201) q[0];
sx q[0];
rz(0.43347484) q[0];
rz(2.6705961) q[1];
sx q[1];
rz(-2.8484671) q[1];
sx q[1];
rz(-0.34117821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6745673) q[0];
sx q[0];
rz(-1.4816947) q[0];
sx q[0];
rz(-3.0568987) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3586832) q[2];
sx q[2];
rz(-2.3723432) q[2];
sx q[2];
rz(1.5616901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2318445) q[1];
sx q[1];
rz(-2.1655786) q[1];
sx q[1];
rz(-2.9308795) q[1];
rz(-pi) q[2];
rz(1.4580995) q[3];
sx q[3];
rz(-0.98033614) q[3];
sx q[3];
rz(-2.0463508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53182536) q[2];
sx q[2];
rz(-2.7957081) q[2];
sx q[2];
rz(-2.7609265) q[2];
rz(-0.62764132) q[3];
sx q[3];
rz(-0.52730477) q[3];
sx q[3];
rz(-0.92831534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8220383) q[0];
sx q[0];
rz(-1.1975937) q[0];
sx q[0];
rz(-2.8354216) q[0];
rz(-0.51433688) q[1];
sx q[1];
rz(-0.75175935) q[1];
sx q[1];
rz(-0.18246442) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771181) q[0];
sx q[0];
rz(-1.6892642) q[0];
sx q[0];
rz(0.76513357) q[0];
rz(-pi) q[1];
rz(2.6128255) q[2];
sx q[2];
rz(-1.6739539) q[2];
sx q[2];
rz(-2.8512252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.565158) q[1];
sx q[1];
rz(-1.7243313) q[1];
sx q[1];
rz(1.9839194) q[1];
rz(-0.45643198) q[3];
sx q[3];
rz(-1.8805537) q[3];
sx q[3];
rz(2.8995958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3719486) q[2];
sx q[2];
rz(-0.19522788) q[2];
sx q[2];
rz(0.39607421) q[2];
rz(-1.1560446) q[3];
sx q[3];
rz(-0.58737415) q[3];
sx q[3];
rz(0.42266947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7523772) q[0];
sx q[0];
rz(-0.14425819) q[0];
sx q[0];
rz(-0.17917646) q[0];
rz(1.3009118) q[1];
sx q[1];
rz(-0.4549028) q[1];
sx q[1];
rz(2.6120344) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0146739) q[0];
sx q[0];
rz(-1.8243941) q[0];
sx q[0];
rz(-2.7932667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0069557351) q[2];
sx q[2];
rz(-1.2901301) q[2];
sx q[2];
rz(2.8988738) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5645819) q[1];
sx q[1];
rz(-1.2839926) q[1];
sx q[1];
rz(-0.69820349) q[1];
rz(-pi) q[2];
rz(1.3934025) q[3];
sx q[3];
rz(-2.1004026) q[3];
sx q[3];
rz(-1.7152928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35310465) q[2];
sx q[2];
rz(-0.64211923) q[2];
sx q[2];
rz(-2.8118964) q[2];
rz(1.6126136) q[3];
sx q[3];
rz(-1.8776882) q[3];
sx q[3];
rz(2.8211856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0223087) q[0];
sx q[0];
rz(-1.6884463) q[0];
sx q[0];
rz(1.6519188) q[0];
rz(-2.655401) q[1];
sx q[1];
rz(-1.9959027) q[1];
sx q[1];
rz(-2.1015658) q[1];
rz(2.9943071) q[2];
sx q[2];
rz(-0.68193692) q[2];
sx q[2];
rz(3.0030967) q[2];
rz(-2.8066951) q[3];
sx q[3];
rz(-2.6819102) q[3];
sx q[3];
rz(0.93536207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
