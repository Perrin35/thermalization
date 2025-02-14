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
rz(1.2394387) q[0];
sx q[0];
rz(4.4702915) q[0];
sx q[0];
rz(9.2128949) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(-0.53242004) q[1];
sx q[1];
rz(-0.37791696) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1013452) q[0];
sx q[0];
rz(-1.57124) q[0];
sx q[0];
rz(1.5803278) q[0];
x q[1];
rz(2.462596) q[2];
sx q[2];
rz(-1.2831935) q[2];
sx q[2];
rz(1.9110514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4204882) q[1];
sx q[1];
rz(-1.2961565) q[1];
sx q[1];
rz(-2.4824449) q[1];
rz(-pi) q[2];
rz(-2.3232949) q[3];
sx q[3];
rz(-0.94486559) q[3];
sx q[3];
rz(-1.675954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8523031) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-3.0618073) q[2];
rz(2.1763109) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6642283) q[0];
sx q[0];
rz(-2.1381162) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(-1.7695919) q[1];
sx q[1];
rz(-1.7141432) q[1];
sx q[1];
rz(1.1044097) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7331284) q[0];
sx q[0];
rz(-0.95703546) q[0];
sx q[0];
rz(-1.852714) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49141617) q[2];
sx q[2];
rz(-2.2960536) q[2];
sx q[2];
rz(-2.8710136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4648393) q[1];
sx q[1];
rz(-1.1285121) q[1];
sx q[1];
rz(-0.16525903) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4711963) q[3];
sx q[3];
rz(-2.1079113) q[3];
sx q[3];
rz(-0.84505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.264512) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(1.4667) q[2];
rz(-3.1125715) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-0.066815289) q[0];
sx q[0];
rz(2.8636041) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(2.7412282) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5059197) q[0];
sx q[0];
rz(-0.95209661) q[0];
sx q[0];
rz(2.1257504) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3157) q[2];
sx q[2];
rz(-2.2410903) q[2];
sx q[2];
rz(0.75302659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4246042) q[1];
sx q[1];
rz(-1.7300055) q[1];
sx q[1];
rz(-0.47674322) q[1];
x q[2];
rz(-0.72680803) q[3];
sx q[3];
rz(-2.6977959) q[3];
sx q[3];
rz(2.3151195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(-2.686783) q[2];
rz(-1.2416035) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95551816) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(0.00051001471) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(0.14437637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3568319) q[0];
sx q[0];
rz(-0.5934754) q[0];
sx q[0];
rz(-1.1952101) q[0];
x q[1];
rz(-1.4207441) q[2];
sx q[2];
rz(-1.320082) q[2];
sx q[2];
rz(-1.5423519) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2541009) q[1];
sx q[1];
rz(-0.93587592) q[1];
sx q[1];
rz(-1.8545521) q[1];
rz(-0.17474971) q[3];
sx q[3];
rz(-2.5876382) q[3];
sx q[3];
rz(0.21077158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1986177) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(2.1790738) q[2];
rz(-1.6019542) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(-2.8008154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2365504) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(1.1849674) q[0];
rz(-0.22625893) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(-1.0386946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7554692) q[0];
sx q[0];
rz(-2.1573108) q[0];
sx q[0];
rz(1.4727794) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4034385) q[2];
sx q[2];
rz(-2.0655491) q[2];
sx q[2];
rz(3.11657) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2782368) q[1];
sx q[1];
rz(-2.0788631) q[1];
sx q[1];
rz(0.79973508) q[1];
rz(-pi) q[2];
rz(-0.065537621) q[3];
sx q[3];
rz(-1.1315232) q[3];
sx q[3];
rz(-2.7807943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7096536) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(-0.31614885) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50452152) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(-2.0380518) q[0];
rz(2.6612813) q[1];
sx q[1];
rz(-2.4791398) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1418614) q[0];
sx q[0];
rz(-1.7787361) q[0];
sx q[0];
rz(2.9491762) q[0];
x q[1];
rz(-0.74540794) q[2];
sx q[2];
rz(-1.3920203) q[2];
sx q[2];
rz(-2.488236) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.054185523) q[1];
sx q[1];
rz(-0.30255908) q[1];
sx q[1];
rz(2.5564479) q[1];
x q[2];
rz(1.8036929) q[3];
sx q[3];
rz(-0.38023708) q[3];
sx q[3];
rz(-1.5986795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20208134) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550605) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(-0.038473815) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-2.9077392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.192894) q[0];
sx q[0];
rz(-1.7561551) q[0];
sx q[0];
rz(-2.9167487) q[0];
x q[1];
rz(1.9373158) q[2];
sx q[2];
rz(-1.2784174) q[2];
sx q[2];
rz(2.2077843) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1894402) q[1];
sx q[1];
rz(-0.92467151) q[1];
sx q[1];
rz(-1.6829856) q[1];
rz(-pi) q[2];
rz(0.41681186) q[3];
sx q[3];
rz(-0.14188611) q[3];
sx q[3];
rz(2.8519423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6528299) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-0.59824198) q[2];
rz(3.0277142) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(-0.84754506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(0.2463499) q[0];
rz(-1.2975289) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(-2.6447703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778567) q[0];
sx q[0];
rz(-1.9851159) q[0];
sx q[0];
rz(-2.4606649) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4250579) q[2];
sx q[2];
rz(-2.2706804) q[2];
sx q[2];
rz(3.0997955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0441974) q[1];
sx q[1];
rz(-0.5438416) q[1];
sx q[1];
rz(-1.103765) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15487352) q[3];
sx q[3];
rz(-1.0836156) q[3];
sx q[3];
rz(0.39184141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.42314998) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(1.9971087) q[2];
rz(1.9874969) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(2.9437039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941187) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(-3.0754454) q[0];
rz(1.1478395) q[1];
sx q[1];
rz(-1.3775974) q[1];
sx q[1];
rz(-0.60417169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9420931) q[0];
sx q[0];
rz(-1.2954933) q[0];
sx q[0];
rz(-1.784174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6148817) q[2];
sx q[2];
rz(-1.4350495) q[2];
sx q[2];
rz(-1.2413687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7614038) q[1];
sx q[1];
rz(-2.1288925) q[1];
sx q[1];
rz(1.3437043) q[1];
x q[2];
rz(2.2630713) q[3];
sx q[3];
rz(-2.1067348) q[3];
sx q[3];
rz(1.5437479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8784647) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(2.7086332) q[2];
rz(1.912502) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.94806725) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(-1.8684335) q[0];
rz(2.676447) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(2.926631) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0419768) q[0];
sx q[0];
rz(-0.93659217) q[0];
sx q[0];
rz(-1.8862104) q[0];
rz(2.0029055) q[2];
sx q[2];
rz(-2.0615163) q[2];
sx q[2];
rz(-0.034772074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7219639) q[1];
sx q[1];
rz(-2.2441494) q[1];
sx q[1];
rz(0.90760214) q[1];
rz(-1.8278024) q[3];
sx q[3];
rz(-1.0238092) q[3];
sx q[3];
rz(-1.5316602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2994069) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(-0.95883933) q[2];
rz(1.9793319) q[3];
sx q[3];
rz(-1.6543038) q[3];
sx q[3];
rz(-1.6963262) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6017629) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(-0.70855793) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(-2.5790527) q[2];
sx q[2];
rz(-1.8481135) q[2];
sx q[2];
rz(-2.6279966) q[2];
rz(-1.3832573) q[3];
sx q[3];
rz(-2.2839727) q[3];
sx q[3];
rz(-1.4061389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
