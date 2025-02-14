OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2321229) q[0];
sx q[0];
rz(-0.98804086) q[0];
sx q[0];
rz(-2.7890132) q[0];
rz(2.4352788) q[1];
sx q[1];
rz(-2.3568454) q[1];
sx q[1];
rz(3.1142601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4135983) q[0];
sx q[0];
rz(-1.9061273) q[0];
sx q[0];
rz(-1.3616614) q[0];
x q[1];
rz(0.17345239) q[2];
sx q[2];
rz(-1.7043742) q[2];
sx q[2];
rz(1.1277778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24522745) q[1];
sx q[1];
rz(-1.4030398) q[1];
sx q[1];
rz(-1.5390525) q[1];
x q[2];
rz(-2.5952847) q[3];
sx q[3];
rz(-0.27537333) q[3];
sx q[3];
rz(-0.34378036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0722384) q[2];
sx q[2];
rz(-1.3212997) q[2];
sx q[2];
rz(1.5864774) q[2];
rz(-0.72748264) q[3];
sx q[3];
rz(-2.1860217) q[3];
sx q[3];
rz(2.6684842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.469406) q[0];
sx q[0];
rz(-1.7962026) q[0];
sx q[0];
rz(-0.94959062) q[0];
rz(-2.8166215) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(-0.51345888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8425861) q[0];
sx q[0];
rz(-2.3761099) q[0];
sx q[0];
rz(2.1730636) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2113234) q[2];
sx q[2];
rz(-0.7604593) q[2];
sx q[2];
rz(2.9414842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9565485) q[1];
sx q[1];
rz(-2.0868851) q[1];
sx q[1];
rz(2.1538877) q[1];
rz(-pi) q[2];
rz(-2.6504032) q[3];
sx q[3];
rz(-2.0955387) q[3];
sx q[3];
rz(0.82651807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.010633858) q[2];
sx q[2];
rz(-1.8124688) q[2];
sx q[2];
rz(1.0884253) q[2];
rz(0.66793495) q[3];
sx q[3];
rz(-2.9974944) q[3];
sx q[3];
rz(-1.7264504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208991) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(1.2574842) q[0];
rz(-0.084818689) q[1];
sx q[1];
rz(-3.0501922) q[1];
sx q[1];
rz(0.65748293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2456928) q[0];
sx q[0];
rz(-2.6728711) q[0];
sx q[0];
rz(2.8492941) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59254139) q[2];
sx q[2];
rz(-2.0503902) q[2];
sx q[2];
rz(1.3035989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.3448101) q[1];
sx q[1];
rz(-1.0940897) q[1];
sx q[1];
rz(-1.3784798) q[1];
rz(-pi) q[2];
x q[2];
rz(2.231953) q[3];
sx q[3];
rz(-1.8170905) q[3];
sx q[3];
rz(1.1681739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76454863) q[2];
sx q[2];
rz(-0.42082861) q[2];
sx q[2];
rz(1.8910889) q[2];
rz(1.4416384) q[3];
sx q[3];
rz(-1.3911894) q[3];
sx q[3];
rz(0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68469754) q[0];
sx q[0];
rz(-2.1853515) q[0];
sx q[0];
rz(-0.60234219) q[0];
rz(0.96386987) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(2.9197555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5103127) q[0];
sx q[0];
rz(-1.5595734) q[0];
sx q[0];
rz(1.5785286) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8922563) q[2];
sx q[2];
rz(-2.0866924) q[2];
sx q[2];
rz(-0.14125401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1067729) q[1];
sx q[1];
rz(-1.9690138) q[1];
sx q[1];
rz(2.5335141) q[1];
x q[2];
rz(0.1476361) q[3];
sx q[3];
rz(-0.45429981) q[3];
sx q[3];
rz(-0.67377485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43231371) q[2];
sx q[2];
rz(-1.9807434) q[2];
sx q[2];
rz(0.088689001) q[2];
rz(2.6426219) q[3];
sx q[3];
rz(-1.6224529) q[3];
sx q[3];
rz(3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(2.2392739) q[0];
rz(0.29014507) q[1];
sx q[1];
rz(-2.2369308) q[1];
sx q[1];
rz(-1.1660928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6571044) q[0];
sx q[0];
rz(-2.1323626) q[0];
sx q[0];
rz(1.7620371) q[0];
rz(0.71641123) q[2];
sx q[2];
rz(-2.7802417) q[2];
sx q[2];
rz(-0.30308576) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2578967) q[1];
sx q[1];
rz(-1.9120815) q[1];
sx q[1];
rz(-0.30053707) q[1];
rz(1.5968634) q[3];
sx q[3];
rz(-1.3192303) q[3];
sx q[3];
rz(-0.1840011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12639283) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(-1.2882721) q[2];
rz(-2.1081693) q[3];
sx q[3];
rz(-2.0234334) q[3];
sx q[3];
rz(1.4748632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425659) q[0];
sx q[0];
rz(-3.0379744) q[0];
sx q[0];
rz(0.15956751) q[0];
rz(-0.60699925) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(-2.0807696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88638611) q[0];
sx q[0];
rz(-1.7632503) q[0];
sx q[0];
rz(2.1783955) q[0];
rz(2.4446746) q[2];
sx q[2];
rz(-0.3437416) q[2];
sx q[2];
rz(-0.70408193) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86058206) q[1];
sx q[1];
rz(-1.7063024) q[1];
sx q[1];
rz(-2.1675088) q[1];
rz(-pi) q[2];
rz(1.6688329) q[3];
sx q[3];
rz(-0.99409249) q[3];
sx q[3];
rz(-2.2876379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1242421) q[2];
sx q[2];
rz(-2.3257747) q[2];
sx q[2];
rz(-0.38786495) q[2];
rz(1.806949) q[3];
sx q[3];
rz(-0.67238656) q[3];
sx q[3];
rz(-1.0268802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(-0.89333308) q[0];
rz(0.54126254) q[1];
sx q[1];
rz(-2.2365139) q[1];
sx q[1];
rz(1.062324) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996595) q[0];
sx q[0];
rz(-0.77882732) q[0];
sx q[0];
rz(2.7657484) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6387894) q[2];
sx q[2];
rz(-2.3927702) q[2];
sx q[2];
rz(2.6523726) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7294638) q[1];
sx q[1];
rz(-2.1264014) q[1];
sx q[1];
rz(1.5985065) q[1];
rz(-pi) q[2];
rz(1.3618062) q[3];
sx q[3];
rz(-0.7721484) q[3];
sx q[3];
rz(-1.0812372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58747753) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(-1.973773) q[2];
rz(-0.81985146) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(2.4988373) q[3];
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
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.576936) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(1.8336953) q[0];
rz(-1.6851743) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(3.0110722) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042386656) q[0];
sx q[0];
rz(-1.644028) q[0];
sx q[0];
rz(1.8797148) q[0];
x q[1];
rz(-2.3442535) q[2];
sx q[2];
rz(-1.2166465) q[2];
sx q[2];
rz(-2.7018099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.83090342) q[1];
sx q[1];
rz(-2.2210418) q[1];
sx q[1];
rz(1.6199153) q[1];
x q[2];
rz(-2.2240488) q[3];
sx q[3];
rz(-0.96644478) q[3];
sx q[3];
rz(2.0487599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1859583) q[2];
sx q[2];
rz(-0.83883494) q[2];
sx q[2];
rz(1.7308509) q[2];
rz(-3.020982) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(-1.6284778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677583) q[0];
sx q[0];
rz(-0.61536106) q[0];
sx q[0];
rz(3.0009785) q[0];
rz(-2.4070542) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(-0.75070423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44245369) q[0];
sx q[0];
rz(-1.2780196) q[0];
sx q[0];
rz(-1.8163766) q[0];
x q[1];
rz(-0.98055697) q[2];
sx q[2];
rz(-0.41008355) q[2];
sx q[2];
rz(-0.29641253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0486601) q[1];
sx q[1];
rz(-1.6322929) q[1];
sx q[1];
rz(-1.3912806) q[1];
x q[2];
rz(-1.6569013) q[3];
sx q[3];
rz(-2.4331581) q[3];
sx q[3];
rz(-1.751386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.004403) q[2];
sx q[2];
rz(-0.69978324) q[2];
sx q[2];
rz(0.75263754) q[2];
rz(0.33958069) q[3];
sx q[3];
rz(-1.3656253) q[3];
sx q[3];
rz(0.021421758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.705377) q[0];
sx q[0];
rz(-0.74420539) q[0];
sx q[0];
rz(0.63802737) q[0];
rz(-0.72884196) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(2.3400838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214695) q[0];
sx q[0];
rz(-0.82372181) q[0];
sx q[0];
rz(2.864564) q[0];
x q[1];
rz(2.3036868) q[2];
sx q[2];
rz(-1.3896754) q[2];
sx q[2];
rz(1.6194026) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70511293) q[1];
sx q[1];
rz(-1.6098238) q[1];
sx q[1];
rz(0.83535933) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16149811) q[3];
sx q[3];
rz(-2.209389) q[3];
sx q[3];
rz(-2.0929008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30404299) q[2];
sx q[2];
rz(-1.7767228) q[2];
sx q[2];
rz(2.020828) q[2];
rz(-0.70901999) q[3];
sx q[3];
rz(-2.6778383) q[3];
sx q[3];
rz(-1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40502248) q[0];
sx q[0];
rz(-2.6800192) q[0];
sx q[0];
rz(-0.17929684) q[0];
rz(-2.2364521) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(0.84409406) q[2];
sx q[2];
rz(-1.8432968) q[2];
sx q[2];
rz(-2.7805614) q[2];
rz(1.0988476) q[3];
sx q[3];
rz(-0.57335486) q[3];
sx q[3];
rz(-2.2344979) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
