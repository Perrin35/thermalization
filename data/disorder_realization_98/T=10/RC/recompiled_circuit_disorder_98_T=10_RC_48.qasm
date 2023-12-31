OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-0.85435003) q[0];
sx q[0];
rz(-2.2358405) q[0];
x q[1];
rz(1.0295463) q[2];
sx q[2];
rz(-1.1272578) q[2];
sx q[2];
rz(2.8461547) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61923164) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(-3.1027604) q[1];
rz(-pi) q[2];
rz(1.2655067) q[3];
sx q[3];
rz(-2.5698834) q[3];
sx q[3];
rz(1.3523462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5271877) q[0];
sx q[0];
rz(-1.5831278) q[0];
sx q[0];
rz(3.1129818) q[0];
rz(-pi) q[1];
rz(-2.2488238) q[2];
sx q[2];
rz(-0.36913482) q[2];
sx q[2];
rz(2.6004651) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3995041) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(-2.8469267) q[1];
rz(-pi) q[2];
rz(0.60450508) q[3];
sx q[3];
rz(-2.1406056) q[3];
sx q[3];
rz(1.2165716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(-1.3495548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2230167) q[0];
sx q[0];
rz(-2.4239459) q[0];
sx q[0];
rz(0.36803228) q[0];
x q[1];
rz(-3.1042728) q[2];
sx q[2];
rz(-1.63675) q[2];
sx q[2];
rz(2.0685591) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.459356) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(-0.45046803) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1447103) q[3];
sx q[3];
rz(-1.5618556) q[3];
sx q[3];
rz(-2.5090891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(0.30113014) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.4105463) q[0];
rz(-0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(3.1052123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43544337) q[0];
sx q[0];
rz(-2.2009146) q[0];
sx q[0];
rz(2.8773727) q[0];
rz(-pi) q[1];
rz(1.245001) q[2];
sx q[2];
rz(-1.0137644) q[2];
sx q[2];
rz(1.4404802) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2991043) q[1];
sx q[1];
rz(-2.4310388) q[1];
sx q[1];
rz(-1.7732265) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23457228) q[3];
sx q[3];
rz(-1.6664701) q[3];
sx q[3];
rz(-2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.9533763) q[2];
rz(2.4711117) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-0.23434815) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7974632) q[0];
sx q[0];
rz(-1.9088233) q[0];
sx q[0];
rz(-2.5812838) q[0];
x q[1];
rz(0.095586153) q[2];
sx q[2];
rz(-2.0628953) q[2];
sx q[2];
rz(1.0620067) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.105141) q[1];
sx q[1];
rz(-2.3804133) q[1];
sx q[1];
rz(0.21703227) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.078951051) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(-2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(2.9210572) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(0.81714001) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(-1.9979427) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691539) q[0];
sx q[0];
rz(-1.5307431) q[0];
sx q[0];
rz(-2.7705631) q[0];
rz(-pi) q[1];
rz(-2.8191889) q[2];
sx q[2];
rz(-2.44256) q[2];
sx q[2];
rz(-2.2088745) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5551344) q[1];
sx q[1];
rz(-2.0124334) q[1];
sx q[1];
rz(-3.0416136) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46623047) q[3];
sx q[3];
rz(-2.5737408) q[3];
sx q[3];
rz(0.0044435244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.75366655) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.8547159) q[0];
rz(-0.6634179) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.9082665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1713283) q[0];
sx q[0];
rz(-1.5818705) q[0];
sx q[0];
rz(1.1867255) q[0];
x q[1];
rz(0.6255409) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(2.2369838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4349605) q[1];
sx q[1];
rz(-0.11206493) q[1];
sx q[1];
rz(-0.12607615) q[1];
rz(-pi) q[2];
rz(-0.4864278) q[3];
sx q[3];
rz(-1.8842116) q[3];
sx q[3];
rz(3.0048971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(1.0664553) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.6960779) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(1.258237) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2381666) q[0];
sx q[0];
rz(-1.4362207) q[0];
sx q[0];
rz(2.0057136) q[0];
rz(-pi) q[1];
rz(2.4969205) q[2];
sx q[2];
rz(-2.1053227) q[2];
sx q[2];
rz(-1.0169741) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80553493) q[1];
sx q[1];
rz(-1.4667257) q[1];
sx q[1];
rz(-0.54688262) q[1];
x q[2];
rz(-0.52950852) q[3];
sx q[3];
rz(-0.92748517) q[3];
sx q[3];
rz(2.1852126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(-2.941926) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(-0.02773157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6513034) q[0];
sx q[0];
rz(-0.9285183) q[0];
sx q[0];
rz(2.519033) q[0];
rz(-pi) q[1];
rz(-0.26720033) q[2];
sx q[2];
rz(-0.88951096) q[2];
sx q[2];
rz(0.53182488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2003277) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(-1.4501249) q[1];
x q[2];
rz(3.1157007) q[3];
sx q[3];
rz(-1.6373487) q[3];
sx q[3];
rz(-0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(-0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(2.2843602) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(2.6877158) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-2.8840816) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2430902) q[0];
sx q[0];
rz(-1.1872963) q[0];
sx q[0];
rz(2.6570508) q[0];
rz(-pi) q[1];
rz(0.72980482) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(-0.44621106) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35044893) q[1];
sx q[1];
rz(-2.2669683) q[1];
sx q[1];
rz(1.5894366) q[1];
rz(-pi) q[2];
rz(-1.0912283) q[3];
sx q[3];
rz(-1.7864831) q[3];
sx q[3];
rz(-1.8857764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3474779) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-0.22656245) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(0.46802855) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(1.0874891) q[3];
sx q[3];
rz(-1.8046422) q[3];
sx q[3];
rz(1.3840152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
