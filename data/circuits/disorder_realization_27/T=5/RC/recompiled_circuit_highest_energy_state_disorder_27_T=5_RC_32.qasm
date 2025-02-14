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
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1013452) q[0];
sx q[0];
rz(-1.57124) q[0];
sx q[0];
rz(1.5612649) q[0];
rz(-2.462596) q[2];
sx q[2];
rz(-1.8583991) q[2];
sx q[2];
rz(-1.2305413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9551785) q[1];
sx q[1];
rz(-2.4354773) q[1];
sx q[1];
rz(-0.43118711) q[1];
x q[2];
rz(-2.3232949) q[3];
sx q[3];
rz(-2.1967271) q[3];
sx q[3];
rz(1.675954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8523031) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-3.0618073) q[2];
rz(0.96528178) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(-0.7307581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7331284) q[0];
sx q[0];
rz(-0.95703546) q[0];
sx q[0];
rz(1.2888786) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0816108) q[2];
sx q[2];
rz(-0.85027611) q[2];
sx q[2];
rz(2.1921076) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8364077) q[1];
sx q[1];
rz(-0.47023222) q[1];
sx q[1];
rz(1.2364619) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9761451) q[3];
sx q[3];
rz(-2.5962127) q[3];
sx q[3];
rz(-0.65217962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.264512) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(-1.4667) q[2];
rz(-3.1125715) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90917176) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(-1.4311283) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(2.7412282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5059197) q[0];
sx q[0];
rz(-0.95209661) q[0];
sx q[0];
rz(1.0158422) q[0];
rz(0.68636559) q[2];
sx q[2];
rz(-1.7698506) q[2];
sx q[2];
rz(0.65716568) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0644933) q[1];
sx q[1];
rz(-1.1005741) q[1];
sx q[1];
rz(-1.749586) q[1];
rz(-0.34137643) q[3];
sx q[3];
rz(-1.2814643) q[3];
sx q[3];
rz(-0.067726243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4572738) q[2];
sx q[2];
rz(-1.6088586) q[2];
sx q[2];
rz(0.45480967) q[2];
rz(1.2416035) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95551816) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(3.1410826) q[0];
rz(-2.5406802) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(-0.14437637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52989456) q[0];
sx q[0];
rz(-1.7774044) q[0];
sx q[0];
rz(-2.1312461) q[0];
x q[1];
rz(-2.8881489) q[2];
sx q[2];
rz(-1.4254693) q[2];
sx q[2];
rz(0.0090473024) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2541009) q[1];
sx q[1];
rz(-0.93587592) q[1];
sx q[1];
rz(1.2870406) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4636643) q[3];
sx q[3];
rz(-2.1153573) q[3];
sx q[3];
rz(0.0060826172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.942975) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(2.1790738) q[2];
rz(1.5396384) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(-0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.9050423) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(1.1849674) q[0];
rz(2.9153337) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(2.102898) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.561956) q[0];
sx q[0];
rz(-2.547894) q[0];
sx q[0];
rz(-2.9953792) q[0];
rz(-pi) q[1];
rz(0.94064297) q[2];
sx q[2];
rz(-2.2046208) q[2];
sx q[2];
rz(1.9535106) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8983137) q[1];
sx q[1];
rz(-0.89363511) q[1];
sx q[1];
rz(0.89666287) q[1];
rz(-pi) q[2];
x q[2];
rz(1.709278) q[3];
sx q[3];
rz(-2.6977728) q[3];
sx q[3];
rz(2.6276772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7096536) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(-1.4656674) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50452152) q[0];
sx q[0];
rz(-2.2383454) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(0.48031131) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(0.59741098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99973122) q[0];
sx q[0];
rz(-1.7787361) q[0];
sx q[0];
rz(0.19241649) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3961847) q[2];
sx q[2];
rz(-1.7495724) q[2];
sx q[2];
rz(2.488236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.589041) q[1];
sx q[1];
rz(-1.8218166) q[1];
sx q[1];
rz(-1.4000721) q[1];
x q[2];
rz(0.091986309) q[3];
sx q[3];
rz(-1.9402656) q[3];
sx q[3];
rz(1.2928158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-2.0651979) q[2];
rz(-1.9988029) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(-0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550605) q[0];
sx q[0];
rz(-1.9831816) q[0];
sx q[0];
rz(0.038473815) q[0];
rz(0.064727457) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-0.23385349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42002871) q[0];
sx q[0];
rz(-1.7917243) q[0];
sx q[0];
rz(-1.7608282) q[0];
x q[1];
rz(1.2042768) q[2];
sx q[2];
rz(-1.2784174) q[2];
sx q[2];
rz(-2.2077843) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.95215248) q[1];
sx q[1];
rz(-2.2169211) q[1];
sx q[1];
rz(1.6829856) q[1];
rz(-0.41681186) q[3];
sx q[3];
rz(-2.9997065) q[3];
sx q[3];
rz(2.8519423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48876277) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(0.59824198) q[2];
rz(-0.11387842) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(2.2940476) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050196) q[0];
sx q[0];
rz(-2.2990655) q[0];
sx q[0];
rz(-2.8952428) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-1.1187226) q[1];
sx q[1];
rz(0.49682239) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637359) q[0];
sx q[0];
rz(-1.1564768) q[0];
sx q[0];
rz(2.4606649) q[0];
rz(2.970817) q[2];
sx q[2];
rz(-2.429212) q[2];
sx q[2];
rz(0.26584372) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0441974) q[1];
sx q[1];
rz(-0.5438416) q[1];
sx q[1];
rz(2.0378276) q[1];
rz(-1.2874576) q[3];
sx q[3];
rz(-2.6322798) q[3];
sx q[3];
rz(2.4278502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7184427) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(1.144484) q[2];
rz(1.1540958) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(0.19788876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941187) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(0.06614729) q[0];
rz(1.1478395) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(-2.537421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8291289) q[0];
sx q[0];
rz(-1.365571) q[0];
sx q[0];
rz(-0.28136307) q[0];
rz(-pi) q[1];
rz(-2.6148817) q[2];
sx q[2];
rz(-1.7065431) q[2];
sx q[2];
rz(-1.2413687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0688404) q[1];
sx q[1];
rz(-1.3786331) q[1];
sx q[1];
rz(-0.56984624) q[1];
rz(-pi) q[2];
rz(-0.87852134) q[3];
sx q[3];
rz(-1.0348579) q[3];
sx q[3];
rz(1.5978447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26312795) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(-0.43295941) q[2];
rz(-1.912502) q[3];
sx q[3];
rz(-1.2058328) q[3];
sx q[3];
rz(1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94806725) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(1.2731592) q[0];
rz(0.46514568) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(-2.926631) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5453542) q[0];
sx q[0];
rz(-2.4430877) q[0];
sx q[0];
rz(-2.7424988) q[0];
rz(-pi) q[1];
rz(-0.5318435) q[2];
sx q[2];
rz(-1.9491247) q[2];
sx q[2];
rz(1.7500306) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7219639) q[1];
sx q[1];
rz(-2.2441494) q[1];
sx q[1];
rz(-0.90760214) q[1];
rz(-pi) q[2];
rz(-0.39542012) q[3];
sx q[3];
rz(-0.59874615) q[3];
sx q[3];
rz(1.0637763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2994069) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(-0.95883933) q[2];
rz(1.1622608) q[3];
sx q[3];
rz(-1.6543038) q[3];
sx q[3];
rz(-1.4452665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(2.4330347) q[1];
sx q[1];
rz(-0.19150145) q[1];
sx q[1];
rz(2.3932744) q[1];
rz(-1.2461927) q[2];
sx q[2];
rz(-1.0321675) q[2];
sx q[2];
rz(2.2553222) q[2];
rz(1.7583354) q[3];
sx q[3];
rz(-2.2839727) q[3];
sx q[3];
rz(-1.4061389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
