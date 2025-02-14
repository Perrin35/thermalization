OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(2.0988965) q[0];
sx q[0];
rz(12.287696) q[0];
rz(4.2884045) q[1];
sx q[1];
rz(0.52462259) q[1];
sx q[1];
rz(8.1537032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3307453) q[0];
sx q[0];
rz(-2.8804104) q[0];
sx q[0];
rz(3.0708341) q[0];
x q[1];
rz(-2.991011) q[2];
sx q[2];
rz(-1.4408334) q[2];
sx q[2];
rz(0.58565631) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9624176) q[1];
sx q[1];
rz(-1.4915052) q[1];
sx q[1];
rz(0.85722629) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60802976) q[3];
sx q[3];
rz(-2.0634867) q[3];
sx q[3];
rz(1.0580491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0311387) q[2];
sx q[2];
rz(-1.9240856) q[2];
sx q[2];
rz(-0.26220599) q[2];
rz(2.0860705) q[3];
sx q[3];
rz(-1.8942098) q[3];
sx q[3];
rz(3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156886) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(-2.9600034) q[0];
rz(-1.4008105) q[1];
sx q[1];
rz(-1.1927651) q[1];
sx q[1];
rz(-2.3006732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4940482) q[0];
sx q[0];
rz(-1.2508014) q[0];
sx q[0];
rz(-1.5522237) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47535954) q[2];
sx q[2];
rz(-0.98612758) q[2];
sx q[2];
rz(-0.66051403) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5121756) q[1];
sx q[1];
rz(-1.4025153) q[1];
sx q[1];
rz(-3.1049033) q[1];
rz(-pi) q[2];
rz(1.0686458) q[3];
sx q[3];
rz(-2.4663894) q[3];
sx q[3];
rz(1.5218376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38222868) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(1.0648897) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.8062785) q[3];
sx q[3];
rz(1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8432124) q[0];
sx q[0];
rz(-1.1792264) q[0];
sx q[0];
rz(2.9425353) q[0];
rz(-2.3152323) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(2.3645649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0039486) q[0];
sx q[0];
rz(-2.9712518) q[0];
sx q[0];
rz(-1.9714703) q[0];
rz(-pi) q[1];
rz(-2.1146745) q[2];
sx q[2];
rz(-1.9444398) q[2];
sx q[2];
rz(-2.0622562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2232845) q[1];
sx q[1];
rz(-1.3121843) q[1];
sx q[1];
rz(2.1373939) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9701403) q[3];
sx q[3];
rz(-2.5094186) q[3];
sx q[3];
rz(-0.20242385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1021154) q[2];
sx q[2];
rz(-1.1481043) q[2];
sx q[2];
rz(2.6285505) q[2];
rz(-2.3593864) q[3];
sx q[3];
rz(-1.2057883) q[3];
sx q[3];
rz(-1.3768844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.9272598) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(1.4904892) q[0];
rz(0.53228846) q[1];
sx q[1];
rz(-2.5106301) q[1];
sx q[1];
rz(-1.5375563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69025264) q[0];
sx q[0];
rz(-1.4079388) q[0];
sx q[0];
rz(1.6785519) q[0];
rz(-pi) q[1];
rz(1.4179211) q[2];
sx q[2];
rz(-2.590153) q[2];
sx q[2];
rz(-1.3496292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5565479) q[1];
sx q[1];
rz(-2.2050207) q[1];
sx q[1];
rz(3.0536431) q[1];
rz(-0.38942899) q[3];
sx q[3];
rz(-2.5797322) q[3];
sx q[3];
rz(-1.6461247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6060467) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(-3.1287076) q[2];
rz(-1.9384725) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(-1.0217246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0032229) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(-1.3997929) q[0];
rz(2.7895582) q[1];
sx q[1];
rz(-1.0540009) q[1];
sx q[1];
rz(0.83784109) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8796757) q[0];
sx q[0];
rz(-0.81955541) q[0];
sx q[0];
rz(2.9330334) q[0];
x q[1];
rz(-0.95372405) q[2];
sx q[2];
rz(-1.0487818) q[2];
sx q[2];
rz(-1.2921289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.074118) q[1];
sx q[1];
rz(-2.5606541) q[1];
sx q[1];
rz(2.9089863) q[1];
x q[2];
rz(0.75148186) q[3];
sx q[3];
rz(-0.22322907) q[3];
sx q[3];
rz(3.1187559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0248523) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.1836729) q[2];
rz(2.8708894) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(-1.6089599) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026641332) q[0];
sx q[0];
rz(-2.5178435) q[0];
sx q[0];
rz(-0.57914105) q[0];
rz(-1.9193513) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(1.1096257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1117119) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(-2.9747308) q[0];
rz(-1.514362) q[2];
sx q[2];
rz(-1.2662958) q[2];
sx q[2];
rz(-0.79848189) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.077495726) q[1];
sx q[1];
rz(-1.239683) q[1];
sx q[1];
rz(1.3084133) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6251276) q[3];
sx q[3];
rz(-1.0871917) q[3];
sx q[3];
rz(-2.2276239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1960725) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(-2.784139) q[2];
rz(-1.1899905) q[3];
sx q[3];
rz(-1.2001195) q[3];
sx q[3];
rz(-1.801871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5424407) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(-0.79208148) q[0];
rz(2.4374841) q[1];
sx q[1];
rz(-2.6857565) q[1];
sx q[1];
rz(-1.5830931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1808727) q[0];
sx q[0];
rz(-1.5678493) q[0];
sx q[0];
rz(-2.6774027) q[0];
x q[1];
rz(2.1941845) q[2];
sx q[2];
rz(-1.0691423) q[2];
sx q[2];
rz(-0.19185184) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5643107) q[1];
sx q[1];
rz(-2.7140712) q[1];
sx q[1];
rz(0.083210817) q[1];
x q[2];
rz(-2.3901863) q[3];
sx q[3];
rz(-1.4045241) q[3];
sx q[3];
rz(-0.64054447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.030297) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(3.1374068) q[2];
rz(-2.8748416) q[3];
sx q[3];
rz(-1.5240074) q[3];
sx q[3];
rz(2.4193616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58077145) q[0];
sx q[0];
rz(-2.4120791) q[0];
sx q[0];
rz(-2.9855477) q[0];
rz(1.5566298) q[1];
sx q[1];
rz(-2.2792351) q[1];
sx q[1];
rz(2.5468266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42389171) q[0];
sx q[0];
rz(-1.0793975) q[0];
sx q[0];
rz(-1.2806847) q[0];
x q[1];
rz(2.4492599) q[2];
sx q[2];
rz(-1.344948) q[2];
sx q[2];
rz(-2.44997) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3799499) q[1];
sx q[1];
rz(-0.98230668) q[1];
sx q[1];
rz(2.2698927) q[1];
rz(-2.4202926) q[3];
sx q[3];
rz(-2.2709322) q[3];
sx q[3];
rz(2.457778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56738371) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(2.1136005) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-0.65516156) q[3];
sx q[3];
rz(1.3348234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780355) q[0];
sx q[0];
rz(-1.4565383) q[0];
sx q[0];
rz(-2.3433319) q[0];
rz(0.80052605) q[1];
sx q[1];
rz(-2.6526178) q[1];
sx q[1];
rz(0.54623234) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52926621) q[0];
sx q[0];
rz(-0.61175013) q[0];
sx q[0];
rz(-1.3806369) q[0];
rz(-pi) q[1];
rz(1.1425708) q[2];
sx q[2];
rz(-2.1167123) q[2];
sx q[2];
rz(-0.8949309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45193397) q[1];
sx q[1];
rz(-2.0469189) q[1];
sx q[1];
rz(2.1228288) q[1];
rz(-pi) q[2];
rz(0.30719884) q[3];
sx q[3];
rz(-0.23715487) q[3];
sx q[3];
rz(3.0320771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3937248) q[2];
sx q[2];
rz(-1.0900898) q[2];
sx q[2];
rz(-3.0785576) q[2];
rz(-2.4437599) q[3];
sx q[3];
rz(-2.883426) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49552712) q[0];
sx q[0];
rz(-0.87764469) q[0];
sx q[0];
rz(-2.7255507) q[0];
rz(1.3620194) q[1];
sx q[1];
rz(-1.1241309) q[1];
sx q[1];
rz(-1.8488041) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.775499) q[0];
sx q[0];
rz(-1.8387357) q[0];
sx q[0];
rz(1.3331205) q[0];
rz(0.45617683) q[2];
sx q[2];
rz(-1.3186698) q[2];
sx q[2];
rz(1.9664362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67329183) q[1];
sx q[1];
rz(-2.3929962) q[1];
sx q[1];
rz(0.8573726) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9927916) q[3];
sx q[3];
rz(-1.8099603) q[3];
sx q[3];
rz(1.07475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68299874) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(-0.90751737) q[2];
rz(1.2609743) q[3];
sx q[3];
rz(-2.5670299) q[3];
sx q[3];
rz(-1.6921836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4954729) q[0];
sx q[0];
rz(-1.6874122) q[0];
sx q[0];
rz(0.73943403) q[0];
rz(-0.44838913) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(-0.35198718) q[2];
sx q[2];
rz(-1.864973) q[2];
sx q[2];
rz(-0.61054338) q[2];
rz(-0.13335107) q[3];
sx q[3];
rz(-0.84682087) q[3];
sx q[3];
rz(-3.1414733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
