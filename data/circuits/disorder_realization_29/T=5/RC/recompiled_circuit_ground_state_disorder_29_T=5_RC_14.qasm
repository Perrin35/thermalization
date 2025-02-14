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
rz(1.1468118) q[1];
sx q[1];
rz(-0.52462259) q[1];
sx q[1];
rz(1.8705179) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88408121) q[0];
sx q[0];
rz(-1.8313098) q[0];
sx q[0];
rz(-1.5519014) q[0];
rz(0.71670436) q[2];
sx q[2];
rz(-2.9430046) q[2];
sx q[2];
rz(-2.8633397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.179175) q[1];
sx q[1];
rz(-1.4915052) q[1];
sx q[1];
rz(-2.2843664) q[1];
rz(-0.60802976) q[3];
sx q[3];
rz(-1.078106) q[3];
sx q[3];
rz(-2.0835435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0311387) q[2];
sx q[2];
rz(-1.9240856) q[2];
sx q[2];
rz(2.8793867) q[2];
rz(-1.0555222) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(-3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5259041) q[0];
sx q[0];
rz(-0.50351244) q[0];
sx q[0];
rz(-0.18158922) q[0];
rz(1.4008105) q[1];
sx q[1];
rz(-1.1927651) q[1];
sx q[1];
rz(-0.84091944) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6475444) q[0];
sx q[0];
rz(-1.2508014) q[0];
sx q[0];
rz(-1.5522237) q[0];
rz(-pi) q[1];
rz(-2.1757751) q[2];
sx q[2];
rz(-0.73558319) q[2];
sx q[2];
rz(-0.090026131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4136728) q[1];
sx q[1];
rz(-0.17219725) q[1];
sx q[1];
rz(1.7834458) q[1];
rz(2.1827994) q[3];
sx q[3];
rz(-1.8763767) q[3];
sx q[3];
rz(-2.6877054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38222868) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(-2.0767029) q[2];
rz(-3.0801638) q[3];
sx q[3];
rz(-1.3353142) q[3];
sx q[3];
rz(1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8432124) q[0];
sx q[0];
rz(-1.1792264) q[0];
sx q[0];
rz(0.19905736) q[0];
rz(2.3152323) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(-2.3645649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26829241) q[0];
sx q[0];
rz(-1.7275294) q[0];
sx q[0];
rz(0.066989338) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42962809) q[2];
sx q[2];
rz(-1.0681392) q[2];
sx q[2];
rz(-2.4328897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18622929) q[1];
sx q[1];
rz(-2.1163773) q[1];
sx q[1];
rz(2.8377691) q[1];
rz(-1.1714524) q[3];
sx q[3];
rz(-0.63217406) q[3];
sx q[3];
rz(-0.20242385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0394773) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(0.51304212) q[2];
rz(0.78220621) q[3];
sx q[3];
rz(-1.2057883) q[3];
sx q[3];
rz(1.7647083) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9272598) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(-1.4904892) q[0];
rz(-0.53228846) q[1];
sx q[1];
rz(-2.5106301) q[1];
sx q[1];
rz(-1.6040364) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.45134) q[0];
sx q[0];
rz(-1.4079388) q[0];
sx q[0];
rz(1.6785519) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1170148) q[2];
sx q[2];
rz(-1.6506631) q[2];
sx q[2];
rz(0.090674222) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73276943) q[1];
sx q[1];
rz(-0.63946001) q[1];
sx q[1];
rz(1.4519522) q[1];
rz(-1.8054078) q[3];
sx q[3];
rz(-2.0861833) q[3];
sx q[3];
rz(-1.0439411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.535546) q[2];
sx q[2];
rz(-2.594785) q[2];
sx q[2];
rz(-0.01288506) q[2];
rz(-1.2031201) q[3];
sx q[3];
rz(-1.67484) q[3];
sx q[3];
rz(2.119868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0032229) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(-1.7417997) q[0];
rz(2.7895582) q[1];
sx q[1];
rz(-1.0540009) q[1];
sx q[1];
rz(0.83784109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45233881) q[0];
sx q[0];
rz(-1.4188915) q[0];
sx q[0];
rz(0.80861966) q[0];
rz(-pi) q[1];
rz(2.3532392) q[2];
sx q[2];
rz(-0.78561312) q[2];
sx q[2];
rz(-2.2503302) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6988236) q[1];
sx q[1];
rz(-1.4439481) q[1];
sx q[1];
rz(2.5731099) q[1];
rz(1.4170333) q[3];
sx q[3];
rz(-1.7332675) q[3];
sx q[3];
rz(-2.3547309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1167404) q[2];
sx q[2];
rz(-2.1768673) q[2];
sx q[2];
rz(1.1836729) q[2];
rz(-2.8708894) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(1.6089599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1149513) q[0];
sx q[0];
rz(-0.6237492) q[0];
sx q[0];
rz(-0.57914105) q[0];
rz(1.2222414) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(1.1096257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1117119) q[0];
sx q[0];
rz(-1.4193901) q[0];
sx q[0];
rz(-0.16686186) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30495651) q[2];
sx q[2];
rz(-1.5169608) q[2];
sx q[2];
rz(2.3862145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0640969) q[1];
sx q[1];
rz(-1.9019097) q[1];
sx q[1];
rz(-1.3084133) q[1];
rz(-pi) q[2];
rz(-0.81619451) q[3];
sx q[3];
rz(-2.4494054) q[3];
sx q[3];
rz(0.029267197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1960725) q[2];
sx q[2];
rz(-1.8438945) q[2];
sx q[2];
rz(-0.35745364) q[2];
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
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59915197) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(0.79208148) q[0];
rz(-2.4374841) q[1];
sx q[1];
rz(-0.45583615) q[1];
sx q[1];
rz(-1.5830931) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9607199) q[0];
sx q[0];
rz(-1.5737434) q[0];
sx q[0];
rz(0.46418993) q[0];
x q[1];
rz(-2.1941845) q[2];
sx q[2];
rz(-1.0691423) q[2];
sx q[2];
rz(-2.9497408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4729135) q[1];
sx q[1];
rz(-1.9967419) q[1];
sx q[1];
rz(1.6086474) q[1];
rz(0.75140636) q[3];
sx q[3];
rz(-1.7370686) q[3];
sx q[3];
rz(0.64054447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11129561) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(0.0041858717) q[2];
rz(2.8748416) q[3];
sx q[3];
rz(-1.6175852) q[3];
sx q[3];
rz(2.4193616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58077145) q[0];
sx q[0];
rz(-2.4120791) q[0];
sx q[0];
rz(0.15604493) q[0];
rz(-1.5566298) q[1];
sx q[1];
rz(-2.2792351) q[1];
sx q[1];
rz(0.59476605) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1536008) q[0];
sx q[0];
rz(-0.56454851) q[0];
sx q[0];
rz(-2.6507244) q[0];
rz(-pi) q[1];
rz(2.4492599) q[2];
sx q[2];
rz(-1.7966447) q[2];
sx q[2];
rz(-0.69162265) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2458098) q[1];
sx q[1];
rz(-1.0059662) q[1];
sx q[1];
rz(0.71706949) q[1];
rz(-0.72130008) q[3];
sx q[3];
rz(-0.87066046) q[3];
sx q[3];
rz(-0.68381468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56738371) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(1.0279921) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-2.4864311) q[3];
sx q[3];
rz(1.8067693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8780355) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(-2.3433319) q[0];
rz(-0.80052605) q[1];
sx q[1];
rz(-2.6526178) q[1];
sx q[1];
rz(-0.54623234) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29833405) q[0];
sx q[0];
rz(-2.1699561) q[0];
sx q[0];
rz(-0.13183044) q[0];
rz(-1.9990218) q[2];
sx q[2];
rz(-2.1167123) q[2];
sx q[2];
rz(2.2466618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.662488) q[1];
sx q[1];
rz(-2.4291385) q[1];
sx q[1];
rz(-2.3478048) q[1];
x q[2];
rz(-0.30719884) q[3];
sx q[3];
rz(-0.23715487) q[3];
sx q[3];
rz(0.10951558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3937248) q[2];
sx q[2];
rz(-2.0515029) q[2];
sx q[2];
rz(-0.063035034) q[2];
rz(-0.6978327) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460655) q[0];
sx q[0];
rz(-2.263948) q[0];
sx q[0];
rz(-2.7255507) q[0];
rz(-1.3620194) q[1];
sx q[1];
rz(-2.0174618) q[1];
sx q[1];
rz(-1.8488041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0344311) q[0];
sx q[0];
rz(-2.7853372) q[0];
sx q[0];
rz(0.70888575) q[0];
rz(-0.52915574) q[2];
sx q[2];
rz(-2.6247027) q[2];
sx q[2];
rz(0.86597499) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-2.0668427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68299874) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(-0.90751737) q[2];
rz(-1.8806184) q[3];
sx q[3];
rz(-2.5670299) q[3];
sx q[3];
rz(1.449409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64611971) q[0];
sx q[0];
rz(-1.4541805) q[0];
sx q[0];
rz(-2.4021586) q[0];
rz(-0.44838913) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(-1.8829968) q[2];
sx q[2];
rz(-1.2345424) q[2];
sx q[2];
rz(0.85415863) q[2];
rz(-3.0082416) q[3];
sx q[3];
rz(-2.2947718) q[3];
sx q[3];
rz(0.00011937192) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
