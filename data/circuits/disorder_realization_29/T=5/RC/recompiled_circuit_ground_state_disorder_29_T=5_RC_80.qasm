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
rz(-1.0426961) q[0];
sx q[0];
rz(0.27867499) q[0];
rz(1.1468118) q[1];
sx q[1];
rz(-0.52462259) q[1];
sx q[1];
rz(1.8705179) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88408121) q[0];
sx q[0];
rz(-1.8313098) q[0];
sx q[0];
rz(-1.5896912) q[0];
rz(-pi) q[1];
rz(2.991011) q[2];
sx q[2];
rz(-1.7007593) q[2];
sx q[2];
rz(0.58565631) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48285741) q[1];
sx q[1];
rz(-2.4244011) q[1];
sx q[1];
rz(-1.6916005) q[1];
rz(0.99156143) q[3];
sx q[3];
rz(-1.0433726) q[3];
sx q[3];
rz(-2.9468731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0311387) q[2];
sx q[2];
rz(-1.217507) q[2];
sx q[2];
rz(-2.8793867) q[2];
rz(-2.0860705) q[3];
sx q[3];
rz(-1.8942098) q[3];
sx q[3];
rz(-3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5259041) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(0.18158922) q[0];
rz(-1.7407821) q[1];
sx q[1];
rz(-1.9488275) q[1];
sx q[1];
rz(-2.3006732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5530295) q[0];
sx q[0];
rz(-0.32051495) q[0];
sx q[0];
rz(-0.055983981) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6662331) q[2];
sx q[2];
rz(-2.1554651) q[2];
sx q[2];
rz(-0.66051403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.629417) q[1];
sx q[1];
rz(-1.7390774) q[1];
sx q[1];
rz(-0.03668935) q[1];
rz(-pi) q[2];
rz(-2.1827994) q[3];
sx q[3];
rz(-1.8763767) q[3];
sx q[3];
rz(-0.45388729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.759364) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(-2.0767029) q[2];
rz(3.0801638) q[3];
sx q[3];
rz(-1.3353142) q[3];
sx q[3];
rz(-1.6399062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29838022) q[0];
sx q[0];
rz(-1.9623663) q[0];
sx q[0];
rz(0.19905736) q[0];
rz(2.3152323) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(0.77702776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0039486) q[0];
sx q[0];
rz(-0.17034082) q[0];
sx q[0];
rz(-1.1701223) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92240693) q[2];
sx q[2];
rz(-0.64903478) q[2];
sx q[2];
rz(-0.051608406) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2232845) q[1];
sx q[1];
rz(-1.3121843) q[1];
sx q[1];
rz(-2.1373939) q[1];
rz(-2.1644209) q[3];
sx q[3];
rz(-1.3389753) q[3];
sx q[3];
rz(2.1013732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0394773) q[2];
sx q[2];
rz(-1.1481043) q[2];
sx q[2];
rz(0.51304212) q[2];
rz(2.3593864) q[3];
sx q[3];
rz(-1.9358044) q[3];
sx q[3];
rz(-1.3768844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9272598) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(-1.6511035) q[0];
rz(-2.6093042) q[1];
sx q[1];
rz(-0.63096255) q[1];
sx q[1];
rz(-1.6040364) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10189942) q[0];
sx q[0];
rz(-2.9465775) q[0];
sx q[0];
rz(-0.57955091) q[0];
rz(0.093393383) q[2];
sx q[2];
rz(-2.1150781) q[2];
sx q[2];
rz(-1.5285847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4088232) q[1];
sx q[1];
rz(-0.63946001) q[1];
sx q[1];
rz(1.6896405) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52738366) q[3];
sx q[3];
rz(-1.77447) q[3];
sx q[3];
rz(0.40959293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6060467) q[2];
sx q[2];
rz(-0.54680768) q[2];
sx q[2];
rz(3.1287076) q[2];
rz(-1.2031201) q[3];
sx q[3];
rz(-1.67484) q[3];
sx q[3];
rz(-1.0217246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0032229) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(1.7417997) q[0];
rz(0.35203448) q[1];
sx q[1];
rz(-2.0875918) q[1];
sx q[1];
rz(-2.3037516) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96124803) q[0];
sx q[0];
rz(-2.367428) q[0];
sx q[0];
rz(1.3526239) q[0];
x q[1];
rz(2.1878686) q[2];
sx q[2];
rz(-1.0487818) q[2];
sx q[2];
rz(1.8494637) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.067474631) q[1];
sx q[1];
rz(-2.5606541) q[1];
sx q[1];
rz(2.9089863) q[1];
rz(-pi) q[2];
rz(0.16437634) q[3];
sx q[3];
rz(-1.419074) q[3];
sx q[3];
rz(0.80899948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0248523) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.9579197) q[2];
rz(-0.27070326) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(1.5326327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026641332) q[0];
sx q[0];
rz(-0.6237492) q[0];
sx q[0];
rz(2.5624516) q[0];
rz(1.2222414) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(-2.0319669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56631382) q[0];
sx q[0];
rz(-1.405861) q[0];
sx q[0];
rz(-1.4172907) q[0];
rz(-pi) q[1];
rz(2.8366361) q[2];
sx q[2];
rz(-1.6246319) q[2];
sx q[2];
rz(-0.75537813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61298227) q[1];
sx q[1];
rz(-0.41944606) q[1];
sx q[1];
rz(-0.64639133) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0274572) q[3];
sx q[3];
rz(-2.0232588) q[3];
sx q[3];
rz(-2.2266091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1960725) q[2];
sx q[2];
rz(-1.8438945) q[2];
sx q[2];
rz(-0.35745364) q[2];
rz(-1.9516021) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(1.3397217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59915197) q[0];
sx q[0];
rz(-2.1245133) q[0];
sx q[0];
rz(2.3495112) q[0];
rz(2.4374841) q[1];
sx q[1];
rz(-2.6857565) q[1];
sx q[1];
rz(1.5584996) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38403758) q[0];
sx q[0];
rz(-2.677394) q[0];
sx q[0];
rz(3.1350101) q[0];
rz(-pi) q[1];
rz(-0.81659813) q[2];
sx q[2];
rz(-0.77864051) q[2];
sx q[2];
rz(-1.1731847) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4729135) q[1];
sx q[1];
rz(-1.9967419) q[1];
sx q[1];
rz(1.6086474) q[1];
rz(1.3450479) q[3];
sx q[3];
rz(-2.3093947) q[3];
sx q[3];
rz(2.3647472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.11129561) q[2];
sx q[2];
rz(-0.78812495) q[2];
sx q[2];
rz(-0.0041858717) q[2];
rz(2.8748416) q[3];
sx q[3];
rz(-1.5240074) q[3];
sx q[3];
rz(-2.4193616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5608212) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(-0.15604493) q[0];
rz(-1.5849628) q[1];
sx q[1];
rz(-0.86235756) q[1];
sx q[1];
rz(-2.5468266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0069591) q[0];
sx q[0];
rz(-1.3158321) q[0];
sx q[0];
rz(2.6322271) q[0];
x q[1];
rz(1.2807219) q[2];
sx q[2];
rz(-2.2422487) q[2];
sx q[2];
rz(-0.69556505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76164272) q[1];
sx q[1];
rz(-2.159286) q[1];
sx q[1];
rz(-0.87169991) q[1];
rz(-0.72799369) q[3];
sx q[3];
rz(-1.0413975) q[3];
sx q[3];
rz(-1.7391143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5742089) q[2];
sx q[2];
rz(-1.5790066) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780355) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(2.3433319) q[0];
rz(-0.80052605) q[1];
sx q[1];
rz(-0.48897484) q[1];
sx q[1];
rz(0.54623234) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8432586) q[0];
sx q[0];
rz(-2.1699561) q[0];
sx q[0];
rz(-0.13183044) q[0];
rz(1.9990218) q[2];
sx q[2];
rz(-2.1167123) q[2];
sx q[2];
rz(-2.2466618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45193397) q[1];
sx q[1];
rz(-1.0946737) q[1];
sx q[1];
rz(2.1228288) q[1];
rz(-pi) q[2];
rz(-1.4978374) q[3];
sx q[3];
rz(-1.7966509) q[3];
sx q[3];
rz(-0.42499229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.74786782) q[2];
sx q[2];
rz(-1.0900898) q[2];
sx q[2];
rz(-0.063035034) q[2];
rz(0.6978327) q[3];
sx q[3];
rz(-2.883426) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(1.2927885) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1071615) q[0];
sx q[0];
rz(-0.35625544) q[0];
sx q[0];
rz(-0.70888575) q[0];
x q[1];
rz(2.6854158) q[2];
sx q[2];
rz(-1.3186698) q[2];
sx q[2];
rz(-1.9664362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4683008) q[1];
sx q[1];
rz(-2.3929962) q[1];
sx q[1];
rz(0.8573726) q[1];
rz(-pi) q[2];
rz(1.9927916) q[3];
sx q[3];
rz(-1.3316324) q[3];
sx q[3];
rz(-2.0668427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68299874) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(-1.8806184) q[3];
sx q[3];
rz(-0.57456273) q[3];
sx q[3];
rz(-1.449409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64611971) q[0];
sx q[0];
rz(-1.4541805) q[0];
sx q[0];
rz(-2.4021586) q[0];
rz(2.6932035) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(2.7896055) q[2];
sx q[2];
rz(-1.864973) q[2];
sx q[2];
rz(-0.61054338) q[2];
rz(3.0082416) q[3];
sx q[3];
rz(-0.84682087) q[3];
sx q[3];
rz(-3.1414733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
