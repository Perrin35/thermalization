OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5689019) q[0];
sx q[0];
rz(-0.98839086) q[0];
sx q[0];
rz(2.6925777) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(0.13154498) q[1];
sx q[1];
rz(8.2934525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1695263) q[0];
sx q[0];
rz(-0.6837877) q[0];
sx q[0];
rz(-2.0440566) q[0];
x q[1];
rz(1.0523316) q[2];
sx q[2];
rz(-2.1396942) q[2];
sx q[2];
rz(2.2232747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3851812) q[1];
sx q[1];
rz(-2.1527421) q[1];
sx q[1];
rz(2.9333098) q[1];
rz(2.1142346) q[3];
sx q[3];
rz(-1.7187331) q[3];
sx q[3];
rz(-1.2191868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.955287) q[2];
sx q[2];
rz(-2.7489642) q[2];
sx q[2];
rz(-0.10360959) q[2];
rz(-0.3668395) q[3];
sx q[3];
rz(-1.47374) q[3];
sx q[3];
rz(1.3175255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1845301) q[0];
sx q[0];
rz(-0.4758895) q[0];
sx q[0];
rz(-2.6153508) q[0];
rz(-1.8980252) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(0.19613656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58914069) q[0];
sx q[0];
rz(-1.6375443) q[0];
sx q[0];
rz(-1.1132973) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19321014) q[2];
sx q[2];
rz(-0.71882283) q[2];
sx q[2];
rz(0.012398243) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6013697) q[1];
sx q[1];
rz(-1.1246846) q[1];
sx q[1];
rz(-2.9838954) q[1];
rz(-pi) q[2];
rz(-1.5603746) q[3];
sx q[3];
rz(-1.5589514) q[3];
sx q[3];
rz(2.7238977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5474995) q[2];
sx q[2];
rz(-0.27442351) q[2];
sx q[2];
rz(-0.37626949) q[2];
rz(-2.2606692) q[3];
sx q[3];
rz(-1.3353525) q[3];
sx q[3];
rz(2.3295565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45102099) q[0];
sx q[0];
rz(-0.32830992) q[0];
sx q[0];
rz(2.1300533) q[0];
rz(-0.66863376) q[1];
sx q[1];
rz(-1.9790383) q[1];
sx q[1];
rz(2.5943601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5315192) q[0];
sx q[0];
rz(-1.3694057) q[0];
sx q[0];
rz(1.5904443) q[0];
rz(-pi) q[1];
x q[1];
rz(2.793712) q[2];
sx q[2];
rz(-0.82150412) q[2];
sx q[2];
rz(-2.2583928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32729748) q[1];
sx q[1];
rz(-1.8117827) q[1];
sx q[1];
rz(0.90118932) q[1];
rz(-2.9487398) q[3];
sx q[3];
rz(-0.82369419) q[3];
sx q[3];
rz(-2.3806662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4554567) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(1.6395462) q[2];
rz(0.57602588) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(-2.7801133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(-1.5107002) q[0];
sx q[0];
rz(-0.80131131) q[0];
sx q[0];
rz(2.8019688) q[0];
rz(-0.54061186) q[1];
sx q[1];
rz(-2.4410591) q[1];
sx q[1];
rz(-0.9300173) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1620063) q[0];
sx q[0];
rz(-1.7135059) q[0];
sx q[0];
rz(3.0703785) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9564863) q[2];
sx q[2];
rz(-2.1197332) q[2];
sx q[2];
rz(-1.576129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.31615694) q[1];
sx q[1];
rz(-0.47940578) q[1];
sx q[1];
rz(1.5569219) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92842438) q[3];
sx q[3];
rz(-1.1517715) q[3];
sx q[3];
rz(-0.85538543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0129619) q[2];
sx q[2];
rz(-1.4350812) q[2];
sx q[2];
rz(-1.5677412) q[2];
rz(-2.3980127) q[3];
sx q[3];
rz(-2.3502374) q[3];
sx q[3];
rz(0.96405205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4670694) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(3.0849482) q[0];
rz(1.4757587) q[1];
sx q[1];
rz(-1.1328127) q[1];
sx q[1];
rz(-1.7830361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4641842) q[0];
sx q[0];
rz(-1.0631936) q[0];
sx q[0];
rz(-2.1305231) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65746324) q[2];
sx q[2];
rz(-1.7320447) q[2];
sx q[2];
rz(2.8828414) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2142519) q[1];
sx q[1];
rz(-2.3433629) q[1];
sx q[1];
rz(2.4386474) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.452946) q[3];
sx q[3];
rz(-2.1292344) q[3];
sx q[3];
rz(1.3056361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8654827) q[2];
sx q[2];
rz(-1.411974) q[2];
sx q[2];
rz(1.126368) q[2];
rz(-1.3364835) q[3];
sx q[3];
rz(-1.0765321) q[3];
sx q[3];
rz(1.0895464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4514076) q[0];
sx q[0];
rz(-2.1162338) q[0];
sx q[0];
rz(-3.0080646) q[0];
rz(2.1639157) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(-2.8547063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8650019) q[0];
sx q[0];
rz(-1.6185068) q[0];
sx q[0];
rz(2.2722831) q[0];
x q[1];
rz(-0.043043625) q[2];
sx q[2];
rz(-1.898694) q[2];
sx q[2];
rz(-2.4509933) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9911954) q[1];
sx q[1];
rz(-2.5270473) q[1];
sx q[1];
rz(-1.1769562) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5941723) q[3];
sx q[3];
rz(-1.35619) q[3];
sx q[3];
rz(1.2494189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4213244) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.6555017) q[2];
rz(-2.7739575) q[3];
sx q[3];
rz(-1.0272107) q[3];
sx q[3];
rz(-1.0953085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7169645) q[0];
sx q[0];
rz(-0.76911887) q[0];
sx q[0];
rz(0.47384438) q[0];
rz(-0.66420707) q[1];
sx q[1];
rz(-1.9240446) q[1];
sx q[1];
rz(1.7582105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63603386) q[0];
sx q[0];
rz(-2.2657388) q[0];
sx q[0];
rz(-0.64993422) q[0];
rz(-pi) q[1];
rz(-0.02772133) q[2];
sx q[2];
rz(-1.9887191) q[2];
sx q[2];
rz(2.6270285) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5210628) q[1];
sx q[1];
rz(-2.8179666) q[1];
sx q[1];
rz(-0.77038295) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3455799) q[3];
sx q[3];
rz(-1.7449505) q[3];
sx q[3];
rz(-0.98291493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15709269) q[2];
sx q[2];
rz(-1.4089156) q[2];
sx q[2];
rz(-2.816693) q[2];
rz(-0.38069185) q[3];
sx q[3];
rz(-0.84563962) q[3];
sx q[3];
rz(2.0438173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070504524) q[0];
sx q[0];
rz(-0.66362137) q[0];
sx q[0];
rz(0.93884236) q[0];
rz(-2.4414869) q[1];
sx q[1];
rz(-0.63799262) q[1];
sx q[1];
rz(0.28900388) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28075179) q[0];
sx q[0];
rz(-1.7993225) q[0];
sx q[0];
rz(3.0620652) q[0];
x q[1];
rz(-2.2174066) q[2];
sx q[2];
rz(-1.980154) q[2];
sx q[2];
rz(1.6381303) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49456143) q[1];
sx q[1];
rz(-0.53003487) q[1];
sx q[1];
rz(-1.3840335) q[1];
x q[2];
rz(2.9634776) q[3];
sx q[3];
rz(-2.7550972) q[3];
sx q[3];
rz(2.9535171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2603904) q[2];
sx q[2];
rz(-1.6198879) q[2];
sx q[2];
rz(2.728906) q[2];
rz(-0.9203426) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(-1.2296366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5145787) q[0];
sx q[0];
rz(-0.98015061) q[0];
sx q[0];
rz(-2.8421616) q[0];
rz(2.0191655) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(1.0312414) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61187295) q[0];
sx q[0];
rz(-2.7752989) q[0];
sx q[0];
rz(-1.6008911) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3716566) q[2];
sx q[2];
rz(-2.0250892) q[2];
sx q[2];
rz(2.9261738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9523176) q[1];
sx q[1];
rz(-0.77042033) q[1];
sx q[1];
rz(-1.1483869) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84800546) q[3];
sx q[3];
rz(-2.3591745) q[3];
sx q[3];
rz(0.54718024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9791744) q[2];
sx q[2];
rz(-1.7981497) q[2];
sx q[2];
rz(-2.8988885) q[2];
rz(-1.8414712) q[3];
sx q[3];
rz(-1.0588812) q[3];
sx q[3];
rz(-0.21568957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31371394) q[0];
sx q[0];
rz(-0.21610459) q[0];
sx q[0];
rz(-1.4208273) q[0];
rz(-2.8315262) q[1];
sx q[1];
rz(-0.87691751) q[1];
sx q[1];
rz(1.5975331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9102073) q[0];
sx q[0];
rz(-0.69150309) q[0];
sx q[0];
rz(-2.4639936) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3256959) q[2];
sx q[2];
rz(-0.35159207) q[2];
sx q[2];
rz(-2.0577672) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3812755) q[1];
sx q[1];
rz(-2.7997428) q[1];
sx q[1];
rz(-0.6980957) q[1];
rz(1.710911) q[3];
sx q[3];
rz(-1.7286674) q[3];
sx q[3];
rz(-2.637459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1342423) q[2];
sx q[2];
rz(-1.4126567) q[2];
sx q[2];
rz(1.3664112) q[2];
rz(0.8775231) q[3];
sx q[3];
rz(-1.6759422) q[3];
sx q[3];
rz(2.2519978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14961814) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(2.2772475) q[1];
sx q[1];
rz(-2.1280011) q[1];
sx q[1];
rz(-0.43307532) q[1];
rz(1.7766118) q[2];
sx q[2];
rz(-1.0792776) q[2];
sx q[2];
rz(0.57224032) q[2];
rz(-1.6771562) q[3];
sx q[3];
rz(-0.48135664) q[3];
sx q[3];
rz(2.4873747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
