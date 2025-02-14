OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8102326) q[0];
sx q[0];
rz(4.3963764) q[0];
sx q[0];
rz(8.7788361) q[0];
rz(-1.4589925) q[1];
sx q[1];
rz(-0.12812935) q[1];
sx q[1];
rz(2.473414) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9409882) q[0];
sx q[0];
rz(-1.683569) q[0];
sx q[0];
rz(-2.8975945) q[0];
x q[1];
rz(-1.1485841) q[2];
sx q[2];
rz(-2.4496884) q[2];
sx q[2];
rz(-0.9949323) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9757802) q[1];
sx q[1];
rz(-0.67953142) q[1];
sx q[1];
rz(-2.2305016) q[1];
rz(2.6553538) q[3];
sx q[3];
rz(-2.231153) q[3];
sx q[3];
rz(0.59729353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6231923) q[2];
sx q[2];
rz(-0.65541583) q[2];
sx q[2];
rz(-2.0584959) q[2];
rz(-2.8862503) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(1.9820836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1622247) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(3.1080143) q[0];
rz(2.0032739) q[1];
sx q[1];
rz(-2.1470224) q[1];
sx q[1];
rz(2.9046362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7743999) q[0];
sx q[0];
rz(-0.26846186) q[0];
sx q[0];
rz(-0.546683) q[0];
rz(-pi) q[1];
rz(1.0944738) q[2];
sx q[2];
rz(-1.4038205) q[2];
sx q[2];
rz(-2.4284438) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8023876) q[1];
sx q[1];
rz(-1.5962114) q[1];
sx q[1];
rz(-1.3783546) q[1];
x q[2];
rz(-2.6941766) q[3];
sx q[3];
rz(-1.0378812) q[3];
sx q[3];
rz(-0.79870236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.529155) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(-3.0212413) q[2];
rz(2.555661) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(1.8732635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.878207) q[0];
sx q[0];
rz(-2.1757941) q[0];
sx q[0];
rz(2.1534488) q[0];
rz(1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(-1.5536701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2399167) q[0];
sx q[0];
rz(-1.5696948) q[0];
sx q[0];
rz(3.1256846) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6695197) q[2];
sx q[2];
rz(-0.1506714) q[2];
sx q[2];
rz(-0.27418384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.37582477) q[1];
sx q[1];
rz(-1.8319523) q[1];
sx q[1];
rz(-3.1118666) q[1];
x q[2];
rz(2.5738945) q[3];
sx q[3];
rz(-1.2097219) q[3];
sx q[3];
rz(-0.89601024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3042018) q[2];
sx q[2];
rz(-1.3434255) q[2];
sx q[2];
rz(-3.0710132) q[2];
rz(-2.6523759) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(-2.9941471) q[3];
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
rz(-1.1864784) q[0];
sx q[0];
rz(-2.4461353) q[0];
sx q[0];
rz(1.7406933) q[0];
rz(-0.1700302) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(1.9084575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9233386) q[0];
sx q[0];
rz(-1.6955351) q[0];
sx q[0];
rz(-1.9943004) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98922421) q[2];
sx q[2];
rz(-2.2007341) q[2];
sx q[2];
rz(-0.0058071227) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.68076949) q[1];
sx q[1];
rz(-1.555838) q[1];
sx q[1];
rz(-1.2911002) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1175851) q[3];
sx q[3];
rz(-2.6769961) q[3];
sx q[3];
rz(-0.28991227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79973811) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(-2.5119761) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-2.5119731) q[3];
sx q[3];
rz(1.7262044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28302309) q[0];
sx q[0];
rz(-2.6333599) q[0];
sx q[0];
rz(-0.10230219) q[0];
rz(1.6434068) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(2.7844875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67960167) q[0];
sx q[0];
rz(-1.7850189) q[0];
sx q[0];
rz(-2.9132183) q[0];
x q[1];
rz(3.018385) q[2];
sx q[2];
rz(-1.5977309) q[2];
sx q[2];
rz(-1.7719442) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60254492) q[1];
sx q[1];
rz(-2.5395406) q[1];
sx q[1];
rz(0.54564387) q[1];
x q[2];
rz(2.6980347) q[3];
sx q[3];
rz(-1.1263811) q[3];
sx q[3];
rz(-2.9068974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47306481) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.8190039) q[2];
rz(-1.708301) q[3];
sx q[3];
rz(-1.3168443) q[3];
sx q[3];
rz(-2.9776261) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54356164) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(1.7033956) q[0];
rz(1.2403129) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(1.7887438) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.529218) q[0];
sx q[0];
rz(-1.1657526) q[0];
sx q[0];
rz(-0.10251001) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2850432) q[2];
sx q[2];
rz(-2.5973598) q[2];
sx q[2];
rz(-2.2533992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6042418) q[1];
sx q[1];
rz(-2.4927995) q[1];
sx q[1];
rz(2.2082445) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8753136) q[3];
sx q[3];
rz(-2.7334573) q[3];
sx q[3];
rz(0.88312393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.214434) q[2];
sx q[2];
rz(-2.7775601) q[2];
sx q[2];
rz(2.7511609) q[2];
rz(-1.3123784) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(-0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696233) q[0];
sx q[0];
rz(-2.2177028) q[0];
sx q[0];
rz(1.4129289) q[0];
rz(1.4631118) q[1];
sx q[1];
rz(-1.5556346) q[1];
sx q[1];
rz(-2.5020592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93776751) q[0];
sx q[0];
rz(-0.97903189) q[0];
sx q[0];
rz(-2.6861565) q[0];
rz(-pi) q[1];
rz(-0.15070559) q[2];
sx q[2];
rz(-1.4960714) q[2];
sx q[2];
rz(0.13241235) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29296103) q[1];
sx q[1];
rz(-2.1994414) q[1];
sx q[1];
rz(-1.9002923) q[1];
rz(-pi) q[2];
rz(-1.1763379) q[3];
sx q[3];
rz(-1.6812075) q[3];
sx q[3];
rz(-2.4019456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9245236) q[2];
sx q[2];
rz(-1.3288493) q[2];
sx q[2];
rz(-0.45219839) q[2];
rz(2.9186115) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(-2.9862459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13101354) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(0.16192326) q[0];
rz(2.7936392) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(2.9071992) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97353712) q[0];
sx q[0];
rz(-0.99322766) q[0];
sx q[0];
rz(-2.4650991) q[0];
rz(-pi) q[1];
rz(-2.4215464) q[2];
sx q[2];
rz(-1.3338425) q[2];
sx q[2];
rz(-2.4958378) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3752611) q[1];
sx q[1];
rz(-2.2041781) q[1];
sx q[1];
rz(-0.23289012) q[1];
rz(0.23080821) q[3];
sx q[3];
rz(-0.2838906) q[3];
sx q[3];
rz(-2.6687255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8552385) q[2];
sx q[2];
rz(-1.0711292) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(-2.3983119) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(-0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8648935) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(1.0869166) q[0];
rz(-1.2326321) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(1.8392275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83643276) q[0];
sx q[0];
rz(-1.6839241) q[0];
sx q[0];
rz(-2.9311084) q[0];
x q[1];
rz(2.1673498) q[2];
sx q[2];
rz(-1.8678209) q[2];
sx q[2];
rz(0.048962083) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.088868695) q[1];
sx q[1];
rz(-1.0430599) q[1];
sx q[1];
rz(2.0305968) q[1];
rz(2.5185561) q[3];
sx q[3];
rz(-1.5554232) q[3];
sx q[3];
rz(-0.29237177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6752601) q[2];
sx q[2];
rz(-1.2559428) q[2];
sx q[2];
rz(2.7272398) q[2];
rz(0.77053344) q[3];
sx q[3];
rz(-1.4221752) q[3];
sx q[3];
rz(0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.82938021) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(-2.6729551) q[0];
rz(1.8679484) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(-0.31148568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92665206) q[0];
sx q[0];
rz(-2.4209341) q[0];
sx q[0];
rz(1.4832627) q[0];
x q[1];
rz(-1.4838672) q[2];
sx q[2];
rz(-1.4540028) q[2];
sx q[2];
rz(2.180336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35986082) q[1];
sx q[1];
rz(-1.6174498) q[1];
sx q[1];
rz(2.0585039) q[1];
rz(-pi) q[2];
rz(2.4021637) q[3];
sx q[3];
rz(-0.66777705) q[3];
sx q[3];
rz(1.8831933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0353388) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(1.3930456) q[2];
rz(-0.92750183) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(2.2813796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504234) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(0.18118478) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(2.9857582) q[2];
sx q[2];
rz(-2.0991785) q[2];
sx q[2];
rz(3.0002158) q[2];
rz(2.2929706) q[3];
sx q[3];
rz(-0.80905882) q[3];
sx q[3];
rz(-0.73425135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
