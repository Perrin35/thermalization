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
rz(0.50245589) q[0];
sx q[0];
rz(4.2491663) q[0];
sx q[0];
rz(7.6585461) q[0];
rz(-3.4604685) q[1];
sx q[1];
rz(4.7825216) q[1];
sx q[1];
rz(8.6934269) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36424822) q[0];
sx q[0];
rz(-2.4349182) q[0];
sx q[0];
rz(1.0663435) q[0];
rz(-pi) q[1];
rz(-0.0053623036) q[2];
sx q[2];
rz(-2.4796133) q[2];
sx q[2];
rz(-2.299451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9125376) q[1];
sx q[1];
rz(-2.3542535) q[1];
sx q[1];
rz(-3.0773735) q[1];
rz(-pi) q[2];
rz(-1.4764686) q[3];
sx q[3];
rz(-1.7324311) q[3];
sx q[3];
rz(1.9548655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88087624) q[2];
sx q[2];
rz(-1.2141576) q[2];
sx q[2];
rz(0.94240776) q[2];
rz(-0.85855329) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.079916) q[0];
sx q[0];
rz(-0.55773568) q[0];
sx q[0];
rz(-3.0829561) q[0];
rz(3.0851641) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-2.0243534) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4898713) q[0];
sx q[0];
rz(-2.3014223) q[0];
sx q[0];
rz(2.1108759) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12877474) q[2];
sx q[2];
rz(-1.2733545) q[2];
sx q[2];
rz(-0.25624881) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8630115) q[1];
sx q[1];
rz(-1.8092691) q[1];
sx q[1];
rz(-0.084785297) q[1];
rz(-pi) q[2];
rz(2.9100203) q[3];
sx q[3];
rz(-0.48677126) q[3];
sx q[3];
rz(-0.51460217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77659082) q[2];
sx q[2];
rz(-2.8114909) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(-2.8255919) q[3];
sx q[3];
rz(-1.6812811) q[3];
sx q[3];
rz(-2.3587904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86421788) q[0];
sx q[0];
rz(-2.6835231) q[0];
sx q[0];
rz(-1.3360485) q[0];
rz(2.863073) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(2.1727402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3144041) q[0];
sx q[0];
rz(-2.1705856) q[0];
sx q[0];
rz(1.1161854) q[0];
x q[1];
rz(0.79040852) q[2];
sx q[2];
rz(-1.1539173) q[2];
sx q[2];
rz(-2.0415155) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90889764) q[1];
sx q[1];
rz(-1.8669584) q[1];
sx q[1];
rz(3.0975268) q[1];
x q[2];
rz(2.8519451) q[3];
sx q[3];
rz(-0.74637369) q[3];
sx q[3];
rz(1.0371072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99732533) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(-0.17511314) q[2];
rz(1.3052321) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(-1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433732) q[0];
sx q[0];
rz(-1.136919) q[0];
sx q[0];
rz(-2.112222) q[0];
rz(-2.3221305) q[1];
sx q[1];
rz(-1.1623323) q[1];
sx q[1];
rz(-2.3447461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5982847) q[0];
sx q[0];
rz(-2.6572022) q[0];
sx q[0];
rz(2.7192857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0060002319) q[2];
sx q[2];
rz(-1.0428535) q[2];
sx q[2];
rz(-1.0061629) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0734892) q[1];
sx q[1];
rz(-1.3655546) q[1];
sx q[1];
rz(-1.2935335) q[1];
rz(-0.66241499) q[3];
sx q[3];
rz(-1.9180505) q[3];
sx q[3];
rz(-1.9062689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0113819) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(1.9127362) q[2];
rz(-0.43904385) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(-1.9108093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586886) q[0];
sx q[0];
rz(-1.5331601) q[0];
sx q[0];
rz(-2.3332692) q[0];
rz(1.7841548) q[1];
sx q[1];
rz(-0.789855) q[1];
sx q[1];
rz(0.70560169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1418512) q[0];
sx q[0];
rz(-1.6054285) q[0];
sx q[0];
rz(2.8409728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23106261) q[2];
sx q[2];
rz(-2.2602644) q[2];
sx q[2];
rz(-1.5901515) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8568154) q[1];
sx q[1];
rz(-1.5681303) q[1];
sx q[1];
rz(-1.8239914) q[1];
x q[2];
rz(-0.79885317) q[3];
sx q[3];
rz(-2.0435413) q[3];
sx q[3];
rz(2.1780739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6610873) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(-2.535848) q[2];
rz(-0.12039603) q[3];
sx q[3];
rz(-1.2984637) q[3];
sx q[3];
rz(-0.81407434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26009387) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(-0.0023728097) q[0];
rz(-0.35898769) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(-2.733309) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15015442) q[0];
sx q[0];
rz(-1.6065734) q[0];
sx q[0];
rz(-1.4442025) q[0];
rz(-pi) q[1];
rz(0.26770182) q[2];
sx q[2];
rz(-2.1609339) q[2];
sx q[2];
rz(-2.6538284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.830891) q[1];
sx q[1];
rz(-2.9017555) q[1];
sx q[1];
rz(-2.8562921) q[1];
rz(-0.95127912) q[3];
sx q[3];
rz(-0.86182098) q[3];
sx q[3];
rz(-3.0326642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6458873) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(0.1296981) q[2];
rz(-0.4194704) q[3];
sx q[3];
rz(-1.9286112) q[3];
sx q[3];
rz(-2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7523338) q[0];
sx q[0];
rz(-2.3663754) q[0];
sx q[0];
rz(0.68354052) q[0];
rz(2.2031802) q[1];
sx q[1];
rz(-1.7529731) q[1];
sx q[1];
rz(-1.2260812) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6594567) q[0];
sx q[0];
rz(-2.1468751) q[0];
sx q[0];
rz(2.0123737) q[0];
rz(-pi) q[1];
rz(0.27023882) q[2];
sx q[2];
rz(-0.61737379) q[2];
sx q[2];
rz(-0.63948217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56537762) q[1];
sx q[1];
rz(-2.5984796) q[1];
sx q[1];
rz(-2.1469223) q[1];
rz(-pi) q[2];
rz(-0.55415537) q[3];
sx q[3];
rz(-2.1423577) q[3];
sx q[3];
rz(-0.038427834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.73416) q[2];
sx q[2];
rz(-1.4340883) q[2];
sx q[2];
rz(-2.4343991) q[2];
rz(-1.6678984) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(-1.7128806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4555175) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(-0.92380512) q[0];
rz(2.0763981) q[1];
sx q[1];
rz(-1.9520452) q[1];
sx q[1];
rz(1.1346029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4489733) q[0];
sx q[0];
rz(-1.53063) q[0];
sx q[0];
rz(-1.19541) q[0];
rz(-pi) q[1];
rz(0.61750268) q[2];
sx q[2];
rz(-0.98101014) q[2];
sx q[2];
rz(1.2794354) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5091619) q[1];
sx q[1];
rz(-1.0405417) q[1];
sx q[1];
rz(-0.061930818) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52319877) q[3];
sx q[3];
rz(-2.4639795) q[3];
sx q[3];
rz(-2.1385156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.137546) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(2.498632) q[2];
rz(2.9709587) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(2.046106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3951025) q[0];
sx q[0];
rz(-0.13245067) q[0];
sx q[0];
rz(2.3353031) q[0];
rz(-0.0099446615) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(0.26201216) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0341971) q[0];
sx q[0];
rz(-0.55343628) q[0];
sx q[0];
rz(-1.8897927) q[0];
rz(-pi) q[1];
rz(-1.7315032) q[2];
sx q[2];
rz(-2.2101058) q[2];
sx q[2];
rz(-0.29603816) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1478473) q[1];
sx q[1];
rz(-0.48300749) q[1];
sx q[1];
rz(-1.1149393) q[1];
rz(-0.32829653) q[3];
sx q[3];
rz(-2.0496164) q[3];
sx q[3];
rz(1.4737859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.749873) q[2];
sx q[2];
rz(-2.5842857) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(0.90397942) q[3];
sx q[3];
rz(-1.6429792) q[3];
sx q[3];
rz(2.2601295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3751462) q[0];
sx q[0];
rz(-1.2297933) q[0];
sx q[0];
rz(-2.4850856) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.1708941) q[1];
sx q[1];
rz(0.62969977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7555576) q[0];
sx q[0];
rz(-0.81869253) q[0];
sx q[0];
rz(0.56171239) q[0];
x q[1];
rz(0.26816396) q[2];
sx q[2];
rz(-2.0879474) q[2];
sx q[2];
rz(-1.535653) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1765285) q[1];
sx q[1];
rz(-1.5986491) q[1];
sx q[1];
rz(-1.615953) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30977003) q[3];
sx q[3];
rz(-2.0783721) q[3];
sx q[3];
rz(-2.6617756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62275824) q[2];
sx q[2];
rz(-0.82509416) q[2];
sx q[2];
rz(-0.68010124) q[2];
rz(-0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(-1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5730561) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(2.9987891) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-1.8129195) q[2];
sx q[2];
rz(-1.6264781) q[2];
sx q[2];
rz(-1.3774058) q[2];
rz(2.18022) q[3];
sx q[3];
rz(-1.6037887) q[3];
sx q[3];
rz(-2.9870839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
