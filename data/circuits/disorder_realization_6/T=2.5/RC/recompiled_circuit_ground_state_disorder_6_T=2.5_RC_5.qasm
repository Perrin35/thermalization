OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4226469) q[0];
sx q[0];
rz(-0.5889686) q[0];
sx q[0];
rz(0.11864057) q[0];
rz(2.150382) q[1];
sx q[1];
rz(-1.041643) q[1];
sx q[1];
rz(0.51377327) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2048747) q[0];
sx q[0];
rz(-0.74130171) q[0];
sx q[0];
rz(0.78420774) q[0];
rz(1.6330326) q[2];
sx q[2];
rz(-2.0449491) q[2];
sx q[2];
rz(-1.4197503) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.895925) q[1];
sx q[1];
rz(-1.2497207) q[1];
sx q[1];
rz(-1.3447869) q[1];
rz(0.051187201) q[3];
sx q[3];
rz(-0.59268206) q[3];
sx q[3];
rz(2.7251303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.674268) q[2];
sx q[2];
rz(-1.5200204) q[2];
sx q[2];
rz(-2.8094214) q[2];
rz(-2.8915571) q[3];
sx q[3];
rz(-1.9146405) q[3];
sx q[3];
rz(-1.2927607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2501204) q[0];
sx q[0];
rz(-1.253506) q[0];
sx q[0];
rz(2.6943595) q[0];
rz(-2.1121292) q[1];
sx q[1];
rz(-2.5599458) q[1];
sx q[1];
rz(3.0404125) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6995938) q[0];
sx q[0];
rz(-1.9706107) q[0];
sx q[0];
rz(-2.6599822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20002983) q[2];
sx q[2];
rz(-1.0009655) q[2];
sx q[2];
rz(-2.9404158) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14462073) q[1];
sx q[1];
rz(-0.50077866) q[1];
sx q[1];
rz(-2.5624496) q[1];
rz(3.0291843) q[3];
sx q[3];
rz(-0.93494697) q[3];
sx q[3];
rz(-2.4684255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0589361) q[2];
sx q[2];
rz(-2.3857748) q[2];
sx q[2];
rz(1.6015046) q[2];
rz(-2.5236409) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(-1.1624973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5471632) q[0];
sx q[0];
rz(-2.1367441) q[0];
sx q[0];
rz(-0.51602236) q[0];
rz(-2.0891321) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(-1.1722391) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7965423) q[0];
sx q[0];
rz(-1.5818514) q[0];
sx q[0];
rz(-0.42406908) q[0];
x q[1];
rz(0.86031057) q[2];
sx q[2];
rz(-1.5154308) q[2];
sx q[2];
rz(-2.4982128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3400628) q[1];
sx q[1];
rz(-2.0593567) q[1];
sx q[1];
rz(-0.21111458) q[1];
x q[2];
rz(-3.1276305) q[3];
sx q[3];
rz(-0.80028557) q[3];
sx q[3];
rz(1.6228907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8664794) q[2];
sx q[2];
rz(-1.8786414) q[2];
sx q[2];
rz(-2.4647554) q[2];
rz(0.46946851) q[3];
sx q[3];
rz(-2.2468086) q[3];
sx q[3];
rz(1.2137871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6770099) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(1.1997892) q[0];
rz(0.59263539) q[1];
sx q[1];
rz(-1.0721782) q[1];
sx q[1];
rz(1.6823654) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30577393) q[0];
sx q[0];
rz(-2.523382) q[0];
sx q[0];
rz(2.799116) q[0];
rz(2.8921669) q[2];
sx q[2];
rz(-0.81455961) q[2];
sx q[2];
rz(-3.0944173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.063733405) q[1];
sx q[1];
rz(-2.1850536) q[1];
sx q[1];
rz(-2.8251404) q[1];
x q[2];
rz(-1.9497037) q[3];
sx q[3];
rz(-0.92273308) q[3];
sx q[3];
rz(-0.60175446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49372855) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(2.4596821) q[2];
rz(-0.41334263) q[3];
sx q[3];
rz(-1.6091434) q[3];
sx q[3];
rz(-1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84363627) q[0];
sx q[0];
rz(-1.4763259) q[0];
sx q[0];
rz(2.8635039) q[0];
rz(-0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(-0.98731891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.322418) q[0];
sx q[0];
rz(-1.8083982) q[0];
sx q[0];
rz(-2.1850719) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5426719) q[2];
sx q[2];
rz(-1.4338981) q[2];
sx q[2];
rz(-1.2281451) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1121171) q[1];
sx q[1];
rz(-1.6320845) q[1];
sx q[1];
rz(-1.960683) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0342399) q[3];
sx q[3];
rz(-1.3349541) q[3];
sx q[3];
rz(0.72808108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.014293369) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(-0.75931749) q[2];
rz(-1.6589818) q[3];
sx q[3];
rz(-1.1967775) q[3];
sx q[3];
rz(-0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229014) q[0];
sx q[0];
rz(-2.3221115) q[0];
sx q[0];
rz(2.4500093) q[0];
rz(0.85451952) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(2.7202594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3499419) q[0];
sx q[0];
rz(-1.5022819) q[0];
sx q[0];
rz(-1.1380026) q[0];
rz(-pi) q[1];
rz(1.2796655) q[2];
sx q[2];
rz(-1.3145245) q[2];
sx q[2];
rz(0.95367431) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5798963) q[1];
sx q[1];
rz(-1.4868127) q[1];
sx q[1];
rz(0.3616228) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.297703) q[3];
sx q[3];
rz(-0.98326245) q[3];
sx q[3];
rz(-1.3466101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2104346) q[2];
sx q[2];
rz(-1.4782108) q[2];
sx q[2];
rz(-2.9388536) q[2];
rz(-1.5984795) q[3];
sx q[3];
rz(-0.32035443) q[3];
sx q[3];
rz(-1.6725484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-2.8572674) q[0];
sx q[0];
rz(-1.4284644) q[0];
sx q[0];
rz(2.7427234) q[0];
rz(-1.9786037) q[1];
sx q[1];
rz(-2.3154924) q[1];
sx q[1];
rz(0.2805447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79518455) q[0];
sx q[0];
rz(-1.1885671) q[0];
sx q[0];
rz(-2.7104177) q[0];
x q[1];
rz(-0.43088308) q[2];
sx q[2];
rz(-0.93965778) q[2];
sx q[2];
rz(0.43090303) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7458492) q[1];
sx q[1];
rz(-1.1458995) q[1];
sx q[1];
rz(0.33558553) q[1];
rz(-pi) q[2];
rz(0.79645282) q[3];
sx q[3];
rz(-2.011125) q[3];
sx q[3];
rz(1.274282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82105381) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(1.4873571) q[2];
rz(-0.94333831) q[3];
sx q[3];
rz(-2.2127559) q[3];
sx q[3];
rz(1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-2.6392537) q[0];
rz(0.54652864) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(-0.46195236) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34794688) q[0];
sx q[0];
rz(-1.0805305) q[0];
sx q[0];
rz(2.8099174) q[0];
x q[1];
rz(1.5080323) q[2];
sx q[2];
rz(-2.4893005) q[2];
sx q[2];
rz(1.7386652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4224902) q[1];
sx q[1];
rz(-1.6797069) q[1];
sx q[1];
rz(-1.9816573) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5580719) q[3];
sx q[3];
rz(-1.0403578) q[3];
sx q[3];
rz(0.96543559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97303566) q[2];
sx q[2];
rz(-0.23739561) q[2];
sx q[2];
rz(0.17042223) q[2];
rz(2.6863344) q[3];
sx q[3];
rz(-1.596343) q[3];
sx q[3];
rz(0.71793238) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1667204) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(0.18044743) q[0];
rz(-0.51668733) q[1];
sx q[1];
rz(-1.5770715) q[1];
sx q[1];
rz(2.7076941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93709263) q[0];
sx q[0];
rz(-2.4334416) q[0];
sx q[0];
rz(0.88237472) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1230311) q[2];
sx q[2];
rz(-1.8862714) q[2];
sx q[2];
rz(1.7231143) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2502807) q[1];
sx q[1];
rz(-1.6906594) q[1];
sx q[1];
rz(0.93617546) q[1];
rz(-0.88017616) q[3];
sx q[3];
rz(-1.343702) q[3];
sx q[3];
rz(-2.2736062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5608998) q[2];
sx q[2];
rz(-2.7036724) q[2];
sx q[2];
rz(-1.8221347) q[2];
rz(2.8699919) q[3];
sx q[3];
rz(-1.3180132) q[3];
sx q[3];
rz(-1.1118838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6962947) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(1.7355504) q[0];
rz(1.4259074) q[1];
sx q[1];
rz(-2.0366663) q[1];
sx q[1];
rz(-1.8941194) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1734764) q[0];
sx q[0];
rz(-2.0531468) q[0];
sx q[0];
rz(2.8456057) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0872074) q[2];
sx q[2];
rz(-1.6882036) q[2];
sx q[2];
rz(-0.94737999) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0669603) q[1];
sx q[1];
rz(-1.6116643) q[1];
sx q[1];
rz(0.34217477) q[1];
x q[2];
rz(-0.62274751) q[3];
sx q[3];
rz(-1.0120262) q[3];
sx q[3];
rz(-1.8468685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71331104) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(0.10078079) q[2];
rz(0.0095602592) q[3];
sx q[3];
rz(-1.5567501) q[3];
sx q[3];
rz(-1.5338219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99209256) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-0.59303444) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(3.0708457) q[2];
sx q[2];
rz(-1.1732709) q[2];
sx q[2];
rz(2.7931961) q[2];
rz(1.5531874) q[3];
sx q[3];
rz(-1.7918158) q[3];
sx q[3];
rz(0.50504897) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
