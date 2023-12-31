OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(-1.9394983) q[0];
sx q[0];
rz(-1.1480968) q[0];
rz(1.2530874) q[1];
sx q[1];
rz(-2.2009067) q[1];
sx q[1];
rz(1.3936477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924661) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(1.5729088) q[0];
rz(2.4601828) q[2];
sx q[2];
rz(-1.0076367) q[2];
sx q[2];
rz(1.5208706) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3323101) q[1];
sx q[1];
rz(-2.6039632) q[1];
sx q[1];
rz(2.5938631) q[1];
x q[2];
rz(-2.0066891) q[3];
sx q[3];
rz(-0.60308686) q[3];
sx q[3];
rz(-3.1325504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(2.9623048) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(0.15549913) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-0.76773983) q[0];
sx q[0];
rz(-3.0088186) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(-2.9002088) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86357147) q[0];
sx q[0];
rz(-1.2153421) q[0];
sx q[0];
rz(1.9832719) q[0];
rz(-2.6163289) q[2];
sx q[2];
rz(-2.004945) q[2];
sx q[2];
rz(0.94793749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1387716) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(-0.9072733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4199735) q[3];
sx q[3];
rz(-1.5449459) q[3];
sx q[3];
rz(-0.54102708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(2.0347118) q[2];
rz(-1.3876623) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.2319516) q[3];
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
rz(-2.7805507) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(-0.45117798) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(-1.0528475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2967865) q[0];
sx q[0];
rz(-2.1150981) q[0];
sx q[0];
rz(1.9980206) q[0];
rz(-0.29849507) q[2];
sx q[2];
rz(-0.91909354) q[2];
sx q[2];
rz(0.088108206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9991418) q[1];
sx q[1];
rz(-1.72292) q[1];
sx q[1];
rz(-1.692148) q[1];
rz(-1.5490407) q[3];
sx q[3];
rz(-0.43020159) q[3];
sx q[3];
rz(-0.83963001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.25049245) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(0.54692522) q[2];
rz(-0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(1.7893715) q[0];
rz(0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(0.23637493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.198092) q[0];
sx q[0];
rz(-0.31790942) q[0];
sx q[0];
rz(2.148118) q[0];
x q[1];
rz(2.0868446) q[2];
sx q[2];
rz(-1.8785254) q[2];
sx q[2];
rz(-2.9204521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35509767) q[1];
sx q[1];
rz(-1.4031283) q[1];
sx q[1];
rz(-3.1349036) q[1];
rz(-pi) q[2];
rz(0.27487367) q[3];
sx q[3];
rz(-2.636424) q[3];
sx q[3];
rz(-1.8531909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(-2.5881055) q[2];
rz(-0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(-0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.8138509) q[0];
sx q[0];
rz(0.70415235) q[0];
rz(-2.0856805) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(-2.5456837) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6248557) q[0];
sx q[0];
rz(-1.7863818) q[0];
sx q[0];
rz(-0.09341021) q[0];
rz(-2.4679568) q[2];
sx q[2];
rz(-1.6792225) q[2];
sx q[2];
rz(-0.078439586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1016741) q[1];
sx q[1];
rz(-0.44821757) q[1];
sx q[1];
rz(0.43804534) q[1];
rz(-pi) q[2];
rz(0.44645198) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.492505) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(-0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(-1.5392039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9838487) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(-2.6668008) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8808522) q[0];
sx q[0];
rz(-1.4319001) q[0];
sx q[0];
rz(-1.9689346) q[0];
x q[1];
rz(-0.31111181) q[2];
sx q[2];
rz(-1.4236441) q[2];
sx q[2];
rz(1.0263024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8413137) q[1];
sx q[1];
rz(-2.0804555) q[1];
sx q[1];
rz(-1.0213486) q[1];
rz(-1.0601677) q[3];
sx q[3];
rz(-1.6276629) q[3];
sx q[3];
rz(-0.97842723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39020145) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(0.93377101) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(1.6850083) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-2.738293) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1093275) q[0];
sx q[0];
rz(-1.9089713) q[0];
sx q[0];
rz(1.1778675) q[0];
x q[1];
rz(-1.2861916) q[2];
sx q[2];
rz(-0.91214123) q[2];
sx q[2];
rz(1.8919485) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.93869081) q[1];
sx q[1];
rz(-0.21015973) q[1];
sx q[1];
rz(1.5797257) q[1];
x q[2];
rz(-1.3808448) q[3];
sx q[3];
rz(-1.5459832) q[3];
sx q[3];
rz(0.97770377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91281259) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0077165724) q[0];
sx q[0];
rz(-3.1354597) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(1.5085295) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(-2.5411434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7136112) q[0];
sx q[0];
rz(-2.5117154) q[0];
sx q[0];
rz(0.71891086) q[0];
rz(-pi) q[1];
rz(0.026293228) q[2];
sx q[2];
rz(-2.0956989) q[2];
sx q[2];
rz(0.80034791) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32158347) q[1];
sx q[1];
rz(-1.4516186) q[1];
sx q[1];
rz(0.18305852) q[1];
rz(-pi) q[2];
rz(1.3622215) q[3];
sx q[3];
rz(-1.2731291) q[3];
sx q[3];
rz(-1.0514333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0042469) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(-2.7159193) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9882934) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(-1.6631888) q[0];
rz(0.93961811) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(2.7240662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077174515) q[0];
sx q[0];
rz(-1.525997) q[0];
sx q[0];
rz(-1.0391462) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9475627) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(2.6434968) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36996499) q[1];
sx q[1];
rz(-1.4039478) q[1];
sx q[1];
rz(0.19373993) q[1];
rz(-0.19374356) q[3];
sx q[3];
rz(-1.994641) q[3];
sx q[3];
rz(1.8684052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(-2.4168329) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(-2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(2.8441692) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(2.0102274) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.92537) q[0];
sx q[0];
rz(-0.67514456) q[0];
sx q[0];
rz(2.569909) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4360043) q[2];
sx q[2];
rz(-2.0863279) q[2];
sx q[2];
rz(-0.74598344) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86913925) q[1];
sx q[1];
rz(-1.9483882) q[1];
sx q[1];
rz(1.8121522) q[1];
x q[2];
rz(-0.60641842) q[3];
sx q[3];
rz(-1.1675721) q[3];
sx q[3];
rz(1.917178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3161105) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-2.4342009) q[2];
rz(-2.3317544) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(2.9530318) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.5384197) q[2];
sx q[2];
rz(-2.0576253) q[2];
sx q[2];
rz(-1.1228592) q[2];
rz(2.5645732) q[3];
sx q[3];
rz(-2.1885625) q[3];
sx q[3];
rz(2.267005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
