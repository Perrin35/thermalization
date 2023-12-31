OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0377169) q[0];
sx q[0];
rz(-1.2020943) q[0];
sx q[0];
rz(-1.9934959) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(-1.3936477) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049126547) q[0];
sx q[0];
rz(-1.8804212) q[0];
sx q[0];
rz(-1.5686839) q[0];
rz(2.2533478) q[2];
sx q[2];
rz(-1.0091072) q[2];
sx q[2];
rz(0.35866666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3323101) q[1];
sx q[1];
rz(-2.6039632) q[1];
sx q[1];
rz(0.54772954) q[1];
x q[2];
rz(-1.0127134) q[3];
sx q[3];
rz(-1.8126243) q[3];
sx q[3];
rz(1.946132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-3.0207108) q[2];
rz(-0.17928784) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(0.13277408) q[0];
rz(-1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(0.24138385) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55610181) q[0];
sx q[0];
rz(-1.9560768) q[0];
sx q[0];
rz(-2.7566064) q[0];
rz(-pi) q[1];
rz(0.74626211) q[2];
sx q[2];
rz(-0.66821874) q[2];
sx q[2];
rz(1.2506739) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0028210359) q[1];
sx q[1];
rz(-1.4797987) q[1];
sx q[1];
rz(-0.9072733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4003795) q[3];
sx q[3];
rz(-0.15300551) q[3];
sx q[3];
rz(-1.9433598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.8227791) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7805507) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(2.0887451) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7171779) q[0];
sx q[0];
rz(-2.4632235) q[0];
sx q[0];
rz(2.5413187) q[0];
rz(-pi) q[1];
rz(-0.29849507) q[2];
sx q[2];
rz(-0.91909354) q[2];
sx q[2];
rz(-3.0534844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82008541) q[1];
sx q[1];
rz(-2.9472889) q[1];
sx q[1];
rz(0.66837515) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1406844) q[3];
sx q[3];
rz(-1.579869) q[3];
sx q[3];
rz(0.7113925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(2.8524103) q[3];
sx q[3];
rz(-2.6047891) q[3];
sx q[3];
rz(-0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1716487) q[0];
sx q[0];
rz(-0.26370731) q[0];
sx q[0];
rz(1.7893715) q[0];
rz(-0.2098473) q[1];
sx q[1];
rz(-2.4317957) q[1];
sx q[1];
rz(2.9052177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.597023) q[0];
sx q[0];
rz(-1.835808) q[0];
sx q[0];
rz(-0.17770627) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7912555) q[2];
sx q[2];
rz(-1.0812034) q[2];
sx q[2];
rz(-1.9621153) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35509767) q[1];
sx q[1];
rz(-1.4031283) q[1];
sx q[1];
rz(-3.1349036) q[1];
rz(-0.27487367) q[3];
sx q[3];
rz(-2.636424) q[3];
sx q[3];
rz(-1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2410879) q[2];
sx q[2];
rz(-2.1630478) q[2];
sx q[2];
rz(-2.5881055) q[2];
rz(-2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(0.64490157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6699162) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(2.4374403) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(2.5456837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6248557) q[0];
sx q[0];
rz(-1.3552109) q[0];
sx q[0];
rz(0.09341021) q[0];
rz(-pi) q[1];
rz(1.709183) q[2];
sx q[2];
rz(-0.90183898) q[2];
sx q[2];
rz(1.578518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1314428) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(-0.4106945) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7704352) q[3];
sx q[3];
rz(-1.3988929) q[3];
sx q[3];
rz(-0.88574823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.492505) q[2];
sx q[2];
rz(-2.0272144) q[2];
sx q[2];
rz(-0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.3144349) q[3];
sx q[3];
rz(1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.15774396) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(-2.1898848) q[0];
rz(-0.47479182) q[1];
sx q[1];
rz(-1.1750849) q[1];
sx q[1];
rz(-0.29528433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2518894) q[0];
sx q[0];
rz(-1.1767052) q[0];
sx q[0];
rz(2.9910812) q[0];
rz(-2.6906602) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(0.97230881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5807647) q[1];
sx q[1];
rz(-1.0974713) q[1];
sx q[1];
rz(2.5614489) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4548364) q[3];
sx q[3];
rz(-0.51350683) q[3];
sx q[3];
rz(-2.6503369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(-1.6823403) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(1.6850083) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(-0.24205762) q[0];
rz(-2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(2.738293) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.032265183) q[0];
sx q[0];
rz(-1.2326213) q[0];
sx q[0];
rz(-1.9637252) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4629668) q[2];
sx q[2];
rz(-1.3468862) q[2];
sx q[2];
rz(-2.9976171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.51822) q[1];
sx q[1];
rz(-1.5726591) q[1];
sx q[1];
rz(1.7809479) q[1];
x q[2];
rz(-3.1163252) q[3];
sx q[3];
rz(-1.380904) q[3];
sx q[3];
rz(-2.5437298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(-2.7071803) q[2];
rz(1.0007535) q[3];
sx q[3];
rz(-2.775511) q[3];
sx q[3];
rz(-2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1338761) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(1.6330632) q[1];
sx q[1];
rz(-1.5155019) q[1];
sx q[1];
rz(2.5411434) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4279815) q[0];
sx q[0];
rz(-0.62987721) q[0];
sx q[0];
rz(0.71891086) q[0];
rz(-1.5254283) q[2];
sx q[2];
rz(-2.6160935) q[2];
sx q[2];
rz(-0.85277992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8200092) q[1];
sx q[1];
rz(-1.4516186) q[1];
sx q[1];
rz(-2.9585341) q[1];
x q[2];
rz(0.30386691) q[3];
sx q[3];
rz(-1.7700723) q[3];
sx q[3];
rz(-0.58135939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(-0.4256734) q[3];
sx q[3];
rz(-2.1843572) q[3];
sx q[3];
rz(11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9882934) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(-2.2019745) q[1];
sx q[1];
rz(-1.8202819) q[1];
sx q[1];
rz(-2.7240662) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077174515) q[0];
sx q[0];
rz(-1.6155956) q[0];
sx q[0];
rz(-1.0391462) q[0];
x q[1];
rz(1.9475627) q[2];
sx q[2];
rz(-2.8586839) q[2];
sx q[2];
rz(2.6434968) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6431943) q[1];
sx q[1];
rz(-0.25499757) q[1];
sx q[1];
rz(2.4229089) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19374356) q[3];
sx q[3];
rz(-1.1469517) q[3];
sx q[3];
rz(1.2731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6614762) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(-0.49368668) q[2];
rz(2.4168329) q[3];
sx q[3];
rz(-1.1586435) q[3];
sx q[3];
rz(2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(0.29742345) q[1];
sx q[1];
rz(-2.90459) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5269276) q[0];
sx q[0];
rz(-2.1242495) q[0];
sx q[0];
rz(-1.1620031) q[0];
rz(-pi) q[1];
rz(0.70558833) q[2];
sx q[2];
rz(-2.0863279) q[2];
sx q[2];
rz(-2.3956092) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28045052) q[1];
sx q[1];
rz(-0.44499731) q[1];
sx q[1];
rz(0.54235561) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60641842) q[3];
sx q[3];
rz(-1.9740205) q[3];
sx q[3];
rz(-1.2244146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3161105) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(0.80983821) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(-0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512882) q[0];
sx q[0];
rz(-2.0177096) q[0];
sx q[0];
rz(2.429005) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(0.75081667) q[2];
sx q[2];
rz(-0.75559794) q[2];
sx q[2];
rz(-2.0970367) q[2];
rz(-2.5645732) q[3];
sx q[3];
rz(-0.95303017) q[3];
sx q[3];
rz(-0.87458761) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
