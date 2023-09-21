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
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0924661) q[0];
sx q[0];
rz(-1.2611715) q[0];
sx q[0];
rz(-1.5686839) q[0];
rz(-pi) q[1];
rz(0.78656466) q[2];
sx q[2];
rz(-2.2872891) q[2];
sx q[2];
rz(2.5094945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4268036) q[1];
sx q[1];
rz(-1.1183294) q[1];
sx q[1];
rz(-1.2697551) q[1];
rz(-pi) q[2];
rz(0.2829708) q[3];
sx q[3];
rz(-1.030778) q[3];
sx q[3];
rz(0.52373826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6538438) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(-0.1208819) q[2];
rz(2.9623048) q[3];
sx q[3];
rz(-2.5458953) q[3];
sx q[3];
rz(2.9860935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091846175) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(-0.13277408) q[0];
rz(1.6800539) q[1];
sx q[1];
rz(-1.5802054) q[1];
sx q[1];
rz(-0.24138385) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7624843) q[0];
sx q[0];
rz(-0.53775162) q[0];
sx q[0];
rz(2.3178029) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52526371) q[2];
sx q[2];
rz(-1.1366476) q[2];
sx q[2];
rz(0.94793749) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6388968) q[1];
sx q[1];
rz(-2.2310889) q[1];
sx q[1];
rz(3.026282) q[1];
x q[2];
rz(-1.4003795) q[3];
sx q[3];
rz(-2.9885871) q[3];
sx q[3];
rz(1.9433598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0138578) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7805507) q[0];
sx q[0];
rz(-2.0401968) q[0];
sx q[0];
rz(2.4011491) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.8809044) q[1];
sx q[1];
rz(-2.0887451) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6360639) q[0];
sx q[0];
rz(-1.9331421) q[0];
sx q[0];
rz(-0.58689582) q[0];
rz(-2.8430976) q[2];
sx q[2];
rz(-2.2224991) q[2];
sx q[2];
rz(0.088108206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82008541) q[1];
sx q[1];
rz(-0.19430375) q[1];
sx q[1];
rz(0.66837515) q[1];
rz(-pi) q[2];
rz(1.1406844) q[3];
sx q[3];
rz(-1.5617237) q[3];
sx q[3];
rz(-0.7113925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8911002) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(2.5946674) q[2];
rz(-2.8524103) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(-0.012399013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(1.7893715) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(0.23637493) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5445697) q[0];
sx q[0];
rz(-1.3057846) q[0];
sx q[0];
rz(-0.17770627) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0868446) q[2];
sx q[2];
rz(-1.8785254) q[2];
sx q[2];
rz(-2.9204521) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.786495) q[1];
sx q[1];
rz(-1.7384643) q[1];
sx q[1];
rz(-0.0066890072) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27487367) q[3];
sx q[3];
rz(-2.636424) q[3];
sx q[3];
rz(1.2884017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2410879) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(0.55348712) q[2];
rz(2.2262946) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6699162) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(0.70415235) q[0];
rz(-1.0559121) q[1];
sx q[1];
rz(-0.45509714) q[1];
sx q[1];
rz(2.5456837) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10396118) q[0];
sx q[0];
rz(-2.9069293) q[0];
sx q[0];
rz(-1.9734567) q[0];
rz(2.4679568) q[2];
sx q[2];
rz(-1.4623702) q[2];
sx q[2];
rz(3.0631531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6223645) q[1];
sx q[1];
rz(-1.16751) q[1];
sx q[1];
rz(-1.3695903) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6951407) q[3];
sx q[3];
rz(-2.7342396) q[3];
sx q[3];
rz(-2.8706467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(0.55111432) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9838487) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(0.95170784) q[0];
rz(0.47479182) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(2.8463083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26074043) q[0];
sx q[0];
rz(-1.4319001) q[0];
sx q[0];
rz(1.9689346) q[0];
rz(-pi) q[1];
rz(-0.45093243) q[2];
sx q[2];
rz(-0.34313289) q[2];
sx q[2];
rz(-0.97230881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5438248) q[1];
sx q[1];
rz(-0.73111594) q[1];
sx q[1];
rz(0.75146971) q[1];
rz(-pi) q[2];
rz(2.0814249) q[3];
sx q[3];
rz(-1.6276629) q[3];
sx q[3];
rz(2.1631654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.3662806) q[2];
sx q[2];
rz(-2.2078216) q[2];
rz(1.4592524) q[3];
sx q[3];
rz(-1.7873584) q[3];
sx q[3];
rz(1.6850083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0657601) q[0];
sx q[0];
rz(-1.1446784) q[0];
sx q[0];
rz(2.899535) q[0];
rz(0.66479713) q[1];
sx q[1];
rz(-1.8533862) q[1];
sx q[1];
rz(-2.738293) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6751854) q[0];
sx q[0];
rz(-1.9403606) q[0];
sx q[0];
rz(0.36375605) q[0];
x q[1];
rz(2.4629668) q[2];
sx q[2];
rz(-1.7947065) q[2];
sx q[2];
rz(-0.14397552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2029018) q[1];
sx q[1];
rz(-2.9314329) q[1];
sx q[1];
rz(1.5797257) q[1];
rz(-pi) q[2];
rz(1.7014916) q[3];
sx q[3];
rz(-0.19154597) q[3];
sx q[3];
rz(-0.46476118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2287801) q[2];
sx q[2];
rz(-1.0711203) q[2];
sx q[2];
rz(2.7071803) q[2];
rz(-1.0007535) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(-2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(2.6469321) q[0];
rz(1.5085295) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52727115) q[0];
sx q[0];
rz(-1.9691756) q[0];
sx q[0];
rz(-2.6398753) q[0];
x q[1];
rz(2.0958488) q[2];
sx q[2];
rz(-1.5480435) q[2];
sx q[2];
rz(-0.75726985) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8200092) q[1];
sx q[1];
rz(-1.6899741) q[1];
sx q[1];
rz(-0.18305852) q[1];
rz(-0.59372254) q[3];
sx q[3];
rz(-2.7799118) q[3];
sx q[3];
rz(2.7152674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0042469) q[2];
sx q[2];
rz(-0.9937976) q[2];
sx q[2];
rz(-0.11631575) q[2];
rz(0.4256734) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(-1*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.9882934) q[0];
sx q[0];
rz(-0.17803742) q[0];
sx q[0];
rz(-1.4784038) q[0];
rz(-0.93961811) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(-2.7240662) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5696213) q[0];
sx q[0];
rz(-2.6082391) q[0];
sx q[0];
rz(-1.6589952) q[0];
rz(-3.0350424) q[2];
sx q[2];
rz(-1.8333734) q[2];
sx q[2];
rz(-3.0343461) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6431943) q[1];
sx q[1];
rz(-0.25499757) q[1];
sx q[1];
rz(-0.71868371) q[1];
rz(-pi) q[2];
rz(-0.19374356) q[3];
sx q[3];
rz(-1.994641) q[3];
sx q[3];
rz(-1.2731874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4801165) q[2];
sx q[2];
rz(-2.5051703) q[2];
sx q[2];
rz(0.49368668) q[2];
rz(-2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(2.8216968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(0.68558145) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(2.0102274) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6146651) q[0];
sx q[0];
rz(-1.0173431) q[0];
sx q[0];
rz(1.1620031) q[0];
x q[1];
rz(0.71815021) q[2];
sx q[2];
rz(-0.8469204) q[2];
sx q[2];
rz(1.3494327) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79216204) q[1];
sx q[1];
rz(-1.3467448) q[1];
sx q[1];
rz(-2.7538124) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64254909) q[3];
sx q[3];
rz(-0.71392871) q[3];
sx q[3];
rz(2.97314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82548213) q[2];
sx q[2];
rz(-1.1967412) q[2];
sx q[2];
rz(0.70739174) q[2];
rz(0.80983821) q[3];
sx q[3];
rz(-0.67088586) q[3];
sx q[3];
rz(0.18856089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(-2.9293625) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(-0.75081667) q[2];
sx q[2];
rz(-2.3859947) q[2];
sx q[2];
rz(1.044556) q[2];
rz(2.2740169) q[3];
sx q[3];
rz(-2.0316364) q[3];
sx q[3];
rz(0.33566725) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
