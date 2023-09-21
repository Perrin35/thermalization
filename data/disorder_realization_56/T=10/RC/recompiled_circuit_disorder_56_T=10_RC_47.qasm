OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1411113) q[0];
sx q[0];
rz(-1.1950462) q[0];
sx q[0];
rz(-1.5953338) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22123863) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(2.0146807) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7386908) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(-0.72523592) q[1];
rz(2.5146033) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(1.0502846) q[2];
rz(1.1132647) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0691836) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(-1.108095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28716921) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(-1.9186583) q[0];
rz(0.42995288) q[2];
sx q[2];
rz(-2.1351095) q[2];
sx q[2];
rz(2.058409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.260028) q[1];
sx q[1];
rz(-1.2177055) q[1];
sx q[1];
rz(-0.41207037) q[1];
rz(1.3012582) q[3];
sx q[3];
rz(-0.88629913) q[3];
sx q[3];
rz(-2.136363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46253282) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(-2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99795139) q[0];
sx q[0];
rz(-2.7086341) q[0];
sx q[0];
rz(0.71822449) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(-2.9253935) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94273401) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(2.0595466) q[1];
rz(-1.987625) q[3];
sx q[3];
rz(-1.7524476) q[3];
sx q[3];
rz(2.3895398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.09459153) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(-2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(0.75727063) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(2.6230985) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-2.8994697) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69232363) q[0];
sx q[0];
rz(-1.7334492) q[0];
sx q[0];
rz(-0.67740324) q[0];
rz(-pi) q[1];
rz(-3.1130586) q[2];
sx q[2];
rz(-0.42978537) q[2];
sx q[2];
rz(-0.98791771) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8072847) q[1];
sx q[1];
rz(-2.0172999) q[1];
sx q[1];
rz(-1.7590894) q[1];
rz(2.7235892) q[3];
sx q[3];
rz(-1.4417218) q[3];
sx q[3];
rz(1.2152745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.018192856) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(-1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(-2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(0.36079303) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-1.1345908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.40588) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(2.1924125) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0005433) q[2];
sx q[2];
rz(-1.2417214) q[2];
sx q[2];
rz(-0.82440257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7736588) q[1];
sx q[1];
rz(-1.415167) q[1];
sx q[1];
rz(-0.44981287) q[1];
x q[2];
rz(-1.6547104) q[3];
sx q[3];
rz(-2.0062371) q[3];
sx q[3];
rz(-0.30121379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(0.57265442) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144433) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.9893601) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5637701) q[0];
sx q[0];
rz(-1.676079) q[0];
sx q[0];
rz(0.012751243) q[0];
rz(-pi) q[1];
rz(0.67201519) q[2];
sx q[2];
rz(-2.5294371) q[2];
sx q[2];
rz(-1.200058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8730948) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(-2.6576659) q[1];
x q[2];
rz(-2.4310962) q[3];
sx q[3];
rz(-2.3264255) q[3];
sx q[3];
rz(1.9424155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6727009) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836442) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-0.3947765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(2.594069) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2497181) q[2];
sx q[2];
rz(-1.2121388) q[2];
sx q[2];
rz(1.1527588) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.863184) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(0.43988887) q[1];
rz(-pi) q[2];
rz(-0.23468252) q[3];
sx q[3];
rz(-2.9838786) q[3];
sx q[3];
rz(1.7705256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(2.8685692) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(1.8008908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.1800982) q[0];
sx q[0];
rz(-0.35238738) q[0];
rz(-1.5199392) q[2];
sx q[2];
rz(-0.87561456) q[2];
sx q[2];
rz(2.4620172) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9400085) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(1.1854978) q[1];
rz(-pi) q[2];
rz(-2.1741381) q[3];
sx q[3];
rz(-2.6199322) q[3];
sx q[3];
rz(2.012901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-2.1179874) q[2];
rz(2.9566531) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-0.83713371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22568208) q[0];
sx q[0];
rz(-1.3689694) q[0];
sx q[0];
rz(2.4542144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8673973) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(-1.1013168) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21048927) q[1];
sx q[1];
rz(-1.466202) q[1];
sx q[1];
rz(1.4180257) q[1];
rz(-2.0765452) q[3];
sx q[3];
rz(-0.55378434) q[3];
sx q[3];
rz(1.6365285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(0.17962757) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76374617) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(2.0843704) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(-2.1616139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36469034) q[0];
sx q[0];
rz(-1.6216462) q[0];
sx q[0];
rz(-1.6965673) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3780031) q[2];
sx q[2];
rz(-2.4739403) q[2];
sx q[2];
rz(1.180336) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.500463) q[1];
sx q[1];
rz(-2.6424721) q[1];
sx q[1];
rz(-1.0541381) q[1];
rz(-1.8944593) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(2.1137962) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733611) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(0.83203075) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(1.0118532) q[2];
sx q[2];
rz(-1.7218628) q[2];
sx q[2];
rz(-2.0414258) q[2];
rz(0.50073033) q[3];
sx q[3];
rz(-1.2320319) q[3];
sx q[3];
rz(-2.4945955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];