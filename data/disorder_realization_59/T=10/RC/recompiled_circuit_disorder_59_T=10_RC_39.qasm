OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479552) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(1.1975343) q[0];
rz(-0.087287993) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(1.0686312) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6359771) q[1];
sx q[1];
rz(-2.639289) q[1];
sx q[1];
rz(-2.99519) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9418342) q[3];
sx q[3];
rz(-1.8098117) q[3];
sx q[3];
rz(-2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(0.74938613) q[2];
rz(-2.1253712) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(2.7367676) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4352903) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(-1.0372112) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(-0.81545365) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6204651) q[0];
sx q[0];
rz(-1.5107811) q[0];
sx q[0];
rz(2.492766) q[0];
rz(-pi) q[1];
rz(0.72950659) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(1.6934998) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.193589) q[1];
sx q[1];
rz(-1.5352328) q[1];
sx q[1];
rz(-0.30102945) q[1];
rz(-pi) q[2];
rz(2.4957982) q[3];
sx q[3];
rz(-1.3919953) q[3];
sx q[3];
rz(1.7088695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(-2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3011424) q[0];
sx q[0];
rz(-1.2121032) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8529352) q[0];
sx q[0];
rz(-1.4212199) q[0];
sx q[0];
rz(-1.268671) q[0];
rz(0.61727662) q[2];
sx q[2];
rz(-1.0059788) q[2];
sx q[2];
rz(-1.9895983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42400186) q[1];
sx q[1];
rz(-1.0832936) q[1];
sx q[1];
rz(1.2815777) q[1];
x q[2];
rz(-3.0619377) q[3];
sx q[3];
rz(-2.0783391) q[3];
sx q[3];
rz(1.5910651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.375305) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(-2.2312009) q[0];
rz(2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7461473) q[0];
sx q[0];
rz(-1.6593885) q[0];
sx q[0];
rz(-3.0625312) q[0];
x q[1];
rz(1.1103815) q[2];
sx q[2];
rz(-1.5693671) q[2];
sx q[2];
rz(-1.187385) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7819314) q[1];
sx q[1];
rz(-2.0205824) q[1];
sx q[1];
rz(0.91314544) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4196017) q[3];
sx q[3];
rz(-0.70221838) q[3];
sx q[3];
rz(-0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3466907) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(2.9798853) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(0.99779469) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.6246187) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.957513) q[0];
sx q[0];
rz(-2.3947869) q[0];
sx q[0];
rz(-2.1225131) q[0];
rz(-pi) q[1];
rz(-0.5337358) q[2];
sx q[2];
rz(-1.2297451) q[2];
sx q[2];
rz(1.595572) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9285674) q[1];
sx q[1];
rz(-1.9145609) q[1];
sx q[1];
rz(-2.7108971) q[1];
x q[2];
rz(2.3715641) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-0.023035223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(1.2639686) q[0];
rz(-0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(2.8009159) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1893038) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(-0.031866372) q[0];
rz(-pi) q[1];
rz(-2.3855626) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(-2.0331969) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6319879) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(-1.8297086) q[1];
rz(-pi) q[2];
rz(0.86990279) q[3];
sx q[3];
rz(-2.5330336) q[3];
sx q[3];
rz(0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-0.9698202) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9724378) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(0.06282839) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(-2.6766434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9879887) q[0];
sx q[0];
rz(-1.1055595) q[0];
sx q[0];
rz(0.21501712) q[0];
rz(2.6480688) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(-1.354419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17864922) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(-2.9220198) q[1];
x q[2];
rz(1.6169448) q[3];
sx q[3];
rz(-1.5449636) q[3];
sx q[3];
rz(1.8125364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(-2.590086) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(3.0632339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4975472) q[0];
sx q[0];
rz(-2.3345778) q[0];
sx q[0];
rz(3.0456196) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8163082) q[2];
sx q[2];
rz(-1.3405372) q[2];
sx q[2];
rz(-2.2895209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52787493) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(-1.2916958) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3571635) q[3];
sx q[3];
rz(-0.16031081) q[3];
sx q[3];
rz(-0.3233288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(-0.25137869) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79779977) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-0.051368512) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(2.267568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6286205) q[0];
sx q[0];
rz(-0.7681094) q[0];
sx q[0];
rz(-3.1290929) q[0];
rz(-2.7810077) q[2];
sx q[2];
rz(-1.9546095) q[2];
sx q[2];
rz(-1.4521445) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8928788) q[1];
sx q[1];
rz(-1.6422179) q[1];
sx q[1];
rz(-1.0700657) q[1];
x q[2];
rz(0.47253982) q[3];
sx q[3];
rz(-0.66415411) q[3];
sx q[3];
rz(3.1411375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(0.61974636) q[2];
rz(1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(-1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.9627409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282384) q[0];
sx q[0];
rz(-1.2984707) q[0];
sx q[0];
rz(2.9288835) q[0];
x q[1];
rz(3.0145698) q[2];
sx q[2];
rz(-1.6633031) q[2];
sx q[2];
rz(-0.18735838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9465543) q[1];
sx q[1];
rz(-2.185501) q[1];
sx q[1];
rz(-2.9123995) q[1];
rz(1.3445271) q[3];
sx q[3];
rz(-1.0591649) q[3];
sx q[3];
rz(0.98480485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62109229) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(3.1065431) q[2];
sx q[2];
rz(-2.0106342) q[2];
sx q[2];
rz(1.0019279) q[2];
rz(2.6562128) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
