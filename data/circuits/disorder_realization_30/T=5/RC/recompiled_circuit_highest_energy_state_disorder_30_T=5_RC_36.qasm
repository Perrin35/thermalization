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
rz(-2.8415866) q[0];
sx q[0];
rz(-2.0553148) q[0];
sx q[0];
rz(0.80817428) q[0];
rz(-1.1440682) q[1];
sx q[1];
rz(-0.77163982) q[1];
sx q[1];
rz(1.5626524) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15641744) q[0];
sx q[0];
rz(-2.1083207) q[0];
sx q[0];
rz(-2.318105) q[0];
rz(-pi) q[1];
rz(-1.3035266) q[2];
sx q[2];
rz(-1.7198945) q[2];
sx q[2];
rz(1.9097415) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6526323) q[1];
sx q[1];
rz(-0.51981407) q[1];
sx q[1];
rz(-1.7955154) q[1];
x q[2];
rz(1.5366301) q[3];
sx q[3];
rz(-1.6164268) q[3];
sx q[3];
rz(-0.51866097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.943104) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(-1.0342106) q[2];
rz(0.69283038) q[3];
sx q[3];
rz(-2.0878017) q[3];
sx q[3];
rz(2.0809765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3028054) q[0];
sx q[0];
rz(-3.1133339) q[0];
sx q[0];
rz(-2.5651108) q[0];
rz(-0.020596404) q[1];
sx q[1];
rz(-0.3978022) q[1];
sx q[1];
rz(-1.0284665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3555455) q[0];
sx q[0];
rz(-1.3554327) q[0];
sx q[0];
rz(0.26305612) q[0];
rz(-pi) q[1];
rz(-1.9443545) q[2];
sx q[2];
rz(-1.2837063) q[2];
sx q[2];
rz(2.8127828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0513269) q[1];
sx q[1];
rz(-1.169186) q[1];
sx q[1];
rz(2.7793867) q[1];
rz(0.50039496) q[3];
sx q[3];
rz(-2.430309) q[3];
sx q[3];
rz(0.82373133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2241406) q[2];
sx q[2];
rz(-0.96106207) q[2];
sx q[2];
rz(1.2646328) q[2];
rz(0.97359109) q[3];
sx q[3];
rz(-1.4507989) q[3];
sx q[3];
rz(0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6384386) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(0.85357443) q[0];
rz(-2.0897934) q[1];
sx q[1];
rz(-0.97426668) q[1];
sx q[1];
rz(-0.015017088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029808345) q[0];
sx q[0];
rz(-1.5795603) q[0];
sx q[0];
rz(0.86017227) q[0];
x q[1];
rz(-1.5921345) q[2];
sx q[2];
rz(-1.4747095) q[2];
sx q[2];
rz(-1.0839407) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5530921) q[1];
sx q[1];
rz(-2.4238844) q[1];
sx q[1];
rz(-3.0293082) q[1];
x q[2];
rz(0.28983966) q[3];
sx q[3];
rz(-1.5095169) q[3];
sx q[3];
rz(2.0743528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1198279) q[2];
sx q[2];
rz(-1.4889577) q[2];
sx q[2];
rz(0.28142288) q[2];
rz(1.4540539) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5467095) q[0];
sx q[0];
rz(-1.895772) q[0];
sx q[0];
rz(-2.967714) q[0];
rz(1.8100544) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(-1.7132267) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7998909) q[0];
sx q[0];
rz(-2.9486994) q[0];
sx q[0];
rz(1.6454205) q[0];
rz(1.6971484) q[2];
sx q[2];
rz(-1.9280806) q[2];
sx q[2];
rz(-2.7972691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54879872) q[1];
sx q[1];
rz(-2.4830635) q[1];
sx q[1];
rz(-3.0786773) q[1];
x q[2];
rz(1.1826199) q[3];
sx q[3];
rz(-1.2283162) q[3];
sx q[3];
rz(1.4701209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4398769) q[2];
sx q[2];
rz(-0.82650799) q[2];
sx q[2];
rz(-2.7244549) q[2];
rz(0.16962984) q[3];
sx q[3];
rz(-3.0483584) q[3];
sx q[3];
rz(0.72181845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193264) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(0.30935031) q[0];
rz(2.3720692) q[1];
sx q[1];
rz(-1.0898432) q[1];
sx q[1];
rz(1.5187029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8474591) q[0];
sx q[0];
rz(-1.3193697) q[0];
sx q[0];
rz(-2.2537116) q[0];
rz(-2.2686636) q[2];
sx q[2];
rz(-2.180122) q[2];
sx q[2];
rz(-0.84963679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7243821) q[1];
sx q[1];
rz(-1.3816178) q[1];
sx q[1];
rz(-3.1114905) q[1];
rz(-1.4368016) q[3];
sx q[3];
rz(-1.4729285) q[3];
sx q[3];
rz(-2.0201473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3994483) q[2];
sx q[2];
rz(-2.1592906) q[2];
sx q[2];
rz(0.99545109) q[2];
rz(0.71850145) q[3];
sx q[3];
rz(-1.813443) q[3];
sx q[3];
rz(-0.79513454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28498483) q[0];
sx q[0];
rz(-1.4549078) q[0];
sx q[0];
rz(-1.5307776) q[0];
rz(-1.4467622) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(-1.7036499) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3510144) q[0];
sx q[0];
rz(-1.5912799) q[0];
sx q[0];
rz(-1.162942) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52891584) q[2];
sx q[2];
rz(-2.4806266) q[2];
sx q[2];
rz(1.5724374) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6561778) q[1];
sx q[1];
rz(-0.5159157) q[1];
sx q[1];
rz(-0.1702711) q[1];
rz(-pi) q[2];
rz(-0.0030869129) q[3];
sx q[3];
rz(-0.75957889) q[3];
sx q[3];
rz(-2.742033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1630254) q[2];
sx q[2];
rz(-1.359553) q[2];
sx q[2];
rz(-1.2320409) q[2];
rz(-1.1466522) q[3];
sx q[3];
rz(-1.3403284) q[3];
sx q[3];
rz(-0.31057096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.6335886) q[0];
sx q[0];
rz(-0.81729832) q[0];
sx q[0];
rz(-2.001413) q[0];
rz(0.89705244) q[1];
sx q[1];
rz(-1.3373809) q[1];
sx q[1];
rz(3.1111029) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.193509) q[0];
sx q[0];
rz(-1.3498257) q[0];
sx q[0];
rz(-0.11917178) q[0];
rz(-pi) q[1];
rz(-0.6426446) q[2];
sx q[2];
rz(-1.8360999) q[2];
sx q[2];
rz(2.93769) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97059879) q[1];
sx q[1];
rz(-2.355651) q[1];
sx q[1];
rz(1.1517609) q[1];
rz(-pi) q[2];
rz(-0.84010853) q[3];
sx q[3];
rz(-1.2222154) q[3];
sx q[3];
rz(-0.19053647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.19411479) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(1.8428295) q[2];
rz(1.6537559) q[3];
sx q[3];
rz(-2.1066809) q[3];
sx q[3];
rz(-1.8208767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.79539913) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(1.7561703) q[0];
rz(-1.2162544) q[1];
sx q[1];
rz(-1.4477891) q[1];
sx q[1];
rz(-0.48929712) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70883553) q[0];
sx q[0];
rz(-2.2422195) q[0];
sx q[0];
rz(-1.5455724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4665746) q[2];
sx q[2];
rz(-2.8036593) q[2];
sx q[2];
rz(-0.80160415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.1810822) q[1];
sx q[1];
rz(-1.1000634) q[1];
sx q[1];
rz(3.0357643) q[1];
rz(-1.1888483) q[3];
sx q[3];
rz(-0.26774613) q[3];
sx q[3];
rz(-1.6588039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74799246) q[2];
sx q[2];
rz(-1.1425428) q[2];
sx q[2];
rz(-0.46621123) q[2];
rz(1.6097869) q[3];
sx q[3];
rz(-1.6377056) q[3];
sx q[3];
rz(-2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4894678) q[0];
sx q[0];
rz(-2.2418699) q[0];
sx q[0];
rz(3.0925282) q[0];
rz(0.24066726) q[1];
sx q[1];
rz(-1.0180232) q[1];
sx q[1];
rz(-3.1191471) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3291859) q[0];
sx q[0];
rz(-2.5617449) q[0];
sx q[0];
rz(1.5270698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2227549) q[2];
sx q[2];
rz(-1.9371667) q[2];
sx q[2];
rz(2.6719246) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6091303) q[1];
sx q[1];
rz(-2.3098619) q[1];
sx q[1];
rz(0.66166454) q[1];
rz(-0.53655728) q[3];
sx q[3];
rz(-0.76013764) q[3];
sx q[3];
rz(0.11790568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13751328) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(2.1732886) q[2];
rz(0.83640313) q[3];
sx q[3];
rz(-2.8439549) q[3];
sx q[3];
rz(-1.77232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7264929) q[0];
sx q[0];
rz(-1.2905755) q[0];
sx q[0];
rz(-0.37877628) q[0];
rz(0.0060161034) q[1];
sx q[1];
rz(-2.6117987) q[1];
sx q[1];
rz(-2.5078497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089687448) q[0];
sx q[0];
rz(-1.9080592) q[0];
sx q[0];
rz(2.428672) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1442398) q[2];
sx q[2];
rz(-1.4359983) q[2];
sx q[2];
rz(-2.1950795) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72530109) q[1];
sx q[1];
rz(-2.0547199) q[1];
sx q[1];
rz(-2.121622) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2692436) q[3];
sx q[3];
rz(-0.81065882) q[3];
sx q[3];
rz(-1.8332421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70574957) q[2];
sx q[2];
rz(-1.5219995) q[2];
sx q[2];
rz(2.5283234) q[2];
rz(0.67522007) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(-0.81775445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262065) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(-2.9551103) q[1];
sx q[1];
rz(-1.3874556) q[1];
sx q[1];
rz(1.5127771) q[1];
rz(-0.45966799) q[2];
sx q[2];
rz(-1.5528233) q[2];
sx q[2];
rz(-0.96832392) q[2];
rz(-2.1929838) q[3];
sx q[3];
rz(-1.6453679) q[3];
sx q[3];
rz(-2.9352376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
