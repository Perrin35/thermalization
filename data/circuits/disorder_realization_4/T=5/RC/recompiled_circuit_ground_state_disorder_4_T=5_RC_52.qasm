OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.088012785) q[0];
sx q[0];
rz(-2.8289284) q[0];
sx q[0];
rz(0.13392681) q[0];
rz(3.1349831) q[1];
sx q[1];
rz(-1.5899038) q[1];
sx q[1];
rz(-0.42491999) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8613409) q[0];
sx q[0];
rz(-1.6854648) q[0];
sx q[0];
rz(0.10114177) q[0];
rz(1.0966461) q[2];
sx q[2];
rz(-0.72410781) q[2];
sx q[2];
rz(-0.5257789) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1462789) q[1];
sx q[1];
rz(-1.2321207) q[1];
sx q[1];
rz(-0.92471735) q[1];
x q[2];
rz(2.4185804) q[3];
sx q[3];
rz(-2.0916846) q[3];
sx q[3];
rz(-0.82244825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6801844) q[2];
sx q[2];
rz(-0.37698656) q[2];
sx q[2];
rz(0.11412966) q[2];
rz(2.7528609) q[3];
sx q[3];
rz(-2.8753493) q[3];
sx q[3];
rz(-1.3346599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0085793063) q[0];
sx q[0];
rz(-2.7393434) q[0];
sx q[0];
rz(0.026799686) q[0];
rz(-1.9642448) q[1];
sx q[1];
rz(-2.4939996) q[1];
sx q[1];
rz(3.1039544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55521727) q[0];
sx q[0];
rz(-2.3669985) q[0];
sx q[0];
rz(-0.019217773) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9178647) q[2];
sx q[2];
rz(-1.7729974) q[2];
sx q[2];
rz(3.1333609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8839542) q[1];
sx q[1];
rz(-2.3515793) q[1];
sx q[1];
rz(-0.86261065) q[1];
rz(1.8438128) q[3];
sx q[3];
rz(-2.3721266) q[3];
sx q[3];
rz(-0.96801341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2521952) q[2];
sx q[2];
rz(-0.9844206) q[2];
sx q[2];
rz(2.9241015) q[2];
rz(-0.43976954) q[3];
sx q[3];
rz(-1.7871126) q[3];
sx q[3];
rz(3.0396438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6169823) q[0];
sx q[0];
rz(-3.1341902) q[0];
sx q[0];
rz(-2.5819085) q[0];
rz(-2.2451378) q[1];
sx q[1];
rz(-0.47210109) q[1];
sx q[1];
rz(-0.40759531) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4993138) q[0];
sx q[0];
rz(-0.98897213) q[0];
sx q[0];
rz(-3.0034742) q[0];
rz(1.4846663) q[2];
sx q[2];
rz(-0.67170947) q[2];
sx q[2];
rz(1.9957146) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.421459) q[1];
sx q[1];
rz(-1.4602666) q[1];
sx q[1];
rz(0.46760749) q[1];
rz(-2.8333307) q[3];
sx q[3];
rz(-2.097766) q[3];
sx q[3];
rz(1.4013724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1233623) q[2];
sx q[2];
rz(-1.9841649) q[2];
sx q[2];
rz(2.1165712) q[2];
rz(1.7068663) q[3];
sx q[3];
rz(-0.44308174) q[3];
sx q[3];
rz(-0.6492492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0822815) q[0];
sx q[0];
rz(-2.870443) q[0];
sx q[0];
rz(-0.49736381) q[0];
rz(1.8479895) q[1];
sx q[1];
rz(-2.135364) q[1];
sx q[1];
rz(-1.2327548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7347438) q[0];
sx q[0];
rz(-0.60601425) q[0];
sx q[0];
rz(0.36263175) q[0];
rz(-1.9903988) q[2];
sx q[2];
rz(-2.2299354) q[2];
sx q[2];
rz(-1.4105547) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9962483) q[1];
sx q[1];
rz(-1.323902) q[1];
sx q[1];
rz(-1.3614142) q[1];
rz(-0.044305459) q[3];
sx q[3];
rz(-2.7409275) q[3];
sx q[3];
rz(1.4728427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.57700145) q[2];
sx q[2];
rz(-0.61508721) q[2];
sx q[2];
rz(2.9992529) q[2];
rz(-1.9650991) q[3];
sx q[3];
rz(-1.0005955) q[3];
sx q[3];
rz(-2.71079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.38554) q[0];
sx q[0];
rz(-0.98243326) q[0];
sx q[0];
rz(-0.2734215) q[0];
rz(-0.39971071) q[1];
sx q[1];
rz(-0.81396657) q[1];
sx q[1];
rz(0.25693691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77784044) q[0];
sx q[0];
rz(-1.8926601) q[0];
sx q[0];
rz(2.1296164) q[0];
rz(0.24946282) q[2];
sx q[2];
rz(-1.9469065) q[2];
sx q[2];
rz(-2.1975236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.6338004) q[1];
sx q[1];
rz(-1.6163428) q[1];
sx q[1];
rz(1.5006658) q[1];
rz(-0.2843379) q[3];
sx q[3];
rz(-0.91129843) q[3];
sx q[3];
rz(-2.4267765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5520681) q[2];
sx q[2];
rz(-2.7858211) q[2];
sx q[2];
rz(0.75999981) q[2];
rz(0.99203569) q[3];
sx q[3];
rz(-1.2090679) q[3];
sx q[3];
rz(-1.1442643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511667) q[0];
sx q[0];
rz(-2.5100584) q[0];
sx q[0];
rz(-2.5354711) q[0];
rz(2.0947314) q[1];
sx q[1];
rz(-1.650834) q[1];
sx q[1];
rz(3.0643588) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55750269) q[0];
sx q[0];
rz(-1.5582005) q[0];
sx q[0];
rz(-1.5303311) q[0];
x q[1];
rz(0.31862835) q[2];
sx q[2];
rz(-1.1201522) q[2];
sx q[2];
rz(0.81540996) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97925767) q[1];
sx q[1];
rz(-1.1282776) q[1];
sx q[1];
rz(-2.6975432) q[1];
rz(-pi) q[2];
rz(-1.3807543) q[3];
sx q[3];
rz(-2.2249208) q[3];
sx q[3];
rz(-0.074617059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91654009) q[2];
sx q[2];
rz(-2.4317135) q[2];
sx q[2];
rz(0.27099657) q[2];
rz(0.58800507) q[3];
sx q[3];
rz(-2.3376412) q[3];
sx q[3];
rz(-0.40905455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005017) q[0];
sx q[0];
rz(-1.5301457) q[0];
sx q[0];
rz(-2.8588168) q[0];
rz(0.68280363) q[1];
sx q[1];
rz(-0.86137259) q[1];
sx q[1];
rz(-0.85321325) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3798987) q[0];
sx q[0];
rz(-2.324484) q[0];
sx q[0];
rz(-1.0495814) q[0];
rz(-pi) q[1];
rz(-0.81099895) q[2];
sx q[2];
rz(-1.2110707) q[2];
sx q[2];
rz(-3.1364417) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0921923) q[1];
sx q[1];
rz(-2.8554648) q[1];
sx q[1];
rz(-0.95502468) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0045347) q[3];
sx q[3];
rz(-0.2314724) q[3];
sx q[3];
rz(1.2616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0679438) q[2];
sx q[2];
rz(-0.51764071) q[2];
sx q[2];
rz(2.850387) q[2];
rz(1.2047042) q[3];
sx q[3];
rz(-2.5681345) q[3];
sx q[3];
rz(-1.5706435) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35029995) q[0];
sx q[0];
rz(-1.5234103) q[0];
sx q[0];
rz(-2.5147901) q[0];
rz(-2.0794012) q[1];
sx q[1];
rz(-2.3419582) q[1];
sx q[1];
rz(-0.83745426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0916679) q[0];
sx q[0];
rz(-0.44418884) q[0];
sx q[0];
rz(-1.8416406) q[0];
rz(-pi) q[1];
rz(-1.8071354) q[2];
sx q[2];
rz(-1.9009942) q[2];
sx q[2];
rz(-2.8446252) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1143226) q[1];
sx q[1];
rz(-2.6114846) q[1];
sx q[1];
rz(1.8370173) q[1];
rz(-pi) q[2];
rz(0.25691311) q[3];
sx q[3];
rz(-0.532386) q[3];
sx q[3];
rz(0.40125602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.31967878) q[2];
sx q[2];
rz(-1.9210812) q[2];
sx q[2];
rz(1.0411881) q[2];
rz(-2.239481) q[3];
sx q[3];
rz(-1.9769042) q[3];
sx q[3];
rz(-0.56727099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5015471) q[0];
sx q[0];
rz(-0.29842672) q[0];
sx q[0];
rz(0.10064594) q[0];
rz(1.3238662) q[1];
sx q[1];
rz(-2.3034425) q[1];
sx q[1];
rz(-2.5460338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32257358) q[0];
sx q[0];
rz(-0.77054502) q[0];
sx q[0];
rz(2.3014803) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8331489) q[2];
sx q[2];
rz(-0.67910087) q[2];
sx q[2];
rz(-0.51788143) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7337949) q[1];
sx q[1];
rz(-2.3721954) q[1];
sx q[1];
rz(2.7831538) q[1];
rz(-pi) q[2];
rz(-0.15138126) q[3];
sx q[3];
rz(-0.80946846) q[3];
sx q[3];
rz(0.71914205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9969534) q[2];
sx q[2];
rz(-0.75012952) q[2];
sx q[2];
rz(1.3422802) q[2];
rz(-0.73623776) q[3];
sx q[3];
rz(-2.8056371) q[3];
sx q[3];
rz(0.098522447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017224273) q[0];
sx q[0];
rz(-2.4898744) q[0];
sx q[0];
rz(-0.29618725) q[0];
rz(1.4172957) q[1];
sx q[1];
rz(-2.2139151) q[1];
sx q[1];
rz(0.16253026) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1862193) q[0];
sx q[0];
rz(-0.027039921) q[0];
sx q[0];
rz(1.3229738) q[0];
rz(-pi) q[1];
rz(2.7835566) q[2];
sx q[2];
rz(-1.0834143) q[2];
sx q[2];
rz(-2.0646937) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0386943) q[1];
sx q[1];
rz(-1.7084645) q[1];
sx q[1];
rz(0.18935151) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82270427) q[3];
sx q[3];
rz(-0.71559042) q[3];
sx q[3];
rz(2.0377318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51809239) q[2];
sx q[2];
rz(-1.2178428) q[2];
sx q[2];
rz(-0.45563844) q[2];
rz(-1.0738922) q[3];
sx q[3];
rz(-2.5116601) q[3];
sx q[3];
rz(-2.5560801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0475273) q[0];
sx q[0];
rz(-1.8210664) q[0];
sx q[0];
rz(-0.6022712) q[0];
rz(2.8216254) q[1];
sx q[1];
rz(-2.0260369) q[1];
sx q[1];
rz(1.8021348) q[1];
rz(0.16945571) q[2];
sx q[2];
rz(-2.8002938) q[2];
sx q[2];
rz(2.4677966) q[2];
rz(-3.0319936) q[3];
sx q[3];
rz(-2.3634956) q[3];
sx q[3];
rz(2.6131438) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
