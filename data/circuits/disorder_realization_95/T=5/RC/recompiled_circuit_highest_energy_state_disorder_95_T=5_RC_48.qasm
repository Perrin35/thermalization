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
rz(-1.7669825) q[0];
sx q[0];
rz(-0.97281015) q[0];
sx q[0];
rz(-1.2306124) q[0];
rz(-1.542701) q[1];
sx q[1];
rz(-0.35294947) q[1];
sx q[1];
rz(-2.7065839) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9821138) q[0];
sx q[0];
rz(-1.4712988) q[0];
sx q[0];
rz(-1.9051108) q[0];
rz(-pi) q[1];
rz(-0.51911609) q[2];
sx q[2];
rz(-2.8808705) q[2];
sx q[2];
rz(0.97445011) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54625816) q[1];
sx q[1];
rz(-1.2915242) q[1];
sx q[1];
rz(-2.3550849) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1633662) q[3];
sx q[3];
rz(-3.1071783) q[3];
sx q[3];
rz(0.84190166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.9445442) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(2.22331) q[2];
rz(-2.7783172) q[3];
sx q[3];
rz(-1.2493635) q[3];
sx q[3];
rz(-1.3242807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0488224) q[0];
sx q[0];
rz(-2.4019882) q[0];
sx q[0];
rz(2.0059465) q[0];
rz(-0.28226918) q[1];
sx q[1];
rz(-1.5156563) q[1];
sx q[1];
rz(2.938882) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9720079) q[0];
sx q[0];
rz(-1.3915194) q[0];
sx q[0];
rz(2.9573134) q[0];
rz(-2.0524998) q[2];
sx q[2];
rz(-1.559245) q[2];
sx q[2];
rz(0.80977075) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7248279) q[1];
sx q[1];
rz(-1.4484805) q[1];
sx q[1];
rz(2.6847647) q[1];
rz(-1.560503) q[3];
sx q[3];
rz(-1.5760311) q[3];
sx q[3];
rz(-2.7017587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.89381605) q[2];
sx q[2];
rz(-2.717369) q[2];
sx q[2];
rz(3.0730548) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.8770542) q[3];
sx q[3];
rz(0.0059277047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0900367) q[0];
sx q[0];
rz(-0.55110252) q[0];
sx q[0];
rz(-2.0447482) q[0];
rz(2.5966273) q[1];
sx q[1];
rz(-2.4577591) q[1];
sx q[1];
rz(2.9473238) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8349154) q[0];
sx q[0];
rz(-0.55473548) q[0];
sx q[0];
rz(2.2771605) q[0];
rz(-pi) q[1];
rz(2.0805919) q[2];
sx q[2];
rz(-0.61697996) q[2];
sx q[2];
rz(1.4088907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.771505) q[1];
sx q[1];
rz(-1.9290826) q[1];
sx q[1];
rz(-0.53046988) q[1];
x q[2];
rz(-1.9220334) q[3];
sx q[3];
rz(-2.7673169) q[3];
sx q[3];
rz(2.5881899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7668309) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(-0.72804803) q[2];
rz(-0.22045615) q[3];
sx q[3];
rz(-2.9900592) q[3];
sx q[3];
rz(-2.5098586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24882889) q[0];
sx q[0];
rz(-1.6680102) q[0];
sx q[0];
rz(1.0695176) q[0];
rz(2.2546774) q[1];
sx q[1];
rz(-2.6353757) q[1];
sx q[1];
rz(-1.3217529) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89267545) q[0];
sx q[0];
rz(-1.5317316) q[0];
sx q[0];
rz(-0.43401735) q[0];
rz(3.07768) q[2];
sx q[2];
rz(-2.9021429) q[2];
sx q[2];
rz(-0.69635812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9001213) q[1];
sx q[1];
rz(-1.4992428) q[1];
sx q[1];
rz(-1.1318737) q[1];
rz(0.5170614) q[3];
sx q[3];
rz(-1.0365465) q[3];
sx q[3];
rz(-1.2599961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6679823) q[2];
sx q[2];
rz(-2.0986291) q[2];
sx q[2];
rz(-1.854151) q[2];
rz(2.1227664) q[3];
sx q[3];
rz(-0.5580709) q[3];
sx q[3];
rz(0.67600018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92893112) q[0];
sx q[0];
rz(-2.9086845) q[0];
sx q[0];
rz(-0.88053298) q[0];
rz(0.13678837) q[1];
sx q[1];
rz(-0.22577481) q[1];
sx q[1];
rz(2.5313964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5545512) q[0];
sx q[0];
rz(-2.4992538) q[0];
sx q[0];
rz(1.5917529) q[0];
x q[1];
rz(-1.4177104) q[2];
sx q[2];
rz(-1.1436074) q[2];
sx q[2];
rz(1.9356021) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1353947) q[1];
sx q[1];
rz(-1.2761226) q[1];
sx q[1];
rz(-1.7735405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.706677) q[3];
sx q[3];
rz(-2.0066094) q[3];
sx q[3];
rz(2.6282981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.60787624) q[2];
sx q[2];
rz(-0.72489649) q[2];
sx q[2];
rz(1.2302715) q[2];
rz(-2.2599334) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(1.3151883) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6864258) q[0];
sx q[0];
rz(-1.5979586) q[0];
sx q[0];
rz(2.0879188) q[0];
rz(1.7651438) q[1];
sx q[1];
rz(-2.0339298) q[1];
sx q[1];
rz(2.460316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72944869) q[0];
sx q[0];
rz(-0.870847) q[0];
sx q[0];
rz(-1.5493718) q[0];
rz(-pi) q[1];
rz(-0.7448911) q[2];
sx q[2];
rz(-1.8305664) q[2];
sx q[2];
rz(1.3448754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5743915) q[1];
sx q[1];
rz(-1.8387477) q[1];
sx q[1];
rz(1.059157) q[1];
x q[2];
rz(1.525649) q[3];
sx q[3];
rz(-0.69250007) q[3];
sx q[3];
rz(3.0255603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.072558746) q[2];
sx q[2];
rz(-2.7636038) q[2];
sx q[2];
rz(-2.6663713) q[2];
rz(-1.2230988) q[3];
sx q[3];
rz(-2.1971072) q[3];
sx q[3];
rz(-0.16330115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2326736) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(-0.21937823) q[0];
rz(-1.0064005) q[1];
sx q[1];
rz(-1.4234797) q[1];
sx q[1];
rz(1.9361608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4426098) q[0];
sx q[0];
rz(-0.97215334) q[0];
sx q[0];
rz(0.026431008) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7458887) q[2];
sx q[2];
rz(-1.8181516) q[2];
sx q[2];
rz(-2.5620714) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.58118966) q[1];
sx q[1];
rz(-2.9415574) q[1];
sx q[1];
rz(-0.086189857) q[1];
x q[2];
rz(1.1262283) q[3];
sx q[3];
rz(-1.6939112) q[3];
sx q[3];
rz(-0.11492226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9203303) q[2];
sx q[2];
rz(-0.39322501) q[2];
sx q[2];
rz(-2.4779251) q[2];
rz(2.4237733) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(-2.9653911) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32888907) q[0];
sx q[0];
rz(-2.1825574) q[0];
sx q[0];
rz(1.9004199) q[0];
rz(-0.42789704) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(1.4235628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20410886) q[0];
sx q[0];
rz(-0.23074575) q[0];
sx q[0];
rz(-0.33270128) q[0];
rz(0.15870416) q[2];
sx q[2];
rz(-2.2825512) q[2];
sx q[2];
rz(-1.2325328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9291722) q[1];
sx q[1];
rz(-2.4732686) q[1];
sx q[1];
rz(-1.6051488) q[1];
rz(-1.6070494) q[3];
sx q[3];
rz(-1.2360308) q[3];
sx q[3];
rz(-2.436155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1443892) q[2];
sx q[2];
rz(-0.68789566) q[2];
sx q[2];
rz(-0.29590657) q[2];
rz(-1.0812673) q[3];
sx q[3];
rz(-1.5238949) q[3];
sx q[3];
rz(2.9937939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57701552) q[0];
sx q[0];
rz(-0.58203375) q[0];
sx q[0];
rz(0.70593315) q[0];
rz(-2.9544746) q[1];
sx q[1];
rz(-0.83693224) q[1];
sx q[1];
rz(-2.6954938) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69639403) q[0];
sx q[0];
rz(-2.9893251) q[0];
sx q[0];
rz(3.1294786) q[0];
x q[1];
rz(-1.0995277) q[2];
sx q[2];
rz(-1.285811) q[2];
sx q[2];
rz(2.1113656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2526568) q[1];
sx q[1];
rz(-1.2351079) q[1];
sx q[1];
rz(-0.97519213) q[1];
rz(-1.3996077) q[3];
sx q[3];
rz(-1.0009196) q[3];
sx q[3];
rz(-2.8448454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17911653) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(-1.7402488) q[2];
rz(-2.5380747) q[3];
sx q[3];
rz(-0.18134914) q[3];
sx q[3];
rz(3.0028499) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1097581) q[0];
sx q[0];
rz(-1.4048445) q[0];
sx q[0];
rz(-3.0314714) q[0];
rz(-1.8203863) q[1];
sx q[1];
rz(-2.2167914) q[1];
sx q[1];
rz(-2.9168911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6455655) q[0];
sx q[0];
rz(-1.9651881) q[0];
sx q[0];
rz(1.0014707) q[0];
rz(-0.16188936) q[2];
sx q[2];
rz(-1.8542552) q[2];
sx q[2];
rz(3.0155011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70485348) q[1];
sx q[1];
rz(-2.3922643) q[1];
sx q[1];
rz(-2.8981461) q[1];
x q[2];
rz(1.4175156) q[3];
sx q[3];
rz(-2.9671921) q[3];
sx q[3];
rz(-2.0896951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8297742) q[2];
sx q[2];
rz(-2.8303787) q[2];
sx q[2];
rz(2.6149926) q[2];
rz(-2.7569568) q[3];
sx q[3];
rz(-1.8943818) q[3];
sx q[3];
rz(2.9068936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93152355) q[0];
sx q[0];
rz(-0.43545224) q[0];
sx q[0];
rz(-0.42238105) q[0];
rz(0.79413636) q[1];
sx q[1];
rz(-1.4310478) q[1];
sx q[1];
rz(-1.4337883) q[1];
rz(0.88243816) q[2];
sx q[2];
rz(-1.0461764) q[2];
sx q[2];
rz(-1.7027693) q[2];
rz(2.4173097) q[3];
sx q[3];
rz(-0.52762994) q[3];
sx q[3];
rz(2.4366196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
