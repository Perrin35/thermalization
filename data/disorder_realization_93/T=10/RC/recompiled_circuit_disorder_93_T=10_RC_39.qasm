OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(-2.2178712) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(-0.3737803) q[0];
rz(-1.853763) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(-1.9622918) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10780653) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(-0.71180196) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63763036) q[3];
sx q[3];
rz(-0.59083592) q[3];
sx q[3];
rz(-1.6319815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2279921) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-2.91495) q[0];
sx q[0];
rz(1.1967999) q[0];
rz(-0.67990047) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.4555567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128132) q[0];
sx q[0];
rz(-1.4038741) q[0];
sx q[0];
rz(2.1524327) q[0];
rz(-pi) q[1];
x q[1];
rz(2.935264) q[2];
sx q[2];
rz(-0.38197877) q[2];
sx q[2];
rz(2.1260335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.96869722) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(-1.5437267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060828408) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(0.5371679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75333726) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(-0.02877409) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.9690537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5867509) q[0];
sx q[0];
rz(-1.8790073) q[0];
sx q[0];
rz(0.59535938) q[0];
x q[1];
rz(0.75540382) q[2];
sx q[2];
rz(-1.3284151) q[2];
sx q[2];
rz(-0.6005477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6781569) q[1];
sx q[1];
rz(-2.0883745) q[1];
sx q[1];
rz(2.8606158) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7964732) q[3];
sx q[3];
rz(-0.25697069) q[3];
sx q[3];
rz(-2.8185609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0671493) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(0.91119901) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-3.0111713) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-2.6180843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2098171) q[0];
sx q[0];
rz(-0.99299252) q[0];
sx q[0];
rz(-2.0156167) q[0];
rz(-0.24387118) q[2];
sx q[2];
rz(-2.802231) q[2];
sx q[2];
rz(-1.7983758) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7586786) q[1];
sx q[1];
rz(-2.3944693) q[1];
sx q[1];
rz(1.1927356) q[1];
rz(-3.0047699) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(2.6329991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(1.4915285) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(2.1496444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0569699) q[0];
sx q[0];
rz(-1.3966494) q[0];
sx q[0];
rz(1.6790381) q[0];
x q[1];
rz(0.0089346272) q[2];
sx q[2];
rz(-1.8736471) q[2];
sx q[2];
rz(2.6645899) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79210287) q[1];
sx q[1];
rz(-1.8794685) q[1];
sx q[1];
rz(-1.5285138) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4335853) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(1.8828132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72702423) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093229175) q[0];
sx q[0];
rz(-1.3472124) q[0];
sx q[0];
rz(-2.6327052) q[0];
rz(-pi) q[1];
rz(1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(2.595682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52275601) q[1];
sx q[1];
rz(-2.6869876) q[1];
sx q[1];
rz(-1.418581) q[1];
rz(-pi) q[2];
rz(0.96958843) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(2.3525402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(-0.56387222) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5086223) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99188995) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(-0.71233149) q[0];
rz(2.2981811) q[2];
sx q[2];
rz(-1.9551829) q[2];
sx q[2];
rz(-0.20993983) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2591178) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(0.74101733) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52845593) q[3];
sx q[3];
rz(-0.65675694) q[3];
sx q[3];
rz(-2.9627851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(2.896893) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.6285508) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(-2.8322463) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.3051422) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37874052) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(2.5255894) q[0];
x q[1];
rz(0.70456409) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(-1.9412083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95272428) q[1];
sx q[1];
rz(-0.75718588) q[1];
sx q[1];
rz(-0.56307478) q[1];
rz(-1.5042138) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(-0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(0.70294356) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2962869) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(-1.9592497) q[0];
rz(2.0800955) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(0.73355567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8150836) q[1];
sx q[1];
rz(-2.106296) q[1];
sx q[1];
rz(0.19170796) q[1];
x q[2];
rz(-0.33396696) q[3];
sx q[3];
rz(-0.98200646) q[3];
sx q[3];
rz(0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(1.8019603) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(-2.6729565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2221453) q[0];
sx q[0];
rz(-1.7788017) q[0];
sx q[0];
rz(2.749445) q[0];
x q[1];
rz(-2.8154545) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(1.067576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74433078) q[1];
sx q[1];
rz(-1.3052193) q[1];
sx q[1];
rz(0.070752146) q[1];
x q[2];
rz(-0.81810276) q[3];
sx q[3];
rz(-2.5460498) q[3];
sx q[3];
rz(-2.9593352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(-1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(2.519683) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(2.2767699) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(-2.2223496) q[3];
sx q[3];
rz(-1.6332492) q[3];
sx q[3];
rz(1.9132683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
