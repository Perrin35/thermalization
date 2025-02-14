OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2562113) q[0];
sx q[0];
rz(-0.71330944) q[0];
sx q[0];
rz(-1.2137086) q[0];
rz(-1.3656536) q[1];
sx q[1];
rz(-2.8627099) q[1];
sx q[1];
rz(-0.20067781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047223481) q[0];
sx q[0];
rz(-2.6728164) q[0];
sx q[0];
rz(1.7550521) q[0];
x q[1];
rz(-0.21733445) q[2];
sx q[2];
rz(-1.4121018) q[2];
sx q[2];
rz(2.8087698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38814784) q[1];
sx q[1];
rz(-0.26784836) q[1];
sx q[1];
rz(-0.74117383) q[1];
rz(-2.8332769) q[3];
sx q[3];
rz(-1.7349958) q[3];
sx q[3];
rz(0.37946821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.063244907) q[2];
sx q[2];
rz(-1.1996317) q[2];
sx q[2];
rz(0.71301785) q[2];
rz(1.2837422) q[3];
sx q[3];
rz(-2.4132437) q[3];
sx q[3];
rz(-0.89948765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17083183) q[0];
sx q[0];
rz(-1.4749227) q[0];
sx q[0];
rz(1.4738039) q[0];
rz(-2.5693192) q[1];
sx q[1];
rz(-1.0621366) q[1];
sx q[1];
rz(-1.7639814) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.361918) q[0];
sx q[0];
rz(-2.5847844) q[0];
sx q[0];
rz(2.1614055) q[0];
rz(-pi) q[1];
rz(2.4430635) q[2];
sx q[2];
rz(-2.4870092) q[2];
sx q[2];
rz(1.6561001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8793844) q[1];
sx q[1];
rz(-1.2405433) q[1];
sx q[1];
rz(0.31008215) q[1];
x q[2];
rz(-0.99239852) q[3];
sx q[3];
rz(-0.36250329) q[3];
sx q[3];
rz(0.9622919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1231692) q[2];
sx q[2];
rz(-1.6320684) q[2];
sx q[2];
rz(1.8040166) q[2];
rz(-2.7959974) q[3];
sx q[3];
rz(-2.5497422) q[3];
sx q[3];
rz(0.93446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086394101) q[0];
sx q[0];
rz(-0.18562695) q[0];
sx q[0];
rz(0.58919543) q[0];
rz(-0.20801726) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(-0.20588188) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97144868) q[0];
sx q[0];
rz(-1.9238233) q[0];
sx q[0];
rz(-2.5021509) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85742204) q[2];
sx q[2];
rz(-1.9586342) q[2];
sx q[2];
rz(-1.1042386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2370473) q[1];
sx q[1];
rz(-0.77213192) q[1];
sx q[1];
rz(0.44293483) q[1];
x q[2];
rz(-3.053756) q[3];
sx q[3];
rz(-1.1655775) q[3];
sx q[3];
rz(-2.3072568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0899293) q[2];
sx q[2];
rz(-2.0970924) q[2];
sx q[2];
rz(0.54452407) q[2];
rz(2.6868467) q[3];
sx q[3];
rz(-0.62362042) q[3];
sx q[3];
rz(0.0039984306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398448) q[0];
sx q[0];
rz(-2.5938617) q[0];
sx q[0];
rz(-2.4416583) q[0];
rz(1.3088538) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(1.4189789) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0467906) q[0];
sx q[0];
rz(-1.7715469) q[0];
sx q[0];
rz(-1.0958452) q[0];
rz(3.1221462) q[2];
sx q[2];
rz(-2.0692973) q[2];
sx q[2];
rz(-2.4011103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8992865) q[1];
sx q[1];
rz(-1.810009) q[1];
sx q[1];
rz(0.71372791) q[1];
rz(2.9955825) q[3];
sx q[3];
rz(-1.2577122) q[3];
sx q[3];
rz(2.2816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2122638) q[2];
sx q[2];
rz(-0.45390359) q[2];
sx q[2];
rz(1.7960499) q[2];
rz(-0.71873194) q[3];
sx q[3];
rz(-2.4001382) q[3];
sx q[3];
rz(-2.45347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6253925) q[0];
sx q[0];
rz(-1.1855519) q[0];
sx q[0];
rz(-0.54310435) q[0];
rz(0.25554666) q[1];
sx q[1];
rz(-2.6592022) q[1];
sx q[1];
rz(-1.3822752) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7092752) q[0];
sx q[0];
rz(-2.0943644) q[0];
sx q[0];
rz(2.0366336) q[0];
rz(-1.7167983) q[2];
sx q[2];
rz(-0.80933648) q[2];
sx q[2];
rz(-0.80235624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.71798793) q[1];
sx q[1];
rz(-1.8306872) q[1];
sx q[1];
rz(2.1433448) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38726728) q[3];
sx q[3];
rz(-0.42902374) q[3];
sx q[3];
rz(-2.9128592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.9690659) q[2];
sx q[2];
rz(-1.4719937) q[2];
sx q[2];
rz(0.36953163) q[2];
rz(-1.3299804) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(1.7867521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1177649) q[0];
sx q[0];
rz(-0.33482877) q[0];
sx q[0];
rz(0.086300015) q[0];
rz(-1.6465126) q[1];
sx q[1];
rz(-1.9638655) q[1];
sx q[1];
rz(0.42974791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54329007) q[0];
sx q[0];
rz(-1.5006646) q[0];
sx q[0];
rz(-1.7599517) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8977857) q[2];
sx q[2];
rz(-2.2026416) q[2];
sx q[2];
rz(1.6666842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4757929) q[1];
sx q[1];
rz(-2.0415039) q[1];
sx q[1];
rz(-1.0757955) q[1];
x q[2];
rz(0.67062258) q[3];
sx q[3];
rz(-1.792893) q[3];
sx q[3];
rz(-2.5012453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0325844) q[2];
sx q[2];
rz(-1.0823559) q[2];
sx q[2];
rz(1.2004131) q[2];
rz(-0.16610185) q[3];
sx q[3];
rz(-2.3717272) q[3];
sx q[3];
rz(0.86336819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85668856) q[0];
sx q[0];
rz(-2.3785474) q[0];
sx q[0];
rz(0.82409182) q[0];
rz(1.2245945) q[1];
sx q[1];
rz(-1.2242182) q[1];
sx q[1];
rz(1.0138938) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4606095) q[0];
sx q[0];
rz(-1.8797726) q[0];
sx q[0];
rz(2.3302484) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61182558) q[2];
sx q[2];
rz(-2.6227747) q[2];
sx q[2];
rz(3.0957565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6142004) q[1];
sx q[1];
rz(-2.1062615) q[1];
sx q[1];
rz(-1.1170858) q[1];
x q[2];
rz(0.39172642) q[3];
sx q[3];
rz(-2.0731032) q[3];
sx q[3];
rz(-1.6592178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8635204) q[2];
sx q[2];
rz(-2.4188228) q[2];
sx q[2];
rz(-1.8966759) q[2];
rz(0.23166367) q[3];
sx q[3];
rz(-1.04117) q[3];
sx q[3];
rz(0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046655) q[0];
sx q[0];
rz(-1.2060839) q[0];
sx q[0];
rz(-0.99223247) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-2.2997586) q[1];
sx q[1];
rz(-1.7281035) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.468576) q[0];
sx q[0];
rz(-1.6454738) q[0];
sx q[0];
rz(-0.59665307) q[0];
x q[1];
rz(-2.0199168) q[2];
sx q[2];
rz(-2.5418021) q[2];
sx q[2];
rz(3.1315294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7602049) q[1];
sx q[1];
rz(-1.5389812) q[1];
sx q[1];
rz(0.29812584) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9416789) q[3];
sx q[3];
rz(-2.6318916) q[3];
sx q[3];
rz(2.5751757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32470545) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(1.7186349) q[2];
rz(-3.0286246) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(0.62937361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9519575) q[0];
sx q[0];
rz(-1.676214) q[0];
sx q[0];
rz(2.5919609) q[0];
rz(-0.59898392) q[1];
sx q[1];
rz(-2.2689029) q[1];
sx q[1];
rz(-1.5113066) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0534759) q[0];
sx q[0];
rz(-1.5742745) q[0];
sx q[0];
rz(0.17668488) q[0];
rz(1.8725996) q[2];
sx q[2];
rz(-1.0827218) q[2];
sx q[2];
rz(-2.7583964) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8558939) q[1];
sx q[1];
rz(-1.0700858) q[1];
sx q[1];
rz(-2.6803174) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8776607) q[3];
sx q[3];
rz(-2.0425386) q[3];
sx q[3];
rz(2.2795171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0281684) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(-2.8070731) q[2];
rz(2.4962375) q[3];
sx q[3];
rz(-1.0048451) q[3];
sx q[3];
rz(-2.9042802) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065011218) q[0];
sx q[0];
rz(-0.33877057) q[0];
sx q[0];
rz(2.5184799) q[0];
rz(0.30162853) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(-1.041144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469676) q[0];
sx q[0];
rz(-2.0026221) q[0];
sx q[0];
rz(0.65460848) q[0];
rz(-pi) q[1];
rz(-2.9567962) q[2];
sx q[2];
rz(-2.6738538) q[2];
sx q[2];
rz(-0.50450215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54339441) q[1];
sx q[1];
rz(-2.3375727) q[1];
sx q[1];
rz(-1.5280185) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0431533) q[3];
sx q[3];
rz(-0.40723793) q[3];
sx q[3];
rz(2.0472297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69680301) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(-0.94318548) q[2];
rz(-1.7275564) q[3];
sx q[3];
rz(-2.2351738) q[3];
sx q[3];
rz(2.3448155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.74054756) q[0];
sx q[0];
rz(-0.39401207) q[0];
sx q[0];
rz(3.0643585) q[0];
rz(-0.41863353) q[1];
sx q[1];
rz(-1.5593465) q[1];
sx q[1];
rz(-2.9895463) q[1];
rz(-0.05337333) q[2];
sx q[2];
rz(-2.5659701) q[2];
sx q[2];
rz(0.26502668) q[2];
rz(2.7100415) q[3];
sx q[3];
rz(-1.3234371) q[3];
sx q[3];
rz(0.32061843) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
