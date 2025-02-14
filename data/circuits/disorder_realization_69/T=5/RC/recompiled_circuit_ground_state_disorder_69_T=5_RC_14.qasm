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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6883301) q[0];
sx q[0];
rz(-1.6536667) q[0];
sx q[0];
rz(2.0327264) q[0];
rz(0.21733445) q[2];
sx q[2];
rz(-1.4121018) q[2];
sx q[2];
rz(0.33282285) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9057718) q[1];
sx q[1];
rz(-1.3911472) q[1];
sx q[1];
rz(-0.19975042) q[1];
rz(-pi) q[2];
rz(-1.7429656) q[3];
sx q[3];
rz(-1.8748314) q[3];
sx q[3];
rz(-1.1393169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0783477) q[2];
sx q[2];
rz(-1.941961) q[2];
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
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17083183) q[0];
sx q[0];
rz(-1.66667) q[0];
sx q[0];
rz(1.4738039) q[0];
rz(-2.5693192) q[1];
sx q[1];
rz(-1.0621366) q[1];
sx q[1];
rz(1.3776113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3085605) q[0];
sx q[0];
rz(-1.8695117) q[0];
sx q[0];
rz(-1.0935944) q[0];
rz(2.6102561) q[2];
sx q[2];
rz(-1.9730933) q[2];
sx q[2];
rz(-0.50237331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26220825) q[1];
sx q[1];
rz(-1.2405433) q[1];
sx q[1];
rz(0.31008215) q[1];
rz(-2.9371524) q[3];
sx q[3];
rz(-1.8722765) q[3];
sx q[3];
rz(-2.7888576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.01842347) q[2];
sx q[2];
rz(-1.5095242) q[2];
sx q[2];
rz(-1.337576) q[2];
rz(0.34559524) q[3];
sx q[3];
rz(-0.59185043) q[3];
sx q[3];
rz(-0.93446294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0551986) q[0];
sx q[0];
rz(-2.9559657) q[0];
sx q[0];
rz(2.5523972) q[0];
rz(2.9335754) q[1];
sx q[1];
rz(-0.55744019) q[1];
sx q[1];
rz(0.20588188) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.170144) q[0];
sx q[0];
rz(-1.9238233) q[0];
sx q[0];
rz(-2.5021509) q[0];
rz(-pi) q[1];
rz(-0.85742204) q[2];
sx q[2];
rz(-1.1829585) q[2];
sx q[2];
rz(-1.1042386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9045453) q[1];
sx q[1];
rz(-0.77213192) q[1];
sx q[1];
rz(-2.6986578) q[1];
rz(-pi) q[2];
rz(1.7725189) q[3];
sx q[3];
rz(-2.7274787) q[3];
sx q[3];
rz(0.61455807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0899293) q[2];
sx q[2];
rz(-2.0970924) q[2];
sx q[2];
rz(2.5970686) q[2];
rz(2.6868467) q[3];
sx q[3];
rz(-0.62362042) q[3];
sx q[3];
rz(-3.1375942) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74314463) q[0];
sx q[0];
rz(-2.5938617) q[0];
sx q[0];
rz(-2.4416583) q[0];
rz(1.3088538) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(-1.7226137) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5154079) q[0];
sx q[0];
rz(-1.1061449) q[0];
sx q[0];
rz(0.22494577) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5350902) q[2];
sx q[2];
rz(-0.49884819) q[2];
sx q[2];
rz(0.6998261) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1261038) q[1];
sx q[1];
rz(-0.88147336) q[1];
sx q[1];
rz(1.2587121) q[1];
rz(-pi) q[2];
rz(-1.1484234) q[3];
sx q[3];
rz(-2.7971533) q[3];
sx q[3];
rz(0.41448739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2122638) q[2];
sx q[2];
rz(-2.6876891) q[2];
sx q[2];
rz(1.7960499) q[2];
rz(-2.4228607) q[3];
sx q[3];
rz(-2.4001382) q[3];
sx q[3];
rz(-0.68812266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5162002) q[0];
sx q[0];
rz(-1.1855519) q[0];
sx q[0];
rz(-2.5984883) q[0];
rz(-2.886046) q[1];
sx q[1];
rz(-2.6592022) q[1];
sx q[1];
rz(-1.3822752) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2204748) q[0];
sx q[0];
rz(-0.68604031) q[0];
sx q[0];
rz(-0.6612079) q[0];
x q[1];
rz(-1.7167983) q[2];
sx q[2];
rz(-0.80933648) q[2];
sx q[2];
rz(2.3392364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68866036) q[1];
sx q[1];
rz(-2.1218461) q[1];
sx q[1];
rz(0.30639415) q[1];
rz(-pi) q[2];
rz(2.7409389) q[3];
sx q[3];
rz(-1.4130428) q[3];
sx q[3];
rz(2.1547013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1725267) q[2];
sx q[2];
rz(-1.669599) q[2];
sx q[2];
rz(-2.772061) q[2];
rz(1.8116123) q[3];
sx q[3];
rz(-0.79052916) q[3];
sx q[3];
rz(1.3548405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1177649) q[0];
sx q[0];
rz(-2.8067639) q[0];
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
rz(2.1275009) q[0];
sx q[0];
rz(-1.7594811) q[0];
sx q[0];
rz(3.0701915) q[0];
rz(-pi) q[1];
rz(1.8977857) q[2];
sx q[2];
rz(-2.2026416) q[2];
sx q[2];
rz(1.4749084) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33510241) q[1];
sx q[1];
rz(-1.1336328) q[1];
sx q[1];
rz(-2.6173068) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4709701) q[3];
sx q[3];
rz(-1.792893) q[3];
sx q[3];
rz(2.5012453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0325844) q[2];
sx q[2];
rz(-1.0823559) q[2];
sx q[2];
rz(1.2004131) q[2];
rz(-2.9754908) q[3];
sx q[3];
rz(-0.76986543) q[3];
sx q[3];
rz(-2.2782245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85668856) q[0];
sx q[0];
rz(-0.76304522) q[0];
sx q[0];
rz(-2.3175008) q[0];
rz(1.2245945) q[1];
sx q[1];
rz(-1.2242182) q[1];
sx q[1];
rz(1.0138938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4606095) q[0];
sx q[0];
rz(-1.8797726) q[0];
sx q[0];
rz(2.3302484) q[0];
rz(-pi) q[1];
rz(2.5297671) q[2];
sx q[2];
rz(-2.6227747) q[2];
sx q[2];
rz(0.045836115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9062324) q[1];
sx q[1];
rz(-0.68720931) q[1];
sx q[1];
rz(-2.5053124) q[1];
x q[2];
rz(-0.9634094) q[3];
sx q[3];
rz(-2.5150891) q[3];
sx q[3];
rz(0.95010786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8635204) q[2];
sx q[2];
rz(-2.4188228) q[2];
sx q[2];
rz(1.2449167) q[2];
rz(2.909929) q[3];
sx q[3];
rz(-1.04117) q[3];
sx q[3];
rz(-0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.046655) q[0];
sx q[0];
rz(-1.9355087) q[0];
sx q[0];
rz(-0.99223247) q[0];
rz(-0.16920432) q[1];
sx q[1];
rz(-2.2997586) q[1];
sx q[1];
rz(-1.7281035) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6730167) q[0];
sx q[0];
rz(-1.6454738) q[0];
sx q[0];
rz(-0.59665307) q[0];
rz(2.0199168) q[2];
sx q[2];
rz(-2.5418021) q[2];
sx q[2];
rz(-3.1315294) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0553356) q[1];
sx q[1];
rz(-0.29976832) q[1];
sx q[1];
rz(-3.0336628) q[1];
rz(-1.9416789) q[3];
sx q[3];
rz(-2.6318916) q[3];
sx q[3];
rz(2.5751757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8168872) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(-1.4229577) q[2];
rz(-3.0286246) q[3];
sx q[3];
rz(-2.8739417) q[3];
sx q[3];
rz(2.512219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1896352) q[0];
sx q[0];
rz(-1.4653787) q[0];
sx q[0];
rz(2.5919609) q[0];
rz(-2.5426087) q[1];
sx q[1];
rz(-0.87268972) q[1];
sx q[1];
rz(-1.5113066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6236512) q[0];
sx q[0];
rz(-1.7474801) q[0];
sx q[0];
rz(-1.5743295) q[0];
x q[1];
rz(-0.51038607) q[2];
sx q[2];
rz(-0.56737075) q[2];
sx q[2];
rz(-0.20287831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0508742) q[1];
sx q[1];
rz(-1.9719187) q[1];
sx q[1];
rz(1.0222597) q[1];
rz(-1.8776607) q[3];
sx q[3];
rz(-2.0425386) q[3];
sx q[3];
rz(0.8620756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1134243) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(0.33451954) q[2];
rz(-2.4962375) q[3];
sx q[3];
rz(-1.0048451) q[3];
sx q[3];
rz(2.9042802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0765814) q[0];
sx q[0];
rz(-2.8028221) q[0];
sx q[0];
rz(2.5184799) q[0];
rz(-0.30162853) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(-2.1004486) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469676) q[0];
sx q[0];
rz(-2.0026221) q[0];
sx q[0];
rz(2.4869842) q[0];
rz(-1.663346) q[2];
sx q[2];
rz(-2.0299533) q[2];
sx q[2];
rz(-0.29806229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5981982) q[1];
sx q[1];
rz(-2.3375727) q[1];
sx q[1];
rz(1.5280185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9275366) q[3];
sx q[3];
rz(-1.7715653) q[3];
sx q[3];
rz(0.014895766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69680301) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(0.94318548) q[2];
rz(1.7275564) q[3];
sx q[3];
rz(-2.2351738) q[3];
sx q[3];
rz(0.79677719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4010451) q[0];
sx q[0];
rz(-0.39401207) q[0];
sx q[0];
rz(3.0643585) q[0];
rz(-2.7229591) q[1];
sx q[1];
rz(-1.5822462) q[1];
sx q[1];
rz(0.15204631) q[1];
rz(-0.05337333) q[2];
sx q[2];
rz(-2.5659701) q[2];
sx q[2];
rz(0.26502668) q[2];
rz(-0.54316212) q[3];
sx q[3];
rz(-0.49351963) q[3];
sx q[3];
rz(-1.7388572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
