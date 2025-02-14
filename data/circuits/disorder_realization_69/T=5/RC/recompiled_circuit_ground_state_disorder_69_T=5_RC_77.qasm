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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0943692) q[0];
sx q[0];
rz(-0.46877623) q[0];
sx q[0];
rz(1.3865406) q[0];
rz(1.7332478) q[2];
sx q[2];
rz(-1.7853569) q[2];
sx q[2];
rz(-1.2728557) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38814784) q[1];
sx q[1];
rz(-0.26784836) q[1];
sx q[1];
rz(2.4004188) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6418125) q[3];
sx q[3];
rz(-2.793514) q[3];
sx q[3];
rz(1.4760555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0783477) q[2];
sx q[2];
rz(-1.941961) q[2];
sx q[2];
rz(-0.71301785) q[2];
rz(-1.8578505) q[3];
sx q[3];
rz(-2.4132437) q[3];
sx q[3];
rz(2.242105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9707608) q[0];
sx q[0];
rz(-1.66667) q[0];
sx q[0];
rz(1.4738039) q[0];
rz(-0.57227349) q[1];
sx q[1];
rz(-2.0794561) q[1];
sx q[1];
rz(-1.7639814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77967465) q[0];
sx q[0];
rz(-2.5847844) q[0];
sx q[0];
rz(2.1614055) q[0];
rz(-pi) q[1];
rz(0.53133659) q[2];
sx q[2];
rz(-1.1684993) q[2];
sx q[2];
rz(2.6392193) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8793844) q[1];
sx q[1];
rz(-1.2405433) q[1];
sx q[1];
rz(-2.8315105) q[1];
rz(-2.1491941) q[3];
sx q[3];
rz(-2.7790894) q[3];
sx q[3];
rz(-2.1793008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1231692) q[2];
sx q[2];
rz(-1.5095242) q[2];
sx q[2];
rz(-1.337576) q[2];
rz(-2.7959974) q[3];
sx q[3];
rz(-2.5497422) q[3];
sx q[3];
rz(-2.2071297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551986) q[0];
sx q[0];
rz(-0.18562695) q[0];
sx q[0];
rz(-2.5523972) q[0];
rz(-2.9335754) q[1];
sx q[1];
rz(-2.5841525) q[1];
sx q[1];
rz(0.20588188) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7939121) q[0];
sx q[0];
rz(-2.1651175) q[0];
sx q[0];
rz(-1.1403313) q[0];
rz(-pi) q[1];
rz(0.49534246) q[2];
sx q[2];
rz(-2.2215507) q[2];
sx q[2];
rz(0.15025727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65218788) q[1];
sx q[1];
rz(-2.2527848) q[1];
sx q[1];
rz(-1.1754065) q[1];
rz(-1.7725189) q[3];
sx q[3];
rz(-0.414114) q[3];
sx q[3];
rz(0.61455807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0899293) q[2];
sx q[2];
rz(-2.0970924) q[2];
sx q[2];
rz(-0.54452407) q[2];
rz(2.6868467) q[3];
sx q[3];
rz(-2.5179722) q[3];
sx q[3];
rz(-0.0039984306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74314463) q[0];
sx q[0];
rz(-0.54773098) q[0];
sx q[0];
rz(-2.4416583) q[0];
rz(1.3088538) q[1];
sx q[1];
rz(-2.0715641) q[1];
sx q[1];
rz(-1.7226137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0467906) q[0];
sx q[0];
rz(-1.3700458) q[0];
sx q[0];
rz(1.0958452) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1221462) q[2];
sx q[2];
rz(-1.0722954) q[2];
sx q[2];
rz(-2.4011103) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8992865) q[1];
sx q[1];
rz(-1.810009) q[1];
sx q[1];
rz(2.4278647) q[1];
rz(-1.9931692) q[3];
sx q[3];
rz(-2.7971533) q[3];
sx q[3];
rz(2.7271053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9293288) q[2];
sx q[2];
rz(-2.6876891) q[2];
sx q[2];
rz(-1.7960499) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6253925) q[0];
sx q[0];
rz(-1.1855519) q[0];
sx q[0];
rz(0.54310435) q[0];
rz(-0.25554666) q[1];
sx q[1];
rz(-2.6592022) q[1];
sx q[1];
rz(-1.7593174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10777625) q[0];
sx q[0];
rz(-1.17125) q[0];
sx q[0];
rz(-2.5679213) q[0];
rz(1.7167983) q[2];
sx q[2];
rz(-0.80933648) q[2];
sx q[2];
rz(-2.3392364) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2321736) q[1];
sx q[1];
rz(-0.62271732) q[1];
sx q[1];
rz(1.1145341) q[1];
x q[2];
rz(2.7409389) q[3];
sx q[3];
rz(-1.4130428) q[3];
sx q[3];
rz(2.1547013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9690659) q[2];
sx q[2];
rz(-1.669599) q[2];
sx q[2];
rz(0.36953163) q[2];
rz(1.8116123) q[3];
sx q[3];
rz(-2.3510635) q[3];
sx q[3];
rz(-1.3548405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0238277) q[0];
sx q[0];
rz(-2.8067639) q[0];
sx q[0];
rz(3.0552926) q[0];
rz(-1.6465126) q[1];
sx q[1];
rz(-1.9638655) q[1];
sx q[1];
rz(0.42974791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7632028) q[0];
sx q[0];
rz(-2.9400005) q[0];
sx q[0];
rz(-1.9283354) q[0];
x q[1];
rz(-2.4835973) q[2];
sx q[2];
rz(-1.832973) q[2];
sx q[2];
rz(-0.10181759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33510241) q[1];
sx q[1];
rz(-2.0079599) q[1];
sx q[1];
rz(-0.52428581) q[1];
rz(-pi) q[2];
rz(-1.2901575) q[3];
sx q[3];
rz(-0.9195111) q[3];
sx q[3];
rz(-0.7574581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1090082) q[2];
sx q[2];
rz(-2.0592368) q[2];
sx q[2];
rz(1.2004131) q[2];
rz(2.9754908) q[3];
sx q[3];
rz(-0.76986543) q[3];
sx q[3];
rz(-0.86336819) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2849041) q[0];
sx q[0];
rz(-0.76304522) q[0];
sx q[0];
rz(0.82409182) q[0];
rz(1.2245945) q[1];
sx q[1];
rz(-1.9173744) q[1];
sx q[1];
rz(-1.0138938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3124634) q[0];
sx q[0];
rz(-2.2861963) q[0];
sx q[0];
rz(0.4146284) q[0];
x q[1];
rz(0.61182558) q[2];
sx q[2];
rz(-2.6227747) q[2];
sx q[2];
rz(-0.045836115) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.287266) q[1];
sx q[1];
rz(-1.9573028) q[1];
sx q[1];
rz(0.58341656) q[1];
rz(-pi) q[2];
rz(2.1070293) q[3];
sx q[3];
rz(-1.9120029) q[3];
sx q[3];
rz(-0.1078913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27807221) q[2];
sx q[2];
rz(-0.72276989) q[2];
sx q[2];
rz(1.8966759) q[2];
rz(2.909929) q[3];
sx q[3];
rz(-2.1004227) q[3];
sx q[3];
rz(0.1483354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046655) q[0];
sx q[0];
rz(-1.9355087) q[0];
sx q[0];
rz(-2.1493602) q[0];
rz(-2.9723883) q[1];
sx q[1];
rz(-0.84183401) q[1];
sx q[1];
rz(-1.7281035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.468576) q[0];
sx q[0];
rz(-1.6454738) q[0];
sx q[0];
rz(0.59665307) q[0];
rz(-0.28861079) q[2];
sx q[2];
rz(-1.0373652) q[2];
sx q[2];
rz(2.6030428) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1991829) q[1];
sx q[1];
rz(-1.272826) q[1];
sx q[1];
rz(-1.5375141) q[1];
rz(-pi) q[2];
rz(1.1999137) q[3];
sx q[3];
rz(-2.6318916) q[3];
sx q[3];
rz(-0.56641691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8168872) q[2];
sx q[2];
rz(-1.6605261) q[2];
sx q[2];
rz(1.4229577) q[2];
rz(0.11296806) q[3];
sx q[3];
rz(-0.26765099) q[3];
sx q[3];
rz(-2.512219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1896352) q[0];
sx q[0];
rz(-1.676214) q[0];
sx q[0];
rz(2.5919609) q[0];
rz(-0.59898392) q[1];
sx q[1];
rz(-2.2689029) q[1];
sx q[1];
rz(-1.5113066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0534759) q[0];
sx q[0];
rz(-1.5742745) q[0];
sx q[0];
rz(0.17668488) q[0];
x q[1];
rz(-2.6312066) q[2];
sx q[2];
rz(-2.5742219) q[2];
sx q[2];
rz(2.9387143) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0508742) q[1];
sx q[1];
rz(-1.9719187) q[1];
sx q[1];
rz(-2.119333) q[1];
rz(-pi) q[2];
rz(1.263932) q[3];
sx q[3];
rz(-2.0425386) q[3];
sx q[3];
rz(0.8620756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0281684) q[2];
sx q[2];
rz(-1.3415965) q[2];
sx q[2];
rz(-2.8070731) q[2];
rz(-2.4962375) q[3];
sx q[3];
rz(-2.1367475) q[3];
sx q[3];
rz(-2.9042802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(0.065011218) q[0];
sx q[0];
rz(-2.8028221) q[0];
sx q[0];
rz(0.6231128) q[0];
rz(2.8399641) q[1];
sx q[1];
rz(-2.5239065) q[1];
sx q[1];
rz(1.041144) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.96466) q[0];
sx q[0];
rz(-2.3752691) q[0];
sx q[0];
rz(0.64789741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6807296) q[2];
sx q[2];
rz(-1.4878556) q[2];
sx q[2];
rz(-1.2316224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99771254) q[1];
sx q[1];
rz(-1.6015983) q[1];
sx q[1];
rz(0.7672337) q[1];
x q[2];
rz(-1.9275366) q[3];
sx q[3];
rz(-1.7715653) q[3];
sx q[3];
rz(3.1266969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.69680301) q[2];
sx q[2];
rz(-2.5371234) q[2];
sx q[2];
rz(-2.1984072) q[2];
rz(-1.7275564) q[3];
sx q[3];
rz(-0.90641886) q[3];
sx q[3];
rz(0.79677719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74054756) q[0];
sx q[0];
rz(-0.39401207) q[0];
sx q[0];
rz(3.0643585) q[0];
rz(-2.7229591) q[1];
sx q[1];
rz(-1.5822462) q[1];
sx q[1];
rz(0.15204631) q[1];
rz(-1.536191) q[2];
sx q[2];
rz(-2.1454951) q[2];
sx q[2];
rz(-2.9401671) q[2];
rz(-2.5984305) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(1.4027355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
