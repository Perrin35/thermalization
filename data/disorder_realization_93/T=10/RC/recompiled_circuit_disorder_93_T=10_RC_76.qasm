OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(-1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0215065) q[0];
sx q[0];
rz(-1.2354295) q[0];
sx q[0];
rz(2.04727) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0338697) q[2];
sx q[2];
rz(-0.31422868) q[2];
sx q[2];
rz(-0.83480922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1957789) q[1];
sx q[1];
rz(-1.073277) q[1];
sx q[1];
rz(-2.0583908) q[1];
rz(-1.190891) q[3];
sx q[3];
rz(-2.0348747) q[3];
sx q[3];
rz(0.78117785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2279921) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(-2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.686036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50549492) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(1.8683744) q[0];
rz(-pi) q[1];
rz(-2.935264) q[2];
sx q[2];
rz(-2.7596139) q[2];
sx q[2];
rz(-1.0155592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.204139) q[1];
sx q[1];
rz(-1.0480282) q[1];
sx q[1];
rz(-0.015603113) q[1];
rz(-3.0807642) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(0.5371679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(-2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75333726) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(0.02877409) q[1];
sx q[1];
rz(-2.0859699) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81329936) q[0];
sx q[0];
rz(-2.1345703) q[0];
sx q[0];
rz(-1.203712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.795479) q[2];
sx q[2];
rz(-0.78595224) q[2];
sx q[2];
rz(-1.2197989) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.064627083) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(-1.1175734) q[1];
rz(-pi) q[2];
rz(-1.7964732) q[3];
sx q[3];
rz(-0.25697069) q[3];
sx q[3];
rz(0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-2.6180843) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1062766) q[0];
sx q[0];
rz(-1.9395394) q[0];
sx q[0];
rz(-2.5160401) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8977215) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(-1.7983758) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.254863) q[1];
sx q[1];
rz(-0.88725315) q[1];
sx q[1];
rz(-2.8121594) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13682271) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(0.50859355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-2.5775487) q[2];
rz(-2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48150912) q[0];
sx q[0];
rz(-0.68843377) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-0.99194828) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49682) q[0];
sx q[0];
rz(-2.9368375) q[0];
sx q[0];
rz(2.5909008) q[0];
rz(-pi) q[1];
rz(-1.2679342) q[2];
sx q[2];
rz(-1.5622683) q[2];
sx q[2];
rz(-2.0451343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3757513) q[1];
sx q[1];
rz(-1.6110794) q[1];
sx q[1];
rz(0.30893107) q[1];
rz(-pi) q[2];
rz(-0.81327849) q[3];
sx q[3];
rz(-1.6653898) q[3];
sx q[3];
rz(0.41155848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(0.27080718) q[2];
rz(-2.9233542) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.4250925) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(-1.3160734) q[0];
rz(-pi) q[1];
rz(1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(2.595682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9565935) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(-1.1207629) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1720042) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(-0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0075334) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(-0.56387222) q[2];
rz(-0.18051906) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-0.40294161) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5086223) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(0.41931835) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99188995) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(0.71233149) q[0];
rz(0.49634883) q[2];
sx q[2];
rz(-2.2349572) q[2];
sx q[2];
rz(-1.6830483) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8824749) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(-2.4005753) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58737289) q[3];
sx q[3];
rz(-1.8837187) q[3];
sx q[3];
rz(-1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(-1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(0.82292557) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(-1.3051422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179498) q[0];
sx q[0];
rz(-0.68730132) q[0];
sx q[0];
rz(2.6108517) q[0];
rz(1.9403946) q[2];
sx q[2];
rz(-2.2410789) q[2];
sx q[2];
rz(0.13424644) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0481938) q[1];
sx q[1];
rz(-1.1953925) q[1];
sx q[1];
rz(-2.4673389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2315024) q[3];
sx q[3];
rz(-1.5486071) q[3];
sx q[3];
rz(-1.0957749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1398853) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(2.4386491) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020333175) q[0];
sx q[0];
rz(-0.54176211) q[0];
sx q[0];
rz(-2.3938177) q[0];
rz(1.0614971) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(2.408037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.798588) q[1];
sx q[1];
rz(-1.4061905) q[1];
sx q[1];
rz(-1.027147) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95548198) q[3];
sx q[3];
rz(-1.2947047) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90074173) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(-0.46863619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2221453) q[0];
sx q[0];
rz(-1.3627909) q[0];
sx q[0];
rz(-2.749445) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87601985) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(1.3181869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.66099) q[1];
sx q[1];
rz(-0.27462474) q[1];
sx q[1];
rz(1.3165228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1115132) q[3];
sx q[3];
rz(-1.9643524) q[3];
sx q[3];
rz(0.72898385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.5293998) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464012) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(0.8016349) q[2];
sx q[2];
rz(-0.87052204) q[2];
sx q[2];
rz(-1.082765) q[2];
rz(1.6735531) q[3];
sx q[3];
rz(-0.65410528) q[3];
sx q[3];
rz(0.42412529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
