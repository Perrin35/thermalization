OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34710994) q[0];
sx q[0];
rz(4.5664909) q[0];
sx q[0];
rz(9.1535887) q[0];
rz(-2.7746692) q[1];
sx q[1];
rz(-0.75777268) q[1];
sx q[1];
rz(-1.4581207) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0163907) q[0];
sx q[0];
rz(-3.1338965) q[0];
sx q[0];
rz(1.9793235) q[0];
rz(-pi) q[1];
rz(-1.0881937) q[2];
sx q[2];
rz(-1.4652243) q[2];
sx q[2];
rz(-2.2029049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.62195146) q[1];
sx q[1];
rz(-0.84712815) q[1];
sx q[1];
rz(0.78459854) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.062964418) q[3];
sx q[3];
rz(-2.8573572) q[3];
sx q[3];
rz(0.96719826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1272507) q[2];
sx q[2];
rz(-2.47561) q[2];
sx q[2];
rz(1.2464397) q[2];
rz(2.1261101) q[3];
sx q[3];
rz(-0.4963488) q[3];
sx q[3];
rz(-0.61771667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366197) q[0];
sx q[0];
rz(-0.07285694) q[0];
sx q[0];
rz(2.4541722) q[0];
rz(2.5303326) q[1];
sx q[1];
rz(-1.6300853) q[1];
sx q[1];
rz(0.90868178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9729823) q[0];
sx q[0];
rz(-1.6618722) q[0];
sx q[0];
rz(0.9206207) q[0];
x q[1];
rz(1.3503752) q[2];
sx q[2];
rz(-2.5155655) q[2];
sx q[2];
rz(1.859425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7707053) q[1];
sx q[1];
rz(-1.0418278) q[1];
sx q[1];
rz(0.7077528) q[1];
x q[2];
rz(0.20199235) q[3];
sx q[3];
rz(-2.575361) q[3];
sx q[3];
rz(-0.30893886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0574135) q[2];
sx q[2];
rz(-0.97687352) q[2];
sx q[2];
rz(-1.9834391) q[2];
rz(2.967945) q[3];
sx q[3];
rz(-1.4519139) q[3];
sx q[3];
rz(2.4524073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.3678906) q[0];
sx q[0];
rz(-0.20562085) q[0];
sx q[0];
rz(2.6389627) q[0];
rz(-1.3357119) q[1];
sx q[1];
rz(-0.10919658) q[1];
sx q[1];
rz(1.7371545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96539474) q[0];
sx q[0];
rz(-1.1317011) q[0];
sx q[0];
rz(-2.0140854) q[0];
rz(-pi) q[1];
rz(-2.6790256) q[2];
sx q[2];
rz(-0.41282648) q[2];
sx q[2];
rz(-3.1011229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81922779) q[1];
sx q[1];
rz(-0.33999264) q[1];
sx q[1];
rz(2.875706) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.037602) q[3];
sx q[3];
rz(-1.0467064) q[3];
sx q[3];
rz(1.0162102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3001083) q[2];
sx q[2];
rz(-2.4418162) q[2];
sx q[2];
rz(-0.67467275) q[2];
rz(2.789433) q[3];
sx q[3];
rz(-1.1511185) q[3];
sx q[3];
rz(1.5153511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34053764) q[0];
sx q[0];
rz(-3.099649) q[0];
sx q[0];
rz(1.9757784) q[0];
rz(-0.55533987) q[1];
sx q[1];
rz(-2.1644939) q[1];
sx q[1];
rz(-1.7305444) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2685617) q[0];
sx q[0];
rz(-1.8202204) q[0];
sx q[0];
rz(-1.5879575) q[0];
rz(1.0434112) q[2];
sx q[2];
rz(-1.4031271) q[2];
sx q[2];
rz(-2.4737918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1852604) q[1];
sx q[1];
rz(-1.0721551) q[1];
sx q[1];
rz(-2.7428328) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8218231) q[3];
sx q[3];
rz(-0.76956144) q[3];
sx q[3];
rz(-1.3057452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1022776) q[2];
sx q[2];
rz(-0.5216051) q[2];
sx q[2];
rz(-2.6585141) q[2];
rz(0.77110243) q[3];
sx q[3];
rz(-1.5604138) q[3];
sx q[3];
rz(-1.7980827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.022472) q[0];
sx q[0];
rz(-0.461687) q[0];
sx q[0];
rz(2.9126677) q[0];
rz(-1.4643033) q[1];
sx q[1];
rz(-2.5720282) q[1];
sx q[1];
rz(0.48008188) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4245473) q[0];
sx q[0];
rz(-1.3120894) q[0];
sx q[0];
rz(-2.8733868) q[0];
rz(-pi) q[1];
rz(0.51182439) q[2];
sx q[2];
rz(-0.49705844) q[2];
sx q[2];
rz(1.9972506) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.11949524) q[1];
sx q[1];
rz(-1.198631) q[1];
sx q[1];
rz(-1.6254025) q[1];
rz(-pi) q[2];
rz(-2.1958417) q[3];
sx q[3];
rz(-1.4212928) q[3];
sx q[3];
rz(-2.3760419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.071712581) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(-1.9055535) q[2];
rz(1.6482327) q[3];
sx q[3];
rz(-2.3107216) q[3];
sx q[3];
rz(-1.4504455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1278095) q[0];
sx q[0];
rz(-2.4007128) q[0];
sx q[0];
rz(-0.86454779) q[0];
rz(-0.8283444) q[1];
sx q[1];
rz(-1.8048077) q[1];
sx q[1];
rz(2.4116662) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975271) q[0];
sx q[0];
rz(-1.1008917) q[0];
sx q[0];
rz(0.92961981) q[0];
rz(-pi) q[1];
rz(-2.8921919) q[2];
sx q[2];
rz(-0.6847544) q[2];
sx q[2];
rz(-0.33011393) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8747507) q[1];
sx q[1];
rz(-2.0952941) q[1];
sx q[1];
rz(-0.28171087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6248405) q[3];
sx q[3];
rz(-2.2934348) q[3];
sx q[3];
rz(2.8703226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56968969) q[2];
sx q[2];
rz(-0.14293417) q[2];
sx q[2];
rz(-1.2804871) q[2];
rz(-1.5661543) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(2.2787826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5363619) q[0];
sx q[0];
rz(-2.1385758) q[0];
sx q[0];
rz(-2.5501116) q[0];
rz(-0.65762562) q[1];
sx q[1];
rz(-1.3720737) q[1];
sx q[1];
rz(0.11810158) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054679019) q[0];
sx q[0];
rz(-0.42092338) q[0];
sx q[0];
rz(2.551159) q[0];
x q[1];
rz(0.12848358) q[2];
sx q[2];
rz(-0.5196577) q[2];
sx q[2];
rz(0.78413397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76712045) q[1];
sx q[1];
rz(-0.53009696) q[1];
sx q[1];
rz(2.9494826) q[1];
x q[2];
rz(2.2198417) q[3];
sx q[3];
rz(-0.27652446) q[3];
sx q[3];
rz(-0.10504237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.20146951) q[2];
sx q[2];
rz(-0.38429364) q[2];
sx q[2];
rz(-1.2598134) q[2];
rz(-1.2093557) q[3];
sx q[3];
rz(-1.7989379) q[3];
sx q[3];
rz(2.2945837) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92846337) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(-1.553836) q[0];
rz(1.8153048) q[1];
sx q[1];
rz(-0.98618788) q[1];
sx q[1];
rz(0.80002588) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5309842) q[0];
sx q[0];
rz(-1.8433231) q[0];
sx q[0];
rz(-0.39604183) q[0];
x q[1];
rz(2.6392512) q[2];
sx q[2];
rz(-2.3141344) q[2];
sx q[2];
rz(-2.7934144) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2721438) q[1];
sx q[1];
rz(-0.71681685) q[1];
sx q[1];
rz(-2.3434756) q[1];
x q[2];
rz(0.38626892) q[3];
sx q[3];
rz(-1.0479386) q[3];
sx q[3];
rz(2.8012365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.942261) q[2];
sx q[2];
rz(-0.76392019) q[2];
sx q[2];
rz(1.2660816) q[2];
rz(-1.0670886) q[3];
sx q[3];
rz(-1.4184003) q[3];
sx q[3];
rz(-0.057723109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1570194) q[0];
sx q[0];
rz(-0.83681256) q[0];
sx q[0];
rz(-2.9465604) q[0];
rz(2.2795279) q[1];
sx q[1];
rz(-1.4007327) q[1];
sx q[1];
rz(-0.21924266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1006323) q[0];
sx q[0];
rz(-1.0182683) q[0];
sx q[0];
rz(2.841903) q[0];
rz(-pi) q[1];
rz(-2.9621808) q[2];
sx q[2];
rz(-1.1288912) q[2];
sx q[2];
rz(-0.68540689) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63263195) q[1];
sx q[1];
rz(-1.0380942) q[1];
sx q[1];
rz(1.6091634) q[1];
rz(2.3211649) q[3];
sx q[3];
rz(-1.9747989) q[3];
sx q[3];
rz(-1.272161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7361043) q[2];
sx q[2];
rz(-0.89724237) q[2];
sx q[2];
rz(-1.415095) q[2];
rz(1.803558) q[3];
sx q[3];
rz(-1.3058563) q[3];
sx q[3];
rz(3.0726748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5445484) q[0];
sx q[0];
rz(-0.91240779) q[0];
sx q[0];
rz(1.190881) q[0];
rz(-1.6944338) q[1];
sx q[1];
rz(-1.8444599) q[1];
sx q[1];
rz(-1.4318633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6042991) q[0];
sx q[0];
rz(-1.3621289) q[0];
sx q[0];
rz(-2.4775726) q[0];
x q[1];
rz(1.9831311) q[2];
sx q[2];
rz(-2.1160876) q[2];
sx q[2];
rz(0.83484989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7427976) q[1];
sx q[1];
rz(-0.94091641) q[1];
sx q[1];
rz(0.50666084) q[1];
rz(-pi) q[2];
x q[2];
rz(0.032037246) q[3];
sx q[3];
rz(-2.3028829) q[3];
sx q[3];
rz(-0.17994954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65149629) q[2];
sx q[2];
rz(-2.8054674) q[2];
sx q[2];
rz(-0.87745848) q[2];
rz(-1.5326001) q[3];
sx q[3];
rz(-1.7209777) q[3];
sx q[3];
rz(-1.0271614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-1.132591) q[0];
sx q[0];
rz(-1.7283716) q[0];
sx q[0];
rz(-0.51474095) q[0];
rz(0.71350907) q[1];
sx q[1];
rz(-2.3206354) q[1];
sx q[1];
rz(-1.6101507) q[1];
rz(1.8590676) q[2];
sx q[2];
rz(-1.9296411) q[2];
sx q[2];
rz(1.267098) q[2];
rz(-1.3446829) q[3];
sx q[3];
rz(-1.3158847) q[3];
sx q[3];
rz(-2.55008) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
