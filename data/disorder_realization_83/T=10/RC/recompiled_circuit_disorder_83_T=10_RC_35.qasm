OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5542334) q[0];
sx q[0];
rz(4.1242546) q[0];
sx q[0];
rz(10.186515) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8964891) q[0];
sx q[0];
rz(-1.6365543) q[0];
sx q[0];
rz(-1.3301646) q[0];
rz(-pi) q[1];
rz(0.6526297) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(-1.8362311) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.57063369) q[1];
sx q[1];
rz(-1.1753517) q[1];
sx q[1];
rz(2.8407437) q[1];
rz(1.9348014) q[3];
sx q[3];
rz(-1.6379106) q[3];
sx q[3];
rz(-1.8458927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649439) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-0.23072492) q[0];
rz(-1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(3.1006295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36386585) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(-1.8198387) q[0];
rz(-2.3054753) q[2];
sx q[2];
rz(-2.7364519) q[2];
sx q[2];
rz(-2.1737614) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4119271) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(2.5018442) q[1];
x q[2];
rz(1.1702234) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.846116) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(-2.4643262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8375741) q[0];
sx q[0];
rz(-1.4025592) q[0];
sx q[0];
rz(-1.6234682) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1352273) q[2];
sx q[2];
rz(-1.5718939) q[2];
sx q[2];
rz(1.416625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.22563572) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(-3.008328) q[1];
x q[2];
rz(-0.67218353) q[3];
sx q[3];
rz(-2.5568092) q[3];
sx q[3];
rz(0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(-2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(-0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-0.051483367) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576142) q[0];
sx q[0];
rz(-1.4371705) q[0];
sx q[0];
rz(2.4806116) q[0];
x q[1];
rz(-3.1191191) q[2];
sx q[2];
rz(-0.60612504) q[2];
sx q[2];
rz(1.4215353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7395775) q[1];
sx q[1];
rz(-2.6699454) q[1];
sx q[1];
rz(0.86400835) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7309534) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(0.48689294) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.8160965) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(2.4868734) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3712758) q[0];
sx q[0];
rz(-1.2219984) q[0];
sx q[0];
rz(-2.6692997) q[0];
x q[1];
rz(-1.7582558) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(2.2802441) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3760738) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(-2.5942624) q[1];
x q[2];
rz(-0.38526411) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(0.011766089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(2.5704685) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(-1.3202745) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90692524) q[0];
sx q[0];
rz(-1.5186131) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-pi) q[1];
rz(-2.5619179) q[2];
sx q[2];
rz(-2.0849166) q[2];
sx q[2];
rz(2.3953953) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0648246) q[1];
sx q[1];
rz(-1.8586564) q[1];
sx q[1];
rz(2.0516146) q[1];
x q[2];
rz(1.8201581) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(2.3251422) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(-0.11418848) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(2.3506929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7366911) q[0];
sx q[0];
rz(-1.4762523) q[0];
sx q[0];
rz(2.1928284) q[0];
rz(-pi) q[1];
rz(0.74366624) q[2];
sx q[2];
rz(-1.8218092) q[2];
sx q[2];
rz(-2.2199092) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8154527) q[1];
sx q[1];
rz(-0.37323144) q[1];
sx q[1];
rz(1.3532072) q[1];
x q[2];
rz(-2.2579455) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2699282) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5131322) q[0];
sx q[0];
rz(-1.5043133) q[0];
sx q[0];
rz(-2.9456338) q[0];
rz(-pi) q[1];
rz(0.0052070219) q[2];
sx q[2];
rz(-1.5723096) q[2];
sx q[2];
rz(3.0031799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2368187) q[1];
sx q[1];
rz(-2.8042951) q[1];
sx q[1];
rz(-0.39223139) q[1];
x q[2];
rz(-1.4905606) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(0.3790006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(0.32661682) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(-1.4183104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083387233) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-0.98186791) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69753416) q[0];
sx q[0];
rz(-0.42675787) q[0];
sx q[0];
rz(-2.8696637) q[0];
rz(-pi) q[1];
rz(0.1083381) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(-0.94891753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.3853449) q[1];
sx q[1];
rz(-0.058538392) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8544159) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(-1.2958131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(3.1072646) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(-0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-0.26836747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4386913) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(-0.86433522) q[0];
rz(-pi) q[1];
rz(-2.1860113) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(-2.5126484) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8908773) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(-0.4472181) q[1];
rz(-pi) q[2];
rz(0.26465613) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(-2.2789126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-2.6230326) q[2];
sx q[2];
rz(-1.8363562) q[2];
sx q[2];
rz(0.73391757) q[2];
rz(0.56223829) q[3];
sx q[3];
rz(-2.7769965) q[3];
sx q[3];
rz(-2.1806352) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
