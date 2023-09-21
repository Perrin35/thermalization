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
rz(-2.1589307) q[0];
sx q[0];
rz(0.76173705) q[0];
rz(-2.3770483) q[1];
sx q[1];
rz(-1.0772871) q[1];
sx q[1];
rz(2.397937) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418158) q[0];
sx q[0];
rz(-1.3306949) q[0];
sx q[0];
rz(-0.067702985) q[0];
x q[1];
rz(1.450199) q[2];
sx q[2];
rz(-0.92168027) q[2];
sx q[2];
rz(2.803034) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2603862) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(1.9828412) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7573962) q[3];
sx q[3];
rz(-2.7717234) q[3];
sx q[3];
rz(-0.44934011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4347697) q[2];
sx q[2];
rz(-0.40199026) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(-0.14262959) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-0.67115152) q[1];
sx q[1];
rz(3.1006295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3687392) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(0.98451891) q[0];
rz(-pi) q[1];
rz(1.2626921) q[2];
sx q[2];
rz(-1.3034046) q[2];
sx q[2];
rz(-1.845713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7663824) q[1];
sx q[1];
rz(-0.76645215) q[1];
sx q[1];
rz(0.68739989) q[1];
rz(0.81539865) q[3];
sx q[3];
rz(-0.5268464) q[3];
sx q[3];
rz(0.53838733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-2.6290821) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(2.4643262) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5325461) q[0];
sx q[0];
rz(-0.17621528) q[0];
sx q[0];
rz(-0.30058582) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5718939) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(2.9874143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22563572) q[1];
sx q[1];
rz(-1.7491873) q[1];
sx q[1];
rz(-0.13326463) q[1];
x q[2];
rz(2.4694091) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(2.4206846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(3.0947321) q[2];
rz(2.3299407) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2617944) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(-1.4021224) q[0];
rz(-pi) q[1];
rz(-1.5863717) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(-1.7473999) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9753032) q[1];
sx q[1];
rz(-1.9235833) q[1];
sx q[1];
rz(0.3198448) q[1];
rz(1.7309534) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(-0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(-1.9481109) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(-3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(-2.5505113) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39010534) q[0];
sx q[0];
rz(-0.57919466) q[0];
sx q[0];
rz(0.67436995) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41886139) q[2];
sx q[2];
rz(-1.7423811) q[2];
sx q[2];
rz(-0.78531839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5419546) q[1];
sx q[1];
rz(-1.1918187) q[1];
sx q[1];
rz(-0.71211262) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22722865) q[3];
sx q[3];
rz(-2.7471746) q[3];
sx q[3];
rz(-1.3487181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(0.067972876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(-1.4978131) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2346674) q[0];
sx q[0];
rz(-1.6229796) q[0];
sx q[0];
rz(-2.232057) q[0];
rz(-pi) q[1];
rz(0.97700714) q[2];
sx q[2];
rz(-1.0736246) q[2];
sx q[2];
rz(2.0055111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.076768) q[1];
sx q[1];
rz(-1.2829363) q[1];
sx q[1];
rz(-1.0899781) q[1];
rz(0.72553708) q[3];
sx q[3];
rz(-2.7751338) q[3];
sx q[3];
rz(3.0612502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(0.027475474) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(2.3506929) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40490155) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(-0.94876429) q[0];
x q[1];
rz(2.3979264) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(0.92168346) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.32614) q[1];
sx q[1];
rz(-0.37323144) q[1];
sx q[1];
rz(1.7883854) q[1];
x q[2];
rz(-0.88364717) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(1.7162232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(2.8875202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2124436) q[0];
sx q[0];
rz(-1.375276) q[0];
sx q[0];
rz(-1.50302) q[0];
rz(-pi) q[1];
rz(-0.28283624) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(1.1495513) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4916363) q[1];
sx q[1];
rz(-1.8815814) q[1];
sx q[1];
rz(1.7040571) q[1];
x q[2];
rz(-3.0244175) q[3];
sx q[3];
rz(-0.60248884) q[3];
sx q[3];
rz(0.52091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6415928) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(0.32661682) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.7232822) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(1.0539508) q[0];
rz(2.9846233) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(2.1597247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40030038) q[0];
sx q[0];
rz(-1.9809082) q[0];
sx q[0];
rz(1.6923231) q[0];
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
rz(pi/2) q[0];
sx q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.3853449) q[1];
sx q[1];
rz(-0.058538392) q[1];
rz(-pi) q[2];
rz(1.8198265) q[3];
sx q[3];
rz(-2.2662275) q[3];
sx q[3];
rz(1.4679366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92423576) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(2.6249264) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(0.26836747) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29753387) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(-0.78155078) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98476121) q[2];
sx q[2];
rz(-1.1981989) q[2];
sx q[2];
rz(-2.7013456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2942218) q[1];
sx q[1];
rz(-2.1547744) q[1];
sx q[1];
rz(1.2482615) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1488232) q[3];
sx q[3];
rz(-1.3475932) q[3];
sx q[3];
rz(0.85152599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(-0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(-1.2673169) q[2];
sx q[2];
rz(-1.0721285) q[2];
sx q[2];
rz(-0.68821651) q[2];
rz(-1.7715122) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
