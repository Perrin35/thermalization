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
rz(0.6495629) q[0];
sx q[0];
rz(-0.49512884) q[0];
sx q[0];
rz(-2.8861217) q[0];
rz(-0.35222346) q[1];
sx q[1];
rz(-1.7429587) q[1];
sx q[1];
rz(1.7359098) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79069239) q[0];
sx q[0];
rz(-3.1381099) q[0];
sx q[0];
rz(-3.0538959) q[0];
rz(-pi) q[1];
rz(2.6256526) q[2];
sx q[2];
rz(-1.4778412) q[2];
sx q[2];
rz(-0.0001212349) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5355994) q[1];
sx q[1];
rz(-0.38315019) q[1];
sx q[1];
rz(1.5799171) q[1];
rz(1.5380404) q[3];
sx q[3];
rz(-2.1385953) q[3];
sx q[3];
rz(-2.7038684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5070255) q[2];
sx q[2];
rz(-3.1105803) q[2];
sx q[2];
rz(2.3090889) q[2];
rz(-0.80820525) q[3];
sx q[3];
rz(-0.015232239) q[3];
sx q[3];
rz(0.25850779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38043624) q[0];
sx q[0];
rz(-0.34270898) q[0];
sx q[0];
rz(-2.9535182) q[0];
rz(-0.07218083) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(1.6292876) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3609027) q[0];
sx q[0];
rz(-1.6854428) q[0];
sx q[0];
rz(0.33495442) q[0];
x q[1];
rz(-3.1109643) q[2];
sx q[2];
rz(-1.5248486) q[2];
sx q[2];
rz(0.46389889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20023274) q[1];
sx q[1];
rz(-1.5116475) q[1];
sx q[1];
rz(-3.0422387) q[1];
x q[2];
rz(-0.80792369) q[3];
sx q[3];
rz(-1.1594541) q[3];
sx q[3];
rz(1.5802671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9709116) q[2];
sx q[2];
rz(-2.4965681) q[2];
sx q[2];
rz(-1.8485273) q[2];
rz(-2.7944148) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(3.0612883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8186571) q[0];
sx q[0];
rz(-1.8523536) q[0];
sx q[0];
rz(0.47054189) q[0];
rz(1.9412387) q[1];
sx q[1];
rz(-2.410694) q[1];
sx q[1];
rz(-1.8847195) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82913933) q[0];
sx q[0];
rz(-2.1020254) q[0];
sx q[0];
rz(0.4108528) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.435661) q[2];
sx q[2];
rz(-1.6622512) q[2];
sx q[2];
rz(-0.29558638) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0378569) q[1];
sx q[1];
rz(-1.5694261) q[1];
sx q[1];
rz(-1.5853154) q[1];
rz(-pi) q[2];
rz(1.1670349) q[3];
sx q[3];
rz(-2.4903273) q[3];
sx q[3];
rz(2.2869956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54258004) q[2];
sx q[2];
rz(-0.97992367) q[2];
sx q[2];
rz(2.2194594) q[2];
rz(-1.3433836) q[3];
sx q[3];
rz(-2.193439) q[3];
sx q[3];
rz(-1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2682997) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(2.701395) q[0];
rz(-1.6104376) q[1];
sx q[1];
rz(-1.6596158) q[1];
sx q[1];
rz(0.24756113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2092106) q[0];
sx q[0];
rz(-0.69057751) q[0];
sx q[0];
rz(-0.85777046) q[0];
x q[1];
rz(2.655022) q[2];
sx q[2];
rz(-0.19127327) q[2];
sx q[2];
rz(-0.74081206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1936613) q[1];
sx q[1];
rz(-0.30821092) q[1];
sx q[1];
rz(-1.8159417) q[1];
rz(-pi) q[2];
rz(-2.4930246) q[3];
sx q[3];
rz(-1.5952791) q[3];
sx q[3];
rz(0.45059965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89746785) q[2];
sx q[2];
rz(-1.030587) q[2];
sx q[2];
rz(-1.6991276) q[2];
rz(0.20081946) q[3];
sx q[3];
rz(-1.9229869) q[3];
sx q[3];
rz(1.1403181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9852801) q[0];
sx q[0];
rz(-2.9887178) q[0];
sx q[0];
rz(-0.54541624) q[0];
rz(0.044005752) q[1];
sx q[1];
rz(-3.1237055) q[1];
sx q[1];
rz(2.6652179) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8524844) q[0];
sx q[0];
rz(-1.327564) q[0];
sx q[0];
rz(1.6575251) q[0];
x q[1];
rz(-0.40619855) q[2];
sx q[2];
rz(-1.1990237) q[2];
sx q[2];
rz(-2.4570454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4112101) q[1];
sx q[1];
rz(-1.3890837) q[1];
sx q[1];
rz(0.86105168) q[1];
x q[2];
rz(-1.3333605) q[3];
sx q[3];
rz(-0.52132208) q[3];
sx q[3];
rz(-2.962474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62006092) q[2];
sx q[2];
rz(-2.4510577) q[2];
sx q[2];
rz(0.36038348) q[2];
rz(2.9195869) q[3];
sx q[3];
rz(-1.6111776) q[3];
sx q[3];
rz(-1.41159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5905269) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(1.5850413) q[0];
rz(-3.0146154) q[1];
sx q[1];
rz(-1.8366837) q[1];
sx q[1];
rz(-3.051905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4596649) q[0];
sx q[0];
rz(-1.0583504) q[0];
sx q[0];
rz(0.95109032) q[0];
rz(2.3562535) q[2];
sx q[2];
rz(-2.0302823) q[2];
sx q[2];
rz(-1.6096514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6818004) q[1];
sx q[1];
rz(-2.5885317) q[1];
sx q[1];
rz(-0.44981594) q[1];
rz(1.2110269) q[3];
sx q[3];
rz(-1.0002975) q[3];
sx q[3];
rz(2.9005172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.108532) q[2];
sx q[2];
rz(-0.33691275) q[2];
sx q[2];
rz(-0.25165558) q[2];
rz(-3.1015977) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(-0.46372908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43117487) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(-0.41496667) q[0];
rz(-1.4893432) q[1];
sx q[1];
rz(-3.0840315) q[1];
sx q[1];
rz(2.8156978) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0123169) q[0];
sx q[0];
rz(-2.0124448) q[0];
sx q[0];
rz(-0.4037767) q[0];
rz(-1.7275179) q[2];
sx q[2];
rz(-0.5236519) q[2];
sx q[2];
rz(2.120842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2532172) q[1];
sx q[1];
rz(-1.0723812) q[1];
sx q[1];
rz(1.2775757) q[1];
x q[2];
rz(-1.332446) q[3];
sx q[3];
rz(-2.8299667) q[3];
sx q[3];
rz(2.1441604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2275647) q[2];
sx q[2];
rz(-1.808337) q[2];
sx q[2];
rz(1.0464767) q[2];
rz(-2.922831) q[3];
sx q[3];
rz(-2.1121787) q[3];
sx q[3];
rz(1.3974765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6662812) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(-2.6783491) q[0];
rz(-2.8766368) q[1];
sx q[1];
rz(-3.1400561) q[1];
sx q[1];
rz(-1.499768) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4508598) q[0];
sx q[0];
rz(-1.5520649) q[0];
sx q[0];
rz(-2.0035726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6959571) q[2];
sx q[2];
rz(-1.5210954) q[2];
sx q[2];
rz(-1.4193064) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3140351) q[1];
sx q[1];
rz(-1.5998987) q[1];
sx q[1];
rz(-1.3489789) q[1];
x q[2];
rz(-1.8231976) q[3];
sx q[3];
rz(-1.5034672) q[3];
sx q[3];
rz(0.52378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90193343) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(-1.8233914) q[2];
rz(-0.0031331172) q[3];
sx q[3];
rz(-1.9414709) q[3];
sx q[3];
rz(1.1656632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5636633) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(-0.71459115) q[0];
rz(-2.928012) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.9307131) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7283413) q[0];
sx q[0];
rz(-1.8759707) q[0];
sx q[0];
rz(-0.11082173) q[0];
x q[1];
rz(-2.4461652) q[2];
sx q[2];
rz(-1.8423242) q[2];
sx q[2];
rz(0.62715392) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9309733) q[1];
sx q[1];
rz(-0.60928357) q[1];
sx q[1];
rz(-2.0740725) q[1];
rz(-pi) q[2];
rz(1.3812495) q[3];
sx q[3];
rz(-1.4154896) q[3];
sx q[3];
rz(1.3819753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1750298) q[2];
sx q[2];
rz(-3.1201456) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(-2.7909732) q[3];
sx q[3];
rz(-1.5703166) q[3];
sx q[3];
rz(-1.143379) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458934) q[0];
sx q[0];
rz(-2.2191255) q[0];
sx q[0];
rz(2.5272227) q[0];
rz(-0.059582926) q[1];
sx q[1];
rz(-3.0501084) q[1];
sx q[1];
rz(1.7170067) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40648288) q[0];
sx q[0];
rz(-0.87672675) q[0];
sx q[0];
rz(-0.31882341) q[0];
x q[1];
rz(0.015083357) q[2];
sx q[2];
rz(-0.13129474) q[2];
sx q[2];
rz(2.7867336) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0265621) q[1];
sx q[1];
rz(-1.6477589) q[1];
sx q[1];
rz(1.413024) q[1];
rz(0.21496694) q[3];
sx q[3];
rz(-1.1720814) q[3];
sx q[3];
rz(-0.48290184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1018579) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(0.52379918) q[2];
rz(-2.0959181) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(2.6773793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47281784) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.7300023) q[1];
sx q[1];
rz(-2.9099606) q[1];
sx q[1];
rz(0.06106613) q[1];
rz(0.3910364) q[2];
sx q[2];
rz(-2.4606065) q[2];
sx q[2];
rz(1.4627964) q[2];
rz(-1.7674592) q[3];
sx q[3];
rz(-2.0296367) q[3];
sx q[3];
rz(2.2214132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
