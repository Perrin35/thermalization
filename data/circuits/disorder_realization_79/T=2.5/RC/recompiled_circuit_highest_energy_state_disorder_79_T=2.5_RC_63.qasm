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
rz(0.25547096) q[0];
rz(2.7893692) q[1];
sx q[1];
rz(-1.398634) q[1];
sx q[1];
rz(1.4056828) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79069239) q[0];
sx q[0];
rz(-3.1381099) q[0];
sx q[0];
rz(-0.087696747) q[0];
x q[1];
rz(0.5159401) q[2];
sx q[2];
rz(-1.4778412) q[2];
sx q[2];
rz(0.0001212349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1683301) q[1];
sx q[1];
rz(-1.5673866) q[1];
sx q[1];
rz(-1.9539321) q[1];
x q[2];
rz(3.0902946) q[3];
sx q[3];
rz(-0.56863943) q[3];
sx q[3];
rz(-2.7647247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63456717) q[2];
sx q[2];
rz(-0.031012379) q[2];
sx q[2];
rz(-2.3090889) q[2];
rz(2.3333874) q[3];
sx q[3];
rz(-3.1263604) q[3];
sx q[3];
rz(-0.25850779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7611564) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(-0.1880745) q[0];
rz(-0.07218083) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(-1.5123051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24969026) q[0];
sx q[0];
rz(-1.9034667) q[0];
sx q[0];
rz(-1.4494677) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6167655) q[2];
sx q[2];
rz(-1.5402003) q[2];
sx q[2];
rz(2.033288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7769216) q[1];
sx q[1];
rz(-1.6699759) q[1];
sx q[1];
rz(1.6302376) q[1];
x q[2];
rz(-2.1339122) q[3];
sx q[3];
rz(-0.84656871) q[3];
sx q[3];
rz(2.7359642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9709116) q[2];
sx q[2];
rz(-0.6450246) q[2];
sx q[2];
rz(-1.8485273) q[2];
rz(2.7944148) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(0.080304317) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3229356) q[0];
sx q[0];
rz(-1.289239) q[0];
sx q[0];
rz(2.6710508) q[0];
rz(-1.9412387) q[1];
sx q[1];
rz(-2.410694) q[1];
sx q[1];
rz(1.8847195) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5393821) q[0];
sx q[0];
rz(-2.4823144) q[0];
sx q[0];
rz(0.97380096) q[0];
rz(3.049301) q[2];
sx q[2];
rz(-1.7053635) q[2];
sx q[2];
rz(-1.8539661) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5804435) q[1];
sx q[1];
rz(-3.1270091) q[1];
sx q[1];
rz(1.6648949) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8506365) q[3];
sx q[3];
rz(-2.1621063) q[3];
sx q[3];
rz(-1.7940429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5990126) q[2];
sx q[2];
rz(-2.161669) q[2];
sx q[2];
rz(-2.2194594) q[2];
rz(1.3433836) q[3];
sx q[3];
rz(-0.94815367) q[3];
sx q[3];
rz(-1.62742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87329292) q[0];
sx q[0];
rz(-2.101185) q[0];
sx q[0];
rz(2.701395) q[0];
rz(-1.6104376) q[1];
sx q[1];
rz(-1.6596158) q[1];
sx q[1];
rz(0.24756113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9151814) q[0];
sx q[0];
rz(-2.0005715) q[0];
sx q[0];
rz(-2.129401) q[0];
rz(-pi) q[1];
rz(0.16952408) q[2];
sx q[2];
rz(-1.6598083) q[2];
sx q[2];
rz(-1.8326056) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94793135) q[1];
sx q[1];
rz(-0.30821092) q[1];
sx q[1];
rz(-1.3256509) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64856802) q[3];
sx q[3];
rz(-1.5463136) q[3];
sx q[3];
rz(2.690993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2441248) q[2];
sx q[2];
rz(-2.1110057) q[2];
sx q[2];
rz(-1.4424651) q[2];
rz(2.9407732) q[3];
sx q[3];
rz(-1.2186058) q[3];
sx q[3];
rz(1.1403181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15631256) q[0];
sx q[0];
rz(-2.9887178) q[0];
sx q[0];
rz(-0.54541624) q[0];
rz(-3.0975869) q[1];
sx q[1];
rz(-0.017887201) q[1];
sx q[1];
rz(-2.6652179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5060318) q[0];
sx q[0];
rz(-0.2579435) q[0];
sx q[0];
rz(0.33584612) q[0];
rz(-pi) q[1];
rz(1.1693994) q[2];
sx q[2];
rz(-1.1938022) q[2];
sx q[2];
rz(0.73125401) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4112101) q[1];
sx q[1];
rz(-1.3890837) q[1];
sx q[1];
rz(2.280541) q[1];
rz(-pi) q[2];
rz(0.13427333) q[3];
sx q[3];
rz(-2.0760787) q[3];
sx q[3];
rz(3.0485632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62006092) q[2];
sx q[2];
rz(-0.69053495) q[2];
sx q[2];
rz(2.7812092) q[2];
rz(2.9195869) q[3];
sx q[3];
rz(-1.5304151) q[3];
sx q[3];
rz(-1.7300026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5510657) q[0];
sx q[0];
rz(-1.5883625) q[0];
sx q[0];
rz(-1.5850413) q[0];
rz(3.0146154) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(-3.051905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.854824) q[0];
sx q[0];
rz(-2.3595605) q[0];
sx q[0];
rz(-2.3402592) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96025567) q[2];
sx q[2];
rz(-0.88443236) q[2];
sx q[2];
rz(-2.6853564) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9759102) q[1];
sx q[1];
rz(-1.0780562) q[1];
sx q[1];
rz(1.3085646) q[1];
rz(-1.2110269) q[3];
sx q[3];
rz(-1.0002975) q[3];
sx q[3];
rz(0.24107547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.108532) q[2];
sx q[2];
rz(-2.8046799) q[2];
sx q[2];
rz(0.25165558) q[2];
rz(3.1015977) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(0.46372908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104178) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(2.726626) q[0];
rz(-1.4893432) q[1];
sx q[1];
rz(-3.0840315) q[1];
sx q[1];
rz(2.8156978) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0123169) q[0];
sx q[0];
rz(-1.1291478) q[0];
sx q[0];
rz(-2.737816) q[0];
x q[1];
rz(-2.0891248) q[2];
sx q[2];
rz(-1.4926693) q[2];
sx q[2];
rz(-2.4555488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.68989) q[1];
sx q[1];
rz(-2.5696515) q[1];
sx q[1];
rz(-2.6533935) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8091466) q[3];
sx q[3];
rz(-2.8299667) q[3];
sx q[3];
rz(0.99743227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2275647) q[2];
sx q[2];
rz(-1.3332557) q[2];
sx q[2];
rz(1.0464767) q[2];
rz(2.922831) q[3];
sx q[3];
rz(-1.0294139) q[3];
sx q[3];
rz(-1.7441162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662812) q[0];
sx q[0];
rz(-0.0175716) q[0];
sx q[0];
rz(2.6783491) q[0];
rz(-0.26495588) q[1];
sx q[1];
rz(-3.1400561) q[1];
sx q[1];
rz(-1.6418246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92058449) q[0];
sx q[0];
rz(-0.43315584) q[0];
sx q[0];
rz(-1.6154352) q[0];
x q[1];
rz(-1.5157264) q[2];
sx q[2];
rz(-2.0158421) q[2];
sx q[2];
rz(2.9663756) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7565175) q[1];
sx q[1];
rz(-0.2236872) q[1];
sx q[1];
rz(1.4392412) q[1];
rz(3.0720674) q[3];
sx q[3];
rz(-1.8226133) q[3];
sx q[3];
rz(-1.0296643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90193343) q[2];
sx q[2];
rz(-2.689211) q[2];
sx q[2];
rz(1.8233914) q[2];
rz(-0.0031331172) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(-1.1656632) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779293) q[0];
sx q[0];
rz(-2.2888373) q[0];
sx q[0];
rz(-0.71459115) q[0];
rz(0.21358061) q[1];
sx q[1];
rz(-0.029284632) q[1];
sx q[1];
rz(-1.2108796) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058556341) q[0];
sx q[0];
rz(-0.3240816) q[0];
sx q[0];
rz(1.9084318) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69542747) q[2];
sx q[2];
rz(-1.2992685) q[2];
sx q[2];
rz(-2.5144387) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9309733) q[1];
sx q[1];
rz(-2.5323091) q[1];
sx q[1];
rz(-2.0740725) q[1];
rz(-pi) q[2];
rz(2.2641422) q[3];
sx q[3];
rz(-2.8971379) q[3];
sx q[3];
rz(-0.86737421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1750298) q[2];
sx q[2];
rz(-3.1201456) q[2];
sx q[2];
rz(0.19801298) q[2];
rz(0.35061947) q[3];
sx q[3];
rz(-1.5703166) q[3];
sx q[3];
rz(-1.143379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458934) q[0];
sx q[0];
rz(-2.2191255) q[0];
sx q[0];
rz(2.5272227) q[0];
rz(0.059582926) q[1];
sx q[1];
rz(-3.0501084) q[1];
sx q[1];
rz(1.4245859) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2587268) q[0];
sx q[0];
rz(-2.388991) q[0];
sx q[0];
rz(1.9310139) q[0];
rz(3.0103127) q[2];
sx q[2];
rz(-1.5688217) q[2];
sx q[2];
rz(1.9406089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0265621) q[1];
sx q[1];
rz(-1.4938338) q[1];
sx q[1];
rz(1.413024) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9779102) q[3];
sx q[3];
rz(-1.7686678) q[3];
sx q[3];
rz(-2.1382633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0397348) q[2];
sx q[2];
rz(-0.27145824) q[2];
sx q[2];
rz(-2.6177935) q[2];
rz(1.0456746) q[3];
sx q[3];
rz(-1.2751251) q[3];
sx q[3];
rz(-2.6773793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47281784) q[0];
sx q[0];
rz(-1.6090809) q[0];
sx q[0];
rz(1.6627298) q[0];
rz(-1.4115903) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(1.8703441) q[2];
sx q[2];
rz(-0.94956492) q[2];
sx q[2];
rz(1.9506394) q[2];
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
