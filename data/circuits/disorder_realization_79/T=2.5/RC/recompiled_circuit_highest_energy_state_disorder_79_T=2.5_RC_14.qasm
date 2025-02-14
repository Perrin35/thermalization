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
rz(-0.35222346) q[1];
sx q[1];
rz(-1.7429587) q[1];
sx q[1];
rz(1.7359098) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79069239) q[0];
sx q[0];
rz(-0.003482799) q[0];
sx q[0];
rz(3.0538959) q[0];
rz(-pi) q[1];
rz(-0.1867577) q[2];
sx q[2];
rz(-2.6180912) q[2];
sx q[2];
rz(-1.40846) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5454332) q[1];
sx q[1];
rz(-1.9539297) q[1];
sx q[1];
rz(3.1379164) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5380404) q[3];
sx q[3];
rz(-2.1385953) q[3];
sx q[3];
rz(-2.7038684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63456717) q[2];
sx q[2];
rz(-3.1105803) q[2];
sx q[2];
rz(2.3090889) q[2];
rz(-0.80820525) q[3];
sx q[3];
rz(-3.1263604) q[3];
sx q[3];
rz(2.8830849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7611564) q[0];
sx q[0];
rz(-2.7988837) q[0];
sx q[0];
rz(-2.9535182) q[0];
rz(0.07218083) q[1];
sx q[1];
rz(-2.1071823) q[1];
sx q[1];
rz(-1.6292876) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24969026) q[0];
sx q[0];
rz(-1.9034667) q[0];
sx q[0];
rz(1.692125) q[0];
rz(-pi) q[1];
rz(-0.98324267) q[2];
sx q[2];
rz(-3.0863783) q[2];
sx q[2];
rz(3.0172341) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9413599) q[1];
sx q[1];
rz(-1.5116475) q[1];
sx q[1];
rz(-3.0422387) q[1];
rz(-0.54297437) q[3];
sx q[3];
rz(-2.256685) q[3];
sx q[3];
rz(0.35561527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.170681) q[2];
sx q[2];
rz(-0.6450246) q[2];
sx q[2];
rz(1.8485273) q[2];
rz(-0.34717789) q[3];
sx q[3];
rz(-0.2694338) q[3];
sx q[3];
rz(-3.0612883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(1.8186571) q[0];
sx q[0];
rz(-1.289239) q[0];
sx q[0];
rz(0.47054189) q[0];
rz(1.9412387) q[1];
sx q[1];
rz(-2.410694) q[1];
sx q[1];
rz(1.2568731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3124533) q[0];
sx q[0];
rz(-2.1020254) q[0];
sx q[0];
rz(0.4108528) q[0];
x q[1];
rz(2.1684709) q[2];
sx q[2];
rz(-2.9785756) q[2];
sx q[2];
rz(-2.4578641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10373579) q[1];
sx q[1];
rz(-1.5694261) q[1];
sx q[1];
rz(-1.5853154) q[1];
rz(-1.9745578) q[3];
sx q[3];
rz(-0.65126538) q[3];
sx q[3];
rz(-2.2869956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5990126) q[2];
sx q[2];
rz(-0.97992367) q[2];
sx q[2];
rz(-2.2194594) q[2];
rz(-1.7982091) q[3];
sx q[3];
rz(-2.193439) q[3];
sx q[3];
rz(-1.5141727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87329292) q[0];
sx q[0];
rz(-1.0404077) q[0];
sx q[0];
rz(-0.44019765) q[0];
rz(1.6104376) q[1];
sx q[1];
rz(-1.6596158) q[1];
sx q[1];
rz(-0.24756113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932382) q[0];
sx q[0];
rz(-2.4510151) q[0];
sx q[0];
rz(2.2838222) q[0];
x q[1];
rz(0.48657067) q[2];
sx q[2];
rz(-2.9503194) q[2];
sx q[2];
rz(-0.74081206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2046856) q[1];
sx q[1];
rz(-1.8695033) q[1];
sx q[1];
rz(-3.0644817) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4930246) q[3];
sx q[3];
rz(-1.5463136) q[3];
sx q[3];
rz(0.45059965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2441248) q[2];
sx q[2];
rz(-2.1110057) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9852801) q[0];
sx q[0];
rz(-0.15287481) q[0];
sx q[0];
rz(2.5961764) q[0];
rz(3.0975869) q[1];
sx q[1];
rz(-3.1237055) q[1];
sx q[1];
rz(-2.6652179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524844) q[0];
sx q[0];
rz(-1.8140287) q[0];
sx q[0];
rz(-1.4840675) q[0];
x q[1];
rz(-2.3628391) q[2];
sx q[2];
rz(-2.5980332) q[2];
sx q[2];
rz(1.5540079) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7303826) q[1];
sx q[1];
rz(-1.752509) q[1];
sx q[1];
rz(0.86105168) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0073193) q[3];
sx q[3];
rz(-1.065514) q[3];
sx q[3];
rz(-3.0485632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.62006092) q[2];
sx q[2];
rz(-2.4510577) q[2];
sx q[2];
rz(2.7812092) q[2];
rz(0.22200577) q[3];
sx q[3];
rz(-1.5304151) q[3];
sx q[3];
rz(1.7300026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5510657) q[0];
sx q[0];
rz(-1.5532302) q[0];
sx q[0];
rz(-1.5850413) q[0];
rz(-0.12697728) q[1];
sx q[1];
rz(-1.304909) q[1];
sx q[1];
rz(-3.051905) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68192775) q[0];
sx q[0];
rz(-1.0583504) q[0];
sx q[0];
rz(-0.95109032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78533919) q[2];
sx q[2];
rz(-1.1113104) q[2];
sx q[2];
rz(-1.6096514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6818004) q[1];
sx q[1];
rz(-0.55306095) q[1];
sx q[1];
rz(2.6917767) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60097127) q[3];
sx q[3];
rz(-1.2699763) q[3];
sx q[3];
rz(-1.6114637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.108532) q[2];
sx q[2];
rz(-2.8046799) q[2];
sx q[2];
rz(-0.25165558) q[2];
rz(-3.1015977) q[3];
sx q[3];
rz(-0.34188855) q[3];
sx q[3];
rz(2.6778636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7104178) q[0];
sx q[0];
rz(-0.091910563) q[0];
sx q[0];
rz(0.41496667) q[0];
rz(-1.4893432) q[1];
sx q[1];
rz(-0.057561189) q[1];
sx q[1];
rz(0.32589486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7975066) q[0];
sx q[0];
rz(-2.5523253) q[0];
sx q[0];
rz(-0.87743585) q[0];
rz(-pi) q[1];
rz(-3.0517111) q[2];
sx q[2];
rz(-2.0873859) q[2];
sx q[2];
rz(0.84026779) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.88837543) q[1];
sx q[1];
rz(-1.0723812) q[1];
sx q[1];
rz(-1.8640169) q[1];
x q[2];
rz(1.2674428) q[3];
sx q[3];
rz(-1.4983431) q[3];
sx q[3];
rz(0.8006351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9140279) q[2];
sx q[2];
rz(-1.3332557) q[2];
sx q[2];
rz(1.0464767) q[2];
rz(-2.922831) q[3];
sx q[3];
rz(-1.0294139) q[3];
sx q[3];
rz(-1.3974765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.47531146) q[0];
sx q[0];
rz(-3.1240211) q[0];
sx q[0];
rz(0.46324357) q[0];
rz(-2.8766368) q[1];
sx q[1];
rz(-0.0015365096) q[1];
sx q[1];
rz(1.499768) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2210082) q[0];
sx q[0];
rz(-2.7084368) q[0];
sx q[0];
rz(-1.6154352) q[0];
rz(-pi) q[1];
rz(3.0266989) q[2];
sx q[2];
rz(-0.44821366) q[2];
sx q[2];
rz(-3.093733) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8275576) q[1];
sx q[1];
rz(-1.5416939) q[1];
sx q[1];
rz(-1.7926137) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0720674) q[3];
sx q[3];
rz(-1.8226133) q[3];
sx q[3];
rz(1.0296643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90193343) q[2];
sx q[2];
rz(-0.45238164) q[2];
sx q[2];
rz(1.8233914) q[2];
rz(0.0031331172) q[3];
sx q[3];
rz(-1.2001218) q[3];
sx q[3];
rz(-1.9759294) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5779293) q[0];
sx q[0];
rz(-0.85275537) q[0];
sx q[0];
rz(-0.71459115) q[0];
rz(-2.928012) q[1];
sx q[1];
rz(-3.112308) q[1];
sx q[1];
rz(-1.9307131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9506257) q[0];
sx q[0];
rz(-1.4651148) q[0];
sx q[0];
rz(-1.8777385) q[0];
rz(-pi) q[1];
rz(-0.69542747) q[2];
sx q[2];
rz(-1.8423242) q[2];
sx q[2];
rz(2.5144387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3397202) q[1];
sx q[1];
rz(-2.0959217) q[1];
sx q[1];
rz(2.8169291) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3812495) q[3];
sx q[3];
rz(-1.726103) q[3];
sx q[3];
rz(-1.7596173) q[3];
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
rz(-2.9435797) q[2];
rz(0.35061947) q[3];
sx q[3];
rz(-1.5712761) q[3];
sx q[3];
rz(-1.9982136) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1853679) q[0];
sx q[0];
rz(-1.3274696) q[0];
sx q[0];
rz(-0.8511935) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0103127) q[2];
sx q[2];
rz(-1.5727709) q[2];
sx q[2];
rz(-1.2009837) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1150306) q[1];
sx q[1];
rz(-1.4938338) q[1];
sx q[1];
rz(-1.7285687) q[1];
x q[2];
rz(-1.102081) q[3];
sx q[3];
rz(-2.6913683) q[3];
sx q[3];
rz(-0.9952105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1018579) q[2];
sx q[2];
rz(-2.8701344) q[2];
sx q[2];
rz(2.6177935) q[2];
rz(1.0456746) q[3];
sx q[3];
rz(-1.8664675) q[3];
sx q[3];
rz(-0.4642134) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6687748) q[0];
sx q[0];
rz(-1.5325118) q[0];
sx q[0];
rz(-1.4788628) q[0];
rz(-1.4115903) q[1];
sx q[1];
rz(-0.23163207) q[1];
sx q[1];
rz(-3.0805265) q[1];
rz(0.642943) q[2];
sx q[2];
rz(-1.8131154) q[2];
sx q[2];
rz(-2.9396069) q[2];
rz(-2.7649391) q[3];
sx q[3];
rz(-2.6451544) q[3];
sx q[3];
rz(-1.3429005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
