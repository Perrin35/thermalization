OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.815925) q[0];
sx q[0];
rz(-0.40580931) q[0];
sx q[0];
rz(1.146846) q[0];
rz(-pi) q[1];
rz(-1.5288562) q[2];
sx q[2];
rz(-2.0176) q[2];
sx q[2];
rz(0.97181335) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5056155) q[1];
sx q[1];
rz(-2.639289) q[1];
sx q[1];
rz(0.14640267) q[1];
rz(-pi) q[2];
rz(-1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(-2.3922065) q[2];
rz(-2.1253712) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(-0.40482503) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7063023) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(2.170927) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(2.326139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0464697) q[0];
sx q[0];
rz(-0.9233343) q[0];
sx q[0];
rz(-1.4955273) q[0];
rz(-pi) q[1];
rz(-2.8380727) q[2];
sx q[2];
rz(-0.75280658) q[2];
sx q[2];
rz(0.34740651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7533469) q[1];
sx q[1];
rz(-1.2699632) q[1];
sx q[1];
rz(1.6080329) q[1];
rz(-2.4957982) q[3];
sx q[3];
rz(-1.3919953) q[3];
sx q[3];
rz(1.4327232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-2.3735235) q[3];
sx q[3];
rz(-2.8320584) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3011424) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(1.8114999) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777269) q[0];
sx q[0];
rz(-0.3361055) q[0];
sx q[0];
rz(2.0396114) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61727662) q[2];
sx q[2];
rz(-1.0059788) q[2];
sx q[2];
rz(-1.9895983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.98993436) q[1];
sx q[1];
rz(-0.5608359) q[1];
sx q[1];
rz(-0.49353091) q[1];
rz(-pi) q[2];
rz(2.0796892) q[3];
sx q[3];
rz(-1.6403927) q[3];
sx q[3];
rz(3.0825465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6040566) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(0.58004722) q[2];
rz(2.3245658) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76628768) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-0.91039175) q[0];
rz(2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3954454) q[0];
sx q[0];
rz(-1.4822042) q[0];
sx q[0];
rz(0.079061411) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0312112) q[2];
sx q[2];
rz(-1.5722256) q[2];
sx q[2];
rz(1.9542076) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6064925) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(0.54775723) q[1];
rz(1.4196017) q[3];
sx q[3];
rz(-0.70221838) q[3];
sx q[3];
rz(-0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.794902) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65790025) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(2.9580341) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(1.516974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.330634) q[0];
sx q[0];
rz(-1.2067544) q[0];
sx q[0];
rz(0.90322687) q[0];
rz(1.9618271) q[2];
sx q[2];
rz(-1.0707676) q[2];
sx q[2];
rz(0.21991877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6301873) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(1.9460815) q[1];
x q[2];
rz(-0.77002854) q[3];
sx q[3];
rz(-2.3838245) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3395485) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-0.3240164) q[2];
rz(-1.8185395) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(-1.6103305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-1.0389675) q[0];
sx q[0];
rz(-1.8776241) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(0.34067672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6402123) q[0];
sx q[0];
rz(-3.0684154) q[0];
sx q[0];
rz(-1.120938) q[0];
rz(-2.8315115) q[2];
sx q[2];
rz(-0.78044621) q[2];
sx q[2];
rz(2.9030637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50960474) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(-1.3118841) q[1];
rz(1.0814704) q[3];
sx q[3];
rz(-1.9483856) q[3];
sx q[3];
rz(1.1946354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(-2.3699956) q[2];
rz(-0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1691549) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(2.7959438) q[0];
rz(-3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-0.46494928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9879887) q[0];
sx q[0];
rz(-1.1055595) q[0];
sx q[0];
rz(-0.21501712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99636997) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(0.043957274) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3073687) q[1];
sx q[1];
rz(-0.24337473) q[1];
sx q[1];
rz(-2.6878396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.081359) q[3];
sx q[3];
rz(-3.0887103) q[3];
sx q[3];
rz(0.75170654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1376301) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(-0.57146227) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.8472291) q[0];
sx q[0];
rz(-2.8572594) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-2.9998144) q[1];
sx q[1];
rz(-0.078358738) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6358444) q[0];
sx q[0];
rz(-2.3730179) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(-pi) q[1];
rz(2.9044754) q[2];
sx q[2];
rz(-1.8097005) q[2];
sx q[2];
rz(-0.6616) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2625418) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(-2.2294728) q[1];
rz(1.3571635) q[3];
sx q[3];
rz(-2.9812818) q[3];
sx q[3];
rz(2.8182639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(1.4096227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79779977) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(0.051368512) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-0.87402469) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4955935) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(1.5828703) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1633515) q[2];
sx q[2];
rz(-1.904084) q[2];
sx q[2];
rz(0.25892205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6898432) q[1];
sx q[1];
rz(-0.50536957) q[1];
sx q[1];
rz(-1.4228574) q[1];
rz(0.47253982) q[3];
sx q[3];
rz(-2.4774385) q[3];
sx q[3];
rz(0.00045517552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.727227) q[2];
sx q[2];
rz(-2.3870769) q[2];
sx q[2];
rz(-0.61974636) q[2];
rz(-1.9571346) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(-1.363389) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0062362) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-0.7243048) q[0];
rz(2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(-1.9627409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9261242) q[0];
sx q[0];
rz(-1.3660396) q[0];
sx q[0];
rz(-1.8490851) q[0];
rz(1.4775425) q[2];
sx q[2];
rz(-1.4443195) q[2];
sx q[2];
rz(-1.7699514) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2420173) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(-0.94355299) q[1];
rz(0.52272777) q[3];
sx q[3];
rz(-1.7676815) q[3];
sx q[3];
rz(0.47375351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.05802352) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(-2.2670238) q[3];
sx q[3];
rz(-0.42566291) q[3];
sx q[3];
rz(-1.6806867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(2.010871) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(2.6562128) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
