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
rz(2.8925722) q[0];
sx q[0];
rz(-2.1828716) q[0];
sx q[0];
rz(0.22554654) q[0];
rz(2.0685937) q[1];
sx q[1];
rz(-0.45773157) q[1];
sx q[1];
rz(0.20888858) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208369) q[0];
sx q[0];
rz(-1.136047) q[0];
sx q[0];
rz(-1.8692383) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63103478) q[2];
sx q[2];
rz(-1.8004144) q[2];
sx q[2];
rz(0.49447021) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35395216) q[1];
sx q[1];
rz(-0.6122511) q[1];
sx q[1];
rz(-2.5933215) q[1];
x q[2];
rz(2.3129102) q[3];
sx q[3];
rz(-0.4503612) q[3];
sx q[3];
rz(2.8652103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1170342) q[2];
sx q[2];
rz(-2.890675) q[2];
sx q[2];
rz(0.92570242) q[2];
rz(-2.3671345) q[3];
sx q[3];
rz(-1.9175994) q[3];
sx q[3];
rz(-1.8584049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9722612) q[0];
sx q[0];
rz(-2.4915578) q[0];
sx q[0];
rz(-1.6768804) q[0];
rz(-0.88893923) q[1];
sx q[1];
rz(-2.2496532) q[1];
sx q[1];
rz(-1.7882141) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52349963) q[0];
sx q[0];
rz(-0.67636469) q[0];
sx q[0];
rz(0.47250611) q[0];
rz(-pi) q[1];
rz(0.36608669) q[2];
sx q[2];
rz(-1.5502068) q[2];
sx q[2];
rz(0.96167694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96129744) q[1];
sx q[1];
rz(-2.5112481) q[1];
sx q[1];
rz(0.35758361) q[1];
x q[2];
rz(2.86496) q[3];
sx q[3];
rz(-2.3440954) q[3];
sx q[3];
rz(-1.5381787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9951524) q[2];
sx q[2];
rz(-2.1285987) q[2];
sx q[2];
rz(0.97274441) q[2];
rz(1.815833) q[3];
sx q[3];
rz(-1.1498068) q[3];
sx q[3];
rz(-2.5241234) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4165118) q[0];
sx q[0];
rz(-1.6575811) q[0];
sx q[0];
rz(0.02956477) q[0];
rz(-2.120453) q[1];
sx q[1];
rz(-0.28426668) q[1];
sx q[1];
rz(2.8254642) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9122304) q[0];
sx q[0];
rz(-1.5490313) q[0];
sx q[0];
rz(0.020177186) q[0];
x q[1];
rz(-2.0029127) q[2];
sx q[2];
rz(-0.2209835) q[2];
sx q[2];
rz(-2.5137171) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21507922) q[1];
sx q[1];
rz(-1.2383019) q[1];
sx q[1];
rz(-2.957587) q[1];
rz(-pi) q[2];
rz(-0.088313266) q[3];
sx q[3];
rz(-2.4188015) q[3];
sx q[3];
rz(-0.73107728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0315087) q[2];
sx q[2];
rz(-1.8434593) q[2];
sx q[2];
rz(-2.617344) q[2];
rz(-0.60018572) q[3];
sx q[3];
rz(-0.46349183) q[3];
sx q[3];
rz(-1.0711099) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4635187) q[0];
sx q[0];
rz(-1.2200032) q[0];
sx q[0];
rz(0.36233166) q[0];
rz(2.433297) q[1];
sx q[1];
rz(-0.34168044) q[1];
sx q[1];
rz(1.1452311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4701242) q[0];
sx q[0];
rz(-0.13397476) q[0];
sx q[0];
rz(-2.2132316) q[0];
rz(-pi) q[1];
rz(2.1378273) q[2];
sx q[2];
rz(-2.3622555) q[2];
sx q[2];
rz(0.61160751) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2140183) q[1];
sx q[1];
rz(-1.8248789) q[1];
sx q[1];
rz(-0.85846445) q[1];
rz(-1.5208552) q[3];
sx q[3];
rz(-1.4068713) q[3];
sx q[3];
rz(-2.0633179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5273253) q[2];
sx q[2];
rz(-2.3136487) q[2];
sx q[2];
rz(3.1329727) q[2];
rz(-0.65420592) q[3];
sx q[3];
rz(-2.4253186) q[3];
sx q[3];
rz(2.8692828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006186) q[0];
sx q[0];
rz(-0.32855496) q[0];
sx q[0];
rz(0.7884489) q[0];
rz(-2.12766) q[1];
sx q[1];
rz(-1.3194111) q[1];
sx q[1];
rz(-3.0533155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28391253) q[0];
sx q[0];
rz(-1.0774892) q[0];
sx q[0];
rz(2.3939214) q[0];
rz(-pi) q[1];
rz(2.9806251) q[2];
sx q[2];
rz(-0.56955284) q[2];
sx q[2];
rz(2.9496824) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83206415) q[1];
sx q[1];
rz(-2.5558379) q[1];
sx q[1];
rz(1.0447211) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5203736) q[3];
sx q[3];
rz(-2.6459624) q[3];
sx q[3];
rz(0.69678604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2225515) q[2];
sx q[2];
rz(-1.7307948) q[2];
sx q[2];
rz(-0.50773531) q[2];
rz(0.12568411) q[3];
sx q[3];
rz(-1.9583743) q[3];
sx q[3];
rz(0.79508933) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8102201) q[0];
sx q[0];
rz(-2.3851244) q[0];
sx q[0];
rz(3.0561225) q[0];
rz(-2.5147009) q[1];
sx q[1];
rz(-1.9848928) q[1];
sx q[1];
rz(1.8524106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7580494) q[0];
sx q[0];
rz(-1.6746759) q[0];
sx q[0];
rz(1.3707177) q[0];
rz(-pi) q[1];
rz(2.3979509) q[2];
sx q[2];
rz(-1.4568475) q[2];
sx q[2];
rz(2.7211962) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8823008) q[1];
sx q[1];
rz(-1.8083824) q[1];
sx q[1];
rz(-0.17452328) q[1];
rz(0.014628476) q[3];
sx q[3];
rz(-1.808015) q[3];
sx q[3];
rz(-0.62495172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.36279303) q[2];
sx q[2];
rz(-0.94393602) q[2];
sx q[2];
rz(-1.8865406) q[2];
rz(-0.058989914) q[3];
sx q[3];
rz(-1.2041644) q[3];
sx q[3];
rz(2.1571933) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4208218) q[0];
sx q[0];
rz(-1.2728007) q[0];
sx q[0];
rz(-0.76501784) q[0];
rz(1.7401241) q[1];
sx q[1];
rz(-0.88686371) q[1];
sx q[1];
rz(0.23745647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.229612) q[0];
sx q[0];
rz(-1.6061826) q[0];
sx q[0];
rz(2.1982212) q[0];
rz(-pi) q[1];
rz(1.4777021) q[2];
sx q[2];
rz(-2.4590465) q[2];
sx q[2];
rz(-1.0034221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73257414) q[1];
sx q[1];
rz(-1.5272101) q[1];
sx q[1];
rz(0.2909046) q[1];
rz(-pi) q[2];
rz(-1.5349489) q[3];
sx q[3];
rz(-2.4447943) q[3];
sx q[3];
rz(-0.076059503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.456363) q[2];
sx q[2];
rz(-2.3380029) q[2];
sx q[2];
rz(-2.2881499) q[2];
rz(1.8026836) q[3];
sx q[3];
rz(-1.572255) q[3];
sx q[3];
rz(0.16551031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45561871) q[0];
sx q[0];
rz(-1.2116665) q[0];
sx q[0];
rz(-2.9562505) q[0];
rz(0.39400426) q[1];
sx q[1];
rz(-1.3961671) q[1];
sx q[1];
rz(-1.0967163) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95581302) q[0];
sx q[0];
rz(-2.5086864) q[0];
sx q[0];
rz(0.88297412) q[0];
x q[1];
rz(-0.11285891) q[2];
sx q[2];
rz(-2.2328909) q[2];
sx q[2];
rz(2.2619132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1450069) q[1];
sx q[1];
rz(-1.2357986) q[1];
sx q[1];
rz(2.5162391) q[1];
rz(-pi) q[2];
rz(0.53727178) q[3];
sx q[3];
rz(-0.66681615) q[3];
sx q[3];
rz(1.7423354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77157053) q[2];
sx q[2];
rz(-0.084150704) q[2];
sx q[2];
rz(2.8328075) q[2];
rz(-1.8874946) q[3];
sx q[3];
rz(-1.8697238) q[3];
sx q[3];
rz(1.1312283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18987385) q[0];
sx q[0];
rz(-1.254344) q[0];
sx q[0];
rz(-2.9217814) q[0];
rz(-1.1687357) q[1];
sx q[1];
rz(-1.3558148) q[1];
sx q[1];
rz(-1.8692325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.397906) q[0];
sx q[0];
rz(-0.96677654) q[0];
sx q[0];
rz(3.1118667) q[0];
rz(0.51080334) q[2];
sx q[2];
rz(-1.496721) q[2];
sx q[2];
rz(-1.6965716) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.90779623) q[1];
sx q[1];
rz(-1.7410454) q[1];
sx q[1];
rz(-2.7336804) q[1];
rz(-0.79325647) q[3];
sx q[3];
rz(-0.66325391) q[3];
sx q[3];
rz(0.59287567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1805588) q[2];
sx q[2];
rz(-2.2982633) q[2];
sx q[2];
rz(2.6832704) q[2];
rz(0.35442963) q[3];
sx q[3];
rz(-1.7731881) q[3];
sx q[3];
rz(1.1752769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4273222) q[0];
sx q[0];
rz(-1.9388119) q[0];
sx q[0];
rz(-0.74139968) q[0];
rz(1.4330014) q[1];
sx q[1];
rz(-2.7129136) q[1];
sx q[1];
rz(3.0523849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2994381) q[0];
sx q[0];
rz(-2.0426867) q[0];
sx q[0];
rz(0.20345511) q[0];
rz(-1.9544425) q[2];
sx q[2];
rz(-1.0133146) q[2];
sx q[2];
rz(-1.5526183) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8073106) q[1];
sx q[1];
rz(-0.34156964) q[1];
sx q[1];
rz(1.4463805) q[1];
rz(-pi) q[2];
rz(-2.4048526) q[3];
sx q[3];
rz(-1.2587446) q[3];
sx q[3];
rz(1.6583575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1409113) q[2];
sx q[2];
rz(-3.0912283) q[2];
sx q[2];
rz(-0.64765206) q[2];
rz(2.7378313) q[3];
sx q[3];
rz(-1.253456) q[3];
sx q[3];
rz(3.0103179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9443001) q[0];
sx q[0];
rz(-2.1680752) q[0];
sx q[0];
rz(-2.0808921) q[0];
rz(0.79509673) q[1];
sx q[1];
rz(-0.92082321) q[1];
sx q[1];
rz(-0.77650741) q[1];
rz(-0.9153824) q[2];
sx q[2];
rz(-1.4861098) q[2];
sx q[2];
rz(-2.8783023) q[2];
rz(2.285801) q[3];
sx q[3];
rz(-0.11920155) q[3];
sx q[3];
rz(-3.0536065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
