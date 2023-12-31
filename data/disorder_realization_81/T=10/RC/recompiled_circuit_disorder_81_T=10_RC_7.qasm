OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(2.2840075) q[1];
sx q[1];
rz(4.6255914) q[1];
sx q[1];
rz(8.9738823) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.742813) q[0];
sx q[0];
rz(-3.0516041) q[0];
sx q[0];
rz(0.16422693) q[0];
rz(2.0499174) q[2];
sx q[2];
rz(-1.4716822) q[2];
sx q[2];
rz(1.6247768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3528459) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(0.43177859) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(2.5205034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(-2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59250295) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-0.26309183) q[0];
rz(-0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037558) q[0];
sx q[0];
rz(-0.058996011) q[0];
sx q[0];
rz(-2.8179413) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19940168) q[2];
sx q[2];
rz(-1.6316895) q[2];
sx q[2];
rz(-2.8859438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7344208) q[1];
sx q[1];
rz(-2.2312575) q[1];
sx q[1];
rz(-2.9829626) q[1];
rz(-pi) q[2];
rz(-1.7632742) q[3];
sx q[3];
rz(-2.2566183) q[3];
sx q[3];
rz(-1.6845077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-2.1388781) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(-2.7705079) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(2.3348715) q[0];
rz(-2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(0.82021964) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41528156) q[0];
sx q[0];
rz(-1.6831241) q[0];
sx q[0];
rz(0.88322722) q[0];
rz(-pi) q[1];
rz(2.9224612) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(2.6205274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.324675) q[1];
sx q[1];
rz(-1.9915238) q[1];
sx q[1];
rz(-2.7540728) q[1];
x q[2];
rz(-1.1060171) q[3];
sx q[3];
rz(-1.4651863) q[3];
sx q[3];
rz(-1.5231903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(-1.320425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60162773) q[0];
sx q[0];
rz(-0.40973445) q[0];
sx q[0];
rz(-1.5967303) q[0];
rz(-1.7772061) q[2];
sx q[2];
rz(-2.0136535) q[2];
sx q[2];
rz(0.70476156) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8561613) q[1];
sx q[1];
rz(-0.74401765) q[1];
sx q[1];
rz(-2.7475949) q[1];
x q[2];
rz(0.96418013) q[3];
sx q[3];
rz(-0.88084953) q[3];
sx q[3];
rz(-2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(2.7992115) q[2];
rz(2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(-2.0006196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(0.3115561) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(1.3006166) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1039935) q[0];
sx q[0];
rz(-2.7866057) q[0];
sx q[0];
rz(-0.07261891) q[0];
rz(-pi) q[1];
x q[1];
rz(2.956203) q[2];
sx q[2];
rz(-2.7692147) q[2];
sx q[2];
rz(-2.7951954) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.865766) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(-0.57002108) q[1];
x q[2];
rz(0.26228321) q[3];
sx q[3];
rz(-1.4296921) q[3];
sx q[3];
rz(-1.2348246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0098003) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(-2.664393) q[2];
rz(-2.9495083) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451097) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(-3.1298424) q[0];
rz(-2.5911962) q[1];
sx q[1];
rz(-1.3563211) q[1];
sx q[1];
rz(-1.5531497) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71205157) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(-2.3343711) q[0];
rz(-pi) q[1];
rz(1.2049098) q[2];
sx q[2];
rz(-1.5503746) q[2];
sx q[2];
rz(2.6423955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1482684) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(3.0118914) q[1];
rz(-2.2736069) q[3];
sx q[3];
rz(-1.5138953) q[3];
sx q[3];
rz(1.9753319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(0.67725956) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(-2.1645434) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2718186) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(-0.13418829) q[0];
rz(-pi) q[1];
rz(1.5708959) q[2];
sx q[2];
rz(-1.4387555) q[2];
sx q[2];
rz(3.0299203) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3358826) q[1];
sx q[1];
rz(-1.6478331) q[1];
sx q[1];
rz(-1.2309993) q[1];
rz(-2.3903923) q[3];
sx q[3];
rz(-2.6874472) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(1.0127257) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(-2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.1241207) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(0.69865984) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69985336) q[0];
sx q[0];
rz(-2.925736) q[0];
sx q[0];
rz(2.9237843) q[0];
rz(-2.2970389) q[2];
sx q[2];
rz(-0.65659467) q[2];
sx q[2];
rz(3.055228) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57751209) q[1];
sx q[1];
rz(-0.022833303) q[1];
sx q[1];
rz(-2.8965685) q[1];
x q[2];
rz(-1.9551679) q[3];
sx q[3];
rz(-2.4723408) q[3];
sx q[3];
rz(-0.8347019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.32968783) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.2930124) q[0];
rz(1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(2.5440149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7062089) q[0];
sx q[0];
rz(-1.7058813) q[0];
sx q[0];
rz(-0.038890966) q[0];
rz(2.3915646) q[2];
sx q[2];
rz(-0.70493297) q[2];
sx q[2];
rz(-2.7761369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55652009) q[1];
sx q[1];
rz(-0.44610281) q[1];
sx q[1];
rz(2.3162566) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5556373) q[3];
sx q[3];
rz(-1.5065985) q[3];
sx q[3];
rz(0.5639329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(-1.2333599) q[2];
rz(-0.90138609) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-0.67217174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925079) q[0];
sx q[0];
rz(-1.5810284) q[0];
sx q[0];
rz(-0.0052878629) q[0];
rz(-1.5800843) q[2];
sx q[2];
rz(-1.0592959) q[2];
sx q[2];
rz(1.7363422) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9384267) q[1];
sx q[1];
rz(-0.90183479) q[1];
sx q[1];
rz(-2.0069564) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.724733) q[3];
sx q[3];
rz(-1.142475) q[3];
sx q[3];
rz(-2.0314914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.9696708) q[2];
sx q[2];
rz(-1.4234067) q[2];
sx q[2];
rz(1.1248551) q[2];
rz(0.29851144) q[3];
sx q[3];
rz(-1.1805503) q[3];
sx q[3];
rz(-1.3211484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
