OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9508096) q[0];
sx q[0];
rz(-1.8291744) q[0];
sx q[0];
rz(1.4825312) q[0];
x q[1];
rz(1.0563645) q[2];
sx q[2];
rz(-1.3008413) q[2];
sx q[2];
rz(-1.1411238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8618776) q[1];
sx q[1];
rz(-1.4289083) q[1];
sx q[1];
rz(3.0008297) q[1];
rz(0.78731491) q[3];
sx q[3];
rz(-1.6978605) q[3];
sx q[3];
rz(0.78391677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.67291659) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(2.2093175) q[2];
rz(-2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(0.96536243) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7070049) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(-2.7804651) q[0];
rz(1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(-0.8180058) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084413962) q[0];
sx q[0];
rz(-1.1857496) q[0];
sx q[0];
rz(1.8684698) q[0];
rz(3.0188574) q[2];
sx q[2];
rz(-1.2018179) q[2];
sx q[2];
rz(0.69532794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4600196) q[1];
sx q[1];
rz(-2.1732554) q[1];
sx q[1];
rz(2.6006992) q[1];
rz(-pi) q[2];
rz(1.2070451) q[3];
sx q[3];
rz(-1.517429) q[3];
sx q[3];
rz(1.766891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8890185) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(-1.3160926) q[3];
sx q[3];
rz(-2.3836453) q[3];
sx q[3];
rz(2.7533598) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7398359) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(2.8443764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39697159) q[0];
sx q[0];
rz(-2.731785) q[0];
sx q[0];
rz(-2.9817392) q[0];
rz(-0.77135135) q[2];
sx q[2];
rz(-2.1043679) q[2];
sx q[2];
rz(0.16650621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.19993648) q[1];
sx q[1];
rz(-1.3841277) q[1];
sx q[1];
rz(3.093064) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2753784) q[3];
sx q[3];
rz(-2.2598887) q[3];
sx q[3];
rz(0.44486526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7685984) q[2];
sx q[2];
rz(-0.82295376) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(1.1887431) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(2.6141613) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(-2.3676681) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(0.68177044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9457152) q[0];
sx q[0];
rz(-1.8316852) q[0];
sx q[0];
rz(2.9257141) q[0];
rz(2.1749928) q[2];
sx q[2];
rz(-1.3912429) q[2];
sx q[2];
rz(2.0008848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75979739) q[1];
sx q[1];
rz(-1.1361546) q[1];
sx q[1];
rz(-0.97310193) q[1];
rz(-pi) q[2];
rz(-2.3247129) q[3];
sx q[3];
rz(-2.6602392) q[3];
sx q[3];
rz(-0.64827418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.008808) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(1.0664252) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(-2.5812896) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(-1.5195742) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9902089) q[0];
sx q[0];
rz(-2.4575893) q[0];
sx q[0];
rz(-2.3955976) q[0];
rz(-0.72603307) q[2];
sx q[2];
rz(-0.75349977) q[2];
sx q[2];
rz(1.3940648) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2625761) q[1];
sx q[1];
rz(-2.3153956) q[1];
sx q[1];
rz(-1.8415585) q[1];
rz(-3.0171379) q[3];
sx q[3];
rz(-2.180047) q[3];
sx q[3];
rz(2.5973158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(-0.2229283) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(-3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3361622) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(-2.0715332) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(1.8575352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239913) q[0];
sx q[0];
rz(-1.0214897) q[0];
sx q[0];
rz(-2.4248289) q[0];
x q[1];
rz(1.9610923) q[2];
sx q[2];
rz(-1.830606) q[2];
sx q[2];
rz(-0.63894546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9027054) q[1];
sx q[1];
rz(-2.0327912) q[1];
sx q[1];
rz(-1.8932896) q[1];
x q[2];
rz(2.0606023) q[3];
sx q[3];
rz(-2.5488857) q[3];
sx q[3];
rz(-0.61471516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47508919) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(0.63759032) q[2];
rz(0.83667886) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(-1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(0.93820757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0465614) q[0];
sx q[0];
rz(-1.1633658) q[0];
sx q[0];
rz(0.33336063) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4110255) q[2];
sx q[2];
rz(-2.0148723) q[2];
sx q[2];
rz(0.72545746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92749121) q[1];
sx q[1];
rz(-0.92988211) q[1];
sx q[1];
rz(-0.44968857) q[1];
rz(1.8754962) q[3];
sx q[3];
rz(-0.88677553) q[3];
sx q[3];
rz(-1.5414433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87970916) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(1.7555457) q[2];
rz(1.322768) q[3];
sx q[3];
rz(-1.991792) q[3];
sx q[3];
rz(3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8742074) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(-1.2738312) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(0.83126718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76038218) q[0];
sx q[0];
rz(-1.7179278) q[0];
sx q[0];
rz(-2.7358664) q[0];
x q[1];
rz(0.12886329) q[2];
sx q[2];
rz(-2.2017751) q[2];
sx q[2];
rz(-1.250759) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2483406) q[1];
sx q[1];
rz(-0.89683956) q[1];
sx q[1];
rz(-2.2158951) q[1];
x q[2];
rz(2.3723699) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(-3.0358918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(0.9643628) q[2];
rz(-1.9780805) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7678541) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(0.97737616) q[0];
rz(1.7550229) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(-1.1057373) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2290105) q[0];
sx q[0];
rz(-1.454103) q[0];
sx q[0];
rz(-2.6936561) q[0];
rz(-pi) q[1];
rz(2.2374723) q[2];
sx q[2];
rz(-1.4881926) q[2];
sx q[2];
rz(0.77043515) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8399664) q[1];
sx q[1];
rz(-1.3857538) q[1];
sx q[1];
rz(1.806083) q[1];
rz(-2.4615199) q[3];
sx q[3];
rz(-2.701093) q[3];
sx q[3];
rz(-2.3624453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-0.14979714) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(1.3214553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6479284) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(-0.1272442) q[0];
rz(1.4808902) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0509335) q[0];
sx q[0];
rz(-1.5899842) q[0];
sx q[0];
rz(0.030254342) q[0];
rz(2.9211505) q[2];
sx q[2];
rz(-2.5219005) q[2];
sx q[2];
rz(0.1534136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2044636) q[1];
sx q[1];
rz(-2.5531054) q[1];
sx q[1];
rz(0.012339331) q[1];
x q[2];
rz(-1.6889257) q[3];
sx q[3];
rz(-1.6932634) q[3];
sx q[3];
rz(0.99692217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1148791) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(-0.090027697) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5933843) q[0];
sx q[0];
rz(-1.1593288) q[0];
sx q[0];
rz(2.819084) q[0];
rz(-0.38800115) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(-2.3315196) q[2];
sx q[2];
rz(-1.0675061) q[2];
sx q[2];
rz(1.6455417) q[2];
rz(-1.2896982) q[3];
sx q[3];
rz(-2.5847808) q[3];
sx q[3];
rz(1.9980711) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
