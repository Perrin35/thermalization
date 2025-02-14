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
rz(0.82290736) q[0];
sx q[0];
rz(-0.35879254) q[0];
sx q[0];
rz(0.86831492) q[0];
rz(-1.4153642) q[1];
sx q[1];
rz(-1.1338898) q[1];
sx q[1];
rz(-0.99376065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39696898) q[0];
sx q[0];
rz(-0.41797371) q[0];
sx q[0];
rz(-0.28095066) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6401071) q[2];
sx q[2];
rz(-1.9619313) q[2];
sx q[2];
rz(-0.30539612) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6958625) q[1];
sx q[1];
rz(-2.5035586) q[1];
sx q[1];
rz(-2.3444434) q[1];
rz(-1.078061) q[3];
sx q[3];
rz(-2.6542814) q[3];
sx q[3];
rz(2.010797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22780861) q[2];
sx q[2];
rz(-1.8081534) q[2];
sx q[2];
rz(0.58161962) q[2];
rz(2.4070814) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(-2.7233126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.2468579) q[0];
sx q[0];
rz(-1.7623836) q[0];
sx q[0];
rz(-0.49767622) q[0];
rz(1.0408164) q[1];
sx q[1];
rz(-2.7383995) q[1];
sx q[1];
rz(-1.9042447) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055641551) q[0];
sx q[0];
rz(-2.4179672) q[0];
sx q[0];
rz(3.0776204) q[0];
x q[1];
rz(-2.9534229) q[2];
sx q[2];
rz(-1.2154748) q[2];
sx q[2];
rz(-2.4057092) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.50857568) q[1];
sx q[1];
rz(-2.687722) q[1];
sx q[1];
rz(1.244844) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26555221) q[3];
sx q[3];
rz(-1.3612116) q[3];
sx q[3];
rz(-0.095106212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70025468) q[2];
sx q[2];
rz(-2.4105218) q[2];
sx q[2];
rz(2.8311484) q[2];
rz(1.9849518) q[3];
sx q[3];
rz(-2.0168596) q[3];
sx q[3];
rz(-1.1667075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0079086) q[0];
sx q[0];
rz(-2.5569361) q[0];
sx q[0];
rz(-2.7957918) q[0];
rz(2.8969104) q[1];
sx q[1];
rz(-0.9318277) q[1];
sx q[1];
rz(1.4366879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61621746) q[0];
sx q[0];
rz(-1.8843009) q[0];
sx q[0];
rz(1.1275191) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1780924) q[2];
sx q[2];
rz(-1.0356734) q[2];
sx q[2];
rz(1.2924043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2776041) q[1];
sx q[1];
rz(-0.79452786) q[1];
sx q[1];
rz(1.2717461) q[1];
x q[2];
rz(1.9448024) q[3];
sx q[3];
rz(-2.0197649) q[3];
sx q[3];
rz(-1.506492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0077385) q[2];
sx q[2];
rz(-1.1275007) q[2];
sx q[2];
rz(0.63068843) q[2];
rz(-0.33356365) q[3];
sx q[3];
rz(-2.1400698) q[3];
sx q[3];
rz(2.1400616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075060464) q[0];
sx q[0];
rz(-0.94828951) q[0];
sx q[0];
rz(-1.7370268) q[0];
rz(-2.7724077) q[1];
sx q[1];
rz(-1.4063947) q[1];
sx q[1];
rz(-0.085748347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4082216) q[0];
sx q[0];
rz(-1.2722172) q[0];
sx q[0];
rz(-2.7135506) q[0];
rz(0.71620415) q[2];
sx q[2];
rz(-2.5599481) q[2];
sx q[2];
rz(0.003613506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.73832522) q[1];
sx q[1];
rz(-1.333263) q[1];
sx q[1];
rz(1.0942142) q[1];
rz(-pi) q[2];
rz(-1.2530009) q[3];
sx q[3];
rz(-1.9488584) q[3];
sx q[3];
rz(-1.8658569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.839445) q[2];
sx q[2];
rz(-1.9044694) q[2];
sx q[2];
rz(-0.5298003) q[2];
rz(3.1032622) q[3];
sx q[3];
rz(-2.4119792) q[3];
sx q[3];
rz(-1.5995601) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2235276) q[0];
sx q[0];
rz(-2.6730838) q[0];
sx q[0];
rz(-1.202762) q[0];
rz(0.27944061) q[1];
sx q[1];
rz(-2.1165106) q[1];
sx q[1];
rz(0.7739982) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5651754) q[0];
sx q[0];
rz(-0.87067184) q[0];
sx q[0];
rz(-0.41482523) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8061637) q[2];
sx q[2];
rz(-2.895439) q[2];
sx q[2];
rz(-2.2858562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3995838) q[1];
sx q[1];
rz(-3.0019433) q[1];
sx q[1];
rz(-1.4781471) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0027917248) q[3];
sx q[3];
rz(-2.1062024) q[3];
sx q[3];
rz(1.1100779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0224132) q[2];
sx q[2];
rz(-0.1410307) q[2];
sx q[2];
rz(2.5775583) q[2];
rz(-0.65308475) q[3];
sx q[3];
rz(-1.1354732) q[3];
sx q[3];
rz(2.691332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7804724) q[0];
sx q[0];
rz(-2.5887964) q[0];
sx q[0];
rz(-0.64315382) q[0];
rz(-1.9505352) q[1];
sx q[1];
rz(-1.4570313) q[1];
sx q[1];
rz(-0.65470421) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768652) q[0];
sx q[0];
rz(-2.961359) q[0];
sx q[0];
rz(0.36156543) q[0];
x q[1];
rz(-0.68393647) q[2];
sx q[2];
rz(-1.6003055) q[2];
sx q[2];
rz(-2.3060407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5712329) q[1];
sx q[1];
rz(-2.1132937) q[1];
sx q[1];
rz(-3.065055) q[1];
rz(-1.1804462) q[3];
sx q[3];
rz(-1.6359513) q[3];
sx q[3];
rz(0.81126838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6341298) q[2];
sx q[2];
rz(-0.37709388) q[2];
sx q[2];
rz(1.4748352) q[2];
rz(-2.6569488) q[3];
sx q[3];
rz(-0.99769297) q[3];
sx q[3];
rz(1.2887597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6396879) q[0];
sx q[0];
rz(-0.26308331) q[0];
sx q[0];
rz(2.4581773) q[0];
rz(-0.13038334) q[1];
sx q[1];
rz(-1.5510635) q[1];
sx q[1];
rz(3.108976) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9630229) q[0];
sx q[0];
rz(-1.3953598) q[0];
sx q[0];
rz(-2.1397207) q[0];
rz(0.24222272) q[2];
sx q[2];
rz(-2.4413779) q[2];
sx q[2];
rz(1.8699353) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2523007) q[1];
sx q[1];
rz(-1.3056763) q[1];
sx q[1];
rz(-0.17269453) q[1];
rz(-pi) q[2];
rz(-2.272761) q[3];
sx q[3];
rz(-1.371939) q[3];
sx q[3];
rz(1.5031888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7183097) q[2];
sx q[2];
rz(-2.0330567) q[2];
sx q[2];
rz(-2.5763467) q[2];
rz(1.7806753) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-0.81418532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3703506) q[0];
sx q[0];
rz(-3.0841565) q[0];
sx q[0];
rz(0.29956079) q[0];
rz(1.4467422) q[1];
sx q[1];
rz(-2.2694777) q[1];
sx q[1];
rz(-0.95132336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0438536) q[0];
sx q[0];
rz(-1.7758177) q[0];
sx q[0];
rz(-0.033173843) q[0];
rz(1.5278242) q[2];
sx q[2];
rz(-2.9298721) q[2];
sx q[2];
rz(-2.4076902) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.044477) q[1];
sx q[1];
rz(-1.591178) q[1];
sx q[1];
rz(-0.42428942) q[1];
rz(-pi) q[2];
rz(-1.3090012) q[3];
sx q[3];
rz(-1.4859395) q[3];
sx q[3];
rz(-0.020430662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1424554) q[2];
sx q[2];
rz(-0.74644011) q[2];
sx q[2];
rz(-0.26774055) q[2];
rz(1.0033876) q[3];
sx q[3];
rz(-0.95971862) q[3];
sx q[3];
rz(-0.80823922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35850152) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(-0.95712334) q[0];
rz(-0.3793017) q[1];
sx q[1];
rz(-0.4117659) q[1];
sx q[1];
rz(-2.9764825) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6388164) q[0];
sx q[0];
rz(-0.90327016) q[0];
sx q[0];
rz(1.1291885) q[0];
rz(-pi) q[1];
rz(0.050641955) q[2];
sx q[2];
rz(-2.0426742) q[2];
sx q[2];
rz(-2.7993921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.86805008) q[1];
sx q[1];
rz(-0.99813491) q[1];
sx q[1];
rz(2.9258201) q[1];
rz(-pi) q[2];
rz(-1.9103138) q[3];
sx q[3];
rz(-0.63106189) q[3];
sx q[3];
rz(-1.5724044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86753201) q[2];
sx q[2];
rz(-1.910285) q[2];
sx q[2];
rz(-2.3629698) q[2];
rz(2.8042931) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(1.5571099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.857665) q[0];
sx q[0];
rz(-2.0511257) q[0];
sx q[0];
rz(2.0528059) q[0];
rz(2.7133443) q[1];
sx q[1];
rz(-2.1183522) q[1];
sx q[1];
rz(2.4931989) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79660611) q[0];
sx q[0];
rz(-1.8341176) q[0];
sx q[0];
rz(-0.4395606) q[0];
x q[1];
rz(-1.4444541) q[2];
sx q[2];
rz(-0.34046945) q[2];
sx q[2];
rz(2.7012417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8377503) q[1];
sx q[1];
rz(-2.492202) q[1];
sx q[1];
rz(-3.0719195) q[1];
rz(-pi) q[2];
rz(-1.4565868) q[3];
sx q[3];
rz(-2.1314959) q[3];
sx q[3];
rz(-1.6778835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26583656) q[2];
sx q[2];
rz(-1.0886322) q[2];
sx q[2];
rz(-2.9618373) q[2];
rz(-1.9836551) q[3];
sx q[3];
rz(-1.4561184) q[3];
sx q[3];
rz(-1.9013083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4614048) q[0];
sx q[0];
rz(-0.15459939) q[0];
sx q[0];
rz(0.79994487) q[0];
rz(1.1663306) q[1];
sx q[1];
rz(-2.1569398) q[1];
sx q[1];
rz(0.91469761) q[1];
rz(1.2733172) q[2];
sx q[2];
rz(-1.3354206) q[2];
sx q[2];
rz(-0.32257783) q[2];
rz(-0.33954444) q[3];
sx q[3];
rz(-1.3673269) q[3];
sx q[3];
rz(0.23585503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
