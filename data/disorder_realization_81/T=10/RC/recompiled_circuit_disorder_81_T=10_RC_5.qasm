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
rz(-1.6575939) q[1];
sx q[1];
rz(-0.4508957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987797) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(2.9773657) q[0];
rz(2.0499174) q[2];
sx q[2];
rz(-1.4716822) q[2];
sx q[2];
rz(1.6247768) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2921819) q[1];
sx q[1];
rz(-2.0690448) q[1];
sx q[1];
rz(-1.317418) q[1];
rz(0.095590683) q[3];
sx q[3];
rz(-2.2523237) q[3];
sx q[3];
rz(2.3694627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(-2.700581) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59250295) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(2.8785008) q[0];
rz(-0.94353765) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3514254) q[0];
sx q[0];
rz(-1.5895491) q[0];
sx q[0];
rz(-0.055939527) q[0];
rz(-0.19940168) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(0.25564889) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1520878) q[1];
sx q[1];
rz(-2.4651335) q[1];
sx q[1];
rz(1.771404) q[1];
x q[2];
rz(-2.9119592) q[3];
sx q[3];
rz(-2.4335055) q[3];
sx q[3];
rz(1.1585483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(1.9821232) q[2];
rz(0.37108478) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(2.8306567) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(-0.80672112) q[0];
rz(-2.9280248) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(-0.82021964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911644) q[0];
sx q[0];
rz(-0.69520742) q[0];
sx q[0];
rz(1.3948963) q[0];
rz(-pi) q[1];
rz(-0.21913146) q[2];
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
rz(1.8169176) q[1];
sx q[1];
rz(-1.9915238) q[1];
sx q[1];
rz(2.7540728) q[1];
x q[2];
rz(1.1060171) q[3];
sx q[3];
rz(-1.4651863) q[3];
sx q[3];
rz(-1.6184023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(0.93079981) q[2];
rz(2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.320425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5399649) q[0];
sx q[0];
rz(-2.7318582) q[0];
sx q[0];
rz(1.5967303) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7337012) q[2];
sx q[2];
rz(-2.6558999) q[2];
sx q[2];
rz(-2.8913468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.79975407) q[1];
sx q[1];
rz(-2.2463887) q[1];
sx q[1];
rz(-1.231133) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5370595) q[3];
sx q[3];
rz(-0.88450888) q[3];
sx q[3];
rz(1.3674919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1057672) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(2.7992115) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8300366) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(1.3006166) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9601701) q[0];
sx q[0];
rz(-1.9248065) q[0];
sx q[0];
rz(1.5976853) q[0];
rz(-1.498921) q[2];
sx q[2];
rz(-1.9364898) q[2];
sx q[2];
rz(-2.9938811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.865766) q[1];
sx q[1];
rz(-1.2424801) q[1];
sx q[1];
rz(0.57002108) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50118581) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(2.3230769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.9921781) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(-0.19208433) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3451097) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(0.011750301) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5884429) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6019183) q[0];
sx q[0];
rz(-1.9214905) q[0];
sx q[0];
rz(-0.35777103) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.021868869) q[2];
sx q[2];
rz(-1.936603) q[2];
sx q[2];
rz(1.0794229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1482684) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(3.0118914) q[1];
rz(2.2736069) q[3];
sx q[3];
rz(-1.5138953) q[3];
sx q[3];
rz(1.1662607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(1.8590415) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-2.0231358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(-2.989785) q[1];
sx q[1];
rz(-1.7671403) q[1];
sx q[1];
rz(0.97704926) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58986321) q[0];
sx q[0];
rz(-1.4954733) q[0];
sx q[0];
rz(-2.548404) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13204079) q[2];
sx q[2];
rz(-1.5706976) q[2];
sx q[2];
rz(1.6824818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44927412) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(1.7983789) q[1];
x q[2];
rz(0.75120039) q[3];
sx q[3];
rz(-2.6874472) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7523505) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(-0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(1.1220804) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(1.2493856) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.218924) q[0];
sx q[0];
rz(-1.7814753) q[0];
sx q[0];
rz(1.6181437) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6685733) q[2];
sx q[2];
rz(-2.0447391) q[2];
sx q[2];
rz(-0.75616403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5640806) q[1];
sx q[1];
rz(-0.022833303) q[1];
sx q[1];
rz(2.8965685) q[1];
rz(-pi) q[2];
rz(-0.28835339) q[3];
sx q[3];
rz(-2.1835612) q[3];
sx q[3];
rz(1.310865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8119048) q[2];
sx q[2];
rz(-2.3554282) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(1.684749) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6417398) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.2930124) q[0];
rz(1.7199843) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(-2.5440149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4353838) q[0];
sx q[0];
rz(-1.4357114) q[0];
sx q[0];
rz(-3.1027017) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75002807) q[2];
sx q[2];
rz(-2.4366597) q[2];
sx q[2];
rz(2.7761369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7087047) q[1];
sx q[1];
rz(-1.8678027) q[1];
sx q[1];
rz(-1.2328641) q[1];
x q[2];
rz(-0.11573128) q[3];
sx q[3];
rz(-2.5525408) q[3];
sx q[3];
rz(2.2310886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.22275816) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(1.9082327) q[2];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047886588) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(0.023660252) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(2.4694209) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8925079) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(0.0052878629) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5615084) q[2];
sx q[2];
rz(-2.0822968) q[2];
sx q[2];
rz(1.4052504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2939261) q[1];
sx q[1];
rz(-2.3617509) q[1];
sx q[1];
rz(2.6508209) q[1];
rz(-pi) q[2];
rz(1.4168596) q[3];
sx q[3];
rz(-1.142475) q[3];
sx q[3];
rz(1.1101013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(2.4178986) q[1];
sx q[1];
rz(-2.1543398) q[1];
sx q[1];
rz(2.2347246) q[1];
rz(1.9696708) q[2];
sx q[2];
rz(-1.7181859) q[2];
sx q[2];
rz(-2.0167375) q[2];
rz(-2.8430812) q[3];
sx q[3];
rz(-1.1805503) q[3];
sx q[3];
rz(-1.3211484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];