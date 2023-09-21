OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4047591) q[0];
sx q[0];
rz(-1.7801378) q[0];
sx q[0];
rz(1.3786432) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(0.4508957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.742813) q[0];
sx q[0];
rz(-0.089988515) q[0];
sx q[0];
rz(-0.16422693) q[0];
x q[1];
rz(-1.0916753) q[2];
sx q[2];
rz(-1.6699104) q[2];
sx q[2];
rz(1.5168158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7887468) q[1];
sx q[1];
rz(-0.55410085) q[1];
sx q[1];
rz(-0.43177859) q[1];
rz(1.6879184) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(0.62108921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(0.84428865) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-2.5355693) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59250295) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(0.26309183) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3514254) q[0];
sx q[0];
rz(-1.5520436) q[0];
sx q[0];
rz(-0.055939527) q[0];
x q[1];
rz(0.19940168) q[2];
sx q[2];
rz(-1.5099031) q[2];
sx q[2];
rz(2.8859438) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26145229) q[1];
sx q[1];
rz(-1.4457236) q[1];
sx q[1];
rz(-2.2373881) q[1];
rz(-pi) q[2];
rz(2.4466189) q[3];
sx q[3];
rz(-1.7193828) q[3];
sx q[3];
rz(-2.9050764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(1.1594695) q[2];
rz(-0.37108478) q[3];
sx q[3];
rz(-1.6371195) q[3];
sx q[3];
rz(2.8306567) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3804669) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(0.80672112) q[0];
rz(2.9280248) q[1];
sx q[1];
rz(-2.6453306) q[1];
sx q[1];
rz(-0.82021964) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263111) q[0];
sx q[0];
rz(-1.6831241) q[0];
sx q[0];
rz(0.88322722) q[0];
rz(2.9224612) q[2];
sx q[2];
rz(-2.0543155) q[2];
sx q[2];
rz(2.6205274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.080938235) q[1];
sx q[1];
rz(-1.2186236) q[1];
sx q[1];
rz(-2.0209795) q[1];
rz(-pi) q[2];
rz(-0.1180325) q[3];
sx q[3];
rz(-1.1088088) q[3];
sx q[3];
rz(-0.0052009728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.31072581) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(0.15549774) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61313066) q[0];
sx q[0];
rz(-2.4202132) q[0];
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
rz(-2.5116918) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(3.1303309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.690372) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(-0.95552432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1530694) q[1];
sx q[1];
rz(-1.8338025) q[1];
sx q[1];
rz(2.4371229) q[1];
rz(0.96418013) q[3];
sx q[3];
rz(-0.88084953) q[3];
sx q[3];
rz(-2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-2.7992115) q[2];
rz(2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(-2.2763021) q[0];
rz(-1.226549) q[1];
sx q[1];
rz(-2.1523235) q[1];
sx q[1];
rz(-1.8409761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375992) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(-3.0689737) q[0];
rz(-pi) q[1];
rz(2.775035) q[2];
sx q[2];
rz(-1.6379116) q[2];
sx q[2];
rz(1.7442489) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.091150065) q[1];
sx q[1];
rz(-1.0346518) q[1];
sx q[1];
rz(-1.9552783) q[1];
rz(-pi) q[2];
rz(0.26228321) q[3];
sx q[3];
rz(-1.7119006) q[3];
sx q[3];
rz(1.2348246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0098003) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(2.9495083) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3451097) q[0];
sx q[0];
rz(-0.61426291) q[0];
sx q[0];
rz(0.011750301) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(-1.5884429) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4295411) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(2.3343711) q[0];
rz(-pi) q[1];
x q[1];
rz(0.021868869) q[2];
sx q[2];
rz(-1.2049897) q[2];
sx q[2];
rz(1.0794229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67571083) q[1];
sx q[1];
rz(-1.4859745) q[1];
sx q[1];
rz(2.431543) q[1];
x q[2];
rz(-2.2736069) q[3];
sx q[3];
rz(-1.5138953) q[3];
sx q[3];
rz(1.9753319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50756303) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.8590415) q[2];
rz(-1.3698618) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-0.67725956) q[0];
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
rz(-2.2718186) q[0];
sx q[0];
rz(-2.5442113) q[0];
sx q[0];
rz(0.13418829) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0095519) q[2];
sx q[2];
rz(-1.570895) q[2];
sx q[2];
rz(1.4591109) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6923185) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(-1.7983789) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7988775) q[3];
sx q[3];
rz(-1.2667155) q[3];
sx q[3];
rz(2.8976687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3892422) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(-2.1288669) q[2];
rz(1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(2.6543806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0174719) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(-0.69865984) q[0];
rz(2.0195122) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.2493856) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4835565) q[0];
sx q[0];
rz(-1.5244966) q[0];
sx q[0];
rz(-0.21090837) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0480568) q[2];
sx q[2];
rz(-1.153423) q[2];
sx q[2];
rz(2.0975031) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3932712) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(-0.022151532) q[1];
rz(-pi) q[2];
rz(0.28835339) q[3];
sx q[3];
rz(-0.95803146) q[3];
sx q[3];
rz(-1.8307277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8119048) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(-1.9630986) q[2];
rz(-1.4568436) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(-2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(-1.8485803) q[0];
rz(-1.4216084) q[1];
sx q[1];
rz(-2.1052108) q[1];
sx q[1];
rz(-2.5440149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7062089) q[0];
sx q[0];
rz(-1.4357114) q[0];
sx q[0];
rz(-3.1027017) q[0];
x q[1];
rz(-2.0963247) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(-1.8906821) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5850726) q[1];
sx q[1];
rz(-2.6954898) q[1];
sx q[1];
rz(0.82533605) q[1];
rz(1.4937917) q[3];
sx q[3];
rz(-0.98620755) q[3];
sx q[3];
rz(2.0921752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9188345) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.2333599) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.6433158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.047886588) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(0.023660252) q[0];
rz(2.1854782) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(0.67217174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8925079) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(3.1363048) q[0];
rz(-0.51151885) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(0.17009232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2031659) q[1];
sx q[1];
rz(-2.2397579) q[1];
sx q[1];
rz(1.1346362) q[1];
rz(-pi) q[2];
x q[2];
rz(0.43283312) q[3];
sx q[3];
rz(-1.4308617) q[3];
sx q[3];
rz(-0.39633745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6293634) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(1.5036748) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
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
