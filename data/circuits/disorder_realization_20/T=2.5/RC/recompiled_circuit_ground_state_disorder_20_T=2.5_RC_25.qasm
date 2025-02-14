OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29937509) q[0];
sx q[0];
rz(-2.8111281) q[0];
sx q[0];
rz(2.0781031) q[0];
rz(-0.039634135) q[1];
sx q[1];
rz(-0.57365817) q[1];
sx q[1];
rz(-1.080245) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1606205) q[0];
sx q[0];
rz(-1.1078796) q[0];
sx q[0];
rz(-3.1391543) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3990381) q[2];
sx q[2];
rz(-1.556201) q[2];
sx q[2];
rz(0.43806048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.049223674) q[1];
sx q[1];
rz(-1.8300841) q[1];
sx q[1];
rz(-0.53569002) q[1];
rz(2.0548212) q[3];
sx q[3];
rz(-0.98376444) q[3];
sx q[3];
rz(-2.0984894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0237191) q[2];
sx q[2];
rz(-1.7261852) q[2];
sx q[2];
rz(-0.37303698) q[2];
rz(-2.5840664) q[3];
sx q[3];
rz(-1.5159461) q[3];
sx q[3];
rz(0.47835866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8525304) q[0];
sx q[0];
rz(-1.7076778) q[0];
sx q[0];
rz(-1.0093932) q[0];
rz(0.32762647) q[1];
sx q[1];
rz(-0.3265003) q[1];
sx q[1];
rz(1.1064628) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0026446) q[0];
sx q[0];
rz(-1.4382285) q[0];
sx q[0];
rz(1.5969921) q[0];
rz(-pi) q[1];
rz(2.4566133) q[2];
sx q[2];
rz(-1.5826591) q[2];
sx q[2];
rz(-1.6983183) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6514837) q[1];
sx q[1];
rz(-0.98020187) q[1];
sx q[1];
rz(0.57140731) q[1];
rz(2.7170466) q[3];
sx q[3];
rz(-0.44711061) q[3];
sx q[3];
rz(3.0101484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24836765) q[2];
sx q[2];
rz(-1.1838341) q[2];
sx q[2];
rz(0.66967213) q[2];
rz(1.1559961) q[3];
sx q[3];
rz(-2.7892734) q[3];
sx q[3];
rz(-0.34935752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.30592331) q[0];
sx q[0];
rz(-2.7140129) q[0];
sx q[0];
rz(1.3315573) q[0];
rz(-1.2003027) q[1];
sx q[1];
rz(-2.8889062) q[1];
sx q[1];
rz(-2.4694064) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6349307) q[0];
sx q[0];
rz(-1.4499393) q[0];
sx q[0];
rz(-2.6664227) q[0];
rz(-0.16803788) q[2];
sx q[2];
rz(-1.1716951) q[2];
sx q[2];
rz(-1.8519515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1609413) q[1];
sx q[1];
rz(-1.13492) q[1];
sx q[1];
rz(-1.3869882) q[1];
x q[2];
rz(2.658031) q[3];
sx q[3];
rz(-0.10933317) q[3];
sx q[3];
rz(-3.1121174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4868698) q[2];
sx q[2];
rz(-1.9630311) q[2];
sx q[2];
rz(-0.85050026) q[2];
rz(0.9507829) q[3];
sx q[3];
rz(-0.50539223) q[3];
sx q[3];
rz(0.90184414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1281857) q[0];
sx q[0];
rz(-0.81939092) q[0];
sx q[0];
rz(-0.73981458) q[0];
rz(-2.9402404) q[1];
sx q[1];
rz(-1.9722152) q[1];
sx q[1];
rz(0.22612017) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27557954) q[0];
sx q[0];
rz(-3.1034307) q[0];
sx q[0];
rz(-1.1665795) q[0];
x q[1];
rz(1.9490384) q[2];
sx q[2];
rz(-1.2653482) q[2];
sx q[2];
rz(-2.7482207) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0836843) q[1];
sx q[1];
rz(-2.1154566) q[1];
sx q[1];
rz(-2.1780531) q[1];
rz(-1.700541) q[3];
sx q[3];
rz(-1.1888388) q[3];
sx q[3];
rz(2.6568535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7390593) q[2];
sx q[2];
rz(-2.2465574) q[2];
sx q[2];
rz(2.1441937) q[2];
rz(-2.7282257) q[3];
sx q[3];
rz(-1.9316542) q[3];
sx q[3];
rz(2.0972142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0053552) q[0];
sx q[0];
rz(-2.4449466) q[0];
sx q[0];
rz(-0.053939017) q[0];
rz(-2.9593762) q[1];
sx q[1];
rz(-0.91582623) q[1];
sx q[1];
rz(2.838476) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16909105) q[0];
sx q[0];
rz(-2.005134) q[0];
sx q[0];
rz(-2.7624755) q[0];
rz(-1.3197754) q[2];
sx q[2];
rz(-0.66749882) q[2];
sx q[2];
rz(-0.72161822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50072155) q[1];
sx q[1];
rz(-0.56001511) q[1];
sx q[1];
rz(-1.9418094) q[1];
x q[2];
rz(0.037542779) q[3];
sx q[3];
rz(-2.3830066) q[3];
sx q[3];
rz(-2.9428218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6984685) q[2];
sx q[2];
rz(-1.8734525) q[2];
sx q[2];
rz(-2.7847086) q[2];
rz(-0.77480626) q[3];
sx q[3];
rz(-1.6322497) q[3];
sx q[3];
rz(-0.43865144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.88724941) q[0];
sx q[0];
rz(-1.9312504) q[0];
sx q[0];
rz(-0.78897011) q[0];
rz(-1.7987159) q[1];
sx q[1];
rz(-1.7306381) q[1];
sx q[1];
rz(-2.3416187) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8784284) q[0];
sx q[0];
rz(-1.557) q[0];
sx q[0];
rz(-0.8850125) q[0];
rz(1.7943013) q[2];
sx q[2];
rz(-2.433521) q[2];
sx q[2];
rz(-1.1184276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9467873) q[1];
sx q[1];
rz(-1.9459263) q[1];
sx q[1];
rz(1.4542143) q[1];
rz(-pi) q[2];
x q[2];
rz(0.81963934) q[3];
sx q[3];
rz(-2.5878518) q[3];
sx q[3];
rz(0.40938974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0903025) q[2];
sx q[2];
rz(-2.3919969) q[2];
sx q[2];
rz(-1.7898111) q[2];
rz(2.6734062) q[3];
sx q[3];
rz(-1.0951833) q[3];
sx q[3];
rz(-2.2075672) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7276329) q[0];
sx q[0];
rz(-2.7233044) q[0];
sx q[0];
rz(0.002451238) q[0];
rz(0.0062423627) q[1];
sx q[1];
rz(-2.7793482) q[1];
sx q[1];
rz(-0.35591602) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.064156) q[0];
sx q[0];
rz(-1.2250803) q[0];
sx q[0];
rz(2.2156209) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8628274) q[2];
sx q[2];
rz(-2.27072) q[2];
sx q[2];
rz(-0.88387671) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1568875) q[1];
sx q[1];
rz(-0.92159373) q[1];
sx q[1];
rz(-2.8125416) q[1];
rz(-pi) q[2];
rz(2.9690817) q[3];
sx q[3];
rz(-0.54076946) q[3];
sx q[3];
rz(1.6921761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1459085) q[2];
sx q[2];
rz(-1.9084787) q[2];
sx q[2];
rz(-1.0710867) q[2];
rz(2.5339825) q[3];
sx q[3];
rz(-1.1039609) q[3];
sx q[3];
rz(1.6149909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.0869658) q[0];
sx q[0];
rz(-0.93769756) q[0];
sx q[0];
rz(1.1338393) q[0];
rz(-1.0519823) q[1];
sx q[1];
rz(-1.6832422) q[1];
sx q[1];
rz(-0.7410616) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0042564226) q[0];
sx q[0];
rz(-1.3110054) q[0];
sx q[0];
rz(-0.57704837) q[0];
rz(-pi) q[1];
rz(-2.1846619) q[2];
sx q[2];
rz(-1.9939878) q[2];
sx q[2];
rz(-2.500071) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.929175) q[1];
sx q[1];
rz(-1.5467531) q[1];
sx q[1];
rz(-0.52530064) q[1];
rz(-pi) q[2];
rz(2.2005733) q[3];
sx q[3];
rz(-1.9004603) q[3];
sx q[3];
rz(0.31808149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1715601) q[2];
sx q[2];
rz(-2.4617608) q[2];
sx q[2];
rz(0.46530923) q[2];
rz(2.4252637) q[3];
sx q[3];
rz(-1.6445487) q[3];
sx q[3];
rz(-1.1819476) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61140907) q[0];
sx q[0];
rz(-2.1937328) q[0];
sx q[0];
rz(-1.9635669) q[0];
rz(2.1317962) q[1];
sx q[1];
rz(-1.806908) q[1];
sx q[1];
rz(-0.10733265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.619894) q[0];
sx q[0];
rz(-1.0891344) q[0];
sx q[0];
rz(-2.480471) q[0];
x q[1];
rz(-0.1520098) q[2];
sx q[2];
rz(-2.3526224) q[2];
sx q[2];
rz(2.0308354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3299249) q[1];
sx q[1];
rz(-2.1582542) q[1];
sx q[1];
rz(-1.8483551) q[1];
rz(1.1501649) q[3];
sx q[3];
rz(-0.64088837) q[3];
sx q[3];
rz(0.057614652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44059077) q[2];
sx q[2];
rz(-0.89459449) q[2];
sx q[2];
rz(-0.44006285) q[2];
rz(-0.20982404) q[3];
sx q[3];
rz(-1.9473636) q[3];
sx q[3];
rz(2.8275209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99301991) q[0];
sx q[0];
rz(-0.58605376) q[0];
sx q[0];
rz(-1.7589737) q[0];
rz(-0.67715174) q[1];
sx q[1];
rz(-1.4644198) q[1];
sx q[1];
rz(2.009353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181732) q[0];
sx q[0];
rz(-1.3480524) q[0];
sx q[0];
rz(1.3370418) q[0];
rz(-pi) q[1];
rz(-2.878849) q[2];
sx q[2];
rz(-2.3296156) q[2];
sx q[2];
rz(0.95974582) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1340227) q[1];
sx q[1];
rz(-1.5948199) q[1];
sx q[1];
rz(2.4244306) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4664425) q[3];
sx q[3];
rz(-1.0680001) q[3];
sx q[3];
rz(0.72489634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9558692) q[2];
sx q[2];
rz(-2.393554) q[2];
sx q[2];
rz(-1.7175187) q[2];
rz(-2.2636223) q[3];
sx q[3];
rz(-0.22017559) q[3];
sx q[3];
rz(-0.018208114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743185) q[0];
sx q[0];
rz(-1.6813288) q[0];
sx q[0];
rz(1.9718476) q[0];
rz(-2.7727903) q[1];
sx q[1];
rz(-1.7316876) q[1];
sx q[1];
rz(0.015451886) q[1];
rz(0.44023505) q[2];
sx q[2];
rz(-1.3095837) q[2];
sx q[2];
rz(-1.5038036) q[2];
rz(-0.5795547) q[3];
sx q[3];
rz(-1.2523635) q[3];
sx q[3];
rz(1.5750332) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
