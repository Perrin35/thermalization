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
rz(2.4616315) q[0];
sx q[0];
rz(7.1890561) q[0];
sx q[0];
rz(11.472975) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(1.9688164) q[1];
sx q[1];
rz(9.8493675) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4571675) q[0];
sx q[0];
rz(-0.20523219) q[0];
sx q[0];
rz(-1.3825018) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0908747) q[2];
sx q[2];
rz(-1.4169622) q[2];
sx q[2];
rz(2.5209629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0478617) q[1];
sx q[1];
rz(-1.3252531) q[1];
sx q[1];
rz(0.081201038) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4154451) q[3];
sx q[3];
rz(-2.4253824) q[3];
sx q[3];
rz(0.20477141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37497416) q[2];
sx q[2];
rz(-1.5593636) q[2];
sx q[2];
rz(1.1589104) q[2];
rz(-2.9186987) q[3];
sx q[3];
rz(-1.4054207) q[3];
sx q[3];
rz(0.29002732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463773) q[0];
sx q[0];
rz(-1.615849) q[0];
sx q[0];
rz(-2.0624397) q[0];
rz(1.7201299) q[1];
sx q[1];
rz(-0.68669569) q[1];
sx q[1];
rz(-3.0770643) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9345029) q[0];
sx q[0];
rz(-0.10182589) q[0];
sx q[0];
rz(-2.1125211) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6653397) q[2];
sx q[2];
rz(-1.5956399) q[2];
sx q[2];
rz(-0.78643878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56386212) q[1];
sx q[1];
rz(-0.81474333) q[1];
sx q[1];
rz(-1.9063063) q[1];
rz(-2.9546176) q[3];
sx q[3];
rz(-2.3114853) q[3];
sx q[3];
rz(2.376698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1530389) q[2];
sx q[2];
rz(-1.3554074) q[2];
sx q[2];
rz(1.1174196) q[2];
rz(-2.2828263) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(0.43706885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7607464) q[0];
sx q[0];
rz(-0.28457156) q[0];
sx q[0];
rz(2.9685156) q[0];
rz(-2.1848047) q[1];
sx q[1];
rz(-1.0280949) q[1];
sx q[1];
rz(-2.079336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0883368) q[0];
sx q[0];
rz(-2.2940141) q[0];
sx q[0];
rz(3.1331314) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64086242) q[2];
sx q[2];
rz(-2.0483565) q[2];
sx q[2];
rz(-1.8030082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30327163) q[1];
sx q[1];
rz(-2.1049989) q[1];
sx q[1];
rz(2.3278585) q[1];
rz(-pi) q[2];
rz(0.84433664) q[3];
sx q[3];
rz(-0.82714836) q[3];
sx q[3];
rz(-1.3177804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7303077) q[2];
sx q[2];
rz(-1.9399119) q[2];
sx q[2];
rz(-2.3968706) q[2];
rz(2.5005285) q[3];
sx q[3];
rz(-1.791626) q[3];
sx q[3];
rz(2.6509269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915801) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(2.2963754) q[0];
rz(2.7382964) q[1];
sx q[1];
rz(-1.5396298) q[1];
sx q[1];
rz(-2.8135615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.302264) q[0];
sx q[0];
rz(-1.7516881) q[0];
sx q[0];
rz(2.9995287) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5045067) q[2];
sx q[2];
rz(-1.3332597) q[2];
sx q[2];
rz(1.4927166) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58807948) q[1];
sx q[1];
rz(-0.57826406) q[1];
sx q[1];
rz(-2.9248613) q[1];
rz(2.0045404) q[3];
sx q[3];
rz(-2.5528926) q[3];
sx q[3];
rz(-0.91982809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24568096) q[2];
sx q[2];
rz(-0.74378219) q[2];
sx q[2];
rz(-2.9823859) q[2];
rz(3.1253452) q[3];
sx q[3];
rz(-1.0280321) q[3];
sx q[3];
rz(-0.73133674) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3943587) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(2.0514945) q[0];
rz(0.12995003) q[1];
sx q[1];
rz(-0.81472412) q[1];
sx q[1];
rz(1.5257588) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14262202) q[0];
sx q[0];
rz(-0.71323133) q[0];
sx q[0];
rz(0.66340982) q[0];
rz(-pi) q[1];
rz(1.4065625) q[2];
sx q[2];
rz(-1.4964074) q[2];
sx q[2];
rz(2.1321572) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33541778) q[1];
sx q[1];
rz(-1.5982143) q[1];
sx q[1];
rz(-2.5078678) q[1];
rz(1.7618175) q[3];
sx q[3];
rz(-0.18409746) q[3];
sx q[3];
rz(0.43162089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3913564) q[2];
sx q[2];
rz(-2.3199234) q[2];
sx q[2];
rz(2.6643122) q[2];
rz(-2.9376049) q[3];
sx q[3];
rz(-0.19652772) q[3];
sx q[3];
rz(0.6240713) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1962965) q[0];
sx q[0];
rz(-0.81552234) q[0];
sx q[0];
rz(-2.3639009) q[0];
rz(0.6174736) q[1];
sx q[1];
rz(-1.6505046) q[1];
sx q[1];
rz(2.2106574) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5701911) q[0];
sx q[0];
rz(-2.5102315) q[0];
sx q[0];
rz(2.3124113) q[0];
rz(2.0495049) q[2];
sx q[2];
rz(-1.9472113) q[2];
sx q[2];
rz(0.37351721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9677862) q[1];
sx q[1];
rz(-1.5266982) q[1];
sx q[1];
rz(-1.2323061) q[1];
rz(-1.0369426) q[3];
sx q[3];
rz(-2.4483238) q[3];
sx q[3];
rz(1.1806837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4882539) q[2];
sx q[2];
rz(-1.3335756) q[2];
sx q[2];
rz(-0.2674357) q[2];
rz(0.93287647) q[3];
sx q[3];
rz(-1.2837912) q[3];
sx q[3];
rz(2.647184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47144181) q[0];
sx q[0];
rz(-1.143456) q[0];
sx q[0];
rz(0.66584051) q[0];
rz(-0.2039856) q[1];
sx q[1];
rz(-2.1603656) q[1];
sx q[1];
rz(-1.9542255) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46118051) q[0];
sx q[0];
rz(-1.102316) q[0];
sx q[0];
rz(1.4673244) q[0];
rz(-pi) q[1];
rz(-2.6186278) q[2];
sx q[2];
rz(-2.209216) q[2];
sx q[2];
rz(1.0467199) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5449808) q[1];
sx q[1];
rz(-2.2447733) q[1];
sx q[1];
rz(-1.2294212) q[1];
rz(1.8245876) q[3];
sx q[3];
rz(-1.521186) q[3];
sx q[3];
rz(-1.8269767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54619706) q[2];
sx q[2];
rz(-0.51837817) q[2];
sx q[2];
rz(2.4714244) q[2];
rz(0.3012805) q[3];
sx q[3];
rz(-1.8471085) q[3];
sx q[3];
rz(2.7231351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.30348521) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(-2.0759034) q[0];
rz(0.65713716) q[1];
sx q[1];
rz(-0.70825759) q[1];
sx q[1];
rz(-1.7117975) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3683726) q[0];
sx q[0];
rz(-0.67267679) q[0];
sx q[0];
rz(1.011959) q[0];
rz(2.6515342) q[2];
sx q[2];
rz(-0.37184162) q[2];
sx q[2];
rz(1.3272367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89420477) q[1];
sx q[1];
rz(-0.46890837) q[1];
sx q[1];
rz(2.6415178) q[1];
rz(-pi) q[2];
rz(-0.18971209) q[3];
sx q[3];
rz(-1.4290775) q[3];
sx q[3];
rz(0.82822463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8810205) q[2];
sx q[2];
rz(-1.842097) q[2];
sx q[2];
rz(2.1232429) q[2];
rz(1.5923422) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(-2.3840267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757979) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(-2.1739668) q[0];
rz(-2.8014917) q[1];
sx q[1];
rz(-0.83931559) q[1];
sx q[1];
rz(-1.0338773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32717413) q[0];
sx q[0];
rz(-2.5670739) q[0];
sx q[0];
rz(-1.3142725) q[0];
rz(-pi) q[1];
rz(-0.56615717) q[2];
sx q[2];
rz(-0.70933178) q[2];
sx q[2];
rz(-2.4475803) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.441923) q[1];
sx q[1];
rz(-1.1811387) q[1];
sx q[1];
rz(-0.81894919) q[1];
rz(-2.5646117) q[3];
sx q[3];
rz(-0.26538244) q[3];
sx q[3];
rz(-0.25091618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10244441) q[2];
sx q[2];
rz(-1.6951963) q[2];
sx q[2];
rz(-3.1026802) q[2];
rz(3.0012722) q[3];
sx q[3];
rz(-3.0571627) q[3];
sx q[3];
rz(1.8261568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4843531) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(-1.8835618) q[0];
rz(-0.21615061) q[1];
sx q[1];
rz(-2.436147) q[1];
sx q[1];
rz(-0.36453882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.011019) q[0];
sx q[0];
rz(-2.6219061) q[0];
sx q[0];
rz(1.0566061) q[0];
rz(-pi) q[1];
rz(-2.3380525) q[2];
sx q[2];
rz(-1.8992956) q[2];
sx q[2];
rz(-2.1568418) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26673022) q[1];
sx q[1];
rz(-1.5577496) q[1];
sx q[1];
rz(3.0860098) q[1];
rz(-pi) q[2];
rz(1.6031426) q[3];
sx q[3];
rz(-1.6286639) q[3];
sx q[3];
rz(-2.3889314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8517427) q[2];
sx q[2];
rz(-1.2052636) q[2];
sx q[2];
rz(0.64129788) q[2];
rz(1.4801721) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(-2.4680468) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4079473) q[0];
sx q[0];
rz(-0.49325627) q[0];
sx q[0];
rz(2.7706964) q[0];
rz(0.98450487) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(-1.7683239) q[2];
sx q[2];
rz(-1.6850204) q[2];
sx q[2];
rz(-1.9584283) q[2];
rz(-2.051864) q[3];
sx q[3];
rz(-1.0952598) q[3];
sx q[3];
rz(2.918603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
