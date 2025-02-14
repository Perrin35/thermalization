OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57269078) q[0];
sx q[0];
rz(-2.1532018) q[0];
sx q[0];
rz(0.44901499) q[0];
rz(-0.16945893) q[1];
sx q[1];
rz(0.13154498) q[1];
sx q[1];
rz(8.2934525) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1695263) q[0];
sx q[0];
rz(-0.6837877) q[0];
sx q[0];
rz(2.0440566) q[0];
rz(-pi) q[1];
rz(-0.63458459) q[2];
sx q[2];
rz(-2.0014844) q[2];
sx q[2];
rz(2.787295) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0180964) q[1];
sx q[1];
rz(-0.61401788) q[1];
sx q[1];
rz(1.2662925) q[1];
x q[2];
rz(-1.851397) q[3];
sx q[3];
rz(-0.56125703) q[3];
sx q[3];
rz(-3.0292976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.955287) q[2];
sx q[2];
rz(-2.7489642) q[2];
sx q[2];
rz(3.0379831) q[2];
rz(-0.3668395) q[3];
sx q[3];
rz(-1.6678526) q[3];
sx q[3];
rz(-1.3175255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1845301) q[0];
sx q[0];
rz(-2.6657031) q[0];
sx q[0];
rz(2.6153508) q[0];
rz(-1.2435675) q[1];
sx q[1];
rz(-1.4105816) q[1];
sx q[1];
rz(2.9454561) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.014482) q[0];
sx q[0];
rz(-2.0271993) q[0];
sx q[0];
rz(3.0672202) q[0];
rz(-2.4320658) q[2];
sx q[2];
rz(-1.4440184) q[2];
sx q[2];
rz(-1.7293872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89333234) q[1];
sx q[1];
rz(-0.47139097) q[1];
sx q[1];
rz(-1.253528) q[1];
rz(-pi) q[2];
rz(-0.72153458) q[3];
sx q[3];
rz(-3.1258158) q[3];
sx q[3];
rz(1.1392913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5940932) q[2];
sx q[2];
rz(-2.8671691) q[2];
sx q[2];
rz(-2.7653232) q[2];
rz(2.2606692) q[3];
sx q[3];
rz(-1.8062402) q[3];
sx q[3];
rz(2.3295565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6905717) q[0];
sx q[0];
rz(-2.8132827) q[0];
sx q[0];
rz(-2.1300533) q[0];
rz(2.4729589) q[1];
sx q[1];
rz(-1.1625544) q[1];
sx q[1];
rz(-2.5943601) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1062463) q[0];
sx q[0];
rz(-1.5515455) q[0];
sx q[0];
rz(-0.20142844) q[0];
x q[1];
rz(0.34788068) q[2];
sx q[2];
rz(-0.82150412) q[2];
sx q[2];
rz(2.2583928) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32729748) q[1];
sx q[1];
rz(-1.3298099) q[1];
sx q[1];
rz(-0.90118932) q[1];
x q[2];
rz(2.3272334) q[3];
sx q[3];
rz(-1.4297155) q[3];
sx q[3];
rz(0.94179487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68613595) q[2];
sx q[2];
rz(-0.28527173) q[2];
sx q[2];
rz(-1.6395462) q[2];
rz(2.5655668) q[3];
sx q[3];
rz(-1.4833769) q[3];
sx q[3];
rz(2.7801133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107002) q[0];
sx q[0];
rz(-0.80131131) q[0];
sx q[0];
rz(0.33962387) q[0];
rz(-0.54061186) q[1];
sx q[1];
rz(-2.4410591) q[1];
sx q[1];
rz(2.2115754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41893533) q[0];
sx q[0];
rz(-1.5003073) q[0];
sx q[0];
rz(-1.427729) q[0];
rz(-1.1851063) q[2];
sx q[2];
rz(-2.1197332) q[2];
sx q[2];
rz(-1.576129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8254357) q[1];
sx q[1];
rz(-0.47940578) q[1];
sx q[1];
rz(-1.5569219) q[1];
rz(-pi) q[2];
rz(0.50765462) q[3];
sx q[3];
rz(-0.99170199) q[3];
sx q[3];
rz(2.7217025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0129619) q[2];
sx q[2];
rz(-1.4350812) q[2];
sx q[2];
rz(-1.5677412) q[2];
rz(-0.74357998) q[3];
sx q[3];
rz(-0.79135528) q[3];
sx q[3];
rz(-2.1775406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6745233) q[0];
sx q[0];
rz(-0.85551298) q[0];
sx q[0];
rz(3.0849482) q[0];
rz(1.4757587) q[1];
sx q[1];
rz(-2.00878) q[1];
sx q[1];
rz(-1.3585565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18902738) q[0];
sx q[0];
rz(-1.0882821) q[0];
sx q[0];
rz(-0.58084647) q[0];
rz(-pi) q[1];
rz(-2.4841294) q[2];
sx q[2];
rz(-1.4095479) q[2];
sx q[2];
rz(0.25875124) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8090933) q[1];
sx q[1];
rz(-2.1488071) q[1];
sx q[1];
rz(-2.1564469) q[1];
rz(-1.6886466) q[3];
sx q[3];
rz(-2.1292344) q[3];
sx q[3];
rz(-1.3056361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.27611) q[2];
sx q[2];
rz(-1.7296187) q[2];
sx q[2];
rz(-2.0152246) q[2];
rz(1.3364835) q[3];
sx q[3];
rz(-2.0650605) q[3];
sx q[3];
rz(-2.0520463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4514076) q[0];
sx q[0];
rz(-1.0253588) q[0];
sx q[0];
rz(3.0080646) q[0];
rz(2.1639157) q[1];
sx q[1];
rz(-0.97594273) q[1];
sx q[1];
rz(0.28688637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8650019) q[0];
sx q[0];
rz(-1.5230858) q[0];
sx q[0];
rz(-0.86930958) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8989766) q[2];
sx q[2];
rz(-1.5300473) q[2];
sx q[2];
rz(-2.2475257) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62090767) q[1];
sx q[1];
rz(-1.0093099) q[1];
sx q[1];
rz(-0.26446277) q[1];
x q[2];
rz(-0.37128504) q[3];
sx q[3];
rz(-2.5142736) q[3];
sx q[3];
rz(-0.016022779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4213244) q[2];
sx q[2];
rz(-1.4807533) q[2];
sx q[2];
rz(1.486091) q[2];
rz(-2.7739575) q[3];
sx q[3];
rz(-2.1143819) q[3];
sx q[3];
rz(1.0953085) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7169645) q[0];
sx q[0];
rz(-0.76911887) q[0];
sx q[0];
rz(2.6677483) q[0];
rz(2.4773856) q[1];
sx q[1];
rz(-1.217548) q[1];
sx q[1];
rz(-1.7582105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48179214) q[0];
sx q[0];
rz(-1.0873902) q[0];
sx q[0];
rz(0.76235911) q[0];
rz(1.9888617) q[2];
sx q[2];
rz(-1.5454614) q[2];
sx q[2];
rz(-1.0449787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72402945) q[1];
sx q[1];
rz(-1.8010407) q[1];
sx q[1];
rz(-1.3413221) q[1];
rz(-pi) q[2];
rz(2.2380736) q[3];
sx q[3];
rz(-0.28378962) q[3];
sx q[3];
rz(1.2354148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9845) q[2];
sx q[2];
rz(-1.4089156) q[2];
sx q[2];
rz(-2.816693) q[2];
rz(2.7609008) q[3];
sx q[3];
rz(-2.295953) q[3];
sx q[3];
rz(-2.0438173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0710881) q[0];
sx q[0];
rz(-0.66362137) q[0];
sx q[0];
rz(-0.93884236) q[0];
rz(-0.70010575) q[1];
sx q[1];
rz(-2.5036) q[1];
sx q[1];
rz(0.28900388) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28075179) q[0];
sx q[0];
rz(-1.7993225) q[0];
sx q[0];
rz(0.079527461) q[0];
x q[1];
rz(2.2174066) q[2];
sx q[2];
rz(-1.980154) q[2];
sx q[2];
rz(1.5034624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6470312) q[1];
sx q[1];
rz(-0.53003487) q[1];
sx q[1];
rz(1.3840335) q[1];
rz(-2.7606332) q[3];
sx q[3];
rz(-1.6376312) q[3];
sx q[3];
rz(1.9240954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2603904) q[2];
sx q[2];
rz(-1.5217047) q[2];
sx q[2];
rz(2.728906) q[2];
rz(-2.2212501) q[3];
sx q[3];
rz(-2.6350382) q[3];
sx q[3];
rz(1.2296366) q[3];
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
rz(1.6270139) q[0];
sx q[0];
rz(-0.98015061) q[0];
sx q[0];
rz(0.29943109) q[0];
rz(1.1224271) q[1];
sx q[1];
rz(-1.2010682) q[1];
sx q[1];
rz(2.1103512) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2107687) q[0];
sx q[0];
rz(-1.5815736) q[0];
sx q[0];
rz(1.9369387) q[0];
rz(-pi) q[1];
rz(-2.6793807) q[2];
sx q[2];
rz(-1.3920857) q[2];
sx q[2];
rz(-1.6978839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6935307) q[1];
sx q[1];
rz(-1.281257) q[1];
sx q[1];
rz(2.2953643) q[1];
rz(-pi) q[2];
rz(-2.5599285) q[3];
sx q[3];
rz(-2.1278893) q[3];
sx q[3];
rz(-0.34639369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1624182) q[2];
sx q[2];
rz(-1.7981497) q[2];
sx q[2];
rz(2.8988885) q[2];
rz(-1.8414712) q[3];
sx q[3];
rz(-2.0827115) q[3];
sx q[3];
rz(0.21568957) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8278787) q[0];
sx q[0];
rz(-2.9254881) q[0];
sx q[0];
rz(-1.4208273) q[0];
rz(-2.8315262) q[1];
sx q[1];
rz(-2.2646751) q[1];
sx q[1];
rz(1.5440595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2472947) q[0];
sx q[0];
rz(-1.159512) q[0];
sx q[0];
rz(-0.572834) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0528131) q[2];
sx q[2];
rz(-1.2301461) q[2];
sx q[2];
rz(1.7973877) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6621993) q[1];
sx q[1];
rz(-1.7879722) q[1];
sx q[1];
rz(-0.26611664) q[1];
rz(-pi) q[2];
rz(-0.15940729) q[3];
sx q[3];
rz(-1.4324354) q[3];
sx q[3];
rz(1.0444928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1342423) q[2];
sx q[2];
rz(-1.728936) q[2];
sx q[2];
rz(-1.3664112) q[2];
rz(-0.8775231) q[3];
sx q[3];
rz(-1.6759422) q[3];
sx q[3];
rz(0.88959488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919745) q[0];
sx q[0];
rz(-1.5852954) q[0];
sx q[0];
rz(-1.4854767) q[0];
rz(-0.86434518) q[1];
sx q[1];
rz(-2.1280011) q[1];
sx q[1];
rz(-0.43307532) q[1];
rz(2.7769185) q[2];
sx q[2];
rz(-0.52959792) q[2];
sx q[2];
rz(-2.9858225) q[2];
rz(-3.0861985) q[3];
sx q[3];
rz(-1.0923891) q[3];
sx q[3];
rz(-0.53434571) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
