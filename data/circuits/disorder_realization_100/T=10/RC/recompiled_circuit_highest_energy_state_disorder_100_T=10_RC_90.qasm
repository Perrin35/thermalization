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
rz(1.1178782) q[0];
sx q[0];
rz(-2.0562545) q[0];
sx q[0];
rz(0.35882741) q[0];
rz(1.2598414) q[1];
sx q[1];
rz(-1.0955732) q[1];
sx q[1];
rz(-1.0395368) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9175076) q[0];
sx q[0];
rz(-1.9244902) q[0];
sx q[0];
rz(0.68266305) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1455457) q[2];
sx q[2];
rz(-3.0910935) q[2];
sx q[2];
rz(2.8157774) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1092069) q[1];
sx q[1];
rz(-2.6646174) q[1];
sx q[1];
rz(2.2677648) q[1];
rz(0.12173481) q[3];
sx q[3];
rz(-2.4678708) q[3];
sx q[3];
rz(0.47704968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.871668) q[2];
sx q[2];
rz(-2.2871127) q[2];
sx q[2];
rz(2.975115) q[2];
rz(-3.0050333) q[3];
sx q[3];
rz(-0.86524335) q[3];
sx q[3];
rz(1.4599919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305552) q[0];
sx q[0];
rz(-2.8189973) q[0];
sx q[0];
rz(0.36001298) q[0];
rz(-2.391138) q[1];
sx q[1];
rz(-0.89025918) q[1];
sx q[1];
rz(0.043936122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.982807) q[0];
sx q[0];
rz(-1.4318529) q[0];
sx q[0];
rz(-1.7278633) q[0];
rz(-1.5055799) q[2];
sx q[2];
rz(-0.64526099) q[2];
sx q[2];
rz(-1.1087904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39831623) q[1];
sx q[1];
rz(-0.77645248) q[1];
sx q[1];
rz(2.8321596) q[1];
rz(-0.93505164) q[3];
sx q[3];
rz(-2.8160444) q[3];
sx q[3];
rz(2.9071484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7875849) q[2];
sx q[2];
rz(-2.6968991) q[2];
sx q[2];
rz(2.1237874) q[2];
rz(2.9338037) q[3];
sx q[3];
rz(-0.60616797) q[3];
sx q[3];
rz(0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4491693) q[0];
sx q[0];
rz(-0.29266161) q[0];
sx q[0];
rz(-1.0798651) q[0];
rz(-1.0049459) q[1];
sx q[1];
rz(-0.69098538) q[1];
sx q[1];
rz(-1.8963922) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8522135) q[0];
sx q[0];
rz(-0.32158467) q[0];
sx q[0];
rz(-2.8711161) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7919134) q[2];
sx q[2];
rz(-0.44860835) q[2];
sx q[2];
rz(-0.3145379) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9745404) q[1];
sx q[1];
rz(-0.61950942) q[1];
sx q[1];
rz(-2.0761247) q[1];
rz(0.16896448) q[3];
sx q[3];
rz(-0.95674435) q[3];
sx q[3];
rz(1.3425296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9006742) q[2];
sx q[2];
rz(-1.2718879) q[2];
sx q[2];
rz(0.34350485) q[2];
rz(-0.014723012) q[3];
sx q[3];
rz(-0.69535178) q[3];
sx q[3];
rz(-1.2408313) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13046509) q[0];
sx q[0];
rz(-1.0437597) q[0];
sx q[0];
rz(0.68487942) q[0];
rz(-0.20826805) q[1];
sx q[1];
rz(-1.3083649) q[1];
sx q[1];
rz(0.79316717) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37453921) q[0];
sx q[0];
rz(-1.4724131) q[0];
sx q[0];
rz(-1.0258425) q[0];
rz(-pi) q[1];
rz(-0.17265138) q[2];
sx q[2];
rz(-1.7829636) q[2];
sx q[2];
rz(1.6546951) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82487166) q[1];
sx q[1];
rz(-2.5301456) q[1];
sx q[1];
rz(1.5427521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21677591) q[3];
sx q[3];
rz(-1.2086676) q[3];
sx q[3];
rz(-1.0877874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76306474) q[2];
sx q[2];
rz(-1.6918162) q[2];
sx q[2];
rz(0.48279631) q[2];
rz(1.5040846) q[3];
sx q[3];
rz(-0.62862325) q[3];
sx q[3];
rz(2.8211236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.095161) q[0];
sx q[0];
rz(-2.910039) q[0];
sx q[0];
rz(3.0401373) q[0];
rz(-0.21508148) q[1];
sx q[1];
rz(-1.9317893) q[1];
sx q[1];
rz(1.2717517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5520855) q[0];
sx q[0];
rz(-2.5906123) q[0];
sx q[0];
rz(2.1349585) q[0];
x q[1];
rz(-2.4319806) q[2];
sx q[2];
rz(-2.4104048) q[2];
sx q[2];
rz(-3.0572011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3569361) q[1];
sx q[1];
rz(-1.1788834) q[1];
sx q[1];
rz(-2.980113) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1121096) q[3];
sx q[3];
rz(-2.7446163) q[3];
sx q[3];
rz(0.1137867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35470206) q[2];
sx q[2];
rz(-2.6933647) q[2];
sx q[2];
rz(2.9917713) q[2];
rz(-0.24770728) q[3];
sx q[3];
rz(-0.37023538) q[3];
sx q[3];
rz(-0.80055922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54514068) q[0];
sx q[0];
rz(-0.62262028) q[0];
sx q[0];
rz(-1.8555634) q[0];
rz(2.3226358) q[1];
sx q[1];
rz(-0.30636925) q[1];
sx q[1];
rz(2.9584296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20907065) q[0];
sx q[0];
rz(-2.5446807) q[0];
sx q[0];
rz(0.36557622) q[0];
rz(-pi) q[1];
rz(-0.18503616) q[2];
sx q[2];
rz(-1.6335979) q[2];
sx q[2];
rz(2.7887995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.10658856) q[1];
sx q[1];
rz(-0.87751203) q[1];
sx q[1];
rz(2.6412852) q[1];
x q[2];
rz(-1.2533893) q[3];
sx q[3];
rz(-2.0485694) q[3];
sx q[3];
rz(2.1368515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5554819) q[2];
sx q[2];
rz(-1.7038944) q[2];
sx q[2];
rz(2.6044031) q[2];
rz(-1.7389343) q[3];
sx q[3];
rz(-1.8776929) q[3];
sx q[3];
rz(-2.5130443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1527767) q[0];
sx q[0];
rz(-2.7429136) q[0];
sx q[0];
rz(-0.11845778) q[0];
rz(-1.1064103) q[1];
sx q[1];
rz(-1.9170008) q[1];
sx q[1];
rz(-2.4647663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7962435) q[0];
sx q[0];
rz(-1.6259832) q[0];
sx q[0];
rz(-0.075166239) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7055461) q[2];
sx q[2];
rz(-2.2206077) q[2];
sx q[2];
rz(1.4092434) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.21444328) q[1];
sx q[1];
rz(-0.24733686) q[1];
sx q[1];
rz(-0.74199583) q[1];
rz(-1.5049797) q[3];
sx q[3];
rz(-0.01465791) q[3];
sx q[3];
rz(-2.3848118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38830882) q[2];
sx q[2];
rz(-0.94677418) q[2];
sx q[2];
rz(-0.050966144) q[2];
rz(-0.1554337) q[3];
sx q[3];
rz(-1.6499465) q[3];
sx q[3];
rz(-2.1211933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.257457) q[0];
sx q[0];
rz(-1.9866885) q[0];
sx q[0];
rz(2.1575523) q[0];
rz(0.46977305) q[1];
sx q[1];
rz(-1.677294) q[1];
sx q[1];
rz(0.22107548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8399682) q[0];
sx q[0];
rz(-1.5519807) q[0];
sx q[0];
rz(3.0473659) q[0];
rz(-pi) q[1];
rz(1.1936155) q[2];
sx q[2];
rz(-2.1862891) q[2];
sx q[2];
rz(-0.35843682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4692336) q[1];
sx q[1];
rz(-1.2738528) q[1];
sx q[1];
rz(0.23632144) q[1];
x q[2];
rz(0.89721591) q[3];
sx q[3];
rz(-0.69475896) q[3];
sx q[3];
rz(-1.0396797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2880963) q[2];
sx q[2];
rz(-2.4595538) q[2];
sx q[2];
rz(-1.811553) q[2];
rz(-0.62411493) q[3];
sx q[3];
rz(-0.95804405) q[3];
sx q[3];
rz(-0.38930711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5617111) q[0];
sx q[0];
rz(-1.6498673) q[0];
sx q[0];
rz(1.0505744) q[0];
rz(-0.9802649) q[1];
sx q[1];
rz(-1.0824208) q[1];
sx q[1];
rz(1.6260653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3526392) q[0];
sx q[0];
rz(-1.727956) q[0];
sx q[0];
rz(-2.609014) q[0];
rz(-0.43584906) q[2];
sx q[2];
rz(-0.51700355) q[2];
sx q[2];
rz(-0.85384254) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61170288) q[1];
sx q[1];
rz(-1.262959) q[1];
sx q[1];
rz(0.3605538) q[1];
x q[2];
rz(-0.6330746) q[3];
sx q[3];
rz(-0.71022034) q[3];
sx q[3];
rz(2.4274277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3773697) q[2];
sx q[2];
rz(-1.2219967) q[2];
sx q[2];
rz(2.4019305) q[2];
rz(2.9861279) q[3];
sx q[3];
rz(-0.37853095) q[3];
sx q[3];
rz(2.9445061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71843201) q[0];
sx q[0];
rz(-0.53373706) q[0];
sx q[0];
rz(-2.0045795) q[0];
rz(2.5088189) q[1];
sx q[1];
rz(-2.4686765) q[1];
sx q[1];
rz(2.4957962) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9912796) q[0];
sx q[0];
rz(-1.2431354) q[0];
sx q[0];
rz(-1.3929266) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8027202) q[2];
sx q[2];
rz(-0.80451369) q[2];
sx q[2];
rz(2.5520012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.551522) q[1];
sx q[1];
rz(-2.3209718) q[1];
sx q[1];
rz(0.094610081) q[1];
rz(-0.31787536) q[3];
sx q[3];
rz(-1.0680305) q[3];
sx q[3];
rz(1.1989145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48468963) q[2];
sx q[2];
rz(-0.36269665) q[2];
sx q[2];
rz(0.038385782) q[2];
rz(-3.0302327) q[3];
sx q[3];
rz(-2.7146118) q[3];
sx q[3];
rz(2.4666983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4511694) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.3896598) q[0];
rz(-1.8545064) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(0.95697516) q[2];
sx q[2];
rz(-1.0697094) q[2];
sx q[2];
rz(-0.30164837) q[2];
rz(-0.18517687) q[3];
sx q[3];
rz(-0.76793302) q[3];
sx q[3];
rz(-1.9938706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
