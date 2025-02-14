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
rz(-2.0237145) q[0];
sx q[0];
rz(-1.0853381) q[0];
sx q[0];
rz(2.7827652) q[0];
rz(1.2598414) q[1];
sx q[1];
rz(2.0460195) q[1];
sx q[1];
rz(10.464315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3921689) q[0];
sx q[0];
rz(-0.75558973) q[0];
sx q[0];
rz(2.6120793) q[0];
rz(-pi) q[1];
rz(3.1141236) q[2];
sx q[2];
rz(-1.6131764) q[2];
sx q[2];
rz(-2.8920763) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4182836) q[1];
sx q[1];
rz(-1.9305349) q[1];
sx q[1];
rz(2.8213001) q[1];
rz(3.0198578) q[3];
sx q[3];
rz(-2.4678708) q[3];
sx q[3];
rz(-0.47704968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.871668) q[2];
sx q[2];
rz(-0.85447997) q[2];
sx q[2];
rz(2.975115) q[2];
rz(3.0050333) q[3];
sx q[3];
rz(-0.86524335) q[3];
sx q[3];
rz(-1.4599919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305552) q[0];
sx q[0];
rz(-0.32259536) q[0];
sx q[0];
rz(2.7815797) q[0];
rz(0.75045466) q[1];
sx q[1];
rz(-2.2513335) q[1];
sx q[1];
rz(-0.043936122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.982807) q[0];
sx q[0];
rz(-1.7097397) q[0];
sx q[0];
rz(1.7278633) q[0];
rz(-pi) q[1];
rz(2.2150351) q[2];
sx q[2];
rz(-1.6100002) q[2];
sx q[2];
rz(2.7317177) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.81961957) q[1];
sx q[1];
rz(-2.3016046) q[1];
sx q[1];
rz(-1.8614443) q[1];
rz(-pi) q[2];
rz(0.93505164) q[3];
sx q[3];
rz(-0.32554829) q[3];
sx q[3];
rz(2.9071484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7875849) q[2];
sx q[2];
rz(-2.6968991) q[2];
sx q[2];
rz(-2.1237874) q[2];
rz(-0.20778896) q[3];
sx q[3];
rz(-2.5354247) q[3];
sx q[3];
rz(2.7868605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4491693) q[0];
sx q[0];
rz(-2.848931) q[0];
sx q[0];
rz(-1.0798651) q[0];
rz(2.1366468) q[1];
sx q[1];
rz(-0.69098538) q[1];
sx q[1];
rz(1.2452004) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2893791) q[0];
sx q[0];
rz(-0.32158467) q[0];
sx q[0];
rz(0.2704765) q[0];
rz(2.7169021) q[2];
sx q[2];
rz(-1.4216558) q[2];
sx q[2];
rz(-1.5678659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1221262) q[1];
sx q[1];
rz(-1.8557185) q[1];
sx q[1];
rz(1.0128922) q[1];
rz(-0.94997378) q[3];
sx q[3];
rz(-1.7086747) q[3];
sx q[3];
rz(2.8153489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9006742) q[2];
sx q[2];
rz(-1.2718879) q[2];
sx q[2];
rz(0.34350485) q[2];
rz(3.1268696) q[3];
sx q[3];
rz(-2.4462409) q[3];
sx q[3];
rz(1.2408313) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13046509) q[0];
sx q[0];
rz(-2.0978329) q[0];
sx q[0];
rz(-0.68487942) q[0];
rz(-2.9333246) q[1];
sx q[1];
rz(-1.8332278) q[1];
sx q[1];
rz(-2.3484255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.255729) q[0];
sx q[0];
rz(-2.1128214) q[0];
sx q[0];
rz(-0.11491187) q[0];
rz(-2.9689413) q[2];
sx q[2];
rz(-1.7829636) q[2];
sx q[2];
rz(-1.6546951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4186331) q[1];
sx q[1];
rz(-1.5868938) q[1];
sx q[1];
rz(2.1820585) q[1];
rz(-pi) q[2];
rz(-1.2007522) q[3];
sx q[3];
rz(-1.368282) q[3];
sx q[3];
rz(-2.5807192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76306474) q[2];
sx q[2];
rz(-1.6918162) q[2];
sx q[2];
rz(0.48279631) q[2];
rz(1.5040846) q[3];
sx q[3];
rz(-2.5129694) q[3];
sx q[3];
rz(-2.8211236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.095161) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(3.0401373) q[0];
rz(-0.21508148) q[1];
sx q[1];
rz(-1.9317893) q[1];
sx q[1];
rz(-1.869841) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5895071) q[0];
sx q[0];
rz(-2.5906123) q[0];
sx q[0];
rz(-1.0066342) q[0];
x q[1];
rz(-1.0418747) q[2];
sx q[2];
rz(-2.1019962) q[2];
sx q[2];
rz(-0.94106969) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3569361) q[1];
sx q[1];
rz(-1.9627092) q[1];
sx q[1];
rz(-2.980113) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21274626) q[3];
sx q[3];
rz(-1.233056) q[3];
sx q[3];
rz(2.4500873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7868906) q[2];
sx q[2];
rz(-2.6933647) q[2];
sx q[2];
rz(-0.1498214) q[2];
rz(-2.8938854) q[3];
sx q[3];
rz(-2.7713573) q[3];
sx q[3];
rz(2.3410334) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.596452) q[0];
sx q[0];
rz(-2.5189724) q[0];
sx q[0];
rz(-1.2860292) q[0];
rz(-0.81895685) q[1];
sx q[1];
rz(-0.30636925) q[1];
sx q[1];
rz(2.9584296) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.055119) q[0];
sx q[0];
rz(-1.3684784) q[0];
sx q[0];
rz(-2.5760465) q[0];
x q[1];
rz(-2.8122453) q[2];
sx q[2];
rz(-0.19528772) q[2];
sx q[2];
rz(1.5415217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3273436) q[1];
sx q[1];
rz(-2.3115989) q[1];
sx q[1];
rz(-1.0472058) q[1];
x q[2];
rz(2.5991393) q[3];
sx q[3];
rz(-2.5748655) q[3];
sx q[3];
rz(2.7572214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5861107) q[2];
sx q[2];
rz(-1.4376983) q[2];
sx q[2];
rz(0.5371896) q[2];
rz(-1.7389343) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(-0.62854832) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9888159) q[0];
sx q[0];
rz(-2.7429136) q[0];
sx q[0];
rz(0.11845778) q[0];
rz(-1.1064103) q[1];
sx q[1];
rz(-1.2245919) q[1];
sx q[1];
rz(2.4647663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9119916) q[0];
sx q[0];
rz(-1.4957447) q[0];
sx q[0];
rz(1.5154535) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9666128) q[2];
sx q[2];
rz(-2.4799441) q[2];
sx q[2];
rz(1.511919) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9271494) q[1];
sx q[1];
rz(-0.24733686) q[1];
sx q[1];
rz(-0.74199583) q[1];
rz(-pi) q[2];
rz(-1.6366129) q[3];
sx q[3];
rz(-3.1269347) q[3];
sx q[3];
rz(-2.3848118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38830882) q[2];
sx q[2];
rz(-0.94677418) q[2];
sx q[2];
rz(0.050966144) q[2];
rz(2.986159) q[3];
sx q[3];
rz(-1.4916462) q[3];
sx q[3];
rz(2.1211933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8841356) q[0];
sx q[0];
rz(-1.9866885) q[0];
sx q[0];
rz(-2.1575523) q[0];
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
rz(1.2673938) q[0];
sx q[0];
rz(-1.4765863) q[0];
sx q[0];
rz(1.5518968) q[0];
rz(-0.48018166) q[2];
sx q[2];
rz(-2.4327185) q[2];
sx q[2];
rz(0.95979662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.672359) q[1];
sx q[1];
rz(-1.2738528) q[1];
sx q[1];
rz(2.9052712) q[1];
rz(-pi) q[2];
rz(-2.662195) q[3];
sx q[3];
rz(-2.0948296) q[3];
sx q[3];
rz(-2.9063922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85349637) q[2];
sx q[2];
rz(-0.6820389) q[2];
sx q[2];
rz(-1.811553) q[2];
rz(0.62411493) q[3];
sx q[3];
rz(-0.95804405) q[3];
sx q[3];
rz(-2.7522855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5617111) q[0];
sx q[0];
rz(-1.4917253) q[0];
sx q[0];
rz(1.0505744) q[0];
rz(2.1613278) q[1];
sx q[1];
rz(-2.0591718) q[1];
sx q[1];
rz(1.5155274) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78895348) q[0];
sx q[0];
rz(-1.4136367) q[0];
sx q[0];
rz(0.53257863) q[0];
rz(-1.3352065) q[2];
sx q[2];
rz(-2.0353999) q[2];
sx q[2];
rz(-2.7795779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5055464) q[1];
sx q[1];
rz(-2.6719173) q[1];
sx q[1];
rz(2.4081025) q[1];
rz(-pi) q[2];
rz(-2.0414166) q[3];
sx q[3];
rz(-2.1242766) q[3];
sx q[3];
rz(0.054892232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3773697) q[2];
sx q[2];
rz(-1.919596) q[2];
sx q[2];
rz(-2.4019305) q[2];
rz(2.9861279) q[3];
sx q[3];
rz(-0.37853095) q[3];
sx q[3];
rz(2.9445061) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71843201) q[0];
sx q[0];
rz(-2.6078556) q[0];
sx q[0];
rz(-1.1370132) q[0];
rz(-0.63277376) q[1];
sx q[1];
rz(-0.67291617) q[1];
sx q[1];
rz(0.64579642) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64090133) q[0];
sx q[0];
rz(-2.7703022) q[0];
sx q[0];
rz(-2.6616606) q[0];
rz(-pi) q[1];
rz(-2.3617451) q[2];
sx q[2];
rz(-1.7371685) q[2];
sx q[2];
rz(-2.322724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.551522) q[1];
sx q[1];
rz(-2.3209718) q[1];
sx q[1];
rz(-3.0469826) q[1];
rz(-pi) q[2];
rz(-2.0876376) q[3];
sx q[3];
rz(-2.5541383) q[3];
sx q[3];
rz(1.7980391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48468963) q[2];
sx q[2];
rz(-2.778896) q[2];
sx q[2];
rz(3.1032069) q[2];
rz(-3.0302327) q[3];
sx q[3];
rz(-0.42698082) q[3];
sx q[3];
rz(-2.4666983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69042324) q[0];
sx q[0];
rz(-1.7043865) q[0];
sx q[0];
rz(-1.7519328) q[0];
rz(-1.2870862) q[1];
sx q[1];
rz(-2.1422287) q[1];
sx q[1];
rz(1.1258463) q[1];
rz(-0.95697516) q[2];
sx q[2];
rz(-2.0718832) q[2];
sx q[2];
rz(2.8399443) q[2];
rz(-1.3948364) q[3];
sx q[3];
rz(-0.8192438) q[3];
sx q[3];
rz(1.4024709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
