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
rz(2.1020558) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0693927) q[0];
sx q[0];
rz(-0.93749267) q[0];
sx q[0];
rz(2.0149489) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5284003) q[2];
sx q[2];
rz(-1.5982407) q[2];
sx q[2];
rz(-1.3224441) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1778587) q[1];
sx q[1];
rz(-1.8699282) q[1];
sx q[1];
rz(1.9480716) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0198578) q[3];
sx q[3];
rz(-0.67372185) q[3];
sx q[3];
rz(0.47704968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.871668) q[2];
sx q[2];
rz(-0.85447997) q[2];
sx q[2];
rz(-0.16647767) q[2];
rz(0.13655937) q[3];
sx q[3];
rz(-0.86524335) q[3];
sx q[3];
rz(1.4599919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305552) q[0];
sx q[0];
rz(-0.32259536) q[0];
sx q[0];
rz(-0.36001298) q[0];
rz(2.391138) q[1];
sx q[1];
rz(-0.89025918) q[1];
sx q[1];
rz(3.0976565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7076516) q[0];
sx q[0];
rz(-1.4152555) q[0];
sx q[0];
rz(0.14065245) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0925748) q[2];
sx q[2];
rz(-2.214458) q[2];
sx q[2];
rz(1.9512392) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94823882) q[1];
sx q[1];
rz(-1.3557503) q[1];
sx q[1];
rz(-0.7521473) q[1];
rz(-pi) q[2];
rz(-0.1978132) q[3];
sx q[3];
rz(-1.3105243) q[3];
sx q[3];
rz(-0.42727268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7875849) q[2];
sx q[2];
rz(-2.6968991) q[2];
sx q[2];
rz(1.0178052) q[2];
rz(-2.9338037) q[3];
sx q[3];
rz(-0.60616797) q[3];
sx q[3];
rz(-0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4491693) q[0];
sx q[0];
rz(-0.29266161) q[0];
sx q[0];
rz(-1.0798651) q[0];
rz(-2.1366468) q[1];
sx q[1];
rz(-0.69098538) q[1];
sx q[1];
rz(1.8963922) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5737138) q[0];
sx q[0];
rz(-1.880293) q[0];
sx q[0];
rz(-1.6595766) q[0];
x q[1];
rz(-2.7919134) q[2];
sx q[2];
rz(-2.6929843) q[2];
sx q[2];
rz(-0.3145379) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1221262) q[1];
sx q[1];
rz(-1.2858741) q[1];
sx q[1];
rz(-1.0128922) q[1];
rz(-2.1916189) q[3];
sx q[3];
rz(-1.7086747) q[3];
sx q[3];
rz(0.32624376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9006742) q[2];
sx q[2];
rz(-1.2718879) q[2];
sx q[2];
rz(2.7980878) q[2];
rz(-3.1268696) q[3];
sx q[3];
rz(-2.4462409) q[3];
sx q[3];
rz(1.9007614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111276) q[0];
sx q[0];
rz(-2.0978329) q[0];
sx q[0];
rz(-2.4567132) q[0];
rz(2.9333246) q[1];
sx q[1];
rz(-1.8332278) q[1];
sx q[1];
rz(2.3484255) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858636) q[0];
sx q[0];
rz(-1.0287713) q[0];
sx q[0];
rz(0.11491187) q[0];
x q[1];
rz(0.17265138) q[2];
sx q[2];
rz(-1.358629) q[2];
sx q[2];
rz(-1.4868975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4186331) q[1];
sx q[1];
rz(-1.5868938) q[1];
sx q[1];
rz(-2.1820585) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0871619) q[3];
sx q[3];
rz(-0.41958365) q[3];
sx q[3];
rz(-0.53158599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3785279) q[2];
sx q[2];
rz(-1.4497764) q[2];
sx q[2];
rz(-0.48279631) q[2];
rz(-1.637508) q[3];
sx q[3];
rz(-0.62862325) q[3];
sx q[3];
rz(-0.32046902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(3.095161) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(-0.10145536) q[0];
rz(0.21508148) q[1];
sx q[1];
rz(-1.2098034) q[1];
sx q[1];
rz(1.2717517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1908784) q[0];
sx q[0];
rz(-1.1125277) q[0];
sx q[0];
rz(-2.8241497) q[0];
x q[1];
rz(-0.70961205) q[2];
sx q[2];
rz(-2.4104048) q[2];
sx q[2];
rz(-0.084391525) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3569361) q[1];
sx q[1];
rz(-1.1788834) q[1];
sx q[1];
rz(-0.16147967) q[1];
rz(-2.1121096) q[3];
sx q[3];
rz(-2.7446163) q[3];
sx q[3];
rz(-0.1137867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7868906) q[2];
sx q[2];
rz(-2.6933647) q[2];
sx q[2];
rz(-2.9917713) q[2];
rz(-2.8938854) q[3];
sx q[3];
rz(-2.7713573) q[3];
sx q[3];
rz(-0.80055922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.596452) q[0];
sx q[0];
rz(-2.5189724) q[0];
sx q[0];
rz(1.8555634) q[0];
rz(-2.3226358) q[1];
sx q[1];
rz(-2.8352234) q[1];
sx q[1];
rz(-0.18316306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20907065) q[0];
sx q[0];
rz(-0.59691197) q[0];
sx q[0];
rz(0.36557622) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6346857) q[2];
sx q[2];
rz(-1.7554635) q[2];
sx q[2];
rz(-1.2062564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10658856) q[1];
sx q[1];
rz(-0.87751203) q[1];
sx q[1];
rz(-0.50030746) q[1];
rz(0.49900238) q[3];
sx q[3];
rz(-1.8516282) q[3];
sx q[3];
rz(2.4256191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5861107) q[2];
sx q[2];
rz(-1.4376983) q[2];
sx q[2];
rz(-0.5371896) q[2];
rz(1.7389343) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(-2.5130443) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1527767) q[0];
sx q[0];
rz(-2.7429136) q[0];
sx q[0];
rz(-3.0231349) q[0];
rz(-1.1064103) q[1];
sx q[1];
rz(-1.9170008) q[1];
sx q[1];
rz(0.67682636) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3453492) q[0];
sx q[0];
rz(-1.6259832) q[0];
sx q[0];
rz(-3.0664264) q[0];
rz(-pi) q[1];
rz(-1.7055461) q[2];
sx q[2];
rz(-2.2206077) q[2];
sx q[2];
rz(1.7323493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5986277) q[1];
sx q[1];
rz(-1.7522546) q[1];
sx q[1];
rz(-1.7398029) q[1];
rz(-pi) q[2];
rz(-0.00096410613) q[3];
sx q[3];
rz(-1.5561702) q[3];
sx q[3];
rz(-2.4506354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7532838) q[2];
sx q[2];
rz(-0.94677418) q[2];
sx q[2];
rz(-3.0906265) q[2];
rz(-0.1554337) q[3];
sx q[3];
rz(-1.4916462) q[3];
sx q[3];
rz(-1.0203993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8841356) q[0];
sx q[0];
rz(-1.1549042) q[0];
sx q[0];
rz(2.1575523) q[0];
rz(-2.6718196) q[1];
sx q[1];
rz(-1.677294) q[1];
sx q[1];
rz(0.22107548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.675908) q[0];
sx q[0];
rz(-3.0455112) q[0];
sx q[0];
rz(0.19739993) q[0];
rz(-pi) q[1];
rz(2.661411) q[2];
sx q[2];
rz(-0.70887414) q[2];
sx q[2];
rz(2.181796) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1578678) q[1];
sx q[1];
rz(-0.37731428) q[1];
sx q[1];
rz(-2.2239216) q[1];
x q[2];
rz(2.1481236) q[3];
sx q[3];
rz(-1.1599891) q[3];
sx q[3];
rz(-1.081117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2880963) q[2];
sx q[2];
rz(-0.6820389) q[2];
sx q[2];
rz(1.3300396) q[2];
rz(0.62411493) q[3];
sx q[3];
rz(-0.95804405) q[3];
sx q[3];
rz(-2.7522855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57988155) q[0];
sx q[0];
rz(-1.6498673) q[0];
sx q[0];
rz(1.0505744) q[0];
rz(-0.9802649) q[1];
sx q[1];
rz(-1.0824208) q[1];
sx q[1];
rz(1.6260653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4517364) q[0];
sx q[0];
rz(-1.045466) q[0];
sx q[0];
rz(1.7527053) q[0];
rz(-pi) q[1];
rz(2.7057436) q[2];
sx q[2];
rz(-0.51700355) q[2];
sx q[2];
rz(-0.85384254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84534553) q[1];
sx q[1];
rz(-1.2279086) q[1];
sx q[1];
rz(1.2432437) q[1];
rz(-0.6330746) q[3];
sx q[3];
rz(-0.71022034) q[3];
sx q[3];
rz(2.4274277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3773697) q[2];
sx q[2];
rz(-1.2219967) q[2];
sx q[2];
rz(-2.4019305) q[2];
rz(2.9861279) q[3];
sx q[3];
rz(-2.7630617) q[3];
sx q[3];
rz(-2.9445061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.71843201) q[0];
sx q[0];
rz(-0.53373706) q[0];
sx q[0];
rz(2.0045795) q[0];
rz(0.63277376) q[1];
sx q[1];
rz(-0.67291617) q[1];
sx q[1];
rz(-0.64579642) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9912796) q[0];
sx q[0];
rz(-1.2431354) q[0];
sx q[0];
rz(-1.748666) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23441961) q[2];
sx q[2];
rz(-2.3478797) q[2];
sx q[2];
rz(2.22375) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72836956) q[1];
sx q[1];
rz(-0.75496208) q[1];
sx q[1];
rz(-1.4697716) q[1];
rz(-pi) q[2];
rz(-2.0955576) q[3];
sx q[3];
rz(-1.8482131) q[3];
sx q[3];
rz(2.926947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48468963) q[2];
sx q[2];
rz(-0.36269665) q[2];
sx q[2];
rz(0.038385782) q[2];
rz(0.11135993) q[3];
sx q[3];
rz(-2.7146118) q[3];
sx q[3];
rz(-0.67489433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69042324) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.3896598) q[0];
rz(1.2870862) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(-2.5512681) q[2];
sx q[2];
rz(-1.0412024) q[2];
sx q[2];
rz(1.5955284) q[2];
rz(-2.3822734) q[3];
sx q[3];
rz(-1.6990468) q[3];
sx q[3];
rz(2.852462) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
