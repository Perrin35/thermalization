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
rz(-2.1844644) q[0];
sx q[0];
rz(-1.6532093) q[0];
sx q[0];
rz(-1.1487577) q[0];
rz(-2.5454638) q[1];
sx q[1];
rz(-2.7355173) q[1];
sx q[1];
rz(-1.6869071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4948954) q[0];
sx q[0];
rz(-2.2410164) q[0];
sx q[0];
rz(2.669932) q[0];
rz(-2.3807736) q[2];
sx q[2];
rz(-1.6502514) q[2];
sx q[2];
rz(-0.75817273) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0599938) q[1];
sx q[1];
rz(-0.76548701) q[1];
sx q[1];
rz(0.19123921) q[1];
rz(2.4835281) q[3];
sx q[3];
rz(-1.4941553) q[3];
sx q[3];
rz(1.5358401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9601606) q[2];
sx q[2];
rz(-0.97603193) q[2];
sx q[2];
rz(1.206548) q[2];
rz(1.1613965) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(-1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96381617) q[0];
sx q[0];
rz(-2.0160567) q[0];
sx q[0];
rz(-0.97208446) q[0];
rz(0.26845911) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(2.6867902) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3773012) q[0];
sx q[0];
rz(-1.2746841) q[0];
sx q[0];
rz(-0.13237093) q[0];
rz(-1.9591397) q[2];
sx q[2];
rz(-1.7162616) q[2];
sx q[2];
rz(-0.59284537) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5059591) q[1];
sx q[1];
rz(-2.035043) q[1];
sx q[1];
rz(1.6281158) q[1];
x q[2];
rz(2.11449) q[3];
sx q[3];
rz(-2.0683401) q[3];
sx q[3];
rz(-0.061391679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1099757) q[2];
sx q[2];
rz(-2.4133108) q[2];
sx q[2];
rz(-2.1419683) q[2];
rz(-3.055618) q[3];
sx q[3];
rz(-1.5857668) q[3];
sx q[3];
rz(-3.1302248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40083945) q[0];
sx q[0];
rz(-1.6735621) q[0];
sx q[0];
rz(2.9174347) q[0];
rz(1.5466746) q[1];
sx q[1];
rz(-1.5983351) q[1];
sx q[1];
rz(-1.8348144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1442566) q[0];
sx q[0];
rz(-1.4486509) q[0];
sx q[0];
rz(-0.55295487) q[0];
rz(-0.73422696) q[2];
sx q[2];
rz(-2.1598744) q[2];
sx q[2];
rz(0.40225077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.063547246) q[1];
sx q[1];
rz(-1.1544466) q[1];
sx q[1];
rz(1.1148808) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8238092) q[3];
sx q[3];
rz(-2.2244277) q[3];
sx q[3];
rz(-0.30994367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33009067) q[2];
sx q[2];
rz(-1.5218488) q[2];
sx q[2];
rz(-2.1211993) q[2];
rz(-1.8118106) q[3];
sx q[3];
rz(-1.1251757) q[3];
sx q[3];
rz(1.6167195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.735585) q[0];
sx q[0];
rz(-2.0841053) q[0];
sx q[0];
rz(1.1757346) q[0];
rz(-2.8658087) q[1];
sx q[1];
rz(-0.83668721) q[1];
sx q[1];
rz(0.63794678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9383134) q[0];
sx q[0];
rz(-1.5601451) q[0];
sx q[0];
rz(0.032708688) q[0];
rz(-pi) q[1];
rz(0.3023382) q[2];
sx q[2];
rz(-0.64322847) q[2];
sx q[2];
rz(-2.0143721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13220683) q[1];
sx q[1];
rz(-1.8957545) q[1];
sx q[1];
rz(-2.4924949) q[1];
rz(-pi) q[2];
rz(0.27880554) q[3];
sx q[3];
rz(-0.33167095) q[3];
sx q[3];
rz(1.7808033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5059169) q[2];
sx q[2];
rz(-1.5636779) q[2];
sx q[2];
rz(1.4471311) q[2];
rz(-0.21008374) q[3];
sx q[3];
rz(-2.688372) q[3];
sx q[3];
rz(0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2337445) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(1.0688758) q[0];
rz(-1.2999889) q[1];
sx q[1];
rz(-1.735382) q[1];
sx q[1];
rz(-1.9357505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17888363) q[0];
sx q[0];
rz(-0.41085748) q[0];
sx q[0];
rz(1.370048) q[0];
x q[1];
rz(2.2236125) q[2];
sx q[2];
rz(-1.5632354) q[2];
sx q[2];
rz(0.078469097) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5609392) q[1];
sx q[1];
rz(-1.4672797) q[1];
sx q[1];
rz(2.2907881) q[1];
x q[2];
rz(0.48128328) q[3];
sx q[3];
rz(-1.4247155) q[3];
sx q[3];
rz(-1.8446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24870366) q[2];
sx q[2];
rz(-0.37962571) q[2];
sx q[2];
rz(-0.068566337) q[2];
rz(-0.93962234) q[3];
sx q[3];
rz(-1.4087804) q[3];
sx q[3];
rz(-1.2287963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762961) q[0];
sx q[0];
rz(-1.3968422) q[0];
sx q[0];
rz(-2.6201541) q[0];
rz(-1.5154845) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(2.5159786) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30197016) q[0];
sx q[0];
rz(-1.5396984) q[0];
sx q[0];
rz(2.1108225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7623903) q[2];
sx q[2];
rz(-1.450945) q[2];
sx q[2];
rz(0.7650991) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3208116) q[1];
sx q[1];
rz(-0.50034517) q[1];
sx q[1];
rz(1.9920177) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7196521) q[3];
sx q[3];
rz(-1.6263537) q[3];
sx q[3];
rz(2.0618604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1336512) q[2];
sx q[2];
rz(-2.5456754) q[2];
sx q[2];
rz(-0.1055183) q[2];
rz(-0.96406913) q[3];
sx q[3];
rz(-1.151842) q[3];
sx q[3];
rz(-2.156637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069020011) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(1.4129114) q[0];
rz(2.7627796) q[1];
sx q[1];
rz(-1.8870995) q[1];
sx q[1];
rz(0.55317318) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39554292) q[0];
sx q[0];
rz(-2.2043214) q[0];
sx q[0];
rz(0.49479158) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4804391) q[2];
sx q[2];
rz(-1.7117662) q[2];
sx q[2];
rz(0.80328926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70771101) q[1];
sx q[1];
rz(-1.6110782) q[1];
sx q[1];
rz(-2.8580185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8314701) q[3];
sx q[3];
rz(-1.7672774) q[3];
sx q[3];
rz(2.6333661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6806014) q[2];
sx q[2];
rz(-2.1059771) q[2];
sx q[2];
rz(2.7294532) q[2];
rz(-0.66017094) q[3];
sx q[3];
rz(-0.5414525) q[3];
sx q[3];
rz(2.6485543) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93765813) q[0];
sx q[0];
rz(-1.0451319) q[0];
sx q[0];
rz(0.20189051) q[0];
rz(2.0837325) q[1];
sx q[1];
rz(-0.36483929) q[1];
sx q[1];
rz(-2.8813664) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3500811) q[0];
sx q[0];
rz(-2.4790194) q[0];
sx q[0];
rz(-2.4988079) q[0];
rz(-2.7041928) q[2];
sx q[2];
rz(-1.0964616) q[2];
sx q[2];
rz(-0.87642821) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.550506) q[1];
sx q[1];
rz(-0.45190736) q[1];
sx q[1];
rz(-0.68303789) q[1];
x q[2];
rz(0.66047858) q[3];
sx q[3];
rz(-1.704543) q[3];
sx q[3];
rz(-1.4123358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1883833) q[2];
sx q[2];
rz(-1.4346432) q[2];
sx q[2];
rz(0.58237135) q[2];
rz(2.6992056) q[3];
sx q[3];
rz(-1.235639) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.22170947) q[0];
sx q[0];
rz(-1.5279122) q[0];
sx q[0];
rz(-1.4269933) q[0];
rz(-1.7859979) q[1];
sx q[1];
rz(-1.4878863) q[1];
sx q[1];
rz(-1.7048763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6528931) q[0];
sx q[0];
rz(-2.821229) q[0];
sx q[0];
rz(-1.6344749) q[0];
rz(-2.6609909) q[2];
sx q[2];
rz(-1.8386823) q[2];
sx q[2];
rz(0.3581697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.88594) q[1];
sx q[1];
rz(-0.24160928) q[1];
sx q[1];
rz(-1.5109946) q[1];
rz(-2.9842078) q[3];
sx q[3];
rz(-1.0495571) q[3];
sx q[3];
rz(-0.97053274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9068678) q[2];
sx q[2];
rz(-1.2805254) q[2];
sx q[2];
rz(1.5391763) q[2];
rz(-0.60798821) q[3];
sx q[3];
rz(-1.7323078) q[3];
sx q[3];
rz(0.68979818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36668396) q[0];
sx q[0];
rz(-2.4803949) q[0];
sx q[0];
rz(-0.82345024) q[0];
rz(0.89556328) q[1];
sx q[1];
rz(-1.7915553) q[1];
sx q[1];
rz(-0.38111883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51682794) q[0];
sx q[0];
rz(-2.4104558) q[0];
sx q[0];
rz(-2.4419489) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7383159) q[2];
sx q[2];
rz(-1.1741271) q[2];
sx q[2];
rz(1.3127017) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1824519) q[1];
sx q[1];
rz(-2.509626) q[1];
sx q[1];
rz(-1.0743293) q[1];
x q[2];
rz(-0.62750979) q[3];
sx q[3];
rz(-1.1063647) q[3];
sx q[3];
rz(-3.0919558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33619189) q[2];
sx q[2];
rz(-0.41173428) q[2];
sx q[2];
rz(-0.98708785) q[2];
rz(1.6030715) q[3];
sx q[3];
rz(-0.58731949) q[3];
sx q[3];
rz(-2.9015598) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9170452) q[0];
sx q[0];
rz(-1.9018835) q[0];
sx q[0];
rz(-1.6214669) q[0];
rz(0.70980258) q[1];
sx q[1];
rz(-2.5820422) q[1];
sx q[1];
rz(-1.6834264) q[1];
rz(-2.8483656) q[2];
sx q[2];
rz(-2.7623889) q[2];
sx q[2];
rz(1.8017906) q[2];
rz(2.5153574) q[3];
sx q[3];
rz(-2.2953485) q[3];
sx q[3];
rz(2.9957537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
