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
rz(0.95712823) q[0];
sx q[0];
rz(4.794802) q[0];
sx q[0];
rz(10.573536) q[0];
rz(-2.5454638) q[1];
sx q[1];
rz(3.547668) q[1];
sx q[1];
rz(14.021056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3341951) q[0];
sx q[0];
rz(-2.3435623) q[0];
sx q[0];
rz(1.0502771) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76081907) q[2];
sx q[2];
rz(-1.4913412) q[2];
sx q[2];
rz(0.75817273) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64950424) q[1];
sx q[1];
rz(-1.7028812) q[1];
sx q[1];
rz(2.3852939) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1249073) q[3];
sx q[3];
rz(-0.66185274) q[3];
sx q[3];
rz(0.13368363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9601606) q[2];
sx q[2];
rz(-2.1655607) q[2];
sx q[2];
rz(1.206548) q[2];
rz(-1.1613965) q[3];
sx q[3];
rz(-2.5738218) q[3];
sx q[3];
rz(1.5706221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1777765) q[0];
sx q[0];
rz(-1.125536) q[0];
sx q[0];
rz(-0.97208446) q[0];
rz(-2.8731335) q[1];
sx q[1];
rz(-1.700289) q[1];
sx q[1];
rz(-0.45480248) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1546611) q[0];
sx q[0];
rz(-1.4442181) q[0];
sx q[0];
rz(1.8693699) q[0];
rz(-pi) q[1];
rz(-1.182453) q[2];
sx q[2];
rz(-1.7162616) q[2];
sx q[2];
rz(-2.5487473) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0510682) q[1];
sx q[1];
rz(-1.5195492) q[1];
sx q[1];
rz(2.6766876) q[1];
rz(-pi) q[2];
rz(0.76105705) q[3];
sx q[3];
rz(-2.421954) q[3];
sx q[3];
rz(2.3005405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1099757) q[2];
sx q[2];
rz(-0.72828186) q[2];
sx q[2];
rz(-2.1419683) q[2];
rz(3.055618) q[3];
sx q[3];
rz(-1.5857668) q[3];
sx q[3];
rz(-0.011367817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7407532) q[0];
sx q[0];
rz(-1.6735621) q[0];
sx q[0];
rz(-2.9174347) q[0];
rz(-1.594918) q[1];
sx q[1];
rz(-1.5983351) q[1];
sx q[1];
rz(1.3067783) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76837158) q[0];
sx q[0];
rz(-2.5766815) q[0];
sx q[0];
rz(2.9119836) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3575338) q[2];
sx q[2];
rz(-2.2360115) q[2];
sx q[2];
rz(2.5248418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94470384) q[1];
sx q[1];
rz(-2.534229) q[1];
sx q[1];
rz(0.783226) q[1];
rz(-pi) q[2];
rz(1.3177835) q[3];
sx q[3];
rz(-0.91716498) q[3];
sx q[3];
rz(-2.831649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33009067) q[2];
sx q[2];
rz(-1.6197438) q[2];
sx q[2];
rz(-1.0203934) q[2];
rz(1.8118106) q[3];
sx q[3];
rz(-2.0164169) q[3];
sx q[3];
rz(-1.5248732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40600768) q[0];
sx q[0];
rz(-1.0574874) q[0];
sx q[0];
rz(1.965858) q[0];
rz(0.27578393) q[1];
sx q[1];
rz(-2.3049054) q[1];
sx q[1];
rz(-0.63794678) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36786554) q[0];
sx q[0];
rz(-1.5380895) q[0];
sx q[0];
rz(1.5601394) q[0];
x q[1];
rz(-2.5204896) q[2];
sx q[2];
rz(-1.7503465) q[2];
sx q[2];
rz(0.68815069) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1012816) q[1];
sx q[1];
rz(-0.71523803) q[1];
sx q[1];
rz(-2.6331227) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6652936) q[3];
sx q[3];
rz(-1.2523942) q[3];
sx q[3];
rz(-1.4868143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63567579) q[2];
sx q[2];
rz(-1.5779147) q[2];
sx q[2];
rz(-1.6944616) q[2];
rz(-2.9315089) q[3];
sx q[3];
rz(-0.45322067) q[3];
sx q[3];
rz(0.92857462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90784812) q[0];
sx q[0];
rz(-2.9892428) q[0];
sx q[0];
rz(2.0727169) q[0];
rz(-1.2999889) q[1];
sx q[1];
rz(-1.4062107) q[1];
sx q[1];
rz(-1.2058421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17888363) q[0];
sx q[0];
rz(-0.41085748) q[0];
sx q[0];
rz(-1.7715447) q[0];
x q[1];
rz(-0.91798012) q[2];
sx q[2];
rz(-1.5632354) q[2];
sx q[2];
rz(-3.0631236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0413549) q[1];
sx q[1];
rz(-2.2861028) q[1];
sx q[1];
rz(-0.13731401) q[1];
rz(-0.48128328) q[3];
sx q[3];
rz(-1.7168772) q[3];
sx q[3];
rz(-1.8446326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.892889) q[2];
sx q[2];
rz(-0.37962571) q[2];
sx q[2];
rz(0.068566337) q[2];
rz(2.2019703) q[3];
sx q[3];
rz(-1.4087804) q[3];
sx q[3];
rz(-1.2287963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.762961) q[0];
sx q[0];
rz(-1.7447504) q[0];
sx q[0];
rz(2.6201541) q[0];
rz(1.6261082) q[1];
sx q[1];
rz(-1.3896959) q[1];
sx q[1];
rz(-0.62561402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30197016) q[0];
sx q[0];
rz(-1.6018943) q[0];
sx q[0];
rz(2.1108225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4418774) q[2];
sx q[2];
rz(-1.9471418) q[2];
sx q[2];
rz(2.3835045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.292899) q[1];
sx q[1];
rz(-2.0239415) q[1];
sx q[1];
rz(-2.921656) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13497495) q[3];
sx q[3];
rz(-0.42536456) q[3];
sx q[3];
rz(0.36799001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0079415) q[2];
sx q[2];
rz(-0.59591728) q[2];
sx q[2];
rz(-3.0360743) q[2];
rz(-0.96406913) q[3];
sx q[3];
rz(-1.9897507) q[3];
sx q[3];
rz(-0.9849557) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0725726) q[0];
sx q[0];
rz(-0.21518406) q[0];
sx q[0];
rz(-1.4129114) q[0];
rz(-0.37881306) q[1];
sx q[1];
rz(-1.2544931) q[1];
sx q[1];
rz(-0.55317318) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2755097) q[0];
sx q[0];
rz(-1.1780773) q[0];
sx q[0];
rz(2.2662972) q[0];
x q[1];
rz(-1.6611536) q[2];
sx q[2];
rz(-1.7117662) q[2];
sx q[2];
rz(0.80328926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70771101) q[1];
sx q[1];
rz(-1.5305145) q[1];
sx q[1];
rz(0.28357419) q[1];
rz(-pi) q[2];
rz(-2.2284248) q[3];
sx q[3];
rz(-2.8165157) q[3];
sx q[3];
rz(2.7108148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46099123) q[2];
sx q[2];
rz(-1.0356156) q[2];
sx q[2];
rz(-2.7294532) q[2];
rz(0.66017094) q[3];
sx q[3];
rz(-2.6001402) q[3];
sx q[3];
rz(-0.49303833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93765813) q[0];
sx q[0];
rz(-2.0964607) q[0];
sx q[0];
rz(-0.20189051) q[0];
rz(2.0837325) q[1];
sx q[1];
rz(-2.7767534) q[1];
sx q[1];
rz(2.8813664) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79151151) q[0];
sx q[0];
rz(-0.66257325) q[0];
sx q[0];
rz(2.4988079) q[0];
rz(-pi) q[1];
rz(-2.7041928) q[2];
sx q[2];
rz(-2.045131) q[2];
sx q[2];
rz(-2.2651644) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34781814) q[1];
sx q[1];
rz(-1.8500237) q[1];
sx q[1];
rz(0.3600959) q[1];
x q[2];
rz(2.4811141) q[3];
sx q[3];
rz(-1.704543) q[3];
sx q[3];
rz(1.4123358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95320931) q[2];
sx q[2];
rz(-1.7069495) q[2];
sx q[2];
rz(2.5592213) q[2];
rz(0.44238704) q[3];
sx q[3];
rz(-1.9059537) q[3];
sx q[3];
rz(-0.83627397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22170947) q[0];
sx q[0];
rz(-1.6136805) q[0];
sx q[0];
rz(-1.4269933) q[0];
rz(-1.3555948) q[1];
sx q[1];
rz(-1.4878863) q[1];
sx q[1];
rz(-1.4367163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6528931) q[0];
sx q[0];
rz(-2.821229) q[0];
sx q[0];
rz(1.6344749) q[0];
x q[1];
rz(-1.8709917) q[2];
sx q[2];
rz(-1.1087024) q[2];
sx q[2];
rz(1.0754881) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.88594) q[1];
sx q[1];
rz(-0.24160928) q[1];
sx q[1];
rz(-1.6305981) q[1];
rz(-pi) q[2];
rz(1.0441699) q[3];
sx q[3];
rz(-1.7071402) q[3];
sx q[3];
rz(2.6201893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9068678) q[2];
sx q[2];
rz(-1.8610672) q[2];
sx q[2];
rz(-1.5391763) q[2];
rz(0.60798821) q[3];
sx q[3];
rz(-1.7323078) q[3];
sx q[3];
rz(2.4517945) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36668396) q[0];
sx q[0];
rz(-0.66119778) q[0];
sx q[0];
rz(-2.3181424) q[0];
rz(2.2460294) q[1];
sx q[1];
rz(-1.7915553) q[1];
sx q[1];
rz(0.38111883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136884) q[0];
sx q[0];
rz(-1.1263337) q[0];
sx q[0];
rz(-0.60143394) q[0];
rz(1.4032768) q[2];
sx q[2];
rz(-1.1741271) q[2];
sx q[2];
rz(-1.3127017) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5504073) q[1];
sx q[1];
rz(-2.1169615) q[1];
sx q[1];
rz(0.3355432) q[1];
x q[2];
rz(-2.1249843) q[3];
sx q[3];
rz(-1.018152) q[3];
sx q[3];
rz(1.2070398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8054008) q[2];
sx q[2];
rz(-2.7298584) q[2];
sx q[2];
rz(0.98708785) q[2];
rz(1.5385212) q[3];
sx q[3];
rz(-2.5542732) q[3];
sx q[3];
rz(0.2400329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9170452) q[0];
sx q[0];
rz(-1.2397091) q[0];
sx q[0];
rz(1.5201257) q[0];
rz(2.4317901) q[1];
sx q[1];
rz(-0.55955049) q[1];
sx q[1];
rz(1.4581663) q[1];
rz(-0.29322704) q[2];
sx q[2];
rz(-0.37920375) q[2];
sx q[2];
rz(-1.3398021) q[2];
rz(-0.74123989) q[3];
sx q[3];
rz(-2.0251353) q[3];
sx q[3];
rz(0.97788772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
