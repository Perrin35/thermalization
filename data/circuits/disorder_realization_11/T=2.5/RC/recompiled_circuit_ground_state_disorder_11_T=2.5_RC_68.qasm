OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1144855) q[0];
sx q[0];
rz(-0.067582421) q[0];
sx q[0];
rz(2.5525868) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(0.20751247) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1197101) q[0];
sx q[0];
rz(-1.3346905) q[0];
sx q[0];
rz(-2.9364999) q[0];
rz(-0.18337266) q[2];
sx q[2];
rz(-1.6499106) q[2];
sx q[2];
rz(2.8734506) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3564975) q[1];
sx q[1];
rz(-1.5743839) q[1];
sx q[1];
rz(-1.5425578) q[1];
rz(-pi) q[2];
rz(-2.9280568) q[3];
sx q[3];
rz(-2.6693404) q[3];
sx q[3];
rz(-2.6708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8091858) q[2];
sx q[2];
rz(-3.1250592) q[2];
sx q[2];
rz(-0.22745505) q[2];
rz(2.4849232) q[3];
sx q[3];
rz(-0.69461099) q[3];
sx q[3];
rz(0.53546661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4559795) q[0];
sx q[0];
rz(-3.1094636) q[0];
sx q[0];
rz(0.44422126) q[0];
rz(-2.0329068) q[1];
sx q[1];
rz(-1.3408835) q[1];
sx q[1];
rz(-1.9734372) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16938528) q[0];
sx q[0];
rz(-0.88645259) q[0];
sx q[0];
rz(1.8207401) q[0];
rz(-1.6158478) q[2];
sx q[2];
rz(-2.0241996) q[2];
sx q[2];
rz(-0.13870961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25623577) q[1];
sx q[1];
rz(-1.1215804) q[1];
sx q[1];
rz(0.083811772) q[1];
rz(0.91932591) q[3];
sx q[3];
rz(-0.72685234) q[3];
sx q[3];
rz(-2.5346699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68219677) q[2];
sx q[2];
rz(-2.0708059) q[2];
sx q[2];
rz(3.0219141) q[2];
rz(-0.97673544) q[3];
sx q[3];
rz(-0.027438199) q[3];
sx q[3];
rz(1.1915092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.479849) q[0];
sx q[0];
rz(-0.16525826) q[0];
sx q[0];
rz(-1.4953493) q[0];
rz(2.7961075) q[1];
sx q[1];
rz(-0.59436878) q[1];
sx q[1];
rz(0.77846175) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0613148) q[0];
sx q[0];
rz(-3.0856371) q[0];
sx q[0];
rz(0.12306889) q[0];
rz(2.5359306) q[2];
sx q[2];
rz(-2.2251943) q[2];
sx q[2];
rz(2.5079648) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.567019) q[1];
sx q[1];
rz(-1.4664093) q[1];
sx q[1];
rz(-1.1280355) q[1];
rz(-pi) q[2];
rz(-0.73466326) q[3];
sx q[3];
rz(-2.0941995) q[3];
sx q[3];
rz(-1.911834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2848795) q[2];
sx q[2];
rz(-3.1046125) q[2];
sx q[2];
rz(1.3899089) q[2];
rz(-0.022627929) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(-2.7271395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49700272) q[0];
sx q[0];
rz(-3.1019326) q[0];
sx q[0];
rz(-0.52763754) q[0];
rz(-0.35609326) q[1];
sx q[1];
rz(-1.5047319) q[1];
sx q[1];
rz(1.5783232) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2882345) q[0];
sx q[0];
rz(-2.4348867) q[0];
sx q[0];
rz(-0.69376365) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4638405) q[2];
sx q[2];
rz(-1.6212109) q[2];
sx q[2];
rz(-2.9843753) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17828748) q[1];
sx q[1];
rz(-1.7366341) q[1];
sx q[1];
rz(2.01009) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08294087) q[3];
sx q[3];
rz(-1.0965875) q[3];
sx q[3];
rz(0.54061962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5712574) q[2];
sx q[2];
rz(-3.1331077) q[2];
sx q[2];
rz(3.0961032) q[2];
rz(-2.3174543) q[3];
sx q[3];
rz(-2.6233311) q[3];
sx q[3];
rz(2.1477108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.6028041) q[0];
sx q[0];
rz(-1.0552009) q[0];
sx q[0];
rz(1.3987199) q[0];
rz(-0.16821965) q[1];
sx q[1];
rz(-1.3361479) q[1];
sx q[1];
rz(0.22163637) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0960064) q[0];
sx q[0];
rz(-1.0500589) q[0];
sx q[0];
rz(2.8686499) q[0];
rz(-pi) q[1];
rz(-1.7847117) q[2];
sx q[2];
rz(-1.6375721) q[2];
sx q[2];
rz(-2.5601088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0757743) q[1];
sx q[1];
rz(-1.5464968) q[1];
sx q[1];
rz(1.1323117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4517253) q[3];
sx q[3];
rz(-0.97056164) q[3];
sx q[3];
rz(-2.8056895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98442709) q[2];
sx q[2];
rz(-0.027313622) q[2];
sx q[2];
rz(1.0013162) q[2];
rz(0.84424132) q[3];
sx q[3];
rz(-3.0735569) q[3];
sx q[3];
rz(-0.74725738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7897414) q[0];
sx q[0];
rz(-0.25775596) q[0];
sx q[0];
rz(1.7716273) q[0];
rz(-0.17579707) q[1];
sx q[1];
rz(-1.5134696) q[1];
sx q[1];
rz(2.083185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31613126) q[0];
sx q[0];
rz(-0.93535813) q[0];
sx q[0];
rz(-2.7881289) q[0];
rz(-pi) q[1];
rz(-1.5878229) q[2];
sx q[2];
rz(-0.58882182) q[2];
sx q[2];
rz(0.39521171) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0936011) q[1];
sx q[1];
rz(-0.8719043) q[1];
sx q[1];
rz(-1.8959787) q[1];
rz(0.36096548) q[3];
sx q[3];
rz(-0.95300337) q[3];
sx q[3];
rz(-1.82919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3093962) q[2];
sx q[2];
rz(-0.0038853566) q[2];
sx q[2];
rz(-2.2957809) q[2];
rz(0.63515615) q[3];
sx q[3];
rz(-0.35717765) q[3];
sx q[3];
rz(2.9966808) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7307067) q[0];
sx q[0];
rz(-0.17179739) q[0];
sx q[0];
rz(2.8806277) q[0];
rz(-1.4313401) q[1];
sx q[1];
rz(-0.15161082) q[1];
sx q[1];
rz(-1.3014911) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430014) q[0];
sx q[0];
rz(-0.81969417) q[0];
sx q[0];
rz(-0.9350594) q[0];
rz(-3.030483) q[2];
sx q[2];
rz(-1.5859563) q[2];
sx q[2];
rz(-2.3738101) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4123685) q[1];
sx q[1];
rz(-0.43392402) q[1];
sx q[1];
rz(-3.0482376) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5673133) q[3];
sx q[3];
rz(-0.26393587) q[3];
sx q[3];
rz(1.6613632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5362376) q[2];
sx q[2];
rz(-1.321512) q[2];
sx q[2];
rz(3.1129254) q[2];
rz(-0.043449314) q[3];
sx q[3];
rz(-0.236792) q[3];
sx q[3];
rz(-1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726149) q[0];
sx q[0];
rz(-2.7942939) q[0];
sx q[0];
rz(0.29796991) q[0];
rz(1.7295674) q[1];
sx q[1];
rz(-2.7822918) q[1];
sx q[1];
rz(-1.5493468) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9095347) q[0];
sx q[0];
rz(-1.2611054) q[0];
sx q[0];
rz(-2.880611) q[0];
rz(-pi) q[1];
rz(1.5875196) q[2];
sx q[2];
rz(-1.564476) q[2];
sx q[2];
rz(0.14594742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1916532) q[1];
sx q[1];
rz(-1.8671163) q[1];
sx q[1];
rz(-0.61298989) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9202254) q[3];
sx q[3];
rz(-1.2949756) q[3];
sx q[3];
rz(0.779895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90184244) q[2];
sx q[2];
rz(-0.032278927) q[2];
sx q[2];
rz(-2.4816404) q[2];
rz(-0.7800855) q[3];
sx q[3];
rz(-1.8476906) q[3];
sx q[3];
rz(-1.9735146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50769794) q[0];
sx q[0];
rz(-0.036490353) q[0];
sx q[0];
rz(-0.80398917) q[0];
rz(2.9102303) q[1];
sx q[1];
rz(-1.7295126) q[1];
sx q[1];
rz(3.0661327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7754525) q[0];
sx q[0];
rz(-2.7260927) q[0];
sx q[0];
rz(2.611931) q[0];
x q[1];
rz(1.6310591) q[2];
sx q[2];
rz(-0.00060877006) q[2];
sx q[2];
rz(3.0920467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8676843) q[1];
sx q[1];
rz(-1.1204506) q[1];
sx q[1];
rz(-0.83012786) q[1];
rz(-1.9583804) q[3];
sx q[3];
rz(-1.6026261) q[3];
sx q[3];
rz(1.4393161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1239473) q[2];
sx q[2];
rz(-2.6804774) q[2];
sx q[2];
rz(2.7258415) q[2];
rz(-1.6425284) q[3];
sx q[3];
rz(-3.1406904) q[3];
sx q[3];
rz(-0.7846964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7163664) q[0];
sx q[0];
rz(-0.056119053) q[0];
sx q[0];
rz(-2.8477493) q[0];
rz(-1.5360688) q[1];
sx q[1];
rz(-2.2686281) q[1];
sx q[1];
rz(1.6764838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241681) q[0];
sx q[0];
rz(-0.10483042) q[0];
sx q[0];
rz(3.0472786) q[0];
rz(-1.15756) q[2];
sx q[2];
rz(-1.0069478) q[2];
sx q[2];
rz(0.68222943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4831381) q[1];
sx q[1];
rz(-2.2171729) q[1];
sx q[1];
rz(1.6049683) q[1];
x q[2];
rz(3.1344101) q[3];
sx q[3];
rz(-1.5752388) q[3];
sx q[3];
rz(-2.9053807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6426223) q[2];
sx q[2];
rz(-0.62752807) q[2];
sx q[2];
rz(1.9210531) q[2];
rz(0.39310655) q[3];
sx q[3];
rz(-3.1246287) q[3];
sx q[3];
rz(0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4569693) q[0];
sx q[0];
rz(-1.4263117) q[0];
sx q[0];
rz(-0.34761467) q[0];
rz(-2.5034703) q[1];
sx q[1];
rz(-0.77824021) q[1];
sx q[1];
rz(2.4660769) q[1];
rz(2.9478922) q[2];
sx q[2];
rz(-1.4932409) q[2];
sx q[2];
rz(1.3328339) q[2];
rz(0.06107851) q[3];
sx q[3];
rz(-0.11299639) q[3];
sx q[3];
rz(-1.1391409) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
