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
rz(1.2122756) q[0];
sx q[0];
rz(4.2733856) q[0];
sx q[0];
rz(10.797664) q[0];
rz(-1.1276487) q[1];
sx q[1];
rz(-1.375066) q[1];
sx q[1];
rz(-2.3553203) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72981156) q[0];
sx q[0];
rz(-1.775911) q[0];
sx q[0];
rz(-1.9461483) q[0];
rz(-pi) q[1];
rz(-0.89531548) q[2];
sx q[2];
rz(-2.3626872) q[2];
sx q[2];
rz(-1.9112183) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3481625) q[1];
sx q[1];
rz(-1.5409267) q[1];
sx q[1];
rz(-1.8601599) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97162515) q[3];
sx q[3];
rz(-1.2279945) q[3];
sx q[3];
rz(2.9620217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11975153) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(-1.2500259) q[2];
rz(2.405808) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(-1.638394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64330548) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(-0.071320891) q[0];
rz(1.9885063) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(-2.6128795) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604989) q[0];
sx q[0];
rz(-2.0369663) q[0];
sx q[0];
rz(2.7866298) q[0];
rz(-pi) q[1];
rz(-2.2430482) q[2];
sx q[2];
rz(-0.29783422) q[2];
sx q[2];
rz(2.1517449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9829927) q[1];
sx q[1];
rz(-2.0865177) q[1];
sx q[1];
rz(-0.92293585) q[1];
rz(-pi) q[2];
rz(2.9170179) q[3];
sx q[3];
rz(-1.16515) q[3];
sx q[3];
rz(-0.39164513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4802287) q[2];
sx q[2];
rz(-0.39863786) q[2];
sx q[2];
rz(-0.028623494) q[2];
rz(-0.35401595) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(0.55907512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70327586) q[0];
sx q[0];
rz(-1.0030614) q[0];
sx q[0];
rz(2.2502374) q[0];
rz(2.640653) q[1];
sx q[1];
rz(-1.5762065) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1315434) q[0];
sx q[0];
rz(-1.5230522) q[0];
sx q[0];
rz(-3.0555807) q[0];
rz(-pi) q[1];
rz(2.8562653) q[2];
sx q[2];
rz(-1.0614191) q[2];
sx q[2];
rz(0.90778186) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.14172983) q[1];
sx q[1];
rz(-1.0838306) q[1];
sx q[1];
rz(2.3906863) q[1];
rz(2.9611582) q[3];
sx q[3];
rz(-1.5984319) q[3];
sx q[3];
rz(-2.7829426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-0.55861837) q[2];
sx q[2];
rz(0.2365665) q[2];
rz(-0.92075721) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(-1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9193566) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(-2.2663569) q[0];
rz(1.5454166) q[1];
sx q[1];
rz(-2.2794006) q[1];
sx q[1];
rz(-0.045305591) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52462477) q[0];
sx q[0];
rz(-0.12147203) q[0];
sx q[0];
rz(1.0018906) q[0];
rz(-pi) q[1];
rz(-2.1854109) q[2];
sx q[2];
rz(-2.3740157) q[2];
sx q[2];
rz(-2.0337348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6958624) q[1];
sx q[1];
rz(-0.67477422) q[1];
sx q[1];
rz(-1.710102) q[1];
x q[2];
rz(1.8604467) q[3];
sx q[3];
rz(-2.5389606) q[3];
sx q[3];
rz(1.2887736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8594325) q[2];
sx q[2];
rz(-0.96264797) q[2];
sx q[2];
rz(1.1939987) q[2];
rz(0.89030877) q[3];
sx q[3];
rz(-1.9186391) q[3];
sx q[3];
rz(2.1931026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0328261) q[0];
sx q[0];
rz(-0.81619167) q[0];
sx q[0];
rz(-0.68242514) q[0];
rz(-1.2884864) q[1];
sx q[1];
rz(-1.5623743) q[1];
sx q[1];
rz(-1.8103745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3159861) q[0];
sx q[0];
rz(-2.2974112) q[0];
sx q[0];
rz(-3.1196703) q[0];
x q[1];
rz(1.9648026) q[2];
sx q[2];
rz(-1.8117684) q[2];
sx q[2];
rz(0.35277982) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2878804) q[1];
sx q[1];
rz(-1.2128218) q[1];
sx q[1];
rz(-1.720721) q[1];
x q[2];
rz(-0.87792895) q[3];
sx q[3];
rz(-2.4932541) q[3];
sx q[3];
rz(0.058365783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0266626) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(-2.4746573) q[2];
rz(0.55495787) q[3];
sx q[3];
rz(-0.68735492) q[3];
sx q[3];
rz(0.2989029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126467) q[0];
sx q[0];
rz(-1.8829367) q[0];
sx q[0];
rz(-0.70372787) q[0];
rz(0.16920371) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(-2.9827548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.279351) q[0];
sx q[0];
rz(-2.9547174) q[0];
sx q[0];
rz(0.86743768) q[0];
rz(-pi) q[1];
rz(-1.6936982) q[2];
sx q[2];
rz(-1.7235061) q[2];
sx q[2];
rz(-1.0433152) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1353246) q[1];
sx q[1];
rz(-0.49990955) q[1];
sx q[1];
rz(0.27854021) q[1];
rz(0.63780906) q[3];
sx q[3];
rz(-1.3862228) q[3];
sx q[3];
rz(1.0705494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3682897) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(-2.2840195) q[2];
rz(-1.9722021) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-0.20839553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(0.17539772) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(0.71863693) q[0];
rz(-0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(2.4729572) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91808337) q[0];
sx q[0];
rz(-1.232172) q[0];
sx q[0];
rz(2.2302367) q[0];
x q[1];
rz(0.35923194) q[2];
sx q[2];
rz(-2.1678535) q[2];
sx q[2];
rz(-0.12166858) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9978511) q[1];
sx q[1];
rz(-2.6392548) q[1];
sx q[1];
rz(-0.84059244) q[1];
x q[2];
rz(0.63124048) q[3];
sx q[3];
rz(-1.5774014) q[3];
sx q[3];
rz(-0.94830482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4947027) q[2];
sx q[2];
rz(-1.1855482) q[2];
sx q[2];
rz(-1.0364214) q[2];
rz(-2.5070665) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(0.79276597) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46102872) q[0];
sx q[0];
rz(-3.0240318) q[0];
sx q[0];
rz(-2.4155937) q[0];
rz(-1.1622102) q[1];
sx q[1];
rz(-1.1770959) q[1];
sx q[1];
rz(1.1964218) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0045753) q[0];
sx q[0];
rz(-1.5992303) q[0];
sx q[0];
rz(-1.397055) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4216719) q[2];
sx q[2];
rz(-0.31055488) q[2];
sx q[2];
rz(2.9322185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6641621) q[1];
sx q[1];
rz(-2.4242085) q[1];
sx q[1];
rz(0.75023164) q[1];
rz(0.67025028) q[3];
sx q[3];
rz(-1.2941417) q[3];
sx q[3];
rz(-1.2069699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2339345) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(-0.93351239) q[2];
rz(-2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(-2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1652949) q[0];
sx q[0];
rz(-0.10460654) q[0];
sx q[0];
rz(3.0883375) q[0];
rz(0.28757295) q[1];
sx q[1];
rz(-2.2122999) q[1];
sx q[1];
rz(2.8823749) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027336425) q[0];
sx q[0];
rz(-1.5467318) q[0];
sx q[0];
rz(-1.5177478) q[0];
rz(-0.13287787) q[2];
sx q[2];
rz(-1.3406375) q[2];
sx q[2];
rz(-1.2990189) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99154918) q[1];
sx q[1];
rz(-0.75356149) q[1];
sx q[1];
rz(-1.4808473) q[1];
x q[2];
rz(0.15786981) q[3];
sx q[3];
rz(-0.77220687) q[3];
sx q[3];
rz(-1.3881573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58664924) q[2];
sx q[2];
rz(-1.8879075) q[2];
sx q[2];
rz(2.9602236) q[2];
rz(-0.3950611) q[3];
sx q[3];
rz(-0.22160465) q[3];
sx q[3];
rz(2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7419389) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(-2.7594866) q[0];
rz(-0.29640472) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(1.872725) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5345602) q[0];
sx q[0];
rz(-1.9968613) q[0];
sx q[0];
rz(-1.0491648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13258719) q[2];
sx q[2];
rz(-0.90182038) q[2];
sx q[2];
rz(-2.7821409) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8812534) q[1];
sx q[1];
rz(-1.2311544) q[1];
sx q[1];
rz(-0.022625523) q[1];
rz(-0.24618547) q[3];
sx q[3];
rz(-2.6757112) q[3];
sx q[3];
rz(-1.4108489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84539139) q[2];
sx q[2];
rz(-2.5808344) q[2];
sx q[2];
rz(-2.7517547) q[2];
rz(-1.2975533) q[3];
sx q[3];
rz(-1.9997948) q[3];
sx q[3];
rz(0.84387422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1525477) q[0];
sx q[0];
rz(-1.7392409) q[0];
sx q[0];
rz(1.637511) q[0];
rz(-1.7851495) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(-0.46089725) q[2];
sx q[2];
rz(-0.97958889) q[2];
sx q[2];
rz(-0.19340672) q[2];
rz(-1.9080347) q[3];
sx q[3];
rz(-2.3261286) q[3];
sx q[3];
rz(0.89589768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
