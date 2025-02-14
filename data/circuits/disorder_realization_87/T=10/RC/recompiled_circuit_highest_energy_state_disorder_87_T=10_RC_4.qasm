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
rz(-1.6753766) q[0];
sx q[0];
rz(-2.197062) q[0];
sx q[0];
rz(2.9401927) q[0];
rz(1.072285) q[1];
sx q[1];
rz(-0.76289248) q[1];
sx q[1];
rz(1.4210757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3823377) q[0];
sx q[0];
rz(-1.6256285) q[0];
sx q[0];
rz(-2.1817529) q[0];
rz(-0.95008738) q[2];
sx q[2];
rz(-2.3664775) q[2];
sx q[2];
rz(2.895854) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3609429) q[1];
sx q[1];
rz(-2.2113178) q[1];
sx q[1];
rz(-2.065413) q[1];
rz(-pi) q[2];
rz(1.9626612) q[3];
sx q[3];
rz(-1.2614935) q[3];
sx q[3];
rz(-2.7397857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1066771) q[2];
sx q[2];
rz(-0.80433977) q[2];
sx q[2];
rz(2.8625281) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(-2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-2.294367) q[0];
sx q[0];
rz(-2.8952059) q[0];
rz(-1.1950182) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(-1.6349207) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5824597) q[0];
sx q[0];
rz(-1.8163067) q[0];
sx q[0];
rz(-2.0364291) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8737239) q[2];
sx q[2];
rz(-1.5298577) q[2];
sx q[2];
rz(0.28266476) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7333201) q[1];
sx q[1];
rz(-2.773914) q[1];
sx q[1];
rz(-1.7321083) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89119567) q[3];
sx q[3];
rz(-2.5053484) q[3];
sx q[3];
rz(1.4435857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(-1.231989) q[2];
rz(-2.039382) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6983637) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(-2.2337636) q[0];
rz(-0.56697956) q[1];
sx q[1];
rz(-0.98273977) q[1];
sx q[1];
rz(-0.10279113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55926052) q[0];
sx q[0];
rz(-2.1301422) q[0];
sx q[0];
rz(1.8190967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6750828) q[2];
sx q[2];
rz(-1.6513746) q[2];
sx q[2];
rz(1.8957418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5788865) q[1];
sx q[1];
rz(-1.9745312) q[1];
sx q[1];
rz(2.5078234) q[1];
x q[2];
rz(2.1249173) q[3];
sx q[3];
rz(-1.8401056) q[3];
sx q[3];
rz(0.96162063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28321442) q[2];
sx q[2];
rz(-0.80867043) q[2];
sx q[2];
rz(-1.4374479) q[2];
rz(-0.39250675) q[3];
sx q[3];
rz(-2.0065353) q[3];
sx q[3];
rz(2.8431622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0482386) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(1.1887953) q[0];
rz(2.8554754) q[1];
sx q[1];
rz(-0.61758271) q[1];
sx q[1];
rz(-1.6292705) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28933576) q[0];
sx q[0];
rz(-1.8484383) q[0];
sx q[0];
rz(1.0278948) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2354765) q[2];
sx q[2];
rz(-2.309707) q[2];
sx q[2];
rz(3.0201833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0864073) q[1];
sx q[1];
rz(-0.98208222) q[1];
sx q[1];
rz(0.22365315) q[1];
rz(-pi) q[2];
rz(2.9305601) q[3];
sx q[3];
rz(-1.0730626) q[3];
sx q[3];
rz(-0.56609234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7901223) q[2];
sx q[2];
rz(-1.0834379) q[2];
sx q[2];
rz(2.4646087) q[2];
rz(1.8949159) q[3];
sx q[3];
rz(-1.8691984) q[3];
sx q[3];
rz(-1.5164794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93382728) q[0];
sx q[0];
rz(-0.25660577) q[0];
sx q[0];
rz(-3.1166792) q[0];
rz(2.086153) q[1];
sx q[1];
rz(-0.98506227) q[1];
sx q[1];
rz(-2.8581462) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.442207) q[0];
sx q[0];
rz(-1.3106579) q[0];
sx q[0];
rz(-2.8752863) q[0];
rz(-pi) q[1];
rz(0.13155664) q[2];
sx q[2];
rz(-1.7649998) q[2];
sx q[2];
rz(2.921791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3230309) q[1];
sx q[1];
rz(-0.50172443) q[1];
sx q[1];
rz(1.6343568) q[1];
x q[2];
rz(-2.2437566) q[3];
sx q[3];
rz(-1.8404318) q[3];
sx q[3];
rz(1.8561038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.31263605) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(-2.8157595) q[2];
rz(-1.8587941) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(-1.5939943) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7376937) q[0];
sx q[0];
rz(-1.6588545) q[0];
sx q[0];
rz(2.7666336) q[0];
rz(-0.47863475) q[1];
sx q[1];
rz(-0.67239434) q[1];
sx q[1];
rz(0.37014827) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70765342) q[0];
sx q[0];
rz(-1.5591199) q[0];
sx q[0];
rz(-0.023109584) q[0];
rz(-pi) q[1];
rz(-2.213573) q[2];
sx q[2];
rz(-0.22980873) q[2];
sx q[2];
rz(1.9191051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9795591) q[1];
sx q[1];
rz(-2.0311072) q[1];
sx q[1];
rz(0.76785894) q[1];
rz(-pi) q[2];
rz(-3.1305976) q[3];
sx q[3];
rz(-0.42783005) q[3];
sx q[3];
rz(1.7181991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7438573) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.3810623) q[2];
rz(2.0928275) q[3];
sx q[3];
rz(-1.3449113) q[3];
sx q[3];
rz(0.63961187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25310707) q[0];
sx q[0];
rz(-2.3024004) q[0];
sx q[0];
rz(2.2174368) q[0];
rz(-2.2249075) q[1];
sx q[1];
rz(-2.075383) q[1];
sx q[1];
rz(-1.7402657) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7766905) q[0];
sx q[0];
rz(-0.49743891) q[0];
sx q[0];
rz(3.0203392) q[0];
rz(-1.3989053) q[2];
sx q[2];
rz(-2.2187833) q[2];
sx q[2];
rz(1.5251708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37937134) q[1];
sx q[1];
rz(-1.4109542) q[1];
sx q[1];
rz(2.400378) q[1];
rz(-0.21715607) q[3];
sx q[3];
rz(-2.0494866) q[3];
sx q[3];
rz(1.1622805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9132729) q[2];
sx q[2];
rz(-1.6797804) q[2];
sx q[2];
rz(1.9165967) q[2];
rz(-3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(-0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5834354) q[0];
sx q[0];
rz(-0.83034101) q[0];
sx q[0];
rz(2.661327) q[0];
rz(-0.14446124) q[1];
sx q[1];
rz(-0.90765777) q[1];
sx q[1];
rz(-0.86437782) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1594161) q[0];
sx q[0];
rz(-1.5242531) q[0];
sx q[0];
rz(-1.7345364) q[0];
x q[1];
rz(-3.046918) q[2];
sx q[2];
rz(-1.6045158) q[2];
sx q[2];
rz(2.0168436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4982953) q[1];
sx q[1];
rz(-1.8113231) q[1];
sx q[1];
rz(2.125419) q[1];
rz(-0.79473991) q[3];
sx q[3];
rz(-0.38215548) q[3];
sx q[3];
rz(-0.69821815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(-0.63507357) q[2];
rz(-2.5687929) q[3];
sx q[3];
rz(-1.3558931) q[3];
sx q[3];
rz(-2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11005814) q[0];
sx q[0];
rz(-0.77339554) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(2.7763413) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(1.431538) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99125049) q[0];
sx q[0];
rz(-1.7479154) q[0];
sx q[0];
rz(0.19512531) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1373491) q[2];
sx q[2];
rz(-2.1449617) q[2];
sx q[2];
rz(3.1176709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.46043432) q[1];
sx q[1];
rz(-1.2660813) q[1];
sx q[1];
rz(2.4439993) q[1];
x q[2];
rz(-1.133119) q[3];
sx q[3];
rz(-1.0854377) q[3];
sx q[3];
rz(0.35972586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.26620904) q[2];
sx q[2];
rz(-0.92140809) q[2];
sx q[2];
rz(-1.1727775) q[2];
rz(-1.4069936) q[3];
sx q[3];
rz(-1.3318136) q[3];
sx q[3];
rz(1.9269358) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239546) q[0];
sx q[0];
rz(-1.1347436) q[0];
sx q[0];
rz(-2.5469575) q[0];
rz(1.0058962) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(2.2524021) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.206736) q[0];
sx q[0];
rz(-1.6214658) q[0];
sx q[0];
rz(1.3502747) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0709549) q[2];
sx q[2];
rz(-1.1114235) q[2];
sx q[2];
rz(0.86462155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31510776) q[1];
sx q[1];
rz(-1.4764305) q[1];
sx q[1];
rz(-3.0815275) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6088283) q[3];
sx q[3];
rz(-1.0547148) q[3];
sx q[3];
rz(2.3446159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95120007) q[2];
sx q[2];
rz(-1.9659646) q[2];
sx q[2];
rz(2.5962043) q[2];
rz(0.96327463) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(-1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6496898) q[0];
sx q[0];
rz(-0.90072537) q[0];
sx q[0];
rz(-1.8186722) q[0];
rz(-2.2979965) q[1];
sx q[1];
rz(-2.0633162) q[1];
sx q[1];
rz(1.8819527) q[1];
rz(0.19914535) q[2];
sx q[2];
rz(-1.1391339) q[2];
sx q[2];
rz(1.7720411) q[2];
rz(-2.4424408) q[3];
sx q[3];
rz(-2.4300162) q[3];
sx q[3];
rz(3.1292849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
