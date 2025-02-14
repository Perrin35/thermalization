OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78703824) q[0];
sx q[0];
rz(-0.7631425) q[0];
sx q[0];
rz(0.95844498) q[0];
rz(-2.7148442) q[1];
sx q[1];
rz(3.5659748) q[1];
sx q[1];
rz(10.397484) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6329968) q[0];
sx q[0];
rz(-0.22211123) q[0];
sx q[0];
rz(2.7048777) q[0];
rz(-2.7673253) q[2];
sx q[2];
rz(-0.18478811) q[2];
sx q[2];
rz(-1.6109465) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0578459) q[1];
sx q[1];
rz(-1.3120756) q[1];
sx q[1];
rz(0.73212917) q[1];
rz(1.7250502) q[3];
sx q[3];
rz(-1.2917545) q[3];
sx q[3];
rz(2.0786301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0613609) q[2];
sx q[2];
rz(-1.9539359) q[2];
sx q[2];
rz(2.3415671) q[2];
rz(-3.0112265) q[3];
sx q[3];
rz(-0.76821199) q[3];
sx q[3];
rz(-0.31726328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80113634) q[0];
sx q[0];
rz(-1.9213333) q[0];
sx q[0];
rz(-2.1521547) q[0];
rz(-2.2438352) q[1];
sx q[1];
rz(-1.0549301) q[1];
sx q[1];
rz(-0.78136939) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5637276) q[0];
sx q[0];
rz(-1.0672309) q[0];
sx q[0];
rz(-0.2204075) q[0];
rz(2.7861681) q[2];
sx q[2];
rz(-0.249042) q[2];
sx q[2];
rz(0.058573478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1220391) q[1];
sx q[1];
rz(-1.0940503) q[1];
sx q[1];
rz(-2.6470584) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7584293) q[3];
sx q[3];
rz(-1.8472611) q[3];
sx q[3];
rz(-0.081646669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.024461688) q[2];
sx q[2];
rz(-1.3776366) q[2];
sx q[2];
rz(-2.3632939) q[2];
rz(-2.4513169) q[3];
sx q[3];
rz(-2.2199151) q[3];
sx q[3];
rz(-1.5604304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86163259) q[0];
sx q[0];
rz(-2.3065688) q[0];
sx q[0];
rz(-2.7296208) q[0];
rz(2.1906134) q[1];
sx q[1];
rz(-2.488766) q[1];
sx q[1];
rz(-1.9297809) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095603) q[0];
sx q[0];
rz(-1.6044084) q[0];
sx q[0];
rz(-1.6092954) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0427709) q[2];
sx q[2];
rz(-1.3490236) q[2];
sx q[2];
rz(-2.2840471) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7344021) q[1];
sx q[1];
rz(-0.71043321) q[1];
sx q[1];
rz(-3.0991952) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7318299) q[3];
sx q[3];
rz(-1.2098243) q[3];
sx q[3];
rz(-1.7019113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.95436207) q[2];
sx q[2];
rz(-0.53202859) q[2];
sx q[2];
rz(-1.5032035) q[2];
rz(3.1014118) q[3];
sx q[3];
rz(-2.2153722) q[3];
sx q[3];
rz(-0.23475501) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16315854) q[0];
sx q[0];
rz(-1.5143159) q[0];
sx q[0];
rz(0.12983313) q[0];
rz(1.8561329) q[1];
sx q[1];
rz(-1.9655656) q[1];
sx q[1];
rz(1.0349549) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87405768) q[0];
sx q[0];
rz(-2.3685936) q[0];
sx q[0];
rz(-0.92138793) q[0];
rz(-pi) q[1];
rz(2.2345654) q[2];
sx q[2];
rz(-1.4849147) q[2];
sx q[2];
rz(-2.2965388) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.6432727) q[1];
sx q[1];
rz(-0.26277143) q[1];
sx q[1];
rz(-2.8064125) q[1];
x q[2];
rz(-2.2620151) q[3];
sx q[3];
rz(-1.1926418) q[3];
sx q[3];
rz(0.86728269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3639565) q[2];
sx q[2];
rz(-1.3866321) q[2];
sx q[2];
rz(-3.0470972) q[2];
rz(-2.7206521) q[3];
sx q[3];
rz(-1.4237483) q[3];
sx q[3];
rz(-1.7054935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41973758) q[0];
sx q[0];
rz(-2.5696745) q[0];
sx q[0];
rz(-3.1255334) q[0];
rz(2.2968966) q[1];
sx q[1];
rz(-0.41929308) q[1];
sx q[1];
rz(-2.2339581) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43953376) q[0];
sx q[0];
rz(-0.82111135) q[0];
sx q[0];
rz(-1.8762212) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0727022) q[2];
sx q[2];
rz(-2.6703983) q[2];
sx q[2];
rz(2.1387177) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4872949) q[1];
sx q[1];
rz(-1.134344) q[1];
sx q[1];
rz(-2.7769293) q[1];
x q[2];
rz(-1.5974371) q[3];
sx q[3];
rz(-1.9147885) q[3];
sx q[3];
rz(2.2709803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5661261) q[2];
sx q[2];
rz(-1.0593654) q[2];
sx q[2];
rz(-3.0779085) q[2];
rz(0.56695402) q[3];
sx q[3];
rz(-2.3072115) q[3];
sx q[3];
rz(2.5440192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-0.31396922) q[0];
sx q[0];
rz(-2.9997928) q[0];
sx q[0];
rz(-1.9837448) q[0];
rz(-1.6572378) q[1];
sx q[1];
rz(-1.6969705) q[1];
sx q[1];
rz(2.8725913) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3753178) q[0];
sx q[0];
rz(-0.14785375) q[0];
sx q[0];
rz(-2.8327441) q[0];
x q[1];
rz(-1.8431191) q[2];
sx q[2];
rz(-1.6629873) q[2];
sx q[2];
rz(2.628679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29998955) q[1];
sx q[1];
rz(-1.4714452) q[1];
sx q[1];
rz(1.7510115) q[1];
rz(2.727521) q[3];
sx q[3];
rz(-1.7709289) q[3];
sx q[3];
rz(-2.4614863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77702648) q[2];
sx q[2];
rz(-2.624161) q[2];
sx q[2];
rz(0.88097921) q[2];
rz(-3.0173054) q[3];
sx q[3];
rz(-1.4275987) q[3];
sx q[3];
rz(0.65829268) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86150259) q[0];
sx q[0];
rz(-0.26743356) q[0];
sx q[0];
rz(2.4524443) q[0];
rz(2.6742477) q[1];
sx q[1];
rz(-0.57599774) q[1];
sx q[1];
rz(2.0435832) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2158581) q[0];
sx q[0];
rz(-1.6491062) q[0];
sx q[0];
rz(2.9758478) q[0];
rz(-pi) q[1];
rz(1.8703837) q[2];
sx q[2];
rz(-1.0429842) q[2];
sx q[2];
rz(0.53003788) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6642521) q[1];
sx q[1];
rz(-2.1108529) q[1];
sx q[1];
rz(-1.5332721) q[1];
rz(-2.2192248) q[3];
sx q[3];
rz(-0.73966129) q[3];
sx q[3];
rz(1.9112223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.565862) q[2];
sx q[2];
rz(-2.2634759) q[2];
sx q[2];
rz(-0.6305387) q[2];
rz(-2.5343043) q[3];
sx q[3];
rz(-2.2346965) q[3];
sx q[3];
rz(-2.0489395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8232515) q[0];
sx q[0];
rz(-0.14597758) q[0];
sx q[0];
rz(-1.7952221) q[0];
rz(1.3660376) q[1];
sx q[1];
rz(-2.4925158) q[1];
sx q[1];
rz(-0.62873658) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7766984) q[0];
sx q[0];
rz(-1.3140251) q[0];
sx q[0];
rz(-2.5904694) q[0];
x q[1];
rz(0.5127617) q[2];
sx q[2];
rz(-2.1760941) q[2];
sx q[2];
rz(-2.5255741) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.47766108) q[1];
sx q[1];
rz(-1.0870502) q[1];
sx q[1];
rz(-0.61772924) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2090861) q[3];
sx q[3];
rz(-0.21245281) q[3];
sx q[3];
rz(0.3750876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2748572) q[2];
sx q[2];
rz(-1.3254415) q[2];
sx q[2];
rz(2.2692915) q[2];
rz(-0.21271475) q[3];
sx q[3];
rz(-1.745696) q[3];
sx q[3];
rz(-2.5581636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.23671959) q[0];
sx q[0];
rz(-1.5308335) q[0];
sx q[0];
rz(-1.8778296) q[0];
rz(1.0172552) q[1];
sx q[1];
rz(-2.2433498) q[1];
sx q[1];
rz(-2.5668626) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1270458) q[0];
sx q[0];
rz(-2.1009863) q[0];
sx q[0];
rz(-1.4958025) q[0];
x q[1];
rz(-0.82860126) q[2];
sx q[2];
rz(-1.7485837) q[2];
sx q[2];
rz(-1.0324024) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8239261) q[1];
sx q[1];
rz(-1.9374018) q[1];
sx q[1];
rz(1.3701141) q[1];
rz(1.0887126) q[3];
sx q[3];
rz(-1.7974554) q[3];
sx q[3];
rz(0.72666336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42402521) q[2];
sx q[2];
rz(-1.8342476) q[2];
sx q[2];
rz(0.032111017) q[2];
rz(-1.1561681) q[3];
sx q[3];
rz(-2.7572032) q[3];
sx q[3];
rz(-1.9305852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1011937) q[0];
sx q[0];
rz(-2.5960584) q[0];
sx q[0];
rz(-2.0174761) q[0];
rz(2.595937) q[1];
sx q[1];
rz(-2.7898495) q[1];
sx q[1];
rz(2.5901332) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41697793) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(-0.12599385) q[0];
x q[1];
rz(1.3842756) q[2];
sx q[2];
rz(-0.2787481) q[2];
sx q[2];
rz(0.84183019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6956967) q[1];
sx q[1];
rz(-2.4425382) q[1];
sx q[1];
rz(2.3732144) q[1];
rz(-pi) q[2];
rz(1.7343246) q[3];
sx q[3];
rz(-1.2683322) q[3];
sx q[3];
rz(1.0176942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5896899) q[2];
sx q[2];
rz(-0.51108131) q[2];
sx q[2];
rz(1.4116633) q[2];
rz(-0.74665135) q[3];
sx q[3];
rz(-1.4614429) q[3];
sx q[3];
rz(2.1667229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.830037) q[0];
sx q[0];
rz(-1.5604326) q[0];
sx q[0];
rz(-1.5695705) q[0];
rz(0.39881067) q[1];
sx q[1];
rz(-1.3347944) q[1];
sx q[1];
rz(1.4904108) q[1];
rz(1.082303) q[2];
sx q[2];
rz(-1.5776967) q[2];
sx q[2];
rz(-2.7803905) q[2];
rz(-0.79979782) q[3];
sx q[3];
rz(-1.3577438) q[3];
sx q[3];
rz(1.8997026) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
