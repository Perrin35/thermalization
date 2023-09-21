OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80823094) q[0];
sx q[0];
rz(-2.1699177) q[0];
sx q[0];
rz(1.4810286) q[0];
rz(-0.12227585) q[1];
sx q[1];
rz(3.0552157) q[1];
sx q[1];
rz(9.4019158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4532783) q[0];
sx q[0];
rz(-0.18916057) q[0];
sx q[0];
rz(-2.2384089) q[0];
rz(-pi) q[1];
rz(2.0055662) q[2];
sx q[2];
rz(-2.044894) q[2];
sx q[2];
rz(-2.9413162) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8926156) q[1];
sx q[1];
rz(-1.2274582) q[1];
sx q[1];
rz(-0.15906072) q[1];
rz(-pi) q[2];
rz(-1.919235) q[3];
sx q[3];
rz(-0.64823965) q[3];
sx q[3];
rz(-0.1227619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4188529) q[2];
sx q[2];
rz(-2.3143694) q[2];
sx q[2];
rz(-2.375405) q[2];
rz(-2.9700759) q[3];
sx q[3];
rz(-0.73232108) q[3];
sx q[3];
rz(2.5550301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8121174) q[0];
sx q[0];
rz(-0.43918878) q[0];
sx q[0];
rz(3.0549333) q[0];
rz(-0.86241972) q[1];
sx q[1];
rz(-2.5984867) q[1];
sx q[1];
rz(-0.08509732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9239685) q[0];
sx q[0];
rz(-2.1930165) q[0];
sx q[0];
rz(1.5505962) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1727824) q[2];
sx q[2];
rz(-1.8245056) q[2];
sx q[2];
rz(1.3041376) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7460691) q[1];
sx q[1];
rz(-1.6568526) q[1];
sx q[1];
rz(-1.6834016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99382932) q[3];
sx q[3];
rz(-0.54518632) q[3];
sx q[3];
rz(2.546437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26248419) q[2];
sx q[2];
rz(-1.5905249) q[2];
sx q[2];
rz(-1.4651728) q[2];
rz(1.9654467) q[3];
sx q[3];
rz(-2.2361103) q[3];
sx q[3];
rz(2.8951077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48433205) q[0];
sx q[0];
rz(-0.13467877) q[0];
sx q[0];
rz(2.235967) q[0];
rz(-0.33272818) q[1];
sx q[1];
rz(-0.85488027) q[1];
sx q[1];
rz(1.3844301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8850088) q[0];
sx q[0];
rz(-0.48261595) q[0];
sx q[0];
rz(0.17871876) q[0];
x q[1];
rz(2.076782) q[2];
sx q[2];
rz(-1.9188606) q[2];
sx q[2];
rz(0.51007523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3570909) q[1];
sx q[1];
rz(-2.38584) q[1];
sx q[1];
rz(0.032553629) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3178187) q[3];
sx q[3];
rz(-2.4274821) q[3];
sx q[3];
rz(-2.2658474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24421346) q[2];
sx q[2];
rz(-2.801557) q[2];
sx q[2];
rz(1.8133694) q[2];
rz(1.3978847) q[3];
sx q[3];
rz(-0.54026794) q[3];
sx q[3];
rz(3.0298997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0760647) q[0];
sx q[0];
rz(-2.9632443) q[0];
sx q[0];
rz(0.67101014) q[0];
rz(1.7975851) q[1];
sx q[1];
rz(-1.9799415) q[1];
sx q[1];
rz(2.9342594) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5642693) q[0];
sx q[0];
rz(-0.042357001) q[0];
sx q[0];
rz(-1.3576515) q[0];
rz(-pi) q[1];
rz(1.3968799) q[2];
sx q[2];
rz(-1.5408278) q[2];
sx q[2];
rz(-2.7001691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3934717) q[1];
sx q[1];
rz(-1.8078184) q[1];
sx q[1];
rz(-0.15740983) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6299423) q[3];
sx q[3];
rz(-1.5582065) q[3];
sx q[3];
rz(2.1714641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.31071445) q[2];
sx q[2];
rz(-2.0534616) q[2];
sx q[2];
rz(-0.80835289) q[2];
rz(0.80777848) q[3];
sx q[3];
rz(-0.72179663) q[3];
sx q[3];
rz(0.19367735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59072524) q[0];
sx q[0];
rz(-0.71722513) q[0];
sx q[0];
rz(0.89548683) q[0];
rz(0.26516178) q[1];
sx q[1];
rz(-2.3064955) q[1];
sx q[1];
rz(-2.6079544) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95414872) q[0];
sx q[0];
rz(-1.2866409) q[0];
sx q[0];
rz(0.28676333) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3445504) q[2];
sx q[2];
rz(-0.81586736) q[2];
sx q[2];
rz(-2.0362542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0698318) q[1];
sx q[1];
rz(-1.3511718) q[1];
sx q[1];
rz(-1.7204309) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2002033) q[3];
sx q[3];
rz(-1.7377059) q[3];
sx q[3];
rz(-3.133579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2395893) q[2];
sx q[2];
rz(-1.5782372) q[2];
sx q[2];
rz(2.5732102) q[2];
rz(-1.2166294) q[3];
sx q[3];
rz(-0.47074461) q[3];
sx q[3];
rz(2.7980428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6333273) q[0];
sx q[0];
rz(-2.2631622) q[0];
sx q[0];
rz(-1.1761965) q[0];
rz(1.5856702) q[1];
sx q[1];
rz(-0.95247477) q[1];
sx q[1];
rz(-1.0046545) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3096136) q[0];
sx q[0];
rz(-1.2011315) q[0];
sx q[0];
rz(1.2020822) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4899646) q[2];
sx q[2];
rz(-1.1751375) q[2];
sx q[2];
rz(2.253502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.172563) q[1];
sx q[1];
rz(-1.8168212) q[1];
sx q[1];
rz(2.5804156) q[1];
rz(-0.011990487) q[3];
sx q[3];
rz(-1.9846091) q[3];
sx q[3];
rz(-0.12366611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7375609) q[2];
sx q[2];
rz(-2.1514251) q[2];
sx q[2];
rz(-0.44815865) q[2];
rz(2.8435977) q[3];
sx q[3];
rz(-0.69152504) q[3];
sx q[3];
rz(0.34860778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.3649243) q[0];
sx q[0];
rz(-1.8402599) q[0];
sx q[0];
rz(-0.65814322) q[0];
rz(-1.8572042) q[1];
sx q[1];
rz(-2.6874266) q[1];
sx q[1];
rz(0.67214322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1561688) q[0];
sx q[0];
rz(-0.16616136) q[0];
sx q[0];
rz(1.6155433) q[0];
rz(-2.0124112) q[2];
sx q[2];
rz(-1.1814983) q[2];
sx q[2];
rz(-1.0600952) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8183648) q[1];
sx q[1];
rz(-2.0548901) q[1];
sx q[1];
rz(1.1546043) q[1];
x q[2];
rz(1.4268731) q[3];
sx q[3];
rz(-0.87464273) q[3];
sx q[3];
rz(-1.5211943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.60823524) q[2];
sx q[2];
rz(-0.63295263) q[2];
sx q[2];
rz(-0.39880025) q[2];
rz(-2.3170025) q[3];
sx q[3];
rz(-1.3922858) q[3];
sx q[3];
rz(2.5430172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44052112) q[0];
sx q[0];
rz(-0.43061391) q[0];
sx q[0];
rz(-0.15821247) q[0];
rz(-1.0868866) q[1];
sx q[1];
rz(-2.465076) q[1];
sx q[1];
rz(-2.3349082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0643443) q[0];
sx q[0];
rz(-1.4532538) q[0];
sx q[0];
rz(2.006152) q[0];
x q[1];
rz(-0.20571795) q[2];
sx q[2];
rz(-1.4652025) q[2];
sx q[2];
rz(2.2149057) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9843236) q[1];
sx q[1];
rz(-1.8586129) q[1];
sx q[1];
rz(-0.70043522) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8716378) q[3];
sx q[3];
rz(-1.5402208) q[3];
sx q[3];
rz(-1.9105063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57758254) q[2];
sx q[2];
rz(-1.492604) q[2];
sx q[2];
rz(-2.1419443) q[2];
rz(-0.10351652) q[3];
sx q[3];
rz(-0.13705702) q[3];
sx q[3];
rz(1.1243533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12061159) q[0];
sx q[0];
rz(-2.3636901) q[0];
sx q[0];
rz(2.4840684) q[0];
rz(-2.855037) q[1];
sx q[1];
rz(-2.2189238) q[1];
sx q[1];
rz(0.34067571) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7332471) q[0];
sx q[0];
rz(-2.354963) q[0];
sx q[0];
rz(-2.0287081) q[0];
rz(-1.2286439) q[2];
sx q[2];
rz(-2.5492382) q[2];
sx q[2];
rz(-2.5387788) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9790736) q[1];
sx q[1];
rz(-2.5538429) q[1];
sx q[1];
rz(-1.2390922) q[1];
x q[2];
rz(-0.70298985) q[3];
sx q[3];
rz(-2.8841416) q[3];
sx q[3];
rz(1.5233745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75491536) q[2];
sx q[2];
rz(-1.9415386) q[2];
sx q[2];
rz(-2.6436451) q[2];
rz(2.8399816) q[3];
sx q[3];
rz(-0.40062723) q[3];
sx q[3];
rz(-2.6224459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-2.4122445) q[0];
sx q[0];
rz(-0.059818581) q[0];
sx q[0];
rz(2.8630032) q[0];
rz(-2.5623698) q[1];
sx q[1];
rz(-2.2021553) q[1];
sx q[1];
rz(0.07671193) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68966507) q[0];
sx q[0];
rz(-1.2566914) q[0];
sx q[0];
rz(0.12977022) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3086583) q[2];
sx q[2];
rz(-0.58912504) q[2];
sx q[2];
rz(-0.72074705) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50347603) q[1];
sx q[1];
rz(-1.5498811) q[1];
sx q[1];
rz(-2.7760844) q[1];
x q[2];
rz(0.46080188) q[3];
sx q[3];
rz(-0.19690234) q[3];
sx q[3];
rz(2.868696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77438337) q[2];
sx q[2];
rz(-0.70720208) q[2];
sx q[2];
rz(-2.6860766) q[2];
rz(0.43141836) q[3];
sx q[3];
rz(-3.016267) q[3];
sx q[3];
rz(-3.07807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9615622) q[0];
sx q[0];
rz(-1.7407692) q[0];
sx q[0];
rz(-2.2041007) q[0];
rz(2.5554399) q[1];
sx q[1];
rz(-1.6393607) q[1];
sx q[1];
rz(1.6929109) q[1];
rz(2.0942681) q[2];
sx q[2];
rz(-2.6101255) q[2];
sx q[2];
rz(2.5892467) q[2];
rz(1.9361817) q[3];
sx q[3];
rz(-1.1689417) q[3];
sx q[3];
rz(0.4369215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];