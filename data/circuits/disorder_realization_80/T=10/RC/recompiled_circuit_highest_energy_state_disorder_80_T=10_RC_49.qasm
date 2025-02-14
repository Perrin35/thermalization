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
rz(-1.9266204) q[0];
sx q[0];
rz(-1.5705234) q[0];
sx q[0];
rz(-2.7745752) q[0];
rz(-0.031232746) q[1];
sx q[1];
rz(-0.90489689) q[1];
sx q[1];
rz(2.4651405) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5117422) q[0];
sx q[0];
rz(-0.29030784) q[0];
sx q[0];
rz(-2.8111893) q[0];
rz(-0.88432856) q[2];
sx q[2];
rz(-1.2791233) q[2];
sx q[2];
rz(2.2560788) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8707968) q[1];
sx q[1];
rz(-0.50179401) q[1];
sx q[1];
rz(1.938799) q[1];
rz(0.2113249) q[3];
sx q[3];
rz(-0.67398807) q[3];
sx q[3];
rz(0.061376288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61554945) q[2];
sx q[2];
rz(-0.56576663) q[2];
sx q[2];
rz(0.88741285) q[2];
rz(-0.081485661) q[3];
sx q[3];
rz(-1.2469651) q[3];
sx q[3];
rz(0.87765774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74251974) q[0];
sx q[0];
rz(-2.6953473) q[0];
sx q[0];
rz(1.9148069) q[0];
rz(0.98608214) q[1];
sx q[1];
rz(-1.6796651) q[1];
sx q[1];
rz(-0.90423924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609294) q[0];
sx q[0];
rz(-0.21337803) q[0];
sx q[0];
rz(-0.55960525) q[0];
rz(-0.23444011) q[2];
sx q[2];
rz(-2.9801705) q[2];
sx q[2];
rz(1.8983253) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0457669) q[1];
sx q[1];
rz(-0.44039044) q[1];
sx q[1];
rz(-2.7604123) q[1];
x q[2];
rz(-2.9504602) q[3];
sx q[3];
rz(-0.88273747) q[3];
sx q[3];
rz(-0.21746527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6282689) q[2];
sx q[2];
rz(-0.29418918) q[2];
sx q[2];
rz(-0.42523709) q[2];
rz(2.6325295) q[3];
sx q[3];
rz(-2.013701) q[3];
sx q[3];
rz(-1.7019255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8486479) q[0];
sx q[0];
rz(-0.078462891) q[0];
sx q[0];
rz(1.5899832) q[0];
rz(-1.4750922) q[1];
sx q[1];
rz(-2.0547129) q[1];
sx q[1];
rz(-2.1045254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35557105) q[0];
sx q[0];
rz(-1.6986587) q[0];
sx q[0];
rz(2.1855445) q[0];
rz(2.6042348) q[2];
sx q[2];
rz(-1.4721057) q[2];
sx q[2];
rz(1.3787998) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4605352) q[1];
sx q[1];
rz(-1.4303416) q[1];
sx q[1];
rz(1.3696425) q[1];
x q[2];
rz(2.6821503) q[3];
sx q[3];
rz(-1.7342576) q[3];
sx q[3];
rz(-1.3889988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40917778) q[2];
sx q[2];
rz(-0.38984177) q[2];
sx q[2];
rz(1.2624435) q[2];
rz(3.0564195) q[3];
sx q[3];
rz(-2.6933935) q[3];
sx q[3];
rz(-2.0195154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.9389116) q[0];
sx q[0];
rz(-2.1124463) q[0];
sx q[0];
rz(-2.9503248) q[0];
rz(2.250504) q[1];
sx q[1];
rz(-2.8137408) q[1];
sx q[1];
rz(-2.0909615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4114496) q[0];
sx q[0];
rz(-0.77783424) q[0];
sx q[0];
rz(1.2744689) q[0];
x q[1];
rz(0.3140788) q[2];
sx q[2];
rz(-0.88523141) q[2];
sx q[2];
rz(0.92957815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5393319) q[1];
sx q[1];
rz(-2.6289399) q[1];
sx q[1];
rz(-0.67474483) q[1];
rz(-pi) q[2];
rz(2.1414143) q[3];
sx q[3];
rz(-1.4550536) q[3];
sx q[3];
rz(1.9952578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27590671) q[2];
sx q[2];
rz(-2.9624717) q[2];
sx q[2];
rz(0.56037819) q[2];
rz(-0.066702453) q[3];
sx q[3];
rz(-1.3542465) q[3];
sx q[3];
rz(-1.0022479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0340843) q[0];
sx q[0];
rz(-1.1488687) q[0];
sx q[0];
rz(2.4565571) q[0];
rz(0.3512474) q[1];
sx q[1];
rz(-0.63397399) q[1];
sx q[1];
rz(1.8018855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9046341) q[0];
sx q[0];
rz(-1.5642421) q[0];
sx q[0];
rz(2.652524) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4701202) q[2];
sx q[2];
rz(-1.9118309) q[2];
sx q[2];
rz(1.070553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1692002) q[1];
sx q[1];
rz(-0.16443096) q[1];
sx q[1];
rz(-1.3878257) q[1];
x q[2];
rz(-2.005607) q[3];
sx q[3];
rz(-1.1204733) q[3];
sx q[3];
rz(-0.45057566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54905218) q[2];
sx q[2];
rz(-1.668674) q[2];
sx q[2];
rz(0.215691) q[2];
rz(-0.14821626) q[3];
sx q[3];
rz(-0.18311466) q[3];
sx q[3];
rz(2.3858002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53629476) q[0];
sx q[0];
rz(-1.7114102) q[0];
sx q[0];
rz(0.17912616) q[0];
rz(-0.99963775) q[1];
sx q[1];
rz(-0.78796402) q[1];
sx q[1];
rz(-2.9771908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59578204) q[0];
sx q[0];
rz(-2.3124957) q[0];
sx q[0];
rz(1.2973644) q[0];
rz(-pi) q[1];
rz(1.4363507) q[2];
sx q[2];
rz(-2.0847581) q[2];
sx q[2];
rz(2.9902637) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8171606) q[1];
sx q[1];
rz(-0.91113976) q[1];
sx q[1];
rz(-1.5615669) q[1];
x q[2];
rz(0.54349676) q[3];
sx q[3];
rz(-2.7201289) q[3];
sx q[3];
rz(-1.8068562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0187692) q[2];
sx q[2];
rz(-1.5888701) q[2];
sx q[2];
rz(-1.4336047) q[2];
rz(-0.99099365) q[3];
sx q[3];
rz(-1.7510479) q[3];
sx q[3];
rz(-1.8979134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7551512) q[0];
sx q[0];
rz(-0.069510892) q[0];
sx q[0];
rz(-2.884927) q[0];
rz(0.060404213) q[1];
sx q[1];
rz(-0.81788617) q[1];
sx q[1];
rz(2.7108257) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2592436) q[0];
sx q[0];
rz(-1.8494864) q[0];
sx q[0];
rz(2.1889127) q[0];
rz(-pi) q[1];
rz(0.97051986) q[2];
sx q[2];
rz(-0.97504967) q[2];
sx q[2];
rz(-2.3733394) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.52331978) q[1];
sx q[1];
rz(-0.91611169) q[1];
sx q[1];
rz(1.0275082) q[1];
rz(-1.7458785) q[3];
sx q[3];
rz(-2.4444067) q[3];
sx q[3];
rz(0.28466636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7595235) q[2];
sx q[2];
rz(-0.52983317) q[2];
sx q[2];
rz(-2.1257373) q[2];
rz(1.2782512) q[3];
sx q[3];
rz(-1.9081554) q[3];
sx q[3];
rz(-0.61346936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.627219) q[0];
sx q[0];
rz(-0.52424279) q[0];
sx q[0];
rz(2.4216477) q[0];
rz(-0.68079692) q[1];
sx q[1];
rz(-2.7578208) q[1];
sx q[1];
rz(2.0489571) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93406308) q[0];
sx q[0];
rz(-1.0813366) q[0];
sx q[0];
rz(-1.5402769) q[0];
rz(-0.60694867) q[2];
sx q[2];
rz(-1.8766093) q[2];
sx q[2];
rz(2.1987777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5515308) q[1];
sx q[1];
rz(-2.2023337) q[1];
sx q[1];
rz(-0.88330357) q[1];
rz(-pi) q[2];
rz(1.6933788) q[3];
sx q[3];
rz(-1.8245909) q[3];
sx q[3];
rz(-3.0021961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1277658) q[2];
sx q[2];
rz(-0.76467815) q[2];
sx q[2];
rz(0.71316767) q[2];
rz(-1.7662883) q[3];
sx q[3];
rz(-1.2262552) q[3];
sx q[3];
rz(1.3885952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16667287) q[0];
sx q[0];
rz(-2.9627934) q[0];
sx q[0];
rz(1.4270225) q[0];
rz(2.7339281) q[1];
sx q[1];
rz(-1.0382321) q[1];
sx q[1];
rz(-2.8020249) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7106066) q[0];
sx q[0];
rz(-0.8961959) q[0];
sx q[0];
rz(2.8047724) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10972523) q[2];
sx q[2];
rz(-0.418677) q[2];
sx q[2];
rz(-0.81594407) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.133278) q[1];
sx q[1];
rz(-2.5247658) q[1];
sx q[1];
rz(2.1471669) q[1];
rz(1.4601213) q[3];
sx q[3];
rz(-0.95451285) q[3];
sx q[3];
rz(-1.2736774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0833544) q[2];
sx q[2];
rz(-2.3149172) q[2];
sx q[2];
rz(0.44462407) q[2];
rz(2.5104163) q[3];
sx q[3];
rz(-1.8034214) q[3];
sx q[3];
rz(-1.1609424) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4525962) q[0];
sx q[0];
rz(-1.3225553) q[0];
sx q[0];
rz(-3.1153862) q[0];
rz(1.7474489) q[1];
sx q[1];
rz(-1.1528287) q[1];
sx q[1];
rz(1.1411512) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0923528) q[0];
sx q[0];
rz(-1.0136407) q[0];
sx q[0];
rz(-2.0505285) q[0];
rz(-pi) q[1];
rz(-0.097392453) q[2];
sx q[2];
rz(-0.44008128) q[2];
sx q[2];
rz(-0.68005622) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2403912) q[1];
sx q[1];
rz(-1.9878042) q[1];
sx q[1];
rz(0.97948018) q[1];
rz(-1.4196542) q[3];
sx q[3];
rz(-1.0060174) q[3];
sx q[3];
rz(-1.4709815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0673087) q[2];
sx q[2];
rz(-1.1003541) q[2];
sx q[2];
rz(1.1624153) q[2];
rz(0.31636604) q[3];
sx q[3];
rz(-1.1917944) q[3];
sx q[3];
rz(1.4828064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3027073) q[0];
sx q[0];
rz(-3.0743297) q[0];
sx q[0];
rz(-0.59826941) q[0];
rz(-3.132598) q[1];
sx q[1];
rz(-0.79575494) q[1];
sx q[1];
rz(1.5473821) q[1];
rz(2.9526906) q[2];
sx q[2];
rz(-2.4707265) q[2];
sx q[2];
rz(2.1553253) q[2];
rz(2.347765) q[3];
sx q[3];
rz(-1.6110653) q[3];
sx q[3];
rz(1.8992784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
