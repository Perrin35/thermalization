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
rz(1.3069557) q[0];
sx q[0];
rz(4.0210273) q[0];
sx q[0];
rz(9.9183912) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077162077) q[0];
sx q[0];
rz(-1.0748942) q[0];
sx q[0];
rz(0.84159633) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2463208) q[2];
sx q[2];
rz(-0.68973422) q[2];
sx q[2];
rz(-1.5472428) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2027758) q[1];
sx q[1];
rz(-0.73828739) q[1];
sx q[1];
rz(2.8144005) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93500626) q[3];
sx q[3];
rz(-1.269644) q[3];
sx q[3];
rz(-3.0060535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1349858) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(-1.1084278) q[2];
rz(-1.8290352) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(2.826622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38158622) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(0.015856892) q[0];
rz(2.1511757) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(-0.86318618) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9445489) q[0];
sx q[0];
rz(-1.6672575) q[0];
sx q[0];
rz(2.405557) q[0];
rz(1.241774) q[2];
sx q[2];
rz(-1.8890427) q[2];
sx q[2];
rz(1.9538188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3962209) q[1];
sx q[1];
rz(-0.52597755) q[1];
sx q[1];
rz(-0.35244103) q[1];
rz(0.84896481) q[3];
sx q[3];
rz(-1.6833653) q[3];
sx q[3];
rz(1.6602885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4414759) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(2.6727943) q[2];
rz(0.68764728) q[3];
sx q[3];
rz(-1.5117437) q[3];
sx q[3];
rz(-0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1402635) q[0];
sx q[0];
rz(-0.92324531) q[0];
sx q[0];
rz(2.4053251) q[0];
rz(-2.7983792) q[1];
sx q[1];
rz(-0.88756573) q[1];
sx q[1];
rz(-0.18377486) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2375419) q[0];
sx q[0];
rz(-1.5045368) q[0];
sx q[0];
rz(2.3973293) q[0];
rz(-pi) q[1];
rz(1.0509346) q[2];
sx q[2];
rz(-0.94280548) q[2];
sx q[2];
rz(-1.3818714) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0597093) q[1];
sx q[1];
rz(-0.78468152) q[1];
sx q[1];
rz(0.42894026) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1322429) q[3];
sx q[3];
rz(-0.61547503) q[3];
sx q[3];
rz(1.4765679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3595714) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(0.51592958) q[2];
rz(-2.0333717) q[3];
sx q[3];
rz(-2.5907232) q[3];
sx q[3];
rz(1.1512604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.16294031) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(-2.6287855) q[0];
rz(2.6817952) q[1];
sx q[1];
rz(-1.5492946) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13614722) q[0];
sx q[0];
rz(-2.9678565) q[0];
sx q[0];
rz(2.4197015) q[0];
x q[1];
rz(0.64900543) q[2];
sx q[2];
rz(-2.634806) q[2];
sx q[2];
rz(-0.014647131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11359195) q[1];
sx q[1];
rz(-1.290844) q[1];
sx q[1];
rz(1.2173613) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34448907) q[3];
sx q[3];
rz(-1.7055344) q[3];
sx q[3];
rz(-1.3400934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1226471) q[2];
sx q[2];
rz(-2.985869) q[2];
sx q[2];
rz(-0.32749185) q[2];
rz(-2.859595) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(1.5871083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.17664385) q[0];
sx q[0];
rz(-0.38589859) q[0];
sx q[0];
rz(-2.1642245) q[0];
rz(0.36879677) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(1.4126973) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4810886) q[0];
sx q[0];
rz(-1.7605283) q[0];
sx q[0];
rz(1.0718179) q[0];
rz(-pi) q[1];
rz(-1.0540038) q[2];
sx q[2];
rz(-1.5828504) q[2];
sx q[2];
rz(-0.64753676) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.818644) q[1];
sx q[1];
rz(-1.6020157) q[1];
sx q[1];
rz(-1.853456) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5630412) q[3];
sx q[3];
rz(-1.2935841) q[3];
sx q[3];
rz(1.2523684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29350027) q[2];
sx q[2];
rz(-0.81587452) q[2];
sx q[2];
rz(2.9917742) q[2];
rz(-0.39854974) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.77001101) q[0];
sx q[0];
rz(-0.12938975) q[0];
sx q[0];
rz(0.55321252) q[0];
rz(2.5999056) q[1];
sx q[1];
rz(-2.0952416) q[1];
sx q[1];
rz(2.6936626) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49419379) q[0];
sx q[0];
rz(-1.9564391) q[0];
sx q[0];
rz(1.4146039) q[0];
rz(-0.80336415) q[2];
sx q[2];
rz(-2.5584872) q[2];
sx q[2];
rz(-1.2606114) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74400157) q[1];
sx q[1];
rz(-2.7251456) q[1];
sx q[1];
rz(-1.9201502) q[1];
rz(0.16334857) q[3];
sx q[3];
rz(-2.8046097) q[3];
sx q[3];
rz(2.8548422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64122671) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(-0.79916239) q[2];
rz(2.7568119) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(2.1414115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36050972) q[0];
sx q[0];
rz(-2.0979083) q[0];
sx q[0];
rz(2.5497896) q[0];
rz(-0.38161713) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(2.9891678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7397933) q[0];
sx q[0];
rz(-1.2619979) q[0];
sx q[0];
rz(-2.1112421) q[0];
x q[1];
rz(0.47758684) q[2];
sx q[2];
rz(-1.7268333) q[2];
sx q[2];
rz(2.0618771) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85411863) q[1];
sx q[1];
rz(-1.9175314) q[1];
sx q[1];
rz(-1.5513913) q[1];
rz(0.8922718) q[3];
sx q[3];
rz(-1.0425869) q[3];
sx q[3];
rz(1.5367374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3064208) q[2];
sx q[2];
rz(-0.99205899) q[2];
sx q[2];
rz(1.6888118) q[2];
rz(0.49550223) q[3];
sx q[3];
rz(-0.82363868) q[3];
sx q[3];
rz(-0.56408322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.38967663) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(-0.16146846) q[0];
rz(0.031919315) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(1.8643103) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93715765) q[0];
sx q[0];
rz(-1.8291408) q[0];
sx q[0];
rz(-2.4169371) q[0];
x q[1];
rz(1.0890245) q[2];
sx q[2];
rz(-2.6725997) q[2];
sx q[2];
rz(-0.8587786) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2252075) q[1];
sx q[1];
rz(-0.3233094) q[1];
sx q[1];
rz(-2.7955452) q[1];
rz(1.444449) q[3];
sx q[3];
rz(-2.845361) q[3];
sx q[3];
rz(-2.031523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45234597) q[2];
sx q[2];
rz(-2.2301058) q[2];
sx q[2];
rz(0.39608836) q[2];
rz(2.7416157) q[3];
sx q[3];
rz(-2.5466099) q[3];
sx q[3];
rz(0.8766492) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24092995) q[0];
sx q[0];
rz(-3.1234968) q[0];
sx q[0];
rz(2.990429) q[0];
rz(0.97686544) q[1];
sx q[1];
rz(-0.38115373) q[1];
sx q[1];
rz(-0.74473286) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3522986) q[0];
sx q[0];
rz(-1.2937102) q[0];
sx q[0];
rz(2.9267071) q[0];
rz(-pi) q[1];
rz(-1.1579723) q[2];
sx q[2];
rz(-1.7016439) q[2];
sx q[2];
rz(1.7826338) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1972629) q[1];
sx q[1];
rz(-1.5722317) q[1];
sx q[1];
rz(-2.4932842) q[1];
rz(-0.35208296) q[3];
sx q[3];
rz(-1.3439281) q[3];
sx q[3];
rz(-2.9405103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.9025977) q[2];
sx q[2];
rz(-2.1935479) q[2];
rz(-2.0840123) q[3];
sx q[3];
rz(-1.4761997) q[3];
sx q[3];
rz(2.0675366) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0483765) q[0];
sx q[0];
rz(-2.8792448) q[0];
sx q[0];
rz(-0.23396215) q[0];
rz(-1.9562862) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(-1.2497466) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.098506) q[0];
sx q[0];
rz(-2.1224501) q[0];
sx q[0];
rz(-1.2292525) q[0];
rz(-pi) q[1];
rz(0.62057497) q[2];
sx q[2];
rz(-1.7399551) q[2];
sx q[2];
rz(0.38693869) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0253623) q[1];
sx q[1];
rz(-2.2409908) q[1];
sx q[1];
rz(-0.7109364) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3951282) q[3];
sx q[3];
rz(-2.2445956) q[3];
sx q[3];
rz(1.8927946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.027111) q[2];
sx q[2];
rz(-1.1610843) q[2];
sx q[2];
rz(-1.7005881) q[2];
rz(-3.0103185) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-0.24391267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048653614) q[0];
sx q[0];
rz(-2.1760512) q[0];
sx q[0];
rz(-2.0078134) q[0];
rz(2.1009905) q[1];
sx q[1];
rz(-2.2941209) q[1];
sx q[1];
rz(2.6170731) q[1];
rz(-0.94338633) q[2];
sx q[2];
rz(-1.3209983) q[2];
sx q[2];
rz(-0.41994195) q[2];
rz(0.20411251) q[3];
sx q[3];
rz(-2.4983046) q[3];
sx q[3];
rz(0.78346969) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
