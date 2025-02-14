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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(-2.057743) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(2.6649063) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81945588) q[0];
sx q[0];
rz(-1.1753723) q[0];
sx q[0];
rz(0.28015341) q[0];
rz(-pi) q[1];
rz(-2.347685) q[2];
sx q[2];
rz(-0.61015266) q[2];
sx q[2];
rz(-2.8086503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.649316) q[1];
sx q[1];
rz(-1.5536123) q[1];
sx q[1];
rz(1.3130929) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45248078) q[3];
sx q[3];
rz(-2.6685331) q[3];
sx q[3];
rz(-0.61247952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1276663) q[2];
sx q[2];
rz(-1.650859) q[2];
sx q[2];
rz(2.9909383) q[2];
rz(1.539218) q[3];
sx q[3];
rz(-0.37709245) q[3];
sx q[3];
rz(-2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056331228) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(-0.37013176) q[0];
rz(-2.6689957) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(-2.1841689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9345624) q[0];
sx q[0];
rz(-2.0683204) q[0];
sx q[0];
rz(3.050404) q[0];
x q[1];
rz(2.022012) q[2];
sx q[2];
rz(-1.7282082) q[2];
sx q[2];
rz(-1.5296641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70499252) q[1];
sx q[1];
rz(-0.7448405) q[1];
sx q[1];
rz(3.0158494) q[1];
x q[2];
rz(2.7368746) q[3];
sx q[3];
rz(-2.4232695) q[3];
sx q[3];
rz(0.61678929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0051673278) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(-2.1180507) q[2];
rz(1.7510022) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(-1.8243779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35109529) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(-0.56280953) q[0];
rz(-1.5142745) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-3.0624342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3252727) q[0];
sx q[0];
rz(-2.5421341) q[0];
sx q[0];
rz(0.58167235) q[0];
rz(-pi) q[1];
rz(-0.076268836) q[2];
sx q[2];
rz(-0.4919211) q[2];
sx q[2];
rz(-2.9894376) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1428433) q[1];
sx q[1];
rz(-0.85341893) q[1];
sx q[1];
rz(-2.1446376) q[1];
x q[2];
rz(-2.4412254) q[3];
sx q[3];
rz(-1.0031962) q[3];
sx q[3];
rz(-1.8710473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2341653) q[2];
sx q[2];
rz(-1.558446) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(-2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(-0.98085421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6818162) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(-0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(-0.2535325) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0363844) q[0];
sx q[0];
rz(-2.0724793) q[0];
sx q[0];
rz(1.2060449) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5229613) q[2];
sx q[2];
rz(-0.19605532) q[2];
sx q[2];
rz(0.679099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56139466) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(-2.8785588) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0757789) q[3];
sx q[3];
rz(-1.2974713) q[3];
sx q[3];
rz(1.2697288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3229052) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(-1.6820071) q[2];
rz(2.5679576) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(-0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-2.0720338) q[0];
sx q[0];
rz(-1.3473508) q[0];
rz(-2.5509293) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(1.7151054) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1388228) q[0];
sx q[0];
rz(-1.9277687) q[0];
sx q[0];
rz(-2.7223552) q[0];
x q[1];
rz(1.1828184) q[2];
sx q[2];
rz(-0.96599865) q[2];
sx q[2];
rz(-1.6590349) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0795143) q[1];
sx q[1];
rz(-0.2582363) q[1];
sx q[1];
rz(-2.9596427) q[1];
rz(-pi) q[2];
rz(-1.6350045) q[3];
sx q[3];
rz(-1.2163463) q[3];
sx q[3];
rz(1.3469157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9023989) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(-1.1725461) q[2];
rz(0.086056195) q[3];
sx q[3];
rz(-2.1890169) q[3];
sx q[3];
rz(-0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758783) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(2.2789047) q[0];
rz(-0.24738303) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(-2.6695796) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0855477) q[0];
sx q[0];
rz(-1.9084832) q[0];
sx q[0];
rz(0.72283904) q[0];
x q[1];
rz(-0.15100592) q[2];
sx q[2];
rz(-1.0698294) q[2];
sx q[2];
rz(-0.46268625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0473049) q[1];
sx q[1];
rz(-1.3940812) q[1];
sx q[1];
rz(0.30728886) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94859) q[3];
sx q[3];
rz(-1.3974117) q[3];
sx q[3];
rz(2.9297048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2315959) q[2];
sx q[2];
rz(-2.1746077) q[2];
sx q[2];
rz(1.3989353) q[2];
rz(1.4553962) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(-2.5808891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-0.48651925) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(0.51862496) q[0];
rz(0.71867603) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(-1.8124883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3603921) q[0];
sx q[0];
rz(-2.5787163) q[0];
sx q[0];
rz(-2.2928574) q[0];
rz(-1.8480214) q[2];
sx q[2];
rz(-1.5223862) q[2];
sx q[2];
rz(-1.2836266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3490263) q[1];
sx q[1];
rz(-1.6036803) q[1];
sx q[1];
rz(-2.3057641) q[1];
rz(-3.1217478) q[3];
sx q[3];
rz(-0.9844508) q[3];
sx q[3];
rz(-1.3278409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79080498) q[2];
sx q[2];
rz(-2.4420276) q[2];
sx q[2];
rz(0.038185509) q[2];
rz(-2.3112678) q[3];
sx q[3];
rz(-1.7540951) q[3];
sx q[3];
rz(0.047688095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4627948) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(1.7713254) q[0];
rz(2.7554152) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(2.4716299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7012335) q[0];
sx q[0];
rz(-1.0770849) q[0];
sx q[0];
rz(-0.53892737) q[0];
rz(-pi) q[1];
rz(-1.6514844) q[2];
sx q[2];
rz(-2.2858857) q[2];
sx q[2];
rz(0.87835588) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8893056) q[1];
sx q[1];
rz(-0.60599594) q[1];
sx q[1];
rz(-1.3713981) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3588328) q[3];
sx q[3];
rz(-1.6073213) q[3];
sx q[3];
rz(-1.9741457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.92114821) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(-2.0849126) q[2];
rz(-1.7250666) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073762745) q[0];
sx q[0];
rz(-1.0404328) q[0];
sx q[0];
rz(-1.006806) q[0];
rz(2.6489068) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(0.83008343) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070979764) q[0];
sx q[0];
rz(-1.604053) q[0];
sx q[0];
rz(0.279056) q[0];
x q[1];
rz(-2.5428204) q[2];
sx q[2];
rz(-1.3922924) q[2];
sx q[2];
rz(2.0395961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0841024) q[1];
sx q[1];
rz(-1.7508645) q[1];
sx q[1];
rz(1.4157285) q[1];
rz(-pi) q[2];
rz(1.1277783) q[3];
sx q[3];
rz(-1.5592056) q[3];
sx q[3];
rz(-3.1168337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.81862187) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-1.0969561) q[2];
rz(-1.8638301) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(-0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85035664) q[0];
sx q[0];
rz(-1.1024029) q[0];
sx q[0];
rz(2.7665603) q[0];
rz(-3.0229783) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(-2.2172701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18436954) q[0];
sx q[0];
rz(-0.65378377) q[0];
sx q[0];
rz(-2.6517646) q[0];
rz(-2.5293328) q[2];
sx q[2];
rz(-0.77185248) q[2];
sx q[2];
rz(1.5764232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0339151) q[1];
sx q[1];
rz(-0.36201358) q[1];
sx q[1];
rz(3.0605143) q[1];
x q[2];
rz(-2.499324) q[3];
sx q[3];
rz(-2.6005201) q[3];
sx q[3];
rz(-2.9716345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2529926) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(-2.2484153) q[2];
rz(-1.1183974) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7669582) q[0];
sx q[0];
rz(-1.081291) q[0];
sx q[0];
rz(-1.3670856) q[0];
rz(-2.959666) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(-0.091808783) q[2];
sx q[2];
rz(-1.257007) q[2];
sx q[2];
rz(-1.9329482) q[2];
rz(1.160865) q[3];
sx q[3];
rz(-2.0515473) q[3];
sx q[3];
rz(-0.69507364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
