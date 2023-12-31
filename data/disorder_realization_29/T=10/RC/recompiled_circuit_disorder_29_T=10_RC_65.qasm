OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(4.6255339) q[0];
sx q[0];
rz(12.892527) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9261949) q[0];
sx q[0];
rz(-2.8370259) q[0];
sx q[0];
rz(-0.79415168) q[0];
rz(-pi) q[1];
rz(-0.60249451) q[2];
sx q[2];
rz(-1.7817111) q[2];
sx q[2];
rz(-2.9162625) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1118288) q[1];
sx q[1];
rz(-1.811869) q[1];
sx q[1];
rz(-1.3328711) q[1];
rz(-1.5845675) q[3];
sx q[3];
rz(-0.78409401) q[3];
sx q[3];
rz(-0.23628326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-0.97936169) q[2];
sx q[2];
rz(0.88511434) q[2];
rz(-0.72201133) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(0.66501578) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(-2.2639993) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.005640429) q[0];
sx q[0];
rz(-0.71605039) q[0];
sx q[0];
rz(0.27766772) q[0];
x q[1];
rz(-1.7052824) q[2];
sx q[2];
rz(-1.209139) q[2];
sx q[2];
rz(-0.59097564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9237325) q[1];
sx q[1];
rz(-1.4638312) q[1];
sx q[1];
rz(0.033543368) q[1];
rz(-pi) q[2];
rz(1.9545994) q[3];
sx q[3];
rz(-1.2591397) q[3];
sx q[3];
rz(2.4812738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39891222) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(1.1304643) q[2];
rz(-1.2997262) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14625064) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(2.8515942) q[0];
rz(2.4747804) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(-3.0677632) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2732669) q[0];
sx q[0];
rz(-0.10859117) q[0];
sx q[0];
rz(-0.6032087) q[0];
rz(-pi) q[1];
rz(-1.2286085) q[2];
sx q[2];
rz(-1.8456736) q[2];
sx q[2];
rz(-2.9800422) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2245582) q[1];
sx q[1];
rz(-0.18637603) q[1];
sx q[1];
rz(-1.264155) q[1];
rz(0.31432398) q[3];
sx q[3];
rz(-2.7735908) q[3];
sx q[3];
rz(0.60955334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6510216) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(-0.64669615) q[2];
rz(-1.1086639) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(-1.1886764) q[0];
rz(-2.1229318) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(-1.7046938) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80285145) q[0];
sx q[0];
rz(-2.2995298) q[0];
sx q[0];
rz(-1.4472423) q[0];
rz(0.75886274) q[2];
sx q[2];
rz(-2.7042411) q[2];
sx q[2];
rz(-0.26668374) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77364591) q[1];
sx q[1];
rz(-1.2229496) q[1];
sx q[1];
rz(2.0639973) q[1];
rz(-1.4401682) q[3];
sx q[3];
rz(-1.4454953) q[3];
sx q[3];
rz(0.13392042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58549762) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(1.1222703) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68200237) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(2.7686152) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.3164828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3554768) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(1.0771846) q[0];
rz(-pi) q[1];
rz(0.76422845) q[2];
sx q[2];
rz(-0.35596213) q[2];
sx q[2];
rz(0.9133577) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.398048) q[1];
sx q[1];
rz(-3.0882356) q[1];
sx q[1];
rz(1.8023026) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76977323) q[3];
sx q[3];
rz(-2.4379745) q[3];
sx q[3];
rz(2.9911656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0405154) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(-0.80580795) q[2];
rz(-2.6082883) q[3];
sx q[3];
rz(-1.1321944) q[3];
sx q[3];
rz(2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7508115) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(3.0506296) q[0];
rz(-2.2816351) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(1.3202753) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108171) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(-2.4369795) q[0];
x q[1];
rz(-0.27310246) q[2];
sx q[2];
rz(-0.58758508) q[2];
sx q[2];
rz(-2.7762129) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7652313) q[1];
sx q[1];
rz(-0.77502854) q[1];
sx q[1];
rz(2.2750957) q[1];
rz(-0.91464197) q[3];
sx q[3];
rz(-1.802889) q[3];
sx q[3];
rz(0.37999145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-2.5406204) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(-1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619974) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(0.55554187) q[0];
rz(3.1069966) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.3909891) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8081252) q[0];
sx q[0];
rz(-1.2398749) q[0];
sx q[0];
rz(2.5740741) q[0];
rz(-1.820307) q[2];
sx q[2];
rz(-1.9267285) q[2];
sx q[2];
rz(-0.78636679) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6185551) q[1];
sx q[1];
rz(-0.68261519) q[1];
sx q[1];
rz(1.7726266) q[1];
x q[2];
rz(1.4963385) q[3];
sx q[3];
rz(-0.7965318) q[3];
sx q[3];
rz(-1.2021241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3283219) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(-2.7461046) q[2];
rz(-1.288712) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(-2.18816) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(-1.4377726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2742845) q[0];
sx q[0];
rz(-2.0831997) q[0];
sx q[0];
rz(-2.9318277) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9713694) q[2];
sx q[2];
rz(-1.4366237) q[2];
sx q[2];
rz(-1.9113049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.21282141) q[1];
sx q[1];
rz(-2.050424) q[1];
sx q[1];
rz(-2.4368083) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1215454) q[3];
sx q[3];
rz(-1.123395) q[3];
sx q[3];
rz(2.945154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(-0.17871857) q[2];
rz(0.86137613) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20539595) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(-2.0478915) q[0];
rz(-2.4049092) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(2.0827983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1595488) q[0];
sx q[0];
rz(-1.0315572) q[0];
sx q[0];
rz(-0.3190785) q[0];
x q[1];
rz(-2.043622) q[2];
sx q[2];
rz(-1.9377922) q[2];
sx q[2];
rz(-0.25640139) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4150548) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(-1.1916257) q[1];
x q[2];
rz(-0.47302795) q[3];
sx q[3];
rz(-1.2086476) q[3];
sx q[3];
rz(0.87272296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(-1.4808562) q[2];
rz(-0.41040928) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695628) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(0.29125443) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(1.4321009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390022) q[0];
sx q[0];
rz(-2.0615091) q[0];
sx q[0];
rz(1.9309994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45551331) q[2];
sx q[2];
rz(-1.3580772) q[2];
sx q[2];
rz(-0.46609512) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1283975) q[1];
sx q[1];
rz(-0.55354512) q[1];
sx q[1];
rz(-3.058606) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0544153) q[3];
sx q[3];
rz(-2.0012337) q[3];
sx q[3];
rz(-3.103053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.46135205) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(-1.6067778) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(-0.46943584) q[3];
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
rz(-1.5466945) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(2.7728511) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(2.4070807) q[2];
sx q[2];
rz(-0.30167689) q[2];
sx q[2];
rz(-2.6405356) q[2];
rz(0.41714824) q[3];
sx q[3];
rz(-2.683831) q[3];
sx q[3];
rz(2.8575069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
