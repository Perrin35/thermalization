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
rz(0.23586805) q[0];
sx q[0];
rz(1.9960825) q[0];
sx q[0];
rz(9.3762015) q[0];
rz(-1.6254758) q[1];
sx q[1];
rz(-1.0928417) q[1];
sx q[1];
rz(0.71576524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19312035) q[0];
sx q[0];
rz(-0.17331757) q[0];
sx q[0];
rz(1.7422218) q[0];
rz(-2.2015436) q[2];
sx q[2];
rz(-0.34206451) q[2];
sx q[2];
rz(-2.2243481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6246031) q[1];
sx q[1];
rz(-1.9698242) q[1];
sx q[1];
rz(-2.3844403) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6244087) q[3];
sx q[3];
rz(-1.2150914) q[3];
sx q[3];
rz(-2.9773447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6526661) q[2];
sx q[2];
rz(-0.55997866) q[2];
sx q[2];
rz(2.0860591) q[2];
rz(-2.7355898) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(-2.3199911) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9838298) q[0];
sx q[0];
rz(-2.1021748) q[0];
sx q[0];
rz(-0.55617547) q[0];
rz(-1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(0.083273085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52516925) q[0];
sx q[0];
rz(-0.94248191) q[0];
sx q[0];
rz(1.9740021) q[0];
rz(-pi) q[1];
x q[1];
rz(0.95319548) q[2];
sx q[2];
rz(-1.0824656) q[2];
sx q[2];
rz(-2.7361779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4652768) q[1];
sx q[1];
rz(-2.0124467) q[1];
sx q[1];
rz(1.8144813) q[1];
x q[2];
rz(0.86822416) q[3];
sx q[3];
rz(-2.4463042) q[3];
sx q[3];
rz(-0.27588683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1978153) q[2];
sx q[2];
rz(-1.3894206) q[2];
sx q[2];
rz(-2.6488292) q[2];
rz(2.1150151) q[3];
sx q[3];
rz(-2.8387098) q[3];
sx q[3];
rz(1.7290285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3781085) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(0.43877959) q[0];
rz(1.4525061) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(0.39295331) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9742115) q[0];
sx q[0];
rz(-2.0618453) q[0];
sx q[0];
rz(-0.031662861) q[0];
rz(-pi) q[1];
rz(0.032756373) q[2];
sx q[2];
rz(-1.0426781) q[2];
sx q[2];
rz(2.817085) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90653268) q[1];
sx q[1];
rz(-1.9882012) q[1];
sx q[1];
rz(0.077666186) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64928253) q[3];
sx q[3];
rz(-1.6341795) q[3];
sx q[3];
rz(-1.6398304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.68985933) q[2];
sx q[2];
rz(-0.84443513) q[2];
sx q[2];
rz(2.8680657) q[2];
rz(-0.61255974) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(-2.2851022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.3278811) q[0];
sx q[0];
rz(-2.5484945) q[0];
sx q[0];
rz(2.3408422) q[0];
rz(0.71209359) q[1];
sx q[1];
rz(-1.1199896) q[1];
sx q[1];
rz(-0.48316479) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8034604) q[0];
sx q[0];
rz(-0.82809859) q[0];
sx q[0];
rz(1.5699734) q[0];
rz(-pi) q[1];
rz(2.0676548) q[2];
sx q[2];
rz(-2.0570847) q[2];
sx q[2];
rz(0.72975791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2702647) q[1];
sx q[1];
rz(-0.92496269) q[1];
sx q[1];
rz(-1.6069018) q[1];
x q[2];
rz(1.9075059) q[3];
sx q[3];
rz(-2.2318865) q[3];
sx q[3];
rz(-2.8291836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.157865) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(1.2845117) q[2];
rz(2.6537248) q[3];
sx q[3];
rz(-1.4372466) q[3];
sx q[3];
rz(1.6871066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3294285) q[0];
sx q[0];
rz(-0.85947961) q[0];
sx q[0];
rz(2.676945) q[0];
rz(2.1931785) q[1];
sx q[1];
rz(-2.3034838) q[1];
sx q[1];
rz(-1.4533739) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8226979) q[0];
sx q[0];
rz(-1.5359211) q[0];
sx q[0];
rz(-0.025392763) q[0];
rz(-pi) q[1];
rz(0.27957966) q[2];
sx q[2];
rz(-2.381392) q[2];
sx q[2];
rz(-0.81499824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0665474) q[1];
sx q[1];
rz(-0.61457026) q[1];
sx q[1];
rz(0.11654186) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0539758) q[3];
sx q[3];
rz(-0.91435104) q[3];
sx q[3];
rz(2.2051728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98443085) q[2];
sx q[2];
rz(-1.0474019) q[2];
sx q[2];
rz(2.7471527) q[2];
rz(-1.1557584) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(-1.8331147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406469) q[0];
sx q[0];
rz(-2.7257958) q[0];
sx q[0];
rz(-1.0619324) q[0];
rz(0.95147079) q[1];
sx q[1];
rz(-0.87363344) q[1];
sx q[1];
rz(1.5207312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5439107) q[0];
sx q[0];
rz(-1.2939699) q[0];
sx q[0];
rz(-0.055276543) q[0];
rz(0.64487793) q[2];
sx q[2];
rz(-2.0487924) q[2];
sx q[2];
rz(2.1410112) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86506063) q[1];
sx q[1];
rz(-2.3569657) q[1];
sx q[1];
rz(-2.9097981) q[1];
rz(-1.3502832) q[3];
sx q[3];
rz(-1.1854544) q[3];
sx q[3];
rz(2.6366608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1157397) q[2];
sx q[2];
rz(-1.9151177) q[2];
sx q[2];
rz(-0.24570492) q[2];
rz(-1.4541516) q[3];
sx q[3];
rz(-1.5116296) q[3];
sx q[3];
rz(-2.3112442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4858953) q[0];
sx q[0];
rz(-0.73807722) q[0];
sx q[0];
rz(0.6947211) q[0];
rz(0.10063902) q[1];
sx q[1];
rz(-1.0357608) q[1];
sx q[1];
rz(2.2763841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14510205) q[0];
sx q[0];
rz(-1.9887513) q[0];
sx q[0];
rz(-1.6161902) q[0];
x q[1];
rz(0.91036441) q[2];
sx q[2];
rz(-1.5857134) q[2];
sx q[2];
rz(-2.2098324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2336967) q[1];
sx q[1];
rz(-1.7204667) q[1];
sx q[1];
rz(2.187633) q[1];
x q[2];
rz(-0.65552222) q[3];
sx q[3];
rz(-2.351774) q[3];
sx q[3];
rz(1.5721377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18409099) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(0.91020477) q[2];
rz(1.4522067) q[3];
sx q[3];
rz(-1.2319177) q[3];
sx q[3];
rz(0.46428251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8978187) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(-2.5183103) q[0];
rz(-0.15577623) q[1];
sx q[1];
rz(-0.77729762) q[1];
sx q[1];
rz(-2.8782841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6197841) q[0];
sx q[0];
rz(-1.5814879) q[0];
sx q[0];
rz(-1.7858206) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.510235) q[2];
sx q[2];
rz(-1.6471631) q[2];
sx q[2];
rz(-1.5983943) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38674445) q[1];
sx q[1];
rz(-1.6476908) q[1];
sx q[1];
rz(-0.57707483) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78531475) q[3];
sx q[3];
rz(-1.974733) q[3];
sx q[3];
rz(0.31111003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2498563) q[2];
sx q[2];
rz(-2.8920434) q[2];
sx q[2];
rz(-2.2535394) q[2];
rz(0.92787162) q[3];
sx q[3];
rz(-1.1279305) q[3];
sx q[3];
rz(-0.97308648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.71858281) q[0];
sx q[0];
rz(-0.50964481) q[0];
sx q[0];
rz(-0.50073671) q[0];
rz(2.0210733) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(0.80148554) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1130509) q[0];
sx q[0];
rz(-3.0705912) q[0];
sx q[0];
rz(3.1138022) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4013585) q[2];
sx q[2];
rz(-2.4317305) q[2];
sx q[2];
rz(2.9346093) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.34757839) q[1];
sx q[1];
rz(-1.9982575) q[1];
sx q[1];
rz(-2.4540837) q[1];
x q[2];
rz(-1.3128223) q[3];
sx q[3];
rz(-1.7593062) q[3];
sx q[3];
rz(2.9303355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72006172) q[2];
sx q[2];
rz(-1.820182) q[2];
sx q[2];
rz(-0.23019543) q[2];
rz(3.1318393) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(2.9856288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4474051) q[0];
sx q[0];
rz(-1.7740086) q[0];
sx q[0];
rz(2.3719924) q[0];
rz(2.1790478) q[1];
sx q[1];
rz(-1.3178408) q[1];
sx q[1];
rz(2.1544971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55961032) q[0];
sx q[0];
rz(-1.6405237) q[0];
sx q[0];
rz(-2.0980673) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9883698) q[2];
sx q[2];
rz(-1.8109057) q[2];
sx q[2];
rz(-2.3599412) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.35008123) q[1];
sx q[1];
rz(-1.529585) q[1];
sx q[1];
rz(-0.033153127) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0760097) q[3];
sx q[3];
rz(-0.85799137) q[3];
sx q[3];
rz(-3.0312169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4276245) q[2];
sx q[2];
rz(-0.93728137) q[2];
sx q[2];
rz(0.06812185) q[2];
rz(2.6918329) q[3];
sx q[3];
rz(-2.7263548) q[3];
sx q[3];
rz(-1.7033887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9925256) q[0];
sx q[0];
rz(-1.608792) q[0];
sx q[0];
rz(1.7444862) q[0];
rz(-0.29009157) q[1];
sx q[1];
rz(-2.7128704) q[1];
sx q[1];
rz(-2.6147978) q[1];
rz(-1.9214966) q[2];
sx q[2];
rz(-2.0439707) q[2];
sx q[2];
rz(1.1743152) q[2];
rz(-0.044779891) q[3];
sx q[3];
rz(-0.3988409) q[3];
sx q[3];
rz(-2.1259784) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
