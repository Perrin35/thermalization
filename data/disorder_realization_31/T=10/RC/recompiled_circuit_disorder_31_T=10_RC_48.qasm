OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(2.5685413) q[0];
sx q[0];
rz(11.723784) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9487171) q[0];
sx q[0];
rz(-0.84515453) q[0];
sx q[0];
rz(-1.8564419) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5904434) q[2];
sx q[2];
rz(-0.95704776) q[2];
sx q[2];
rz(-0.27054271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2394489) q[1];
sx q[1];
rz(-1.4414756) q[1];
sx q[1];
rz(2.7620478) q[1];
rz(-pi) q[2];
rz(1.6928715) q[3];
sx q[3];
rz(-2.1616518) q[3];
sx q[3];
rz(-1.6216244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6136916) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(-2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801385) q[0];
sx q[0];
rz(-2.2606235) q[0];
sx q[0];
rz(2.8161312) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017529537) q[0];
sx q[0];
rz(-1.5520099) q[0];
sx q[0];
rz(1.5879052) q[0];
x q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.7413505) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7696015) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(3.0197057) q[1];
rz(-pi) q[2];
rz(1.8335908) q[3];
sx q[3];
rz(-1.7495219) q[3];
sx q[3];
rz(2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(-0.47131053) q[3];
sx q[3];
rz(-1.703197) q[3];
sx q[3];
rz(-0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8283591) q[0];
sx q[0];
rz(-1.4947083) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(-0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-2.0498958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1611623) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(0.79361332) q[0];
rz(2.2268779) q[2];
sx q[2];
rz(-1.2351742) q[2];
sx q[2];
rz(-0.4682954) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.97869067) q[1];
sx q[1];
rz(-0.86135094) q[1];
sx q[1];
rz(2.6497926) q[1];
rz(-pi) q[2];
rz(-1.5063498) q[3];
sx q[3];
rz(-1.6985053) q[3];
sx q[3];
rz(0.32999048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(1.3736003) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(2.7048892) q[0];
rz(-2.9084335) q[1];
sx q[1];
rz(-1.8893087) q[1];
sx q[1];
rz(-2.8312347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6757641) q[0];
sx q[0];
rz(-0.40256631) q[0];
sx q[0];
rz(-2.7990544) q[0];
rz(-pi) q[1];
rz(0.15375806) q[2];
sx q[2];
rz(-2.4506844) q[2];
sx q[2];
rz(1.549987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0879285) q[1];
sx q[1];
rz(-1.9287319) q[1];
sx q[1];
rz(3.1217561) q[1];
rz(1.849732) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.9524277) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(0.24965832) q[0];
rz(-1.5769618) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91709671) q[0];
sx q[0];
rz(-1.774569) q[0];
sx q[0];
rz(2.8208371) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(0.86415926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1340027) q[1];
sx q[1];
rz(-2.3569199) q[1];
sx q[1];
rz(0.44962928) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8364041) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.27328086) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.3195066) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34981397) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(2.8836024) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.649883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0598037) q[0];
sx q[0];
rz(-0.44133082) q[0];
sx q[0];
rz(0.23131891) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67955534) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(0.92781767) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.7548435) q[1];
sx q[1];
rz(-2.1347087) q[1];
sx q[1];
rz(2.2350603) q[1];
rz(-1.8317354) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(0.84364676) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(2.2475524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-0.31232467) q[1];
sx q[1];
rz(3.1076028) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2798529) q[0];
sx q[0];
rz(-1.4089157) q[0];
sx q[0];
rz(0.78761657) q[0];
rz(-pi) q[1];
rz(-2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(0.63461441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9719203) q[1];
sx q[1];
rz(-1.6213413) q[1];
sx q[1];
rz(-0.89269841) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9339318) q[3];
sx q[3];
rz(-0.62519473) q[3];
sx q[3];
rz(-3.0561662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.9411496) q[3];
sx q[3];
rz(0.14311895) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(1.4454909) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96955339) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(-2.0041549) q[0];
x q[1];
rz(-2.3169575) q[2];
sx q[2];
rz(-2.3687009) q[2];
sx q[2];
rz(-0.075721272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98390025) q[1];
sx q[1];
rz(-2.1196113) q[1];
sx q[1];
rz(1.7556612) q[1];
x q[2];
rz(-2.3950855) q[3];
sx q[3];
rz(-0.81614796) q[3];
sx q[3];
rz(-1.5619123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(0.40714804) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(2.8919162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(0.30977419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6471841) q[0];
sx q[0];
rz(-1.697288) q[0];
sx q[0];
rz(2.6148318) q[0];
x q[1];
rz(1.3938815) q[2];
sx q[2];
rz(-1.5801016) q[2];
sx q[2];
rz(-2.6975346) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0032776912) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(-3.0823207) q[1];
rz(-pi) q[2];
rz(0.36848948) q[3];
sx q[3];
rz(-2.512305) q[3];
sx q[3];
rz(-3.1143509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(-1.2072198) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(-0.65565482) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(-1.4996128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3033894) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(-0.10586664) q[0];
rz(-3.0707804) q[2];
sx q[2];
rz(-2.376308) q[2];
sx q[2];
rz(-0.27809696) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26350281) q[1];
sx q[1];
rz(-2.6260758) q[1];
sx q[1];
rz(-2.0104736) q[1];
rz(-1.7435944) q[3];
sx q[3];
rz(-1.5683335) q[3];
sx q[3];
rz(0.60009225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(2.1949027) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832227) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(3.070667) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(2.0680188) q[2];
sx q[2];
rz(-1.9855269) q[2];
sx q[2];
rz(-1.2863458) q[2];
rz(2.5881913) q[3];
sx q[3];
rz(-0.80080606) q[3];
sx q[3];
rz(-1.3047119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
