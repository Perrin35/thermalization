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
rz(0.69005203) q[0];
sx q[0];
rz(8.5741841) q[0];
sx q[0];
rz(9.2246715) q[0];
rz(2.9432358) q[1];
sx q[1];
rz(-2.618572) q[1];
sx q[1];
rz(-1.3318292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057077335) q[0];
sx q[0];
rz(-0.50327089) q[0];
sx q[0];
rz(-2.0664735) q[0];
x q[1];
rz(2.4780689) q[2];
sx q[2];
rz(-2.4276456) q[2];
sx q[2];
rz(2.6235142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3568452) q[1];
sx q[1];
rz(-0.6797528) q[1];
sx q[1];
rz(0.076514449) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0588714) q[3];
sx q[3];
rz(-2.045407) q[3];
sx q[3];
rz(3.0623013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.072210463) q[2];
sx q[2];
rz(-0.36082265) q[2];
sx q[2];
rz(1.8270095) q[2];
rz(0.9032816) q[3];
sx q[3];
rz(-1.7777781) q[3];
sx q[3];
rz(-0.34590736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8697934) q[0];
sx q[0];
rz(-1.6558187) q[0];
sx q[0];
rz(0.85579175) q[0];
rz(2.3675512) q[1];
sx q[1];
rz(-1.9776521) q[1];
sx q[1];
rz(-0.78786293) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5408527) q[0];
sx q[0];
rz(-0.2517646) q[0];
sx q[0];
rz(1.6600079) q[0];
rz(-pi) q[1];
x q[1];
rz(2.79736) q[2];
sx q[2];
rz(-1.5308799) q[2];
sx q[2];
rz(-1.0913864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3968421) q[1];
sx q[1];
rz(-1.9991268) q[1];
sx q[1];
rz(-2.7672) q[1];
x q[2];
rz(-1.3608906) q[3];
sx q[3];
rz(-0.91448254) q[3];
sx q[3];
rz(2.8950952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.065980109) q[2];
sx q[2];
rz(-1.9115261) q[2];
sx q[2];
rz(-0.39609972) q[2];
rz(1.5710477) q[3];
sx q[3];
rz(-1.5038265) q[3];
sx q[3];
rz(1.3554696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9539255) q[0];
sx q[0];
rz(-1.6463771) q[0];
sx q[0];
rz(3.0809825) q[0];
rz(-2.4568779) q[1];
sx q[1];
rz(-1.4002607) q[1];
sx q[1];
rz(2.8025467) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.402917) q[0];
sx q[0];
rz(-2.5581723) q[0];
sx q[0];
rz(-3.1039799) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57444797) q[2];
sx q[2];
rz(-0.22413218) q[2];
sx q[2];
rz(1.4719065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10968929) q[1];
sx q[1];
rz(-1.8432143) q[1];
sx q[1];
rz(1.8930356) q[1];
x q[2];
rz(-1.5590844) q[3];
sx q[3];
rz(-2.3475966) q[3];
sx q[3];
rz(2.7744966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.48245779) q[2];
sx q[2];
rz(-1.277801) q[2];
sx q[2];
rz(2.6711312) q[2];
rz(0.86236924) q[3];
sx q[3];
rz(-2.6926398) q[3];
sx q[3];
rz(1.5615777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793133) q[0];
sx q[0];
rz(-1.1184432) q[0];
sx q[0];
rz(-0.40801868) q[0];
rz(-1.3511924) q[1];
sx q[1];
rz(-1.2974757) q[1];
sx q[1];
rz(1.8203576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119638) q[0];
sx q[0];
rz(-0.59998673) q[0];
sx q[0];
rz(2.5866051) q[0];
x q[1];
rz(-2.2259813) q[2];
sx q[2];
rz(-1.5024937) q[2];
sx q[2];
rz(3.0223522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3258354) q[1];
sx q[1];
rz(-1.6920631) q[1];
sx q[1];
rz(-1.6440637) q[1];
x q[2];
rz(-1.9923237) q[3];
sx q[3];
rz(-1.9644004) q[3];
sx q[3];
rz(2.7019661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9185751) q[2];
sx q[2];
rz(-1.431798) q[2];
sx q[2];
rz(-2.4187386) q[2];
rz(0.84960788) q[3];
sx q[3];
rz(-1.9690211) q[3];
sx q[3];
rz(-2.1037219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866813) q[0];
sx q[0];
rz(-0.51690042) q[0];
sx q[0];
rz(0.55996672) q[0];
rz(-0.2050744) q[1];
sx q[1];
rz(-1.2566902) q[1];
sx q[1];
rz(-2.8505039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2084853) q[0];
sx q[0];
rz(-1.4015084) q[0];
sx q[0];
rz(-0.27219682) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0801513) q[2];
sx q[2];
rz(-0.80020088) q[2];
sx q[2];
rz(2.7903008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6711802) q[1];
sx q[1];
rz(-0.7268097) q[1];
sx q[1];
rz(1.2695168) q[1];
rz(-pi) q[2];
rz(0.17125968) q[3];
sx q[3];
rz(-2.3823839) q[3];
sx q[3];
rz(-1.9943135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2228955) q[2];
sx q[2];
rz(-1.4783858) q[2];
sx q[2];
rz(0.29829868) q[2];
rz(-1.1693303) q[3];
sx q[3];
rz(-2.1746641) q[3];
sx q[3];
rz(-1.5112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8339612) q[0];
sx q[0];
rz(-1.4606322) q[0];
sx q[0];
rz(0.6024012) q[0];
rz(0.3715474) q[1];
sx q[1];
rz(-1.5448152) q[1];
sx q[1];
rz(-1.0317624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69799549) q[0];
sx q[0];
rz(-1.5549608) q[0];
sx q[0];
rz(-1.7450733) q[0];
rz(-2.6844031) q[2];
sx q[2];
rz(-2.6962099) q[2];
sx q[2];
rz(3.074844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25615869) q[1];
sx q[1];
rz(-2.6226461) q[1];
sx q[1];
rz(-2.7536254) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57547456) q[3];
sx q[3];
rz(-2.1361975) q[3];
sx q[3];
rz(-2.7421212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3370257) q[2];
sx q[2];
rz(-2.0402543) q[2];
sx q[2];
rz(-1.5186914) q[2];
rz(0.8484146) q[3];
sx q[3];
rz(-2.266326) q[3];
sx q[3];
rz(-3.1019822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1161716) q[0];
sx q[0];
rz(-0.051055901) q[0];
sx q[0];
rz(-2.1436932) q[0];
rz(-0.2746703) q[1];
sx q[1];
rz(-2.2614567) q[1];
sx q[1];
rz(-1.3708699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3034015) q[0];
sx q[0];
rz(-2.5104021) q[0];
sx q[0];
rz(3.1188008) q[0];
rz(1.2148803) q[2];
sx q[2];
rz(-1.4265665) q[2];
sx q[2];
rz(3.1145417) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8805367) q[1];
sx q[1];
rz(-0.66268259) q[1];
sx q[1];
rz(2.7438394) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3418535) q[3];
sx q[3];
rz(-2.1465786) q[3];
sx q[3];
rz(2.6588001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85713282) q[2];
sx q[2];
rz(-2.1993115) q[2];
sx q[2];
rz(-0.029646309) q[2];
rz(0.61521411) q[3];
sx q[3];
rz(-1.3968202) q[3];
sx q[3];
rz(2.7992547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0520332) q[0];
sx q[0];
rz(-2.679306) q[0];
sx q[0];
rz(0.7830559) q[0];
rz(-2.0217333) q[1];
sx q[1];
rz(-2.4744787) q[1];
sx q[1];
rz(-1.3979744) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2180041) q[0];
sx q[0];
rz(-2.5177885) q[0];
sx q[0];
rz(-0.82360928) q[0];
rz(-pi) q[1];
rz(-1.8552165) q[2];
sx q[2];
rz(-1.7952953) q[2];
sx q[2];
rz(-2.7530991) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82247558) q[1];
sx q[1];
rz(-2.2386595) q[1];
sx q[1];
rz(-3.130359) q[1];
x q[2];
rz(0.21381883) q[3];
sx q[3];
rz(-2.345746) q[3];
sx q[3];
rz(2.4406432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8457501) q[2];
sx q[2];
rz(-0.58640277) q[2];
sx q[2];
rz(-1.3207377) q[2];
rz(-1.0197506) q[3];
sx q[3];
rz(-1.7176065) q[3];
sx q[3];
rz(-1.7216871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418697) q[0];
sx q[0];
rz(-2.0482735) q[0];
sx q[0];
rz(2.1871908) q[0];
rz(1.9873387) q[1];
sx q[1];
rz(-2.4955165) q[1];
sx q[1];
rz(1.0844213) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99710871) q[0];
sx q[0];
rz(-0.86291828) q[0];
sx q[0];
rz(0.055521418) q[0];
rz(-pi) q[1];
rz(1.8516638) q[2];
sx q[2];
rz(-1.3251588) q[2];
sx q[2];
rz(-2.7410067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3855648) q[1];
sx q[1];
rz(-1.4984908) q[1];
sx q[1];
rz(-2.0438309) q[1];
rz(0.19799216) q[3];
sx q[3];
rz(-0.75847236) q[3];
sx q[3];
rz(0.93659479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.083100975) q[2];
sx q[2];
rz(-2.9497171) q[2];
sx q[2];
rz(-2.7044738) q[2];
rz(1.7982177) q[3];
sx q[3];
rz(-1.3513887) q[3];
sx q[3];
rz(3.0962211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0992391) q[0];
sx q[0];
rz(-1.0004685) q[0];
sx q[0];
rz(1.7310671) q[0];
rz(2.9208185) q[1];
sx q[1];
rz(-0.53123728) q[1];
sx q[1];
rz(2.208362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30770424) q[0];
sx q[0];
rz(-1.5571463) q[0];
sx q[0];
rz(-0.74539124) q[0];
x q[1];
rz(0.48671203) q[2];
sx q[2];
rz(-0.67811869) q[2];
sx q[2];
rz(-2.202452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5367914) q[1];
sx q[1];
rz(-1.4385421) q[1];
sx q[1];
rz(-0.79562809) q[1];
rz(-1.4187063) q[3];
sx q[3];
rz(-2.503431) q[3];
sx q[3];
rz(1.9877951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.71020469) q[2];
sx q[2];
rz(-0.34222558) q[2];
sx q[2];
rz(-2.5101275) q[2];
rz(2.4225875) q[3];
sx q[3];
rz(-1.3471666) q[3];
sx q[3];
rz(-2.9001111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030293) q[0];
sx q[0];
rz(-1.6521492) q[0];
sx q[0];
rz(-1.8304992) q[0];
rz(-0.44166625) q[1];
sx q[1];
rz(-1.3413981) q[1];
sx q[1];
rz(-1.9955019) q[1];
rz(-0.031671192) q[2];
sx q[2];
rz(-1.1902255) q[2];
sx q[2];
rz(-2.5524216) q[2];
rz(2.5720015) q[3];
sx q[3];
rz(-1.6574331) q[3];
sx q[3];
rz(2.1059753) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
