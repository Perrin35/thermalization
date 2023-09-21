OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(2.8154362) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(1.3500554) q[1];
sx q[1];
rz(4.6842484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21539772) q[0];
sx q[0];
rz(-0.30456671) q[0];
sx q[0];
rz(2.347441) q[0];
rz(-pi) q[1];
rz(-2.5390981) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(-2.9162625) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5427248) q[1];
sx q[1];
rz(-1.3398783) q[1];
sx q[1];
rz(-2.8938107) q[1];
x q[2];
rz(-0.78674973) q[3];
sx q[3];
rz(-1.5610715) q[3];
sx q[3];
rz(1.7973289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8618384) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(2.2564783) q[2];
rz(-2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9279813) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(2.0425178) q[0];
rz(-0.66501578) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(2.2639993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7859902) q[0];
sx q[0];
rz(-0.88760932) q[0];
sx q[0];
rz(-1.3366633) q[0];
x q[1];
rz(-0.36466937) q[2];
sx q[2];
rz(-1.696535) q[2];
sx q[2];
rz(-2.113935) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3493537) q[1];
sx q[1];
rz(-1.6041479) q[1];
sx q[1];
rz(1.4637714) q[1];
rz(-pi) q[2];
rz(-1.9545994) q[3];
sx q[3];
rz(-1.2591397) q[3];
sx q[3];
rz(0.66031885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.39891222) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(-2.0111283) q[2];
rz(-1.8418664) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14625064) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(0.28999844) q[0];
rz(2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(3.0677632) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2732669) q[0];
sx q[0];
rz(-3.0330015) q[0];
sx q[0];
rz(-2.5383839) q[0];
x q[1];
rz(-1.9129842) q[2];
sx q[2];
rz(-1.295919) q[2];
sx q[2];
rz(-2.9800422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9170345) q[1];
sx q[1];
rz(-2.9552166) q[1];
sx q[1];
rz(1.264155) q[1];
x q[2];
rz(-2.79014) q[3];
sx q[3];
rz(-1.4593399) q[3];
sx q[3];
rz(0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6510216) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(-2.4948965) q[2];
rz(-2.0329287) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(-2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17764238) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(1.0186609) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(1.7046938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80285145) q[0];
sx q[0];
rz(-0.84206284) q[0];
sx q[0];
rz(-1.4472423) q[0];
rz(2.8145153) q[2];
sx q[2];
rz(-1.8665258) q[2];
sx q[2];
rz(-2.0138274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9096133) q[1];
sx q[1];
rz(-0.59514272) q[1];
sx q[1];
rz(-0.91722782) q[1];
rz(-3.0152263) q[3];
sx q[3];
rz(-1.4411981) q[3];
sx q[3];
rz(1.453293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.556095) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(1.1222703) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(-0.99075738) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595903) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(0.37297747) q[0];
rz(-0.22398082) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(1.8251098) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1236498) q[0];
sx q[0];
rz(-2.0262449) q[0];
sx q[0];
rz(0.42670593) q[0];
rz(-0.26222783) q[2];
sx q[2];
rz(-1.3272459) q[2];
sx q[2];
rz(1.3893931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.398048) q[1];
sx q[1];
rz(-0.053357031) q[1];
sx q[1];
rz(1.3392901) q[1];
x q[2];
rz(0.54721197) q[3];
sx q[3];
rz(-2.0378761) q[3];
sx q[3];
rz(2.3576749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0405154) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(-2.3357847) q[2];
rz(-0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(-1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39078113) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(3.0506296) q[0];
rz(0.85995752) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.8213173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1307756) q[0];
sx q[0];
rz(-0.59033075) q[0];
sx q[0];
rz(0.70461313) q[0];
rz(-pi) q[1];
rz(-2.8684902) q[2];
sx q[2];
rz(-2.5540076) q[2];
sx q[2];
rz(0.36537974) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7399866) q[1];
sx q[1];
rz(-1.1005797) q[1];
sx q[1];
rz(2.2120038) q[1];
rz(0.28989132) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(2.1260726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(-1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(-0.55554187) q[0];
rz(0.034596054) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.7506036) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3334675) q[0];
sx q[0];
rz(-1.9017178) q[0];
sx q[0];
rz(2.5740741) q[0];
x q[1];
rz(1.820307) q[2];
sx q[2];
rz(-1.9267285) q[2];
sx q[2];
rz(0.78636679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2512867) q[1];
sx q[1];
rz(-1.444) q[1];
sx q[1];
rz(-0.89819737) q[1];
rz(-0.075918003) q[3];
sx q[3];
rz(-0.77709353) q[3];
sx q[3];
rz(1.0958375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3283219) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.288712) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(0.95343268) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(1.4377726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19251084) q[0];
sx q[0];
rz(-1.7532945) q[0];
sx q[0];
rz(1.048868) q[0];
rz(-1.4346801) q[2];
sx q[2];
rz(-1.7394749) q[2];
sx q[2];
rz(0.31751925) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4096262) q[1];
sx q[1];
rz(-2.1831174) q[1];
sx q[1];
rz(2.1698976) q[1];
x q[2];
rz(1.612547) q[3];
sx q[3];
rz(-2.6937727) q[3];
sx q[3];
rz(2.8988422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(-2.9628741) q[2];
rz(-2.2802165) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20539595) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(1.0937011) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-1.0587943) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98204389) q[0];
sx q[0];
rz(-1.0315572) q[0];
sx q[0];
rz(-0.3190785) q[0];
x q[1];
rz(-2.2718614) q[2];
sx q[2];
rz(-2.5517002) q[2];
sx q[2];
rz(-1.215495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72653786) q[1];
sx q[1];
rz(-1.448436) q[1];
sx q[1];
rz(-1.9499669) q[1];
rz(-1.1684253) q[3];
sx q[3];
rz(-1.1306922) q[3];
sx q[3];
rz(-0.518706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.4808562) q[2];
rz(-0.41040928) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720298) q[0];
sx q[0];
rz(-0.27150387) q[0];
sx q[0];
rz(0.29125443) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(-1.7094918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3647389) q[0];
sx q[0];
rz(-0.59989444) q[0];
sx q[0];
rz(-2.5584496) q[0];
x q[1];
rz(-1.8068238) q[2];
sx q[2];
rz(-2.0152976) q[2];
sx q[2];
rz(-1.0016463) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11066477) q[1];
sx q[1];
rz(-2.1222161) q[1];
sx q[1];
rz(1.5196147) q[1];
rz(0.087177353) q[3];
sx q[3];
rz(-1.140359) q[3];
sx q[3];
rz(3.103053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46135205) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(-2.0001901) q[2];
rz(-1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(0.46943584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-2.4070807) q[2];
sx q[2];
rz(-2.8399158) q[2];
sx q[2];
rz(0.5010571) q[2];
rz(0.42320078) q[3];
sx q[3];
rz(-1.7508218) q[3];
sx q[3];
rz(-1.4765061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];