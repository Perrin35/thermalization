OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7967427) q[0];
sx q[0];
rz(-2.8673708) q[0];
sx q[0];
rz(2.5728777) q[0];
rz(-5.0721726) q[1];
sx q[1];
rz(0.99744263) q[1];
sx q[1];
rz(9.6923516) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9833004) q[0];
sx q[0];
rz(-2.112297) q[0];
sx q[0];
rz(-1.3356317) q[0];
rz(1.2368343) q[2];
sx q[2];
rz(-1.4638454) q[2];
sx q[2];
rz(-3.1283875) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.888511) q[1];
sx q[1];
rz(-0.0086697658) q[1];
sx q[1];
rz(1.289283) q[1];
rz(-1.9851793) q[3];
sx q[3];
rz(-0.15358812) q[3];
sx q[3];
rz(-1.6382662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8895175) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(1.1313103) q[2];
rz(0.45025292) q[3];
sx q[3];
rz(-0.69142747) q[3];
sx q[3];
rz(-0.94436193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(1.407628) q[0];
rz(2.6990926) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(2.5462467) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2059518) q[0];
sx q[0];
rz(-1.8186339) q[0];
sx q[0];
rz(1.3252844) q[0];
x q[1];
rz(-2.2240665) q[2];
sx q[2];
rz(-1.1703614) q[2];
sx q[2];
rz(1.8446856) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7390531) q[1];
sx q[1];
rz(-1.1418704) q[1];
sx q[1];
rz(1.5651171) q[1];
x q[2];
rz(-0.54259681) q[3];
sx q[3];
rz(-0.68477453) q[3];
sx q[3];
rz(-0.34943504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2596316) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(-0.76914966) q[2];
rz(-1.361557) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(2.5462525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24460569) q[0];
sx q[0];
rz(-2.2266882) q[0];
sx q[0];
rz(2.8047674) q[0];
rz(-1.4312076) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(0.097188458) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1203994) q[0];
sx q[0];
rz(-1.0836856) q[0];
sx q[0];
rz(2.9117286) q[0];
rz(-pi) q[1];
rz(1.4732185) q[2];
sx q[2];
rz(-1.9505902) q[2];
sx q[2];
rz(1.9623836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0984368) q[1];
sx q[1];
rz(-1.7088026) q[1];
sx q[1];
rz(0.34855493) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6146856) q[3];
sx q[3];
rz(-2.1873173) q[3];
sx q[3];
rz(-1.34672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97239697) q[2];
sx q[2];
rz(-1.376386) q[2];
sx q[2];
rz(0.81673679) q[2];
rz(-0.9225325) q[3];
sx q[3];
rz(-0.43729344) q[3];
sx q[3];
rz(-1.1923265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794491) q[0];
sx q[0];
rz(-1.3379931) q[0];
sx q[0];
rz(2.2747967) q[0];
rz(1.2358933) q[1];
sx q[1];
rz(-2.0265323) q[1];
sx q[1];
rz(-2.8996276) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075027) q[0];
sx q[0];
rz(-2.3781812) q[0];
sx q[0];
rz(-1.1795189) q[0];
x q[1];
rz(0.073512065) q[2];
sx q[2];
rz(-1.7318372) q[2];
sx q[2];
rz(1.1390151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45631187) q[1];
sx q[1];
rz(-1.1667098) q[1];
sx q[1];
rz(-1.7194952) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5210152) q[3];
sx q[3];
rz(-2.0672653) q[3];
sx q[3];
rz(-2.1703326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26291004) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(2.6061457) q[2];
rz(2.6483436) q[3];
sx q[3];
rz(-0.89087629) q[3];
sx q[3];
rz(-2.6935327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0754452) q[0];
sx q[0];
rz(-0.15114052) q[0];
sx q[0];
rz(-0.84392631) q[0];
rz(1.4752202) q[1];
sx q[1];
rz(-1.4862783) q[1];
sx q[1];
rz(-2.7484238) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1895904) q[0];
sx q[0];
rz(-1.8247427) q[0];
sx q[0];
rz(0.80811951) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8507502) q[2];
sx q[2];
rz(-2.484349) q[2];
sx q[2];
rz(0.90279451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1816493) q[1];
sx q[1];
rz(-1.2928597) q[1];
sx q[1];
rz(2.6647869) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6541071) q[3];
sx q[3];
rz(-0.84934399) q[3];
sx q[3];
rz(-0.3075499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4142485) q[2];
sx q[2];
rz(-2.9372637) q[2];
sx q[2];
rz(-2.7434529) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-2.3530493) q[3];
sx q[3];
rz(-1.3032234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141465) q[0];
sx q[0];
rz(-1.0108203) q[0];
sx q[0];
rz(-0.8557125) q[0];
rz(0.69264597) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(0.2690424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1097668) q[0];
sx q[0];
rz(-0.19309154) q[0];
sx q[0];
rz(1.4433799) q[0];
x q[1];
rz(-0.034560692) q[2];
sx q[2];
rz(-1.0644352) q[2];
sx q[2];
rz(2.0848372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4384296) q[1];
sx q[1];
rz(-2.3932578) q[1];
sx q[1];
rz(-1.8412043) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9773922) q[3];
sx q[3];
rz(-1.3815666) q[3];
sx q[3];
rz(0.38092074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1143703) q[2];
sx q[2];
rz(-1.8222858) q[2];
sx q[2];
rz(-3.031292) q[2];
rz(2.2731764) q[3];
sx q[3];
rz(-1.1766368) q[3];
sx q[3];
rz(1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480943) q[0];
sx q[0];
rz(-0.49999923) q[0];
sx q[0];
rz(0.86135832) q[0];
rz(-1.6294847) q[1];
sx q[1];
rz(-1.6810828) q[1];
sx q[1];
rz(0.82383627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31164661) q[0];
sx q[0];
rz(-0.59893805) q[0];
sx q[0];
rz(-1.3715368) q[0];
rz(-pi) q[1];
rz(-0.74987133) q[2];
sx q[2];
rz(-1.8786542) q[2];
sx q[2];
rz(-1.4726382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41125248) q[1];
sx q[1];
rz(-1.3578051) q[1];
sx q[1];
rz(-0.61720444) q[1];
x q[2];
rz(1.2560315) q[3];
sx q[3];
rz(-0.94292313) q[3];
sx q[3];
rz(1.0726269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6350101) q[2];
sx q[2];
rz(-2.625605) q[2];
sx q[2];
rz(0.79279509) q[2];
rz(-3.1363764) q[3];
sx q[3];
rz(-0.79013932) q[3];
sx q[3];
rz(-1.4504356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.97013) q[0];
sx q[0];
rz(-1.1928394) q[0];
sx q[0];
rz(2.9520853) q[0];
rz(2.3686523) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(-2.5453087) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0590116) q[0];
sx q[0];
rz(-2.4333911) q[0];
sx q[0];
rz(-0.26074683) q[0];
rz(-pi) q[1];
rz(0.79769602) q[2];
sx q[2];
rz(-2.0828649) q[2];
sx q[2];
rz(1.860581) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2345727) q[1];
sx q[1];
rz(-1.4661745) q[1];
sx q[1];
rz(-1.6638882) q[1];
x q[2];
rz(2.6176386) q[3];
sx q[3];
rz(-2.2125681) q[3];
sx q[3];
rz(0.061836035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4160055) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(-1.928398) q[2];
rz(-0.98617918) q[3];
sx q[3];
rz(-1.4158019) q[3];
sx q[3];
rz(-1.8825611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0202494) q[0];
sx q[0];
rz(-1.5556524) q[0];
sx q[0];
rz(-1.1100618) q[0];
rz(2.8129261) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(-1.2967671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.219562) q[0];
sx q[0];
rz(-2.2264105) q[0];
sx q[0];
rz(0.11830413) q[0];
x q[1];
rz(-1.4234366) q[2];
sx q[2];
rz(-0.40032101) q[2];
sx q[2];
rz(1.6860698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4420524) q[1];
sx q[1];
rz(-1.5083543) q[1];
sx q[1];
rz(2.6477835) q[1];
x q[2];
rz(-2.8696257) q[3];
sx q[3];
rz(-2.4578531) q[3];
sx q[3];
rz(-2.7955987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0686331) q[2];
sx q[2];
rz(-0.98560846) q[2];
sx q[2];
rz(3.0375321) q[2];
rz(2.0067298) q[3];
sx q[3];
rz(-1.7770504) q[3];
sx q[3];
rz(-0.22741905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4858911) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(-2.4110598) q[0];
rz(-0.55745521) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(2.7117859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308425) q[0];
sx q[0];
rz(-1.5702899) q[0];
sx q[0];
rz(2.3798126) q[0];
x q[1];
rz(-2.7676959) q[2];
sx q[2];
rz(-1.3329525) q[2];
sx q[2];
rz(1.4542945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81372637) q[1];
sx q[1];
rz(-1.9414492) q[1];
sx q[1];
rz(-3.0560545) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12390512) q[3];
sx q[3];
rz(-1.2031816) q[3];
sx q[3];
rz(-1.0197365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3320015) q[2];
sx q[2];
rz(-2.1198699) q[2];
sx q[2];
rz(0.29279718) q[2];
rz(0.14033595) q[3];
sx q[3];
rz(-2.1122746) q[3];
sx q[3];
rz(0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.77107) q[0];
sx q[0];
rz(-1.4650383) q[0];
sx q[0];
rz(-2.9248206) q[0];
rz(-2.4222005) q[1];
sx q[1];
rz(-1.2750625) q[1];
sx q[1];
rz(0.19663179) q[1];
rz(-2.0784081) q[2];
sx q[2];
rz(-1.4127991) q[2];
sx q[2];
rz(0.71411919) q[2];
rz(0.1340014) q[3];
sx q[3];
rz(-1.25105) q[3];
sx q[3];
rz(2.7950263) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
