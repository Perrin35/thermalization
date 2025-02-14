OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7944827) q[0];
sx q[0];
rz(-1.4248983) q[0];
sx q[0];
rz(-2.8704034) q[0];
rz(-2.7746692) q[1];
sx q[1];
rz(-0.75777268) q[1];
sx q[1];
rz(1.6834719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.125202) q[0];
sx q[0];
rz(-3.1338965) q[0];
sx q[0];
rz(-1.9793235) q[0];
rz(-1.0881937) q[2];
sx q[2];
rz(-1.6763684) q[2];
sx q[2];
rz(2.2029049) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5196412) q[1];
sx q[1];
rz(-0.84712815) q[1];
sx q[1];
rz(-0.78459854) q[1];
rz(1.5524158) q[3];
sx q[3];
rz(-1.2871398) q[3];
sx q[3];
rz(-2.2399834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.014341982) q[2];
sx q[2];
rz(-2.47561) q[2];
sx q[2];
rz(1.2464397) q[2];
rz(2.1261101) q[3];
sx q[3];
rz(-0.4963488) q[3];
sx q[3];
rz(-0.61771667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.604973) q[0];
sx q[0];
rz(-3.0687357) q[0];
sx q[0];
rz(2.4541722) q[0];
rz(2.5303326) q[1];
sx q[1];
rz(-1.5115073) q[1];
sx q[1];
rz(2.2329109) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8584367) q[0];
sx q[0];
rz(-0.65560616) q[0];
sx q[0];
rz(-1.7205419) q[0];
rz(-pi) q[1];
rz(2.9848027) q[2];
sx q[2];
rz(-2.1794381) q[2];
sx q[2];
rz(1.551924) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3492891) q[1];
sx q[1];
rz(-2.1667203) q[1];
sx q[1];
rz(0.915058) q[1];
rz(1.4439529) q[3];
sx q[3];
rz(-2.1241564) q[3];
sx q[3];
rz(3.0707095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0574135) q[2];
sx q[2];
rz(-0.97687352) q[2];
sx q[2];
rz(-1.9834391) q[2];
rz(-0.17364764) q[3];
sx q[3];
rz(-1.6896788) q[3];
sx q[3];
rz(-2.4524073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3678906) q[0];
sx q[0];
rz(-2.9359718) q[0];
sx q[0];
rz(-0.50262991) q[0];
rz(1.8058808) q[1];
sx q[1];
rz(-0.10919658) q[1];
sx q[1];
rz(1.7371545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40622845) q[0];
sx q[0];
rz(-1.9695008) q[0];
sx q[0];
rz(-0.47946341) q[0];
x q[1];
rz(0.4625671) q[2];
sx q[2];
rz(-0.41282648) q[2];
sx q[2];
rz(0.040469801) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5002513) q[1];
sx q[1];
rz(-1.6585357) q[1];
sx q[1];
rz(-0.32890202) q[1];
rz(-pi) q[2];
rz(-0.59111528) q[3];
sx q[3];
rz(-2.0264668) q[3];
sx q[3];
rz(-2.299813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8414843) q[2];
sx q[2];
rz(-0.69977641) q[2];
sx q[2];
rz(-2.4669199) q[2];
rz(-2.789433) q[3];
sx q[3];
rz(-1.1511185) q[3];
sx q[3];
rz(-1.5153511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.34053764) q[0];
sx q[0];
rz(-0.041943701) q[0];
sx q[0];
rz(-1.9757784) q[0];
rz(-2.5862528) q[1];
sx q[1];
rz(-2.1644939) q[1];
sx q[1];
rz(1.7305444) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87303091) q[0];
sx q[0];
rz(-1.3213722) q[0];
sx q[0];
rz(1.5879575) q[0];
x q[1];
rz(-1.0434112) q[2];
sx q[2];
rz(-1.7384656) q[2];
sx q[2];
rz(0.66780082) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6785612) q[1];
sx q[1];
rz(-2.5138066) q[1];
sx q[1];
rz(2.1902172) q[1];
rz(-pi) q[2];
rz(-1.3197696) q[3];
sx q[3];
rz(-2.3720312) q[3];
sx q[3];
rz(-1.3057452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1022776) q[2];
sx q[2];
rz(-0.5216051) q[2];
sx q[2];
rz(-0.48307854) q[2];
rz(-2.3704902) q[3];
sx q[3];
rz(-1.5811788) q[3];
sx q[3];
rz(-1.34351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.022472) q[0];
sx q[0];
rz(-0.461687) q[0];
sx q[0];
rz(0.22892496) q[0];
rz(1.6772894) q[1];
sx q[1];
rz(-0.56956446) q[1];
sx q[1];
rz(-0.48008188) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0374566) q[0];
sx q[0];
rz(-2.7711282) q[0];
sx q[0];
rz(2.3568826) q[0];
rz(-0.44179113) q[2];
sx q[2];
rz(-1.8065154) q[2];
sx q[2];
rz(-0.032203151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.029705392) q[1];
sx q[1];
rz(-2.7656274) q[1];
sx q[1];
rz(0.13891797) q[1];
rz(2.1958417) q[3];
sx q[3];
rz(-1.7202999) q[3];
sx q[3];
rz(0.76555071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0698801) q[2];
sx q[2];
rz(-2.0166848) q[2];
sx q[2];
rz(1.2360392) q[2];
rz(1.6482327) q[3];
sx q[3];
rz(-0.83087102) q[3];
sx q[3];
rz(1.4504455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1278095) q[0];
sx q[0];
rz(-2.4007128) q[0];
sx q[0];
rz(-0.86454779) q[0];
rz(-0.8283444) q[1];
sx q[1];
rz(-1.336785) q[1];
sx q[1];
rz(0.72992647) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1407699) q[0];
sx q[0];
rz(-1.0082908) q[0];
sx q[0];
rz(2.5767586) q[0];
x q[1];
rz(0.66940825) q[2];
sx q[2];
rz(-1.727549) q[2];
sx q[2];
rz(2.0956958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3507281) q[1];
sx q[1];
rz(-2.5525188) q[1];
sx q[1];
rz(1.1228611) q[1];
rz(3.0804068) q[3];
sx q[3];
rz(-2.4172999) q[3];
sx q[3];
rz(0.35288399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.56968969) q[2];
sx q[2];
rz(-0.14293417) q[2];
sx q[2];
rz(-1.8611056) q[2];
rz(-1.5754383) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(-2.2787826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60523072) q[0];
sx q[0];
rz(-2.1385758) q[0];
sx q[0];
rz(-0.59148106) q[0];
rz(-0.65762562) q[1];
sx q[1];
rz(-1.3720737) q[1];
sx q[1];
rz(-3.0234911) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0765083) q[0];
sx q[0];
rz(-1.3413091) q[0];
sx q[0];
rz(-2.7855532) q[0];
rz(0.51609765) q[2];
sx q[2];
rz(-1.6344667) q[2];
sx q[2];
rz(-2.2432567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76712045) q[1];
sx q[1];
rz(-2.6114957) q[1];
sx q[1];
rz(0.19211001) q[1];
rz(-pi) q[2];
rz(-2.2198417) q[3];
sx q[3];
rz(-2.8650682) q[3];
sx q[3];
rz(3.0365503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9401231) q[2];
sx q[2];
rz(-0.38429364) q[2];
sx q[2];
rz(1.8817792) q[2];
rz(-1.2093557) q[3];
sx q[3];
rz(-1.3426547) q[3];
sx q[3];
rz(0.847009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2131293) q[0];
sx q[0];
rz(-0.49638143) q[0];
sx q[0];
rz(1.553836) q[0];
rz(1.3262879) q[1];
sx q[1];
rz(-0.98618788) q[1];
sx q[1];
rz(2.3415668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84811178) q[0];
sx q[0];
rz(-1.9514553) q[0];
sx q[0];
rz(1.8649376) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3800092) q[2];
sx q[2];
rz(-1.933145) q[2];
sx q[2];
rz(-1.5630388) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0985089) q[1];
sx q[1];
rz(-1.0810163) q[1];
sx q[1];
rz(0.54651641) q[1];
x q[2];
rz(1.0142087) q[3];
sx q[3];
rz(-1.9032904) q[3];
sx q[3];
rz(-1.0300762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19933166) q[2];
sx q[2];
rz(-2.3776725) q[2];
sx q[2];
rz(1.8755111) q[2];
rz(1.0670886) q[3];
sx q[3];
rz(-1.7231924) q[3];
sx q[3];
rz(3.0838695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98457321) q[0];
sx q[0];
rz(-2.3047801) q[0];
sx q[0];
rz(-0.19503221) q[0];
rz(-0.86206478) q[1];
sx q[1];
rz(-1.74086) q[1];
sx q[1];
rz(0.21924266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5685312) q[0];
sx q[0];
rz(-0.6210621) q[0];
sx q[0];
rz(2.0173418) q[0];
rz(-pi) q[1];
rz(-1.1226038) q[2];
sx q[2];
rz(-1.4087798) q[2];
sx q[2];
rz(2.3336099) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2229206) q[1];
sx q[1];
rz(-1.6038451) q[1];
sx q[1];
rz(-2.6085684) q[1];
rz(-pi) q[2];
rz(-0.82042779) q[3];
sx q[3];
rz(-1.1667937) q[3];
sx q[3];
rz(1.272161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7361043) q[2];
sx q[2];
rz(-2.2443503) q[2];
sx q[2];
rz(-1.7264977) q[2];
rz(1.3380346) q[3];
sx q[3];
rz(-1.8357364) q[3];
sx q[3];
rz(3.0726748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59704429) q[0];
sx q[0];
rz(-2.2291849) q[0];
sx q[0];
rz(-1.190881) q[0];
rz(1.6944338) q[1];
sx q[1];
rz(-1.2971327) q[1];
sx q[1];
rz(1.7097293) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9473744) q[0];
sx q[0];
rz(-2.2179513) q[0];
sx q[0];
rz(-1.8334652) q[0];
rz(-pi) q[1];
rz(1.9831311) q[2];
sx q[2];
rz(-2.1160876) q[2];
sx q[2];
rz(0.83484989) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7427976) q[1];
sx q[1];
rz(-0.94091641) q[1];
sx q[1];
rz(-0.50666084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.032037246) q[3];
sx q[3];
rz(-2.3028829) q[3];
sx q[3];
rz(-2.9616431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.65149629) q[2];
sx q[2];
rz(-2.8054674) q[2];
sx q[2];
rz(0.87745848) q[2];
rz(1.5326001) q[3];
sx q[3];
rz(-1.7209777) q[3];
sx q[3];
rz(-2.1144313) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.132591) q[0];
sx q[0];
rz(-1.4132211) q[0];
sx q[0];
rz(2.6268517) q[0];
rz(2.4280836) q[1];
sx q[1];
rz(-0.82095725) q[1];
sx q[1];
rz(1.531442) q[1];
rz(-2.4930231) q[2];
sx q[2];
rz(-0.45634698) q[2];
sx q[2];
rz(1.9683051) q[2];
rz(2.8803181) q[3];
sx q[3];
rz(-1.3521104) q[3];
sx q[3];
rz(2.1043652) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
