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
rz(-1.4581207) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0163907) q[0];
sx q[0];
rz(-0.0076961829) q[0];
sx q[0];
rz(1.9793235) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7952806) q[2];
sx q[2];
rz(-2.6484688) q[2];
sx q[2];
rz(0.43362793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5344421) q[1];
sx q[1];
rz(-1.0118224) q[1];
sx q[1];
rz(0.89627756) q[1];
x q[2];
rz(-3.0786282) q[3];
sx q[3];
rz(-0.28423542) q[3];
sx q[3];
rz(0.96719826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.014341982) q[2];
sx q[2];
rz(-0.6659826) q[2];
sx q[2];
rz(-1.8951529) q[2];
rz(2.1261101) q[3];
sx q[3];
rz(-0.4963488) q[3];
sx q[3];
rz(-0.61771667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.604973) q[0];
sx q[0];
rz(-0.07285694) q[0];
sx q[0];
rz(2.4541722) q[0];
rz(-0.61126002) q[1];
sx q[1];
rz(-1.5115073) q[1];
sx q[1];
rz(-0.90868178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9729823) q[0];
sx q[0];
rz(-1.4797204) q[0];
sx q[0];
rz(-2.220972) q[0];
rz(-pi) q[1];
rz(1.3503752) q[2];
sx q[2];
rz(-0.62602717) q[2];
sx q[2];
rz(-1.859425) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7327888) q[1];
sx q[1];
rz(-2.286274) q[1];
sx q[1];
rz(-2.4092731) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4439529) q[3];
sx q[3];
rz(-1.0174362) q[3];
sx q[3];
rz(-0.07088319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0574135) q[2];
sx q[2];
rz(-2.1647191) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77370206) q[0];
sx q[0];
rz(-0.20562085) q[0];
sx q[0];
rz(-0.50262991) q[0];
rz(1.8058808) q[1];
sx q[1];
rz(-3.0323961) q[1];
sx q[1];
rz(1.4044382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.805946) q[0];
sx q[0];
rz(-2.5281161) q[0];
sx q[0];
rz(-2.4015266) q[0];
rz(2.6790256) q[2];
sx q[2];
rz(-0.41282648) q[2];
sx q[2];
rz(-0.040469801) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.100443) q[1];
sx q[1];
rz(-1.2432069) q[1];
sx q[1];
rz(-1.4781152) q[1];
rz(2.5504774) q[3];
sx q[3];
rz(-1.1151259) q[3];
sx q[3];
rz(2.299813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8414843) q[2];
sx q[2];
rz(-0.69977641) q[2];
sx q[2];
rz(0.67467275) q[2];
rz(-0.35215968) q[3];
sx q[3];
rz(-1.1511185) q[3];
sx q[3];
rz(-1.6262416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34053764) q[0];
sx q[0];
rz(-3.099649) q[0];
sx q[0];
rz(1.1658143) q[0];
rz(-0.55533987) q[1];
sx q[1];
rz(-2.1644939) q[1];
sx q[1];
rz(-1.7305444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80361405) q[0];
sx q[0];
rz(-0.25000152) q[0];
sx q[0];
rz(-3.0743272) q[0];
rz(-pi) q[1];
rz(1.8952266) q[2];
sx q[2];
rz(-0.55098767) q[2];
sx q[2];
rz(1.18206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6785612) q[1];
sx q[1];
rz(-0.62778607) q[1];
sx q[1];
rz(2.1902172) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23616236) q[3];
sx q[3];
rz(-2.3104226) q[3];
sx q[3];
rz(2.1787852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0393151) q[2];
sx q[2];
rz(-2.6199876) q[2];
sx q[2];
rz(0.48307854) q[2];
rz(0.77110243) q[3];
sx q[3];
rz(-1.5604138) q[3];
sx q[3];
rz(1.34351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.1191206) q[0];
sx q[0];
rz(-2.6799057) q[0];
sx q[0];
rz(2.9126677) q[0];
rz(-1.6772894) q[1];
sx q[1];
rz(-2.5720282) q[1];
sx q[1];
rz(-0.48008188) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71704536) q[0];
sx q[0];
rz(-1.8295032) q[0];
sx q[0];
rz(0.26820582) q[0];
rz(2.6998015) q[2];
sx q[2];
rz(-1.3350772) q[2];
sx q[2];
rz(-3.1093895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1118873) q[1];
sx q[1];
rz(-2.7656274) q[1];
sx q[1];
rz(-0.13891797) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94575092) q[3];
sx q[3];
rz(-1.7202999) q[3];
sx q[3];
rz(2.3760419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0698801) q[2];
sx q[2];
rz(-1.1249079) q[2];
sx q[2];
rz(-1.9055535) q[2];
rz(1.49336) q[3];
sx q[3];
rz(-0.83087102) q[3];
sx q[3];
rz(1.6911471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1278095) q[0];
sx q[0];
rz(-2.4007128) q[0];
sx q[0];
rz(2.2770449) q[0];
rz(2.3132482) q[1];
sx q[1];
rz(-1.336785) q[1];
sx q[1];
rz(-2.4116662) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8720388) q[0];
sx q[0];
rz(-2.3667496) q[0];
sx q[0];
rz(0.86685769) q[0];
rz(-pi) q[1];
rz(-1.7696769) q[2];
sx q[2];
rz(-0.91105295) q[2];
sx q[2];
rz(-0.64780647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3507281) q[1];
sx q[1];
rz(-0.58907382) q[1];
sx q[1];
rz(-1.1228611) q[1];
rz(-pi) q[2];
rz(1.6248405) q[3];
sx q[3];
rz(-0.84815787) q[3];
sx q[3];
rz(2.8703226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56968969) q[2];
sx q[2];
rz(-2.9986585) q[2];
sx q[2];
rz(-1.2804871) q[2];
rz(1.5754383) q[3];
sx q[3];
rz(-2.1011293) q[3];
sx q[3];
rz(2.2787826) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5363619) q[0];
sx q[0];
rz(-1.0030168) q[0];
sx q[0];
rz(2.5501116) q[0];
rz(0.65762562) q[1];
sx q[1];
rz(-1.769519) q[1];
sx q[1];
rz(-3.0234911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5629092) q[0];
sx q[0];
rz(-1.9171035) q[0];
sx q[0];
rz(1.8150532) q[0];
rz(-pi) q[1];
rz(1.4976229) q[2];
sx q[2];
rz(-1.055848) q[2];
sx q[2];
rz(0.63637444) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.98886314) q[1];
sx q[1];
rz(-2.0901457) q[1];
sx q[1];
rz(-1.4593655) q[1];
rz(-pi) q[2];
rz(-2.2198417) q[3];
sx q[3];
rz(-0.27652446) q[3];
sx q[3];
rz(0.10504237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20146951) q[2];
sx q[2];
rz(-2.757299) q[2];
sx q[2];
rz(1.2598134) q[2];
rz(1.2093557) q[3];
sx q[3];
rz(-1.3426547) q[3];
sx q[3];
rz(-0.847009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92846337) q[0];
sx q[0];
rz(-2.6452112) q[0];
sx q[0];
rz(1.553836) q[0];
rz(1.8153048) q[1];
sx q[1];
rz(-0.98618788) q[1];
sx q[1];
rz(-2.3415668) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84811178) q[0];
sx q[0];
rz(-1.9514553) q[0];
sx q[0];
rz(-1.8649376) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7615835) q[2];
sx q[2];
rz(-1.933145) q[2];
sx q[2];
rz(1.5630388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86944888) q[1];
sx q[1];
rz(-0.71681685) q[1];
sx q[1];
rz(0.79811704) q[1];
rz(-pi) q[2];
rz(-0.38626892) q[3];
sx q[3];
rz(-2.093654) q[3];
sx q[3];
rz(2.8012365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.942261) q[2];
sx q[2];
rz(-2.3776725) q[2];
sx q[2];
rz(1.2660816) q[2];
rz(-1.0670886) q[3];
sx q[3];
rz(-1.7231924) q[3];
sx q[3];
rz(-3.0838695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98457321) q[0];
sx q[0];
rz(-0.83681256) q[0];
sx q[0];
rz(-2.9465604) q[0];
rz(-2.2795279) q[1];
sx q[1];
rz(-1.74086) q[1];
sx q[1];
rz(-0.21924266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57306141) q[0];
sx q[0];
rz(-0.6210621) q[0];
sx q[0];
rz(2.0173418) q[0];
x q[1];
rz(2.9621808) q[2];
sx q[2];
rz(-1.1288912) q[2];
sx q[2];
rz(-2.4561858) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2229206) q[1];
sx q[1];
rz(-1.5377475) q[1];
sx q[1];
rz(2.6085684) q[1];
rz(-pi) q[2];
rz(-2.6126543) q[3];
sx q[3];
rz(-2.248507) q[3];
sx q[3];
rz(2.4916745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4054883) q[2];
sx q[2];
rz(-0.89724237) q[2];
sx q[2];
rz(1.7264977) q[2];
rz(1.3380346) q[3];
sx q[3];
rz(-1.8357364) q[3];
sx q[3];
rz(3.0726748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59704429) q[0];
sx q[0];
rz(-0.91240779) q[0];
sx q[0];
rz(1.9507116) q[0];
rz(-1.4471588) q[1];
sx q[1];
rz(-1.2971327) q[1];
sx q[1];
rz(1.7097293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53729355) q[0];
sx q[0];
rz(-1.3621289) q[0];
sx q[0];
rz(-0.66402004) q[0];
rz(-0.5848627) q[2];
sx q[2];
rz(-1.9205127) q[2];
sx q[2];
rz(2.62874) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9976342) q[1];
sx q[1];
rz(-1.1678417) q[1];
sx q[1];
rz(2.2657299) q[1];
rz(-2.3031381) q[3];
sx q[3];
rz(-1.5946232) q[3];
sx q[3];
rz(1.4122653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4900964) q[2];
sx q[2];
rz(-0.33612529) q[2];
sx q[2];
rz(2.2641342) q[2];
rz(1.6089926) q[3];
sx q[3];
rz(-1.7209777) q[3];
sx q[3];
rz(-1.0271614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.132591) q[0];
sx q[0];
rz(-1.4132211) q[0];
sx q[0];
rz(2.6268517) q[0];
rz(-2.4280836) q[1];
sx q[1];
rz(-2.3206354) q[1];
sx q[1];
rz(-1.6101507) q[1];
rz(1.8590676) q[2];
sx q[2];
rz(-1.9296411) q[2];
sx q[2];
rz(1.267098) q[2];
rz(1.3446829) q[3];
sx q[3];
rz(-1.825708) q[3];
sx q[3];
rz(0.59151266) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
