OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8146347) q[0];
sx q[0];
rz(3.1231413) q[0];
sx q[0];
rz(9.2508247) q[0];
rz(-1.3375018) q[1];
sx q[1];
rz(-0.87556535) q[1];
sx q[1];
rz(-2.8101885) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0259779) q[0];
sx q[0];
rz(-1.4237119) q[0];
sx q[0];
rz(0.30814104) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69281225) q[2];
sx q[2];
rz(-1.1651784) q[2];
sx q[2];
rz(1.931153) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.64664272) q[1];
sx q[1];
rz(-0.70716049) q[1];
sx q[1];
rz(2.8487514) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9865206) q[3];
sx q[3];
rz(-1.5382681) q[3];
sx q[3];
rz(-0.28964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6203674) q[2];
sx q[2];
rz(-1.8769033) q[2];
sx q[2];
rz(-0.59086665) q[2];
rz(-1.8477731) q[3];
sx q[3];
rz(-1.7655244) q[3];
sx q[3];
rz(-2.6025492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6227601) q[0];
sx q[0];
rz(-0.22587124) q[0];
sx q[0];
rz(-0.87509218) q[0];
rz(0.092763364) q[1];
sx q[1];
rz(-1.0567254) q[1];
sx q[1];
rz(1.0088395) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.09966775) q[0];
sx q[0];
rz(-2.0001051) q[0];
sx q[0];
rz(-0.03089182) q[0];
x q[1];
rz(-0.039142056) q[2];
sx q[2];
rz(-0.58696568) q[2];
sx q[2];
rz(-2.0330737) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0752901) q[1];
sx q[1];
rz(-1.7199893) q[1];
sx q[1];
rz(2.2239524) q[1];
rz(-0.013539516) q[3];
sx q[3];
rz(-1.5126745) q[3];
sx q[3];
rz(0.45765578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86866093) q[2];
sx q[2];
rz(-1.2596143) q[2];
sx q[2];
rz(-1.3966365) q[2];
rz(2.7331288) q[3];
sx q[3];
rz(-2.4558081) q[3];
sx q[3];
rz(-2.7062611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246966) q[0];
sx q[0];
rz(-0.28852391) q[0];
sx q[0];
rz(2.1911279) q[0];
rz(1.2735927) q[1];
sx q[1];
rz(-1.344695) q[1];
sx q[1];
rz(-1.7705932) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4038611) q[0];
sx q[0];
rz(-1.4329684) q[0];
sx q[0];
rz(2.4859684) q[0];
rz(-pi) q[1];
rz(-2.7080405) q[2];
sx q[2];
rz(-1.4656275) q[2];
sx q[2];
rz(0.97135949) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85670722) q[1];
sx q[1];
rz(-2.6213985) q[1];
sx q[1];
rz(-1.3150645) q[1];
x q[2];
rz(2.9955735) q[3];
sx q[3];
rz(-0.64490025) q[3];
sx q[3];
rz(-0.19949808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3889968) q[2];
sx q[2];
rz(-1.7277371) q[2];
sx q[2];
rz(-1.625212) q[2];
rz(0.73836941) q[3];
sx q[3];
rz(-1.6298031) q[3];
sx q[3];
rz(0.71877688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7826638) q[0];
sx q[0];
rz(-1.5072701) q[0];
sx q[0];
rz(1.2231109) q[0];
rz(0.75543985) q[1];
sx q[1];
rz(-2.0016045) q[1];
sx q[1];
rz(3.0208407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8275195) q[0];
sx q[0];
rz(-1.4085007) q[0];
sx q[0];
rz(1.7364962) q[0];
rz(-1.3322387) q[2];
sx q[2];
rz(-1.0285796) q[2];
sx q[2];
rz(1.2921367) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11343458) q[1];
sx q[1];
rz(-1.0453117) q[1];
sx q[1];
rz(-2.3281086) q[1];
rz(-pi) q[2];
rz(3.122943) q[3];
sx q[3];
rz(-1.6183162) q[3];
sx q[3];
rz(1.7619676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52189031) q[2];
sx q[2];
rz(-1.4554687) q[2];
sx q[2];
rz(1.5187368) q[2];
rz(1.8146727) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(2.0324619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0066852) q[0];
sx q[0];
rz(-1.6970072) q[0];
sx q[0];
rz(-2.4476442) q[0];
rz(-0.42849439) q[1];
sx q[1];
rz(-0.65281147) q[1];
sx q[1];
rz(1.8983715) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0611785) q[0];
sx q[0];
rz(-1.616487) q[0];
sx q[0];
rz(-1.7663598) q[0];
rz(-pi) q[1];
rz(0.92840009) q[2];
sx q[2];
rz(-2.8619741) q[2];
sx q[2];
rz(-1.9595264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.32520884) q[1];
sx q[1];
rz(-1.7933791) q[1];
sx q[1];
rz(1.1735664) q[1];
x q[2];
rz(1.6749773) q[3];
sx q[3];
rz(-1.0910463) q[3];
sx q[3];
rz(1.0436077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1474233) q[2];
sx q[2];
rz(-1.3012412) q[2];
sx q[2];
rz(2.6541138) q[2];
rz(-0.068988919) q[3];
sx q[3];
rz(-0.65620771) q[3];
sx q[3];
rz(0.60855734) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92031205) q[0];
sx q[0];
rz(-2.1423036) q[0];
sx q[0];
rz(-2.6699303) q[0];
rz(0.77700514) q[1];
sx q[1];
rz(-1.124137) q[1];
sx q[1];
rz(-1.5583001) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651833) q[0];
sx q[0];
rz(-1.5411955) q[0];
sx q[0];
rz(1.5628003) q[0];
rz(-pi) q[1];
x q[1];
rz(2.110811) q[2];
sx q[2];
rz(-1.2031297) q[2];
sx q[2];
rz(0.45258507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4397227) q[1];
sx q[1];
rz(-1.1367203) q[1];
sx q[1];
rz(-1.3324655) q[1];
rz(-1.8040001) q[3];
sx q[3];
rz(-0.85790715) q[3];
sx q[3];
rz(1.9645129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9787057) q[2];
sx q[2];
rz(-1.4894314) q[2];
sx q[2];
rz(1.9784652) q[2];
rz(-0.85824054) q[3];
sx q[3];
rz(-1.6909928) q[3];
sx q[3];
rz(-2.3024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6100886) q[0];
sx q[0];
rz(-1.4634345) q[0];
sx q[0];
rz(0.32291821) q[0];
rz(-1.3195751) q[1];
sx q[1];
rz(-2.0009305) q[1];
sx q[1];
rz(1.742977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.271287) q[0];
sx q[0];
rz(-1.7338949) q[0];
sx q[0];
rz(2.968279) q[0];
rz(-pi) q[1];
rz(-2.182102) q[2];
sx q[2];
rz(-1.2388907) q[2];
sx q[2];
rz(1.4997945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9515032) q[1];
sx q[1];
rz(-1.5384337) q[1];
sx q[1];
rz(0.75143674) q[1];
x q[2];
rz(0.12945505) q[3];
sx q[3];
rz(-1.6627747) q[3];
sx q[3];
rz(1.5935073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6814331) q[2];
sx q[2];
rz(-2.5761649) q[2];
sx q[2];
rz(2.4088755) q[2];
rz(0.021746246) q[3];
sx q[3];
rz(-1.5214058) q[3];
sx q[3];
rz(-2.963613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1499504) q[0];
sx q[0];
rz(-0.1217148) q[0];
sx q[0];
rz(-0.44418401) q[0];
rz(-1.3581879) q[1];
sx q[1];
rz(-1.5779747) q[1];
sx q[1];
rz(1.1202728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0585384) q[0];
sx q[0];
rz(-2.3069366) q[0];
sx q[0];
rz(3.0025981) q[0];
x q[1];
rz(-1.1728196) q[2];
sx q[2];
rz(-1.1160154) q[2];
sx q[2];
rz(-2.9000521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4746589) q[1];
sx q[1];
rz(-0.28945112) q[1];
sx q[1];
rz(0.54929026) q[1];
rz(2.4230749) q[3];
sx q[3];
rz(-2.5542298) q[3];
sx q[3];
rz(-1.6247251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84919763) q[2];
sx q[2];
rz(-0.85453832) q[2];
sx q[2];
rz(1.9067859) q[2];
rz(2.1644927) q[3];
sx q[3];
rz(-1.8465202) q[3];
sx q[3];
rz(-2.5375514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6848171) q[0];
sx q[0];
rz(-0.026739459) q[0];
sx q[0];
rz(2.7944414) q[0];
rz(-2.8390362) q[1];
sx q[1];
rz(-2.0358993) q[1];
sx q[1];
rz(-3.1355296) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6183229) q[0];
sx q[0];
rz(-1.5875641) q[0];
sx q[0];
rz(1.7527333) q[0];
rz(-pi) q[1];
rz(-1.9337186) q[2];
sx q[2];
rz(-2.4293682) q[2];
sx q[2];
rz(2.1729529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2711941) q[1];
sx q[1];
rz(-0.70576233) q[1];
sx q[1];
rz(-1.0949542) q[1];
x q[2];
rz(1.6177032) q[3];
sx q[3];
rz(-0.79104086) q[3];
sx q[3];
rz(-0.49680809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2431295) q[2];
sx q[2];
rz(-2.1738238) q[2];
sx q[2];
rz(-1.3904307) q[2];
rz(-9/(1*pi)) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(2.0738585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4047563) q[0];
sx q[0];
rz(-2.3911349) q[0];
sx q[0];
rz(2.2682037) q[0];
rz(2.5923173) q[1];
sx q[1];
rz(-2.4875689) q[1];
sx q[1];
rz(-0.80680791) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20414509) q[0];
sx q[0];
rz(-0.68430005) q[0];
sx q[0];
rz(-2.7214976) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15676683) q[2];
sx q[2];
rz(-1.6878551) q[2];
sx q[2];
rz(-2.5396404) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0650821) q[1];
sx q[1];
rz(-1.8190093) q[1];
sx q[1];
rz(-0.11947554) q[1];
x q[2];
rz(-1.1519506) q[3];
sx q[3];
rz(-0.65336334) q[3];
sx q[3];
rz(1.9532596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6640545) q[2];
sx q[2];
rz(-1.8287903) q[2];
sx q[2];
rz(2.14892) q[2];
rz(0.70332518) q[3];
sx q[3];
rz(-1.1329009) q[3];
sx q[3];
rz(-2.4542184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0032046) q[0];
sx q[0];
rz(-0.78835543) q[0];
sx q[0];
rz(2.0302571) q[0];
rz(0.55029021) q[1];
sx q[1];
rz(-2.7688409) q[1];
sx q[1];
rz(0.54584835) q[1];
rz(1.4081804) q[2];
sx q[2];
rz(-2.2615216) q[2];
sx q[2];
rz(0.79590454) q[2];
rz(-2.417682) q[3];
sx q[3];
rz(-2.5528298) q[3];
sx q[3];
rz(-0.34686723) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
