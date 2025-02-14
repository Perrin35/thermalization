OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.326958) q[0];
sx q[0];
rz(-3.1231413) q[0];
sx q[0];
rz(2.9676394) q[0];
rz(1.8040909) q[1];
sx q[1];
rz(-2.2660273) q[1];
sx q[1];
rz(2.8101885) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1179292) q[0];
sx q[0];
rz(-2.8011596) q[0];
sx q[0];
rz(-2.6871919) q[0];
rz(-pi) q[1];
rz(-2.4487804) q[2];
sx q[2];
rz(-1.1651784) q[2];
sx q[2];
rz(-1.2104397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64664272) q[1];
sx q[1];
rz(-2.4344322) q[1];
sx q[1];
rz(-2.8487514) q[1];
rz(1.5378733) q[3];
sx q[3];
rz(-1.7257856) q[3];
sx q[3];
rz(1.2862358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5212253) q[2];
sx q[2];
rz(-1.2646893) q[2];
sx q[2];
rz(2.550726) q[2];
rz(-1.8477731) q[3];
sx q[3];
rz(-1.7655244) q[3];
sx q[3];
rz(-2.6025492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5188326) q[0];
sx q[0];
rz(-0.22587124) q[0];
sx q[0];
rz(2.2665005) q[0];
rz(-3.0488293) q[1];
sx q[1];
rz(-2.0848672) q[1];
sx q[1];
rz(2.1327532) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4839904) q[0];
sx q[0];
rz(-1.5988841) q[0];
sx q[0];
rz(1.1413069) q[0];
rz(0.58661239) q[2];
sx q[2];
rz(-1.5924708) q[2];
sx q[2];
rz(2.7119111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61774594) q[1];
sx q[1];
rz(-0.92611317) q[1];
sx q[1];
rz(-2.9545386) q[1];
x q[2];
rz(-3.1280531) q[3];
sx q[3];
rz(-1.6289181) q[3];
sx q[3];
rz(0.45765578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.86866093) q[2];
sx q[2];
rz(-1.8819784) q[2];
sx q[2];
rz(-1.3966365) q[2];
rz(-0.40846387) q[3];
sx q[3];
rz(-2.4558081) q[3];
sx q[3];
rz(-2.7062611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8946266) q[0];
sx q[0];
rz(-2.8530687) q[0];
sx q[0];
rz(-0.95046473) q[0];
rz(-1.8679999) q[1];
sx q[1];
rz(-1.7968977) q[1];
sx q[1];
rz(1.7705932) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34370658) q[0];
sx q[0];
rz(-2.4737353) q[0];
sx q[0];
rz(2.9178828) q[0];
rz(-1.6865963) q[2];
sx q[2];
rz(-1.1398004) q[2];
sx q[2];
rz(2.4935982) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1493686) q[1];
sx q[1];
rz(-1.0691339) q[1];
sx q[1];
rz(-0.14389529) q[1];
rz(-1.6798065) q[3];
sx q[3];
rz(-0.9338769) q[3];
sx q[3];
rz(-3.1240841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7525959) q[2];
sx q[2];
rz(-1.4138556) q[2];
sx q[2];
rz(-1.5163806) q[2];
rz(2.4032232) q[3];
sx q[3];
rz(-1.5117896) q[3];
sx q[3];
rz(-2.4228158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.3589288) q[0];
sx q[0];
rz(-1.6343225) q[0];
sx q[0];
rz(-1.9184817) q[0];
rz(-0.75543985) q[1];
sx q[1];
rz(-1.1399882) q[1];
sx q[1];
rz(-0.12075195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31407315) q[0];
sx q[0];
rz(-1.733092) q[0];
sx q[0];
rz(1.4050964) q[0];
rz(1.809354) q[2];
sx q[2];
rz(-1.0285796) q[2];
sx q[2];
rz(1.2921367) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1270406) q[1];
sx q[1];
rz(-0.93441957) q[1];
sx q[1];
rz(-2.4680952) q[1];
rz(-1.5232682) q[3];
sx q[3];
rz(-1.589425) q[3];
sx q[3];
rz(-2.9513074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6197023) q[2];
sx q[2];
rz(-1.4554687) q[2];
sx q[2];
rz(1.5187368) q[2];
rz(1.8146727) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(-1.1091308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0066852) q[0];
sx q[0];
rz(-1.6970072) q[0];
sx q[0];
rz(2.4476442) q[0];
rz(0.42849439) q[1];
sx q[1];
rz(-0.65281147) q[1];
sx q[1];
rz(1.2432212) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4994298) q[0];
sx q[0];
rz(-1.7661531) q[0];
sx q[0];
rz(-3.0950154) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9712293) q[2];
sx q[2];
rz(-1.3479832) q[2];
sx q[2];
rz(0.52056584) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.32520884) q[1];
sx q[1];
rz(-1.3482136) q[1];
sx q[1];
rz(1.9680262) q[1];
x q[2];
rz(2.9443191) q[3];
sx q[3];
rz(-2.6515238) q[3];
sx q[3];
rz(2.3207612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99416939) q[2];
sx q[2];
rz(-1.3012412) q[2];
sx q[2];
rz(-2.6541138) q[2];
rz(0.068988919) q[3];
sx q[3];
rz(-2.4853849) q[3];
sx q[3];
rz(-2.5330353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92031205) q[0];
sx q[0];
rz(-0.99928904) q[0];
sx q[0];
rz(2.6699303) q[0];
rz(-2.3645875) q[1];
sx q[1];
rz(-2.0174556) q[1];
sx q[1];
rz(1.5583001) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5013103) q[0];
sx q[0];
rz(-3.1109312) q[0];
sx q[0];
rz(2.8778381) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0307817) q[2];
sx q[2];
rz(-1.938463) q[2];
sx q[2];
rz(0.45258507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9635545) q[1];
sx q[1];
rz(-0.49152546) q[1];
sx q[1];
rz(0.47103931) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26119097) q[3];
sx q[3];
rz(-2.3979366) q[3];
sx q[3];
rz(-2.3128775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16288699) q[2];
sx q[2];
rz(-1.6521613) q[2];
sx q[2];
rz(1.1631274) q[2];
rz(0.85824054) q[3];
sx q[3];
rz(-1.6909928) q[3];
sx q[3];
rz(2.3024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.531504) q[0];
sx q[0];
rz(-1.6781582) q[0];
sx q[0];
rz(-0.32291821) q[0];
rz(-1.8220176) q[1];
sx q[1];
rz(-2.0009305) q[1];
sx q[1];
rz(1.3986157) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4482898) q[0];
sx q[0];
rz(-0.23742659) q[0];
sx q[0];
rz(0.7619995) q[0];
rz(-pi) q[1];
rz(0.9594907) q[2];
sx q[2];
rz(-1.2388907) q[2];
sx q[2];
rz(1.4997945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4109282) q[1];
sx q[1];
rz(-0.81984869) q[1];
sx q[1];
rz(1.6150722) q[1];
rz(-1.4780462) q[3];
sx q[3];
rz(-1.6997011) q[3];
sx q[3];
rz(-0.034667548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4601595) q[2];
sx q[2];
rz(-0.56542772) q[2];
sx q[2];
rz(0.73271712) q[2];
rz(3.1198464) q[3];
sx q[3];
rz(-1.5214058) q[3];
sx q[3];
rz(-0.17797962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9916423) q[0];
sx q[0];
rz(-3.0198779) q[0];
sx q[0];
rz(-0.44418401) q[0];
rz(1.7834047) q[1];
sx q[1];
rz(-1.5779747) q[1];
sx q[1];
rz(-2.0213199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5356784) q[0];
sx q[0];
rz(-1.6736503) q[0];
sx q[0];
rz(-0.82983183) q[0];
rz(-pi) q[1];
rz(0.6702126) q[2];
sx q[2];
rz(-0.59500098) q[2];
sx q[2];
rz(-2.6196644) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43441916) q[1];
sx q[1];
rz(-1.4212233) q[1];
sx q[1];
rz(-2.892848) q[1];
rz(-pi) q[2];
rz(-0.71851774) q[3];
sx q[3];
rz(-2.5542298) q[3];
sx q[3];
rz(-1.6247251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.292395) q[2];
sx q[2];
rz(-0.85453832) q[2];
sx q[2];
rz(-1.9067859) q[2];
rz(-0.97709996) q[3];
sx q[3];
rz(-1.2950725) q[3];
sx q[3];
rz(-0.60404122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6848171) q[0];
sx q[0];
rz(-0.026739459) q[0];
sx q[0];
rz(0.34715125) q[0];
rz(0.30255643) q[1];
sx q[1];
rz(-1.1056933) q[1];
sx q[1];
rz(-0.0060630719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232698) q[0];
sx q[0];
rz(-1.5540285) q[0];
sx q[0];
rz(-1.7527333) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2498908) q[2];
sx q[2];
rz(-1.804934) q[2];
sx q[2];
rz(-0.32227031) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2711941) q[1];
sx q[1];
rz(-0.70576233) q[1];
sx q[1];
rz(-2.0466385) q[1];
x q[2];
rz(-0.047386332) q[3];
sx q[3];
rz(-2.3607254) q[3];
sx q[3];
rz(0.56347195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2431295) q[2];
sx q[2];
rz(-2.1738238) q[2];
sx q[2];
rz(-1.3904307) q[2];
rz(9/(1*pi)) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(1.0677342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4047563) q[0];
sx q[0];
rz(-0.75045776) q[0];
sx q[0];
rz(0.87338895) q[0];
rz(0.5492754) q[1];
sx q[1];
rz(-0.6540238) q[1];
sx q[1];
rz(2.3347847) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20414509) q[0];
sx q[0];
rz(-0.68430005) q[0];
sx q[0];
rz(-0.42009507) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64555706) q[2];
sx q[2];
rz(-0.19536138) q[2];
sx q[2];
rz(-2.8091009) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0765105) q[1];
sx q[1];
rz(-1.3225833) q[1];
sx q[1];
rz(3.0221171) q[1];
x q[2];
rz(-0.30183123) q[3];
sx q[3];
rz(-2.1595397) q[3];
sx q[3];
rz(1.6993239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6640545) q[2];
sx q[2];
rz(-1.8287903) q[2];
sx q[2];
rz(2.14892) q[2];
rz(2.4382675) q[3];
sx q[3];
rz(-1.1329009) q[3];
sx q[3];
rz(-0.68737427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.9481683) q[2];
sx q[2];
rz(-2.4350567) q[2];
sx q[2];
rz(-2.5977313) q[2];
rz(-2.6777018) q[3];
sx q[3];
rz(-1.9474467) q[3];
sx q[3];
rz(1.8578702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
