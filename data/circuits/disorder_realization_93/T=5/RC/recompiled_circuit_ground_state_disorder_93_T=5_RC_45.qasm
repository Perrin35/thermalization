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
rz(-0.0184514) q[0];
sx q[0];
rz(-2.9676394) q[0];
rz(1.8040909) q[1];
sx q[1];
rz(-2.2660273) q[1];
sx q[1];
rz(2.8101885) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023663432) q[0];
sx q[0];
rz(-2.8011596) q[0];
sx q[0];
rz(-0.4544008) q[0];
rz(-2.549667) q[2];
sx q[2];
rz(-2.3560696) q[2];
sx q[2];
rz(-0.80410081) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4949499) q[1];
sx q[1];
rz(-0.70716049) q[1];
sx q[1];
rz(-2.8487514) q[1];
rz(-pi) q[2];
rz(2.9865206) q[3];
sx q[3];
rz(-1.6033246) q[3];
sx q[3];
rz(0.28964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5212253) q[2];
sx q[2];
rz(-1.2646893) q[2];
sx q[2];
rz(-0.59086665) q[2];
rz(-1.2938195) q[3];
sx q[3];
rz(-1.3760682) q[3];
sx q[3];
rz(-2.6025492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6227601) q[0];
sx q[0];
rz(-2.9157214) q[0];
sx q[0];
rz(0.87509218) q[0];
rz(3.0488293) q[1];
sx q[1];
rz(-1.0567254) q[1];
sx q[1];
rz(-1.0088395) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025564145) q[0];
sx q[0];
rz(-0.43034986) q[0];
sx q[0];
rz(-1.6381646) q[0];
rz(-pi) q[1];
rz(0.039142056) q[2];
sx q[2];
rz(-0.58696568) q[2];
sx q[2];
rz(2.0330737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0663026) q[1];
sx q[1];
rz(-1.7199893) q[1];
sx q[1];
rz(-2.2239524) q[1];
x q[2];
rz(1.6289234) q[3];
sx q[3];
rz(-1.5572797) q[3];
sx q[3];
rz(-1.112354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.86866093) q[2];
sx q[2];
rz(-1.8819784) q[2];
sx q[2];
rz(1.3966365) q[2];
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
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8946266) q[0];
sx q[0];
rz(-0.28852391) q[0];
sx q[0];
rz(-0.95046473) q[0];
rz(1.8679999) q[1];
sx q[1];
rz(-1.344695) q[1];
sx q[1];
rz(-1.3709995) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7978861) q[0];
sx q[0];
rz(-2.4737353) q[0];
sx q[0];
rz(2.9178828) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7080405) q[2];
sx q[2];
rz(-1.4656275) q[2];
sx q[2];
rz(2.1702332) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6505988) q[1];
sx q[1];
rz(-1.6968602) q[1];
sx q[1];
rz(-2.0768512) q[1];
rz(-pi) q[2];
rz(0.14601918) q[3];
sx q[3];
rz(-0.64490025) q[3];
sx q[3];
rz(-2.9420946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3889968) q[2];
sx q[2];
rz(-1.4138556) q[2];
sx q[2];
rz(-1.625212) q[2];
rz(-0.73836941) q[3];
sx q[3];
rz(-1.5117896) q[3];
sx q[3];
rz(-2.4228158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3589288) q[0];
sx q[0];
rz(-1.6343225) q[0];
sx q[0];
rz(1.2231109) q[0];
rz(-2.3861528) q[1];
sx q[1];
rz(-1.1399882) q[1];
sx q[1];
rz(-3.0208407) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6530891) q[0];
sx q[0];
rz(-2.9101744) q[0];
sx q[0];
rz(2.3525224) q[0];
rz(0.37379548) q[2];
sx q[2];
rz(-2.5540562) q[2];
sx q[2];
rz(0.85174207) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.014552) q[1];
sx q[1];
rz(-2.2071731) q[1];
sx q[1];
rz(-2.4680952) q[1];
x q[2];
rz(1.6183245) q[3];
sx q[3];
rz(-1.589425) q[3];
sx q[3];
rz(-2.9513074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6197023) q[2];
sx q[2];
rz(-1.4554687) q[2];
sx q[2];
rz(1.6228559) q[2];
rz(-1.3269199) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(-1.1091308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0066852) q[0];
sx q[0];
rz(-1.4445855) q[0];
sx q[0];
rz(-0.69394845) q[0];
rz(-0.42849439) q[1];
sx q[1];
rz(-0.65281147) q[1];
sx q[1];
rz(1.8983715) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.87782) q[0];
sx q[0];
rz(-0.20076361) q[0];
sx q[0];
rz(-1.8018886) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2131926) q[2];
sx q[2];
rz(-0.27961857) q[2];
sx q[2];
rz(-1.9595264) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1532395) q[1];
sx q[1];
rz(-1.1838938) q[1];
sx q[1];
rz(0.24068618) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4666154) q[3];
sx q[3];
rz(-2.0505464) q[3];
sx q[3];
rz(1.0436077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1474233) q[2];
sx q[2];
rz(-1.8403515) q[2];
sx q[2];
rz(-2.6541138) q[2];
rz(-0.068988919) q[3];
sx q[3];
rz(-0.65620771) q[3];
sx q[3];
rz(-2.5330353) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92031205) q[0];
sx q[0];
rz(-2.1423036) q[0];
sx q[0];
rz(2.6699303) q[0];
rz(0.77700514) q[1];
sx q[1];
rz(-2.0174556) q[1];
sx q[1];
rz(-1.5832925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37640939) q[0];
sx q[0];
rz(-1.5411955) q[0];
sx q[0];
rz(1.5628003) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92783096) q[2];
sx q[2];
rz(-0.6428679) q[2];
sx q[2];
rz(0.57807482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.70187) q[1];
sx q[1];
rz(-2.0048723) q[1];
sx q[1];
rz(1.3324655) q[1];
rz(-1.3375925) q[3];
sx q[3];
rz(-0.85790715) q[3];
sx q[3];
rz(1.1770798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9787057) q[2];
sx q[2];
rz(-1.4894314) q[2];
sx q[2];
rz(-1.1631274) q[2];
rz(-2.2833521) q[3];
sx q[3];
rz(-1.6909928) q[3];
sx q[3];
rz(2.3024018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.531504) q[0];
sx q[0];
rz(-1.6781582) q[0];
sx q[0];
rz(2.8186744) q[0];
rz(-1.8220176) q[1];
sx q[1];
rz(-1.1406621) q[1];
sx q[1];
rz(-1.3986157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4482898) q[0];
sx q[0];
rz(-0.23742659) q[0];
sx q[0];
rz(-0.7619995) q[0];
x q[1];
rz(2.1115913) q[2];
sx q[2];
rz(-0.68533603) q[2];
sx q[2];
rz(-2.6353996) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7955048) q[1];
sx q[1];
rz(-2.3895956) q[1];
sx q[1];
rz(3.094207) q[1];
x q[2];
rz(-0.12945505) q[3];
sx q[3];
rz(-1.4788179) q[3];
sx q[3];
rz(1.5935073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4601595) q[2];
sx q[2];
rz(-2.5761649) q[2];
sx q[2];
rz(-2.4088755) q[2];
rz(-3.1198464) q[3];
sx q[3];
rz(-1.5214058) q[3];
sx q[3];
rz(0.17797962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1499504) q[0];
sx q[0];
rz(-3.0198779) q[0];
sx q[0];
rz(-0.44418401) q[0];
rz(-1.3581879) q[1];
sx q[1];
rz(-1.5636179) q[1];
sx q[1];
rz(2.0213199) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60591423) q[0];
sx q[0];
rz(-1.6736503) q[0];
sx q[0];
rz(0.82983183) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6539069) q[2];
sx q[2];
rz(-1.2151866) q[2];
sx q[2];
rz(-1.5118701) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0430482) q[1];
sx q[1];
rz(-1.8167059) q[1];
sx q[1];
rz(-1.7250452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4230749) q[3];
sx q[3];
rz(-0.58736283) q[3];
sx q[3];
rz(-1.6247251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84919763) q[2];
sx q[2];
rz(-0.85453832) q[2];
sx q[2];
rz(-1.9067859) q[2];
rz(2.1644927) q[3];
sx q[3];
rz(-1.8465202) q[3];
sx q[3];
rz(-2.5375514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(3.1355296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13841471) q[0];
sx q[0];
rz(-0.18269953) q[0];
sx q[0];
rz(-1.6632141) q[0];
x q[1];
rz(0.89170189) q[2];
sx q[2];
rz(-1.3366586) q[2];
sx q[2];
rz(2.8193223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67601678) q[1];
sx q[1];
rz(-2.1853059) q[1];
sx q[1];
rz(-0.37218233) q[1];
rz(-0.047386332) q[3];
sx q[3];
rz(-2.3607254) q[3];
sx q[3];
rz(-2.5781207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89846316) q[2];
sx q[2];
rz(-2.1738238) q[2];
sx q[2];
rz(1.3904307) q[2];
rz(-0.27680382) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(1.0677342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7368363) q[0];
sx q[0];
rz(-2.3911349) q[0];
sx q[0];
rz(-2.2682037) q[0];
rz(-0.5492754) q[1];
sx q[1];
rz(-0.6540238) q[1];
sx q[1];
rz(-2.3347847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0334369) q[0];
sx q[0];
rz(-1.3100385) q[0];
sx q[0];
rz(-2.5013862) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64555706) q[2];
sx q[2];
rz(-2.9462313) q[2];
sx q[2];
rz(-0.33249172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0650821) q[1];
sx q[1];
rz(-1.3225833) q[1];
sx q[1];
rz(3.0221171) q[1];
rz(0.96050519) q[3];
sx q[3];
rz(-1.3209855) q[3];
sx q[3];
rz(-0.042674289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6640545) q[2];
sx q[2];
rz(-1.8287903) q[2];
sx q[2];
rz(-0.99267268) q[2];
rz(2.4382675) q[3];
sx q[3];
rz(-2.0086918) q[3];
sx q[3];
rz(-2.4542184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-2.5913024) q[1];
sx q[1];
rz(-2.7688409) q[1];
sx q[1];
rz(0.54584835) q[1];
rz(2.9481683) q[2];
sx q[2];
rz(-0.70653595) q[2];
sx q[2];
rz(0.54386137) q[2];
rz(1.1543858) q[3];
sx q[3];
rz(-1.9999123) q[3];
sx q[3];
rz(-2.6725389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
