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
rz(-1.2920657) q[0];
sx q[0];
rz(-2.3176365) q[0];
sx q[0];
rz(0.92744654) q[0];
rz(2.0916341) q[1];
sx q[1];
rz(-0.97133049) q[1];
sx q[1];
rz(-0.8144905) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096734418) q[0];
sx q[0];
rz(-1.1219026) q[0];
sx q[0];
rz(0.4273703) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61887069) q[2];
sx q[2];
rz(-2.2437988) q[2];
sx q[2];
rz(-0.25970632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8100639) q[1];
sx q[1];
rz(-2.0671856) q[1];
sx q[1];
rz(-2.1138825) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4708704) q[3];
sx q[3];
rz(-2.1938087) q[3];
sx q[3];
rz(-2.5581806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.77259511) q[2];
sx q[2];
rz(-1.7721704) q[2];
sx q[2];
rz(-2.7102846) q[2];
rz(1.4946233) q[3];
sx q[3];
rz(-1.5458115) q[3];
sx q[3];
rz(-2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876708) q[0];
sx q[0];
rz(-0.19667721) q[0];
sx q[0];
rz(-1.824463) q[0];
rz(-1.0493086) q[1];
sx q[1];
rz(-2.6056555) q[1];
sx q[1];
rz(-0.61872331) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65370377) q[0];
sx q[0];
rz(-1.426322) q[0];
sx q[0];
rz(-2.5313951) q[0];
x q[1];
rz(-2.2631386) q[2];
sx q[2];
rz(-1.8900423) q[2];
sx q[2];
rz(-2.648166) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2545027) q[1];
sx q[1];
rz(-1.2225593) q[1];
sx q[1];
rz(-3.00249) q[1];
rz(0.57024184) q[3];
sx q[3];
rz(-0.61755731) q[3];
sx q[3];
rz(2.9887426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.416136) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(0.49187342) q[2];
rz(-1.5623931) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(0.57945848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1270776) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(2.8614817) q[0];
rz(3.0666871) q[1];
sx q[1];
rz(-2.7066051) q[1];
sx q[1];
rz(-1.7710955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6579191) q[0];
sx q[0];
rz(-0.83951211) q[0];
sx q[0];
rz(0.28261225) q[0];
rz(-2.5016194) q[2];
sx q[2];
rz(-2.1378064) q[2];
sx q[2];
rz(2.8550573) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31077267) q[1];
sx q[1];
rz(-1.8337063) q[1];
sx q[1];
rz(0.57057686) q[1];
x q[2];
rz(0.13223572) q[3];
sx q[3];
rz(-1.5521129) q[3];
sx q[3];
rz(-1.4052237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8568153) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(2.3397297) q[2];
rz(-0.91184584) q[3];
sx q[3];
rz(-0.82404476) q[3];
sx q[3];
rz(2.1626933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3367679) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(2.5208852) q[0];
rz(2.8630818) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(-1.0034358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706948) q[0];
sx q[0];
rz(-1.2549632) q[0];
sx q[0];
rz(-1.210808) q[0];
rz(-1.7264257) q[2];
sx q[2];
rz(-2.7162726) q[2];
sx q[2];
rz(2.164054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4651854) q[1];
sx q[1];
rz(-0.27449755) q[1];
sx q[1];
rz(-2.9390915) q[1];
rz(-2.5430319) q[3];
sx q[3];
rz(-1.374647) q[3];
sx q[3];
rz(0.3071827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8371381) q[2];
sx q[2];
rz(-2.0292323) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(-2.0153913) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(-0.57536212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758301) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(-0.78731147) q[0];
rz(0.85834223) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(-2.0599005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3560548) q[0];
sx q[0];
rz(-0.83844664) q[0];
sx q[0];
rz(2.4466956) q[0];
x q[1];
rz(2.0400286) q[2];
sx q[2];
rz(-1.1936371) q[2];
sx q[2];
rz(0.18845972) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8199205) q[1];
sx q[1];
rz(-1.5770438) q[1];
sx q[1];
rz(-2.856918) q[1];
x q[2];
rz(2.1882964) q[3];
sx q[3];
rz(-2.3205608) q[3];
sx q[3];
rz(1.2987069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(2.2990189) q[2];
rz(2.8652371) q[3];
sx q[3];
rz(-2.1601951) q[3];
sx q[3];
rz(1.3493376) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2862947) q[0];
sx q[0];
rz(-1.8305625) q[0];
sx q[0];
rz(-1.9516161) q[0];
rz(1.2225993) q[1];
sx q[1];
rz(-1.2346376) q[1];
sx q[1];
rz(-0.58293265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806429) q[0];
sx q[0];
rz(-1.1034729) q[0];
sx q[0];
rz(1.1790465) q[0];
x q[1];
rz(-1.5491484) q[2];
sx q[2];
rz(-0.57367838) q[2];
sx q[2];
rz(1.4230185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5546023) q[1];
sx q[1];
rz(-1.490056) q[1];
sx q[1];
rz(-2.9029773) q[1];
rz(-0.56598466) q[3];
sx q[3];
rz(-0.71403394) q[3];
sx q[3];
rz(2.4772754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0588093) q[2];
sx q[2];
rz(-1.8913816) q[2];
sx q[2];
rz(-2.1018551) q[2];
rz(2.5522363) q[3];
sx q[3];
rz(-1.6683234) q[3];
sx q[3];
rz(-2.1350433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7650918) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.804922) q[0];
rz(2.916015) q[1];
sx q[1];
rz(-2.070319) q[1];
sx q[1];
rz(-0.76990661) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5107489) q[0];
sx q[0];
rz(-1.1026369) q[0];
sx q[0];
rz(1.3965142) q[0];
rz(-pi) q[1];
rz(0.82775292) q[2];
sx q[2];
rz(-1.079664) q[2];
sx q[2];
rz(-2.998005) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6698031) q[1];
sx q[1];
rz(-0.21179971) q[1];
sx q[1];
rz(0.36630156) q[1];
rz(0.60524551) q[3];
sx q[3];
rz(-1.5063933) q[3];
sx q[3];
rz(-1.6689672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5169107) q[2];
sx q[2];
rz(-1.7485488) q[2];
sx q[2];
rz(-2.1262271) q[2];
rz(2.9376302) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(1.7502194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90149752) q[0];
sx q[0];
rz(-2.5101341) q[0];
sx q[0];
rz(-2.748306) q[0];
rz(-0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(3.111305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05370985) q[0];
sx q[0];
rz(-1.6018512) q[0];
sx q[0];
rz(0.24760084) q[0];
rz(-pi) q[1];
x q[1];
rz(1.302839) q[2];
sx q[2];
rz(-0.7795802) q[2];
sx q[2];
rz(-1.5991581) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5133678) q[1];
sx q[1];
rz(-2.1163883) q[1];
sx q[1];
rz(1.5978769) q[1];
x q[2];
rz(0.26408402) q[3];
sx q[3];
rz(-1.4435766) q[3];
sx q[3];
rz(0.41632392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0344737) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(2.0195473) q[2];
rz(-2.0864887) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(-1.25157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17290641) q[0];
sx q[0];
rz(-0.98772573) q[0];
sx q[0];
rz(-0.77447844) q[0];
rz(0.6340181) q[1];
sx q[1];
rz(-1.0301544) q[1];
sx q[1];
rz(-1.0672306) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2174847) q[0];
sx q[0];
rz(-1.587968) q[0];
sx q[0];
rz(0.73826684) q[0];
rz(-pi) q[1];
rz(-1.2634694) q[2];
sx q[2];
rz(-1.6174223) q[2];
sx q[2];
rz(1.4263375) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9561056) q[1];
sx q[1];
rz(-1.2602578) q[1];
sx q[1];
rz(-2.8794672) q[1];
x q[2];
rz(-2.0499633) q[3];
sx q[3];
rz(-1.7288343) q[3];
sx q[3];
rz(0.8534067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.60712236) q[2];
sx q[2];
rz(-0.12132135) q[2];
sx q[2];
rz(1.2657451) q[2];
rz(-0.040146116) q[3];
sx q[3];
rz(-2.7669192) q[3];
sx q[3];
rz(0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48113111) q[0];
sx q[0];
rz(-2.2534695) q[0];
sx q[0];
rz(-1.2598502) q[0];
rz(-1.4866359) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(0.0025509603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5887506) q[0];
sx q[0];
rz(-1.6903983) q[0];
sx q[0];
rz(-1.4445067) q[0];
rz(-pi) q[1];
rz(1.6348636) q[2];
sx q[2];
rz(-1.5803711) q[2];
sx q[2];
rz(-2.9970882) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8067183) q[1];
sx q[1];
rz(-2.2069227) q[1];
sx q[1];
rz(1.0735372) q[1];
rz(-pi) q[2];
rz(-3.1151089) q[3];
sx q[3];
rz(-1.1364352) q[3];
sx q[3];
rz(-2.6514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.586414) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(-2.4671538) q[2];
rz(-1.1174348) q[3];
sx q[3];
rz(-1.0267461) q[3];
sx q[3];
rz(-2.8286772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085070327) q[0];
sx q[0];
rz(-1.1310348) q[0];
sx q[0];
rz(-0.52666589) q[0];
rz(0.3479192) q[1];
sx q[1];
rz(-1.9031453) q[1];
sx q[1];
rz(-1.3487945) q[1];
rz(-0.34571503) q[2];
sx q[2];
rz(-1.7713265) q[2];
sx q[2];
rz(-2.3878018) q[2];
rz(-1.4244637) q[3];
sx q[3];
rz(-1.5531333) q[3];
sx q[3];
rz(0.2266758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
