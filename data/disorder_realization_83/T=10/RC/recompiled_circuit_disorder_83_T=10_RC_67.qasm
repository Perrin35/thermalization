OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(2.3798556) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640472) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(1.8403948) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6526297) q[2];
sx q[2];
rz(-1.666781) q[2];
sx q[2];
rz(-1.3053615) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.57063369) q[1];
sx q[1];
rz(-1.1753517) q[1];
sx q[1];
rz(2.8407437) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9348014) q[3];
sx q[3];
rz(-1.6379106) q[3];
sx q[3];
rz(1.8458927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(0.14262959) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(-2.9108677) q[0];
rz(1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-0.040963106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36386585) q[0];
sx q[0];
rz(-1.7330609) q[0];
sx q[0];
rz(-1.3217539) q[0];
x q[1];
rz(2.8616521) q[2];
sx q[2];
rz(-1.8676114) q[2];
sx q[2];
rz(2.9505626) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37521024) q[1];
sx q[1];
rz(-2.3751405) q[1];
sx q[1];
rz(0.68739989) q[1];
rz(-pi) q[2];
rz(1.1702234) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(-2.7924201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.42999643) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-0.45021737) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-2.4643262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8375741) q[0];
sx q[0];
rz(-1.7390334) q[0];
sx q[0];
rz(-1.6234682) q[0];
x q[1];
rz(1.5718939) q[2];
sx q[2];
rz(-1.5771616) q[2];
sx q[2];
rz(2.9874143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2689506) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(-2.2058669) q[1];
rz(-2.4694091) q[3];
sx q[3];
rz(-0.58478343) q[3];
sx q[3];
rz(0.72090805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(-0.046860524) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(-2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1470404) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(1.1725918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(-1.4021224) q[0];
x q[1];
rz(1.5863717) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(1.7473999) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51845156) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(1.9407942) q[1];
rz(-pi) q[2];
x q[2];
rz(1.138932) q[3];
sx q[3];
rz(-1.6383088) q[3];
sx q[3];
rz(2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(2.6546997) q[2];
rz(1.9481109) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(-1.8160965) q[0];
rz(-0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(2.4868734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7692208) q[0];
sx q[0];
rz(-2.0125611) q[0];
sx q[0];
rz(1.1830933) q[0];
rz(-pi) q[1];
rz(1.3833369) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(2.2802441) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7655189) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(-0.54733025) q[1];
rz(-2.7563285) q[3];
sx q[3];
rz(-1.6574727) q[3];
sx q[3];
rz(3.1298266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(1.3202745) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.6437795) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70443557) q[0];
sx q[0];
rz(-0.91059443) q[0];
sx q[0];
rz(0.066083834) q[0];
rz(0.97700714) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-2.0055111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.3469997) q[1];
sx q[1];
rz(-1.1113249) q[1];
sx q[1];
rz(-0.32230349) q[1];
rz(-1.3214345) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41247535) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(-2.3506929) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7366911) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(0.94876429) q[0];
x q[1];
rz(-2.3979264) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(-0.92168346) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.32614) q[1];
sx q[1];
rz(-2.7683612) q[1];
sx q[1];
rz(1.3532072) q[1];
x q[2];
rz(-2.2579455) q[3];
sx q[3];
rz(-1.2847319) q[3];
sx q[3];
rz(-1.7162232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2699282) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-0.86501914) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4928116) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(0.53721792) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8763037) q[0];
sx q[0];
rz(-0.2067925) q[0];
sx q[0];
rz(0.32949038) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28283624) q[2];
sx q[2];
rz(-0.0054224646) q[2];
sx q[2];
rz(-1.9920414) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4916363) q[1];
sx q[1];
rz(-1.8815814) q[1];
sx q[1];
rz(1.4375356) q[1];
x q[2];
rz(2.5423126) q[3];
sx q[3];
rz(-1.6370956) q[3];
sx q[3];
rz(-1.146572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(2.8059778) q[2];
rz(0.32661682) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(2.0876419) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(0.98186791) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4440585) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(2.8696637) q[0];
x q[1];
rz(0.1083381) q[2];
sx q[2];
rz(-2.1370558) q[2];
sx q[2];
rz(-2.1926751) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8874444) q[1];
sx q[1];
rz(-1.6283298) q[1];
sx q[1];
rz(1.7565586) q[1];
rz(-pi) q[2];
rz(-2.8544159) q[3];
sx q[3];
rz(-2.4099726) q[3];
sx q[3];
rz(-1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(0.99739933) q[2];
rz(2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2845594) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(2.6249264) q[0];
rz(-2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(2.8732252) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29753387) q[0];
sx q[0];
rz(-0.8807655) q[0];
sx q[0];
rz(0.78155078) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95558138) q[2];
sx q[2];
rz(-2.4590543) q[2];
sx q[2];
rz(-2.5126484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4587935) q[1];
sx q[1];
rz(-1.8384215) q[1];
sx q[1];
rz(-2.5330179) q[1];
rz(-2.1488232) q[3];
sx q[3];
rz(-1.7939995) q[3];
sx q[3];
rz(-2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-2.3582017) q[2];
sx q[2];
rz(-0.56383413) q[2];
rz(-1.1901723) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47678369) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(-2.6396991) q[2];
sx q[2];
rz(-2.5645651) q[2];
sx q[2];
rz(-1.2679451) q[2];
rz(-1.7715122) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
