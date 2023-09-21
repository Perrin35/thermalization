OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(4.1424799) q[0];
sx q[0];
rz(9.6371798) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(-1.2815055) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0467487) q[0];
sx q[0];
rz(-1.9108512) q[0];
sx q[0];
rz(0.21324555) q[0];
rz(-3.1060495) q[2];
sx q[2];
rz(-2.4487552) q[2];
sx q[2];
rz(-0.37301829) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9077397) q[1];
sx q[1];
rz(-0.45957652) q[1];
sx q[1];
rz(-1.185226) q[1];
x q[2];
rz(1.699502) q[3];
sx q[3];
rz(-1.5081815) q[3];
sx q[3];
rz(2.6821729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.32221258) q[2];
sx q[2];
rz(-3.0266422) q[2];
sx q[2];
rz(-1.6248576) q[2];
rz(-0.24762282) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6973998) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(-0.63931757) q[1];
sx q[1];
rz(-2.1289861) q[1];
sx q[1];
rz(1.4888391) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.981502) q[0];
sx q[0];
rz(-1.0233876) q[0];
sx q[0];
rz(1.3840958) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4376051) q[2];
sx q[2];
rz(-0.96103243) q[2];
sx q[2];
rz(-3.0733382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.78755806) q[1];
sx q[1];
rz(-2.2984142) q[1];
sx q[1];
rz(-2.3393199) q[1];
x q[2];
rz(1.5403455) q[3];
sx q[3];
rz(-1.4342562) q[3];
sx q[3];
rz(0.88324916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8504101) q[2];
sx q[2];
rz(-0.098878421) q[2];
sx q[2];
rz(0.29671159) q[2];
rz(-2.8091649) q[3];
sx q[3];
rz(-1.6231096) q[3];
sx q[3];
rz(1.9077574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333703) q[0];
sx q[0];
rz(-1.2587661) q[0];
sx q[0];
rz(2.6112774) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-1.1214316) q[1];
sx q[1];
rz(-1.189032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9246763) q[0];
sx q[0];
rz(-2.8042256) q[0];
sx q[0];
rz(-0.49305537) q[0];
rz(2.8934115) q[2];
sx q[2];
rz(-1.6096874) q[2];
sx q[2];
rz(1.789202) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6991068) q[1];
sx q[1];
rz(-2.2190296) q[1];
sx q[1];
rz(2.4590765) q[1];
x q[2];
rz(-2.4431908) q[3];
sx q[3];
rz(-1.9365371) q[3];
sx q[3];
rz(1.6549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.6685852) q[2];
sx q[2];
rz(1.2515602) q[2];
rz(2.6873612) q[3];
sx q[3];
rz(-0.43332228) q[3];
sx q[3];
rz(-2.7864578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.3635451) q[0];
sx q[0];
rz(-2.8319034) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(-1.9569424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458762) q[0];
sx q[0];
rz(-1.3564975) q[0];
sx q[0];
rz(1.7488259) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2934154) q[2];
sx q[2];
rz(-2.4978673) q[2];
sx q[2];
rz(1.4959469) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8682755) q[1];
sx q[1];
rz(-1.454121) q[1];
sx q[1];
rz(-0.97881808) q[1];
rz(2.437856) q[3];
sx q[3];
rz(-2.6061213) q[3];
sx q[3];
rz(2.717201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9081395) q[2];
sx q[2];
rz(-2.2557857) q[2];
sx q[2];
rz(-0.63956368) q[2];
rz(0.8574287) q[3];
sx q[3];
rz(-1.144751) q[3];
sx q[3];
rz(-0.48979428) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63067591) q[0];
sx q[0];
rz(-2.5998721) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(-0.10294542) q[1];
sx q[1];
rz(-1.0504477) q[1];
sx q[1];
rz(-0.72881126) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136912) q[0];
sx q[0];
rz(-1.3662845) q[0];
sx q[0];
rz(1.6708899) q[0];
rz(0.16291933) q[2];
sx q[2];
rz(-2.2883121) q[2];
sx q[2];
rz(1.3416854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0923857) q[1];
sx q[1];
rz(-1.8930463) q[1];
sx q[1];
rz(0.21578034) q[1];
rz(-pi) q[2];
rz(-2.1843188) q[3];
sx q[3];
rz(-2.0271218) q[3];
sx q[3];
rz(2.7077655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0430498) q[2];
sx q[2];
rz(-0.35991943) q[2];
sx q[2];
rz(0.70971242) q[2];
rz(2.5085311) q[3];
sx q[3];
rz(-1.9934306) q[3];
sx q[3];
rz(-3.0098797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1028041) q[0];
sx q[0];
rz(-0.10237256) q[0];
sx q[0];
rz(0.88371712) q[0];
rz(-0.9206413) q[1];
sx q[1];
rz(-2.0795152) q[1];
sx q[1];
rz(0.19764915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716498) q[0];
sx q[0];
rz(-0.93402223) q[0];
sx q[0];
rz(2.3417579) q[0];
rz(-pi) q[1];
rz(1.6502041) q[2];
sx q[2];
rz(-1.4646155) q[2];
sx q[2];
rz(0.90027819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3870961) q[1];
sx q[1];
rz(-1.5281786) q[1];
sx q[1];
rz(0.45511647) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75943767) q[3];
sx q[3];
rz(-1.0348088) q[3];
sx q[3];
rz(-2.939455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.343507) q[2];
sx q[2];
rz(-0.075604288) q[2];
rz(1.4893701) q[3];
sx q[3];
rz(-0.39947709) q[3];
sx q[3];
rz(-2.3125125) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7224715) q[0];
sx q[0];
rz(-2.2561181) q[0];
sx q[0];
rz(0.042073123) q[0];
rz(1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.9721608) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5116611) q[0];
sx q[0];
rz(-2.9637664) q[0];
sx q[0];
rz(1.9566262) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61755007) q[2];
sx q[2];
rz(-0.46717656) q[2];
sx q[2];
rz(3.0385366) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52160727) q[1];
sx q[1];
rz(-1.0343401) q[1];
sx q[1];
rz(1.9636088) q[1];
x q[2];
rz(-0.38798214) q[3];
sx q[3];
rz(-1.6098607) q[3];
sx q[3];
rz(-2.8475873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3183257) q[2];
sx q[2];
rz(-2.0264758) q[2];
sx q[2];
rz(-3.1271093) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.2967916) q[3];
sx q[3];
rz(0.55317944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13977215) q[0];
sx q[0];
rz(-0.35035366) q[0];
sx q[0];
rz(2.5792504) q[0];
rz(1.4315804) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6477752) q[0];
sx q[0];
rz(-1.6543979) q[0];
sx q[0];
rz(-1.093822) q[0];
rz(-pi) q[1];
x q[1];
rz(0.55191314) q[2];
sx q[2];
rz(-2.1898824) q[2];
sx q[2];
rz(2.1053095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3265877) q[1];
sx q[1];
rz(-1.5851138) q[1];
sx q[1];
rz(-2.307595) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3433414) q[3];
sx q[3];
rz(-1.0700873) q[3];
sx q[3];
rz(-2.8727369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3056425) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(-0.28253728) q[2];
rz(2.1841168) q[3];
sx q[3];
rz(-1.3922393) q[3];
sx q[3];
rz(-2.6628475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95865059) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.4022934) q[0];
rz(2.4422586) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(1.8797849) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9937447) q[0];
sx q[0];
rz(-1.9980668) q[0];
sx q[0];
rz(-1.2742313) q[0];
x q[1];
rz(-0.26913531) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(-3.105643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3952179) q[1];
sx q[1];
rz(-2.4002541) q[1];
sx q[1];
rz(-2.7625601) q[1];
rz(2.6120139) q[3];
sx q[3];
rz(-1.7486608) q[3];
sx q[3];
rz(1.1416658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4385779) q[2];
sx q[2];
rz(-1.6377178) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.7137182) q[3];
sx q[3];
rz(-0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2465729) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.1599468) q[0];
rz(-3.0874522) q[1];
sx q[1];
rz(-1.4777947) q[1];
sx q[1];
rz(2.0711526) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69660891) q[0];
sx q[0];
rz(-2.3009926) q[0];
sx q[0];
rz(-0.22459774) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3355867) q[2];
sx q[2];
rz(-1.5431879) q[2];
sx q[2];
rz(-2.1105786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.479446) q[1];
sx q[1];
rz(-1.3227191) q[1];
sx q[1];
rz(-0.76920385) q[1];
rz(-pi) q[2];
rz(-1.7461734) q[3];
sx q[3];
rz(-1.326582) q[3];
sx q[3];
rz(2.0004686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3728309) q[2];
sx q[2];
rz(-2.5925345) q[2];
sx q[2];
rz(-1.6147511) q[2];
rz(-2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(-2.4894864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6282745) q[0];
sx q[0];
rz(-1.6074629) q[0];
sx q[0];
rz(0.71832023) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(1.6127916) q[2];
sx q[2];
rz(-2.4091346) q[2];
sx q[2];
rz(-0.11382881) q[2];
rz(-2.0680239) q[3];
sx q[3];
rz(-2.0484925) q[3];
sx q[3];
rz(-1.9852553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
