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
rz(-1.0128385) q[0];
sx q[0];
rz(-1.9586118) q[0];
sx q[0];
rz(2.8943789) q[0];
rz(-1.5850868) q[1];
sx q[1];
rz(-2.1272008) q[1];
sx q[1];
rz(-2.1870764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92468372) q[0];
sx q[0];
rz(-1.9449073) q[0];
sx q[0];
rz(3.1103136) q[0];
rz(-0.45329161) q[2];
sx q[2];
rz(-1.5201836) q[2];
sx q[2];
rz(-2.5091189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7407511) q[1];
sx q[1];
rz(-2.2993326) q[1];
sx q[1];
rz(0.42899112) q[1];
rz(-pi) q[2];
rz(-1.263721) q[3];
sx q[3];
rz(-0.82160866) q[3];
sx q[3];
rz(1.6916569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9556094) q[2];
sx q[2];
rz(-0.60638967) q[2];
sx q[2];
rz(0.28140226) q[2];
rz(1.9143117) q[3];
sx q[3];
rz(-1.809092) q[3];
sx q[3];
rz(-1.9277771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1948552) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(0.52014822) q[1];
sx q[1];
rz(-1.9303493) q[1];
sx q[1];
rz(-2.928226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0555455) q[0];
sx q[0];
rz(-2.9790857) q[0];
sx q[0];
rz(1.1799728) q[0];
rz(-pi) q[1];
rz(2.2198898) q[2];
sx q[2];
rz(-1.1496759) q[2];
sx q[2];
rz(0.77563167) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0736463) q[1];
sx q[1];
rz(-1.6151307) q[1];
sx q[1];
rz(1.4774051) q[1];
rz(-pi) q[2];
rz(1.0454093) q[3];
sx q[3];
rz(-1.3297218) q[3];
sx q[3];
rz(1.2607393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0352036) q[2];
sx q[2];
rz(-2.8530402) q[2];
sx q[2];
rz(-0.2717379) q[2];
rz(-2.4954097) q[3];
sx q[3];
rz(-1.313442) q[3];
sx q[3];
rz(1.2077695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81994098) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(-0.2680378) q[0];
rz(2.9053814) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(1.4150298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38430308) q[0];
sx q[0];
rz(-1.4688604) q[0];
sx q[0];
rz(1.524964) q[0];
x q[1];
rz(0.36249749) q[2];
sx q[2];
rz(-1.7300743) q[2];
sx q[2];
rz(-0.25937072) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27576646) q[1];
sx q[1];
rz(-1.2938465) q[1];
sx q[1];
rz(-2.6867178) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67309849) q[3];
sx q[3];
rz(-1.7894701) q[3];
sx q[3];
rz(2.297087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1872824) q[2];
sx q[2];
rz(-2.8128746) q[2];
sx q[2];
rz(-1.7598565) q[2];
rz(2.7749744) q[3];
sx q[3];
rz(-1.4724052) q[3];
sx q[3];
rz(3.0439175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919375) q[0];
sx q[0];
rz(-0.20579919) q[0];
sx q[0];
rz(0.016121443) q[0];
rz(-1.6751809) q[1];
sx q[1];
rz(-1.6203531) q[1];
sx q[1];
rz(-0.88502562) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9217236) q[0];
sx q[0];
rz(-1.1678639) q[0];
sx q[0];
rz(-1.3833439) q[0];
rz(-0.22815223) q[2];
sx q[2];
rz(-1.2127664) q[2];
sx q[2];
rz(-1.4756853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6975034) q[1];
sx q[1];
rz(-0.61629811) q[1];
sx q[1];
rz(0.69774903) q[1];
x q[2];
rz(1.0543941) q[3];
sx q[3];
rz(-1.9266911) q[3];
sx q[3];
rz(-3.0404224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57899388) q[2];
sx q[2];
rz(-2.6851974) q[2];
sx q[2];
rz(-1.5835416) q[2];
rz(-1.7700178) q[3];
sx q[3];
rz(-1.2746425) q[3];
sx q[3];
rz(-2.9461327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850995) q[0];
sx q[0];
rz(-1.766196) q[0];
sx q[0];
rz(2.9920355) q[0];
rz(2.0984446) q[1];
sx q[1];
rz(-2.1688192) q[1];
sx q[1];
rz(-0.8358039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59398932) q[0];
sx q[0];
rz(-0.41446742) q[0];
sx q[0];
rz(-2.2311121) q[0];
x q[1];
rz(2.2611375) q[2];
sx q[2];
rz(-1.5635075) q[2];
sx q[2];
rz(1.7178423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52302915) q[1];
sx q[1];
rz(-1.5715741) q[1];
sx q[1];
rz(-1.7788299) q[1];
rz(-0.94344027) q[3];
sx q[3];
rz(-2.8356878) q[3];
sx q[3];
rz(2.4750575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.122637) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(2.7254851) q[2];
rz(0.082402669) q[3];
sx q[3];
rz(-0.96937886) q[3];
sx q[3];
rz(-2.940322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0788954) q[0];
sx q[0];
rz(-2.1128928) q[0];
sx q[0];
rz(-1.0908352) q[0];
rz(1.1385607) q[1];
sx q[1];
rz(-1.1999612) q[1];
sx q[1];
rz(0.60637766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17601062) q[0];
sx q[0];
rz(-2.2041498) q[0];
sx q[0];
rz(1.2562344) q[0];
x q[1];
rz(-0.44754812) q[2];
sx q[2];
rz(-2.2879763) q[2];
sx q[2];
rz(-2.5765975) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55058897) q[1];
sx q[1];
rz(-0.85500756) q[1];
sx q[1];
rz(-0.71859931) q[1];
rz(-1.6929469) q[3];
sx q[3];
rz(-0.85616099) q[3];
sx q[3];
rz(-0.84343101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31327569) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(-3.0307148) q[2];
rz(-1.5634792) q[3];
sx q[3];
rz(-2.9302247) q[3];
sx q[3];
rz(1.8160688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.204708) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(-2.4844266) q[0];
rz(-1.8183297) q[1];
sx q[1];
rz(-2.3902049) q[1];
sx q[1];
rz(2.0164067) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96699636) q[0];
sx q[0];
rz(-1.5558934) q[0];
sx q[0];
rz(3.1257939) q[0];
rz(-1.6680587) q[2];
sx q[2];
rz(-1.8832369) q[2];
sx q[2];
rz(0.51431235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70874309) q[1];
sx q[1];
rz(-1.9545991) q[1];
sx q[1];
rz(-2.6118181) q[1];
rz(-pi) q[2];
rz(2.7636307) q[3];
sx q[3];
rz(-2.4306647) q[3];
sx q[3];
rz(0.42676778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40921673) q[2];
sx q[2];
rz(-1.1820619) q[2];
sx q[2];
rz(-0.96119514) q[2];
rz(0.0088648908) q[3];
sx q[3];
rz(-0.88518393) q[3];
sx q[3];
rz(-1.9302906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0601785) q[0];
sx q[0];
rz(-0.1135122) q[0];
sx q[0];
rz(-2.7139582) q[0];
rz(-2.0291406) q[1];
sx q[1];
rz(-1.8266269) q[1];
sx q[1];
rz(-2.1449259) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8450981) q[0];
sx q[0];
rz(-2.0581492) q[0];
sx q[0];
rz(-2.4250406) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8585343) q[2];
sx q[2];
rz(-0.78327228) q[2];
sx q[2];
rz(-2.9976792) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.06293077) q[1];
sx q[1];
rz(-1.0109631) q[1];
sx q[1];
rz(2.8275851) q[1];
x q[2];
rz(2.8248722) q[3];
sx q[3];
rz(-2.5473928) q[3];
sx q[3];
rz(0.34689981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6796278) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(0.648554) q[2];
rz(2.8280761) q[3];
sx q[3];
rz(-2.2532012) q[3];
sx q[3];
rz(-0.052481767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.40470966) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(-0.96694651) q[0];
rz(2.1838358) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(1.0831833) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082302467) q[0];
sx q[0];
rz(-1.3246312) q[0];
sx q[0];
rz(2.7013426) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71523209) q[2];
sx q[2];
rz(-0.87022793) q[2];
sx q[2];
rz(-0.45203766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31551595) q[1];
sx q[1];
rz(-0.89102645) q[1];
sx q[1];
rz(-0.51425528) q[1];
rz(-pi) q[2];
rz(-3.1414491) q[3];
sx q[3];
rz(-0.4793491) q[3];
sx q[3];
rz(2.0574613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12072418) q[2];
sx q[2];
rz(-0.82896295) q[2];
sx q[2];
rz(-1.4618358) q[2];
rz(-0.36137897) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(-2.1593275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35230377) q[0];
sx q[0];
rz(-0.73065773) q[0];
sx q[0];
rz(0.58255449) q[0];
rz(-1.7763058) q[1];
sx q[1];
rz(-1.0368232) q[1];
sx q[1];
rz(0.22886151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.992975) q[0];
sx q[0];
rz(-1.1753191) q[0];
sx q[0];
rz(1.3735885) q[0];
rz(-2.9624356) q[2];
sx q[2];
rz(-1.3115885) q[2];
sx q[2];
rz(1.1863866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4881674) q[1];
sx q[1];
rz(-1.5060802) q[1];
sx q[1];
rz(1.4117227) q[1];
x q[2];
rz(0.39005847) q[3];
sx q[3];
rz(-1.9771246) q[3];
sx q[3];
rz(-1.4618946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87054306) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(1.7331227) q[2];
rz(0.72077858) q[3];
sx q[3];
rz(-1.9920789) q[3];
sx q[3];
rz(-1.6225947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0558753) q[0];
sx q[0];
rz(-1.4376823) q[0];
sx q[0];
rz(1.8327376) q[0];
rz(1.5695288) q[1];
sx q[1];
rz(-2.5253898) q[1];
sx q[1];
rz(0.22402221) q[1];
rz(1.189255) q[2];
sx q[2];
rz(-2.5425662) q[2];
sx q[2];
rz(-2.7903916) q[2];
rz(-1.4103945) q[3];
sx q[3];
rz(-1.7216202) q[3];
sx q[3];
rz(-1.3826821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
