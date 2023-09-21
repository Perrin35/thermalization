OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(6.531125) q[0];
sx q[0];
rz(8.6046435) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(3.9305384) q[1];
sx q[1];
rz(9.5126704) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9633006) q[0];
sx q[0];
rz(-2.0110197) q[0];
sx q[0];
rz(-0.51142366) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7152039) q[2];
sx q[2];
rz(-2.0585367) q[2];
sx q[2];
rz(-2.3418155) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3327515) q[1];
sx q[1];
rz(-1.0665227) q[1];
sx q[1];
rz(1.5732952) q[1];
x q[2];
rz(1.5486693) q[3];
sx q[3];
rz(-2.8828354) q[3];
sx q[3];
rz(0.72507897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(-1.5287483) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(-1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(0.55364451) q[0];
rz(-1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.9083317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.935826) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(0.68769023) q[0];
rz(-pi) q[1];
rz(-0.1594752) q[2];
sx q[2];
rz(-0.6808241) q[2];
sx q[2];
rz(-2.8603539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.711449) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(-1.364691) q[1];
rz(-pi) q[2];
rz(-0.66744653) q[3];
sx q[3];
rz(-1.3818041) q[3];
sx q[3];
rz(2.8611956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0040940293) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(3.1393576) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8925791) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(2.2900443) q[0];
rz(-0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.203513) q[0];
sx q[0];
rz(-0.31103125) q[0];
sx q[0];
rz(0.87840338) q[0];
rz(-pi) q[1];
rz(-1.115828) q[2];
sx q[2];
rz(-1.3340545) q[2];
sx q[2];
rz(0.67982212) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5350266) q[1];
sx q[1];
rz(-0.92381322) q[1];
sx q[1];
rz(1.4024629) q[1];
rz(-2.381388) q[3];
sx q[3];
rz(-2.1225727) q[3];
sx q[3];
rz(2.9387568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-2.6181347) q[2];
rz(2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(-2.3160034) q[0];
rz(1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.694214) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71890812) q[0];
sx q[0];
rz(-1.8005162) q[0];
sx q[0];
rz(0.16750383) q[0];
x q[1];
rz(-1.6214192) q[2];
sx q[2];
rz(-2.9630337) q[2];
sx q[2];
rz(0.52390097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64926681) q[1];
sx q[1];
rz(-0.71486231) q[1];
sx q[1];
rz(-1.3267172) q[1];
rz(-pi) q[2];
rz(-2.7778266) q[3];
sx q[3];
rz(-1.3437265) q[3];
sx q[3];
rz(2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(3.1269126) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99575627) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(0.68403912) q[0];
rz(1.1902635) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(1.9285944) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0048464674) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(-0.61704163) q[0];
rz(-0.80981363) q[2];
sx q[2];
rz(-0.81293101) q[2];
sx q[2];
rz(0.44250689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.1274384) q[1];
sx q[1];
rz(-1.452983) q[1];
sx q[1];
rz(-0.25728667) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9641987) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(2.859476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(-0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260139) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(-2.8552326) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(-2.5568331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8557381) q[0];
sx q[0];
rz(-1.2585088) q[0];
sx q[0];
rz(-0.62494846) q[0];
x q[1];
rz(-0.61442394) q[2];
sx q[2];
rz(-1.5922058) q[2];
sx q[2];
rz(1.9911839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11440052) q[1];
sx q[1];
rz(-2.235734) q[1];
sx q[1];
rz(2.7726638) q[1];
x q[2];
rz(0.29295178) q[3];
sx q[3];
rz(-0.17360273) q[3];
sx q[3];
rz(2.0041182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.5131081) q[2];
rz(2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(-3.0896297) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(2.5351977) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0604027) q[0];
sx q[0];
rz(-0.58433825) q[0];
sx q[0];
rz(-1.1330182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7658471) q[2];
sx q[2];
rz(-0.94259113) q[2];
sx q[2];
rz(1.1496161) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28601521) q[1];
sx q[1];
rz(-0.77076036) q[1];
sx q[1];
rz(1.7325749) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9803195) q[3];
sx q[3];
rz(-1.2425353) q[3];
sx q[3];
rz(-1.8144153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2018044) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(1.4846876) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(-2.1915961) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(-0.22769134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75922155) q[0];
sx q[0];
rz(-1.641944) q[0];
sx q[0];
rz(-1.1573777) q[0];
x q[1];
rz(1.20485) q[2];
sx q[2];
rz(-1.1507251) q[2];
sx q[2];
rz(-0.36047381) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3701054) q[1];
sx q[1];
rz(-2.5574554) q[1];
sx q[1];
rz(-2.7258956) q[1];
x q[2];
rz(-0.065683059) q[3];
sx q[3];
rz(-2.5816397) q[3];
sx q[3];
rz(0.63510676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0630539) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7446328) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(-2.3378085) q[0];
rz(-2.0869758) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(-1.1484336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.046612) q[0];
sx q[0];
rz(-1.5946832) q[0];
sx q[0];
rz(-1.4987962) q[0];
rz(-0.36275136) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(-0.19778684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8879537) q[1];
sx q[1];
rz(-1.9983074) q[1];
sx q[1];
rz(2.0931787) q[1];
rz(-2.7896342) q[3];
sx q[3];
rz(-1.6037914) q[3];
sx q[3];
rz(1.9381423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-1.0174948) q[2];
sx q[2];
rz(-2.3633374) q[2];
rz(2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8675999) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(0.69828066) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90726501) q[0];
sx q[0];
rz(-0.31120473) q[0];
sx q[0];
rz(-2.9051203) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8928705) q[2];
sx q[2];
rz(-3.0634355) q[2];
sx q[2];
rz(0.31101481) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.56863927) q[1];
sx q[1];
rz(-2.0268974) q[1];
sx q[1];
rz(-0.24104636) q[1];
x q[2];
rz(-0.24153696) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(-2.9709771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0239821) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(0.84021604) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8828076) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(0.029126833) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(0.60224709) q[2];
sx q[2];
rz(-1.2802274) q[2];
sx q[2];
rz(1.0614492) q[2];
rz(-2.9914231) q[3];
sx q[3];
rz(-2.2554845) q[3];
sx q[3];
rz(3.1156202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
