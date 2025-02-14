OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(8.2052054) q[0];
rz(3.976877) q[1];
sx q[1];
rz(4.4983954) q[1];
sx q[1];
rz(8.5228336) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57563462) q[0];
sx q[0];
rz(-2.1778641) q[0];
sx q[0];
rz(-1.3835039) q[0];
rz(-pi) q[1];
rz(0.38036687) q[2];
sx q[2];
rz(-1.3970301) q[2];
sx q[2];
rz(-1.997239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.343051) q[1];
sx q[1];
rz(-1.7447116) q[1];
sx q[1];
rz(1.2760722) q[1];
rz(-pi) q[2];
rz(2.2252623) q[3];
sx q[3];
rz(-2.188465) q[3];
sx q[3];
rz(-0.82181069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96282643) q[2];
sx q[2];
rz(-2.1361394) q[2];
sx q[2];
rz(-2.315305) q[2];
rz(1.4055584) q[3];
sx q[3];
rz(-0.30705753) q[3];
sx q[3];
rz(-2.3981986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590782) q[0];
sx q[0];
rz(-0.40575108) q[0];
sx q[0];
rz(1.7412809) q[0];
rz(1.1236313) q[1];
sx q[1];
rz(-2.7476937) q[1];
sx q[1];
rz(-2.0434911) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0880249) q[0];
sx q[0];
rz(-2.9505749) q[0];
sx q[0];
rz(2.5294526) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3674521) q[2];
sx q[2];
rz(-0.4565276) q[2];
sx q[2];
rz(2.4696333) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2929165) q[1];
sx q[1];
rz(-1.6458938) q[1];
sx q[1];
rz(3.0203222) q[1];
rz(-1.1693177) q[3];
sx q[3];
rz(-2.7352834) q[3];
sx q[3];
rz(-3.0412948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0565679) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(1.8918096) q[2];
rz(-1.362484) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(-1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5054841) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(2.8662477) q[0];
rz(2.6823726) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(3.088248) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5864075) q[0];
sx q[0];
rz(-1.4925028) q[0];
sx q[0];
rz(2.9281062) q[0];
rz(-pi) q[1];
x q[1];
rz(1.533639) q[2];
sx q[2];
rz(-0.46122631) q[2];
sx q[2];
rz(-0.30468582) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5596603) q[1];
sx q[1];
rz(-1.5301063) q[1];
sx q[1];
rz(-0.55198021) q[1];
rz(0.20823222) q[3];
sx q[3];
rz(-1.2484115) q[3];
sx q[3];
rz(-2.4875708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3205545) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(1.6374755) q[2];
rz(1.5004246) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(2.8659081) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42030537) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(-0.48459184) q[0];
rz(1.2424319) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(0.34908435) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8139127) q[0];
sx q[0];
rz(-2.591344) q[0];
sx q[0];
rz(-1.489086) q[0];
rz(-pi) q[1];
rz(1.3626839) q[2];
sx q[2];
rz(-0.35018626) q[2];
sx q[2];
rz(1.1454358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5308456) q[1];
sx q[1];
rz(-1.3813003) q[1];
sx q[1];
rz(-2.9086472) q[1];
rz(-pi) q[2];
rz(2.7233299) q[3];
sx q[3];
rz(-1.0111448) q[3];
sx q[3];
rz(-0.5681611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35859534) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-2.6915468) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(-2.94256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435476) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(-0.099844649) q[0];
rz(1.3205344) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(-2.5340396) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71220416) q[0];
sx q[0];
rz(-0.11183248) q[0];
sx q[0];
rz(2.4080316) q[0];
x q[1];
rz(-0.54536087) q[2];
sx q[2];
rz(-1.4358476) q[2];
sx q[2];
rz(-0.0098740059) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6441356) q[1];
sx q[1];
rz(-1.0007528) q[1];
sx q[1];
rz(-1.312478) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2149548) q[3];
sx q[3];
rz(-0.27493335) q[3];
sx q[3];
rz(-2.7235018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1943835) q[2];
sx q[2];
rz(-0.76346976) q[2];
sx q[2];
rz(-2.6109931) q[2];
rz(2.7001906) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(-2.1921564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43668231) q[0];
sx q[0];
rz(-1.3446151) q[0];
sx q[0];
rz(-2.5675024) q[0];
rz(-2.2615945) q[1];
sx q[1];
rz(-2.0206385) q[1];
sx q[1];
rz(-1.7990187) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18595565) q[0];
sx q[0];
rz(-1.7711207) q[0];
sx q[0];
rz(-1.093051) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9832575) q[2];
sx q[2];
rz(-1.165373) q[2];
sx q[2];
rz(-0.5328005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7751392) q[1];
sx q[1];
rz(-0.43700686) q[1];
sx q[1];
rz(-0.98458146) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4785512) q[3];
sx q[3];
rz(-1.6812857) q[3];
sx q[3];
rz(-2.7821531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-0.28417045) q[2];
sx q[2];
rz(-2.6944842) q[2];
rz(1.0304662) q[3];
sx q[3];
rz(-1.6298529) q[3];
sx q[3];
rz(0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36987385) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(-2.7287927) q[0];
rz(-1.075047) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(0.80361754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8014288) q[0];
sx q[0];
rz(-1.0699843) q[0];
sx q[0];
rz(2.1636971) q[0];
x q[1];
rz(-0.92757954) q[2];
sx q[2];
rz(-2.0394435) q[2];
sx q[2];
rz(-1.2609294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4839125) q[1];
sx q[1];
rz(-1.5849304) q[1];
sx q[1];
rz(-2.3573228) q[1];
rz(1.647023) q[3];
sx q[3];
rz(-1.7496568) q[3];
sx q[3];
rz(3.0146527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86218086) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(1.6406406) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-2.3707135) q[3];
sx q[3];
rz(-1.4890495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.415446) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(-0.55794445) q[0];
rz(-2.5803512) q[1];
sx q[1];
rz(-1.7702425) q[1];
sx q[1];
rz(1.4161313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87420207) q[0];
sx q[0];
rz(-0.54055291) q[0];
sx q[0];
rz(2.1442599) q[0];
rz(-pi) q[1];
rz(-2.8475262) q[2];
sx q[2];
rz(-1.9379788) q[2];
sx q[2];
rz(-0.73193708) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76991612) q[1];
sx q[1];
rz(-1.2425819) q[1];
sx q[1];
rz(3.0382243) q[1];
rz(0.21702311) q[3];
sx q[3];
rz(-0.76730928) q[3];
sx q[3];
rz(-1.2485152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8884362) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-1.0003132) q[2];
rz(3.0194164) q[3];
sx q[3];
rz(-1.5001985) q[3];
sx q[3];
rz(-2.9848671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67824739) q[0];
sx q[0];
rz(-1.954701) q[0];
sx q[0];
rz(3.1372702) q[0];
rz(2.9777572) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.3660627) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9095847) q[0];
sx q[0];
rz(-0.9238832) q[0];
sx q[0];
rz(-0.034897403) q[0];
rz(-pi) q[1];
rz(2.4539656) q[2];
sx q[2];
rz(-2.673871) q[2];
sx q[2];
rz(0.91910884) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7228702) q[1];
sx q[1];
rz(-1.8218177) q[1];
sx q[1];
rz(-1.6081564) q[1];
rz(2.2774599) q[3];
sx q[3];
rz(-0.88418761) q[3];
sx q[3];
rz(0.1089801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(-2.5223993) q[2];
rz(-2.5891417) q[3];
sx q[3];
rz(-1.3055472) q[3];
sx q[3];
rz(2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7904952) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(-0.47875324) q[0];
rz(-1.152773) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(-2.1943888) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9626434) q[0];
sx q[0];
rz(-1.6048204) q[0];
sx q[0];
rz(-1.7555139) q[0];
rz(-pi) q[1];
rz(2.1586015) q[2];
sx q[2];
rz(-1.3653367) q[2];
sx q[2];
rz(-3.0794155) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6852764) q[1];
sx q[1];
rz(-2.3446977) q[1];
sx q[1];
rz(-2.0048098) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7889616) q[3];
sx q[3];
rz(-1.336634) q[3];
sx q[3];
rz(-1.682483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(2.8749386) q[2];
rz(-2.0742778) q[3];
sx q[3];
rz(-1.2131194) q[3];
sx q[3];
rz(0.96323577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(3.0530479) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(2.0883941) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(-3.13022) q[2];
sx q[2];
rz(-0.24200242) q[2];
sx q[2];
rz(0.92462362) q[2];
rz(-0.57907818) q[3];
sx q[3];
rz(-1.3828207) q[3];
sx q[3];
rz(-2.4170392) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
