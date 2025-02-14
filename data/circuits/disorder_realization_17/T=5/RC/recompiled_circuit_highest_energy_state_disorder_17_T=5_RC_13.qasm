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
rz(2.0789335) q[0];
sx q[0];
rz(3.9570424) q[0];
sx q[0];
rz(11.856356) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(11.234565) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3478617) q[0];
sx q[0];
rz(-1.2958741) q[0];
sx q[0];
rz(-1.8040276) q[0];
rz(-pi) q[1];
rz(-2.5153036) q[2];
sx q[2];
rz(-2.7262027) q[2];
sx q[2];
rz(0.20520575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9985245) q[1];
sx q[1];
rz(-1.3245163) q[1];
sx q[1];
rz(2.1889436) q[1];
rz(2.2043742) q[3];
sx q[3];
rz(-1.211418) q[3];
sx q[3];
rz(2.0546953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4189202) q[2];
sx q[2];
rz(-1.2166497) q[2];
sx q[2];
rz(-2.326272) q[2];
rz(-0.80960649) q[3];
sx q[3];
rz(-1.4987192) q[3];
sx q[3];
rz(-2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.078449) q[0];
sx q[0];
rz(-1.0493295) q[0];
sx q[0];
rz(3.0310042) q[0];
rz(-2.1965006) q[1];
sx q[1];
rz(-2.7022305) q[1];
sx q[1];
rz(-1.5623215) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7720582) q[0];
sx q[0];
rz(-0.86449776) q[0];
sx q[0];
rz(0.95916551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0192102) q[2];
sx q[2];
rz(-1.1875936) q[2];
sx q[2];
rz(-2.2926083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89797276) q[1];
sx q[1];
rz(-2.0531516) q[1];
sx q[1];
rz(1.4528758) q[1];
rz(2.6591127) q[3];
sx q[3];
rz(-2.2594995) q[3];
sx q[3];
rz(0.83079332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9940146) q[2];
sx q[2];
rz(-2.9958604) q[2];
sx q[2];
rz(-2.6914524) q[2];
rz(-1.133793) q[3];
sx q[3];
rz(-1.6885875) q[3];
sx q[3];
rz(-0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728773) q[0];
sx q[0];
rz(-2.1935538) q[0];
sx q[0];
rz(0.28133389) q[0];
rz(1.3866792) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(-2.5414355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0735787) q[0];
sx q[0];
rz(-1.8610483) q[0];
sx q[0];
rz(-1.5223548) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1776206) q[2];
sx q[2];
rz(-1.6168878) q[2];
sx q[2];
rz(0.54296903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.197613) q[1];
sx q[1];
rz(-2.463732) q[1];
sx q[1];
rz(-0.53450905) q[1];
rz(0.26140499) q[3];
sx q[3];
rz(-2.5525744) q[3];
sx q[3];
rz(-2.1847069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0967789) q[2];
sx q[2];
rz(-1.1324977) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(2.7438296) q[3];
sx q[3];
rz(-2.4989276) q[3];
sx q[3];
rz(-0.1325632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9699049) q[0];
sx q[0];
rz(-0.02709087) q[0];
sx q[0];
rz(-2.0890253) q[0];
rz(-1.196208) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(-0.41950163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1284258) q[0];
sx q[0];
rz(-1.8343529) q[0];
sx q[0];
rz(2.5530784) q[0];
x q[1];
rz(2.9151462) q[2];
sx q[2];
rz(-2.1927877) q[2];
sx q[2];
rz(0.43606191) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7345404) q[1];
sx q[1];
rz(-2.190935) q[1];
sx q[1];
rz(-0.46874115) q[1];
rz(1.2798843) q[3];
sx q[3];
rz(-2.0613163) q[3];
sx q[3];
rz(-1.377493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9516051) q[2];
sx q[2];
rz(-2.0230468) q[2];
sx q[2];
rz(1.0901964) q[2];
rz(-2.6103141) q[3];
sx q[3];
rz(-2.4254906) q[3];
sx q[3];
rz(-1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215609) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(1.3979727) q[0];
rz(-0.13294237) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(0.001210777) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.087861) q[0];
sx q[0];
rz(-2.7606761) q[0];
sx q[0];
rz(0.043405966) q[0];
rz(-2.9820061) q[2];
sx q[2];
rz(-1.8529466) q[2];
sx q[2];
rz(-2.7196376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1493133) q[1];
sx q[1];
rz(-0.94301254) q[1];
sx q[1];
rz(1.1924442) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3962383) q[3];
sx q[3];
rz(-0.15197411) q[3];
sx q[3];
rz(-1.1308972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91263897) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(-2.9310628) q[2];
rz(-2.1052965) q[3];
sx q[3];
rz(-0.784289) q[3];
sx q[3];
rz(1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9482816) q[0];
sx q[0];
rz(-1.9183777) q[0];
sx q[0];
rz(-0.40570983) q[0];
rz(-1.3544719) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(0.92232651) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852823) q[0];
sx q[0];
rz(-1.6772224) q[0];
sx q[0];
rz(-1.2805416) q[0];
x q[1];
rz(-1.7013266) q[2];
sx q[2];
rz(-2.4707327) q[2];
sx q[2];
rz(-2.6595569) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.90627024) q[1];
sx q[1];
rz(-0.962469) q[1];
sx q[1];
rz(-3.0023605) q[1];
rz(0.57862307) q[3];
sx q[3];
rz(-2.4994183) q[3];
sx q[3];
rz(-2.7901621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.9269383) q[2];
sx q[2];
rz(2.7396438) q[2];
rz(-1.2881783) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(-1.8624381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.5227018) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(-3.0772305) q[0];
rz(2.3035658) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(1.8109969) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9143599) q[0];
sx q[0];
rz(-1.6961251) q[0];
sx q[0];
rz(-0.57526875) q[0];
rz(-pi) q[1];
rz(-2.4766972) q[2];
sx q[2];
rz(-2.4192296) q[2];
sx q[2];
rz(1.6160053) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76112642) q[1];
sx q[1];
rz(-2.5674324) q[1];
sx q[1];
rz(-0.4131743) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95860062) q[3];
sx q[3];
rz(-1.4293839) q[3];
sx q[3];
rz(2.8116016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4055206) q[2];
sx q[2];
rz(-1.5084927) q[2];
sx q[2];
rz(0.68354496) q[2];
rz(-0.73379597) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(-0.84764135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(1.1731359) q[0];
rz(0.37551156) q[1];
sx q[1];
rz(-2.651732) q[1];
sx q[1];
rz(1.0367905) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11060729) q[0];
sx q[0];
rz(-1.5944423) q[0];
sx q[0];
rz(1.551669) q[0];
x q[1];
rz(-1.8470959) q[2];
sx q[2];
rz(-1.7299998) q[2];
sx q[2];
rz(-0.72180787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4040096) q[1];
sx q[1];
rz(-0.77464691) q[1];
sx q[1];
rz(-0.20139106) q[1];
x q[2];
rz(2.6508207) q[3];
sx q[3];
rz(-1.0939071) q[3];
sx q[3];
rz(-2.8618438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5658687) q[2];
sx q[2];
rz(-0.91369358) q[2];
sx q[2];
rz(2.4110528) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-0.52339619) q[3];
sx q[3];
rz(-1.8714347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4850979) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(3.094161) q[0];
rz(0.221953) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(0.91112959) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8442366) q[0];
sx q[0];
rz(-0.043023303) q[0];
sx q[0];
rz(-2.4639315) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5510727) q[2];
sx q[2];
rz(-0.64177401) q[2];
sx q[2];
rz(2.3373147) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4621689) q[1];
sx q[1];
rz(-2.7738681) q[1];
sx q[1];
rz(-0.035178758) q[1];
rz(1.0044596) q[3];
sx q[3];
rz(-1.7438423) q[3];
sx q[3];
rz(0.53849788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86965108) q[2];
sx q[2];
rz(-2.7615774) q[2];
sx q[2];
rz(-0.16858777) q[2];
rz(0.21909675) q[3];
sx q[3];
rz(-2.2897661) q[3];
sx q[3];
rz(1.3271837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1472226) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(-2.6841573) q[0];
rz(1.2262729) q[1];
sx q[1];
rz(-0.33718449) q[1];
sx q[1];
rz(-1.7785243) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1049688) q[0];
sx q[0];
rz(-1.6311121) q[0];
sx q[0];
rz(0.6078267) q[0];
rz(-2.5117433) q[2];
sx q[2];
rz(-1.3721264) q[2];
sx q[2];
rz(-2.3185286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41238989) q[1];
sx q[1];
rz(-0.41657428) q[1];
sx q[1];
rz(2.7135486) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0567597) q[3];
sx q[3];
rz(-2.2515577) q[3];
sx q[3];
rz(2.5670402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3739796) q[2];
sx q[2];
rz(-0.55459905) q[2];
sx q[2];
rz(-0.45200959) q[2];
rz(0.40397817) q[3];
sx q[3];
rz(-1.2510866) q[3];
sx q[3];
rz(2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44611888) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(2.6192464) q[1];
sx q[1];
rz(-0.53032395) q[1];
sx q[1];
rz(1.2293336) q[1];
rz(-2.859237) q[2];
sx q[2];
rz(-1.3236125) q[2];
sx q[2];
rz(-0.57409928) q[2];
rz(1.4215076) q[3];
sx q[3];
rz(-1.310077) q[3];
sx q[3];
rz(2.6826912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
