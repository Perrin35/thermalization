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
rz(0.63737386) q[0];
sx q[0];
rz(-2.0847991) q[0];
sx q[0];
rz(2.3699397) q[0];
rz(2.9906988) q[1];
sx q[1];
rz(-2.8291191) q[1];
sx q[1];
rz(1.6627275) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61687311) q[0];
sx q[0];
rz(-1.2390854) q[0];
sx q[0];
rz(1.2613368) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0885299) q[2];
sx q[2];
rz(-1.3446756) q[2];
sx q[2];
rz(2.256611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3841074) q[1];
sx q[1];
rz(-0.84534422) q[1];
sx q[1];
rz(1.2164335) q[1];
rz(-pi) q[2];
rz(-2.462384) q[3];
sx q[3];
rz(-1.7082885) q[3];
sx q[3];
rz(-2.7161571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5375157) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(-0.83211952) q[2];
rz(2.4845947) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(0.35201296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81545365) q[0];
sx q[0];
rz(-2.5322545) q[0];
sx q[0];
rz(0.64999181) q[0];
rz(-3.0187712) q[1];
sx q[1];
rz(-0.42418066) q[1];
sx q[1];
rz(2.4688683) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8685659) q[0];
sx q[0];
rz(-1.582524) q[0];
sx q[0];
rz(1.5290029) q[0];
x q[1];
rz(-2.4240875) q[2];
sx q[2];
rz(-0.40552545) q[2];
sx q[2];
rz(-0.44446352) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9422226) q[1];
sx q[1];
rz(-1.7666715) q[1];
sx q[1];
rz(-2.1327399) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1150713) q[3];
sx q[3];
rz(-1.0511729) q[3];
sx q[3];
rz(-2.5260203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0279205) q[2];
sx q[2];
rz(-1.5936667) q[2];
sx q[2];
rz(2.7552354) q[2];
rz(1.9733285) q[3];
sx q[3];
rz(-1.2315653) q[3];
sx q[3];
rz(1.108235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.89762178) q[0];
sx q[0];
rz(-1.4742999) q[0];
sx q[0];
rz(3.0834055) q[0];
rz(-2.8367786) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(-0.85495368) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8182504) q[0];
sx q[0];
rz(-1.5691681) q[0];
sx q[0];
rz(3.073284) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.57010713) q[2];
sx q[2];
rz(-2.3666581) q[2];
sx q[2];
rz(-0.75054815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7185229) q[1];
sx q[1];
rz(-0.81941019) q[1];
sx q[1];
rz(-1.8526444) q[1];
x q[2];
rz(0.26400395) q[3];
sx q[3];
rz(-1.0314634) q[3];
sx q[3];
rz(-1.7738938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0893687) q[2];
sx q[2];
rz(-1.3244018) q[2];
sx q[2];
rz(1.8734056) q[2];
rz(-2.7005699) q[3];
sx q[3];
rz(-0.62097582) q[3];
sx q[3];
rz(2.2997901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12255254) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(0.71504354) q[0];
rz(2.7252281) q[1];
sx q[1];
rz(-0.48960296) q[1];
sx q[1];
rz(-0.29410902) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.153115) q[0];
sx q[0];
rz(-1.5128969) q[0];
sx q[0];
rz(-2.2672976) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8744743) q[2];
sx q[2];
rz(-0.21409245) q[2];
sx q[2];
rz(-2.10432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3015705) q[1];
sx q[1];
rz(-1.4403655) q[1];
sx q[1];
rz(-2.954847) q[1];
rz(-pi) q[2];
rz(0.77463051) q[3];
sx q[3];
rz(-2.3916158) q[3];
sx q[3];
rz(0.74912723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9010345) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(-2.412839) q[2];
rz(-0.48480222) q[3];
sx q[3];
rz(-2.1383643) q[3];
sx q[3];
rz(2.9639967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4575551) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(-0.27873248) q[0];
rz(0.067525603) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(2.6356437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8794969) q[0];
sx q[0];
rz(-0.56096062) q[0];
sx q[0];
rz(-3.0039588) q[0];
x q[1];
rz(0.47053703) q[2];
sx q[2];
rz(-0.22122771) q[2];
sx q[2];
rz(-0.063723353) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5488779) q[1];
sx q[1];
rz(-1.028182) q[1];
sx q[1];
rz(-0.9701258) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26624532) q[3];
sx q[3];
rz(-2.2266386) q[3];
sx q[3];
rz(-0.4474934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2013596) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(-0.10227164) q[2];
rz(-0.30397948) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(2.6612018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85222307) q[0];
sx q[0];
rz(-2.3134573) q[0];
sx q[0];
rz(2.4612259) q[0];
rz(-0.96087372) q[1];
sx q[1];
rz(-1.5870794) q[1];
sx q[1];
rz(0.56112498) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2645123) q[0];
sx q[0];
rz(-1.7753557) q[0];
sx q[0];
rz(1.2007942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9343453) q[2];
sx q[2];
rz(-0.85734136) q[2];
sx q[2];
rz(2.0240473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1316072) q[1];
sx q[1];
rz(-1.5161884) q[1];
sx q[1];
rz(-2.8427678) q[1];
x q[2];
rz(-2.0797727) q[3];
sx q[3];
rz(-0.99409311) q[3];
sx q[3];
rz(0.95147607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23201021) q[2];
sx q[2];
rz(-0.97753111) q[2];
sx q[2];
rz(1.5634465) q[2];
rz(2.1252508) q[3];
sx q[3];
rz(-1.4278102) q[3];
sx q[3];
rz(-2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.72961724) q[0];
sx q[0];
rz(-2.602674) q[0];
sx q[0];
rz(2.6475661) q[0];
rz(1.5771075) q[1];
sx q[1];
rz(-2.1421075) q[1];
sx q[1];
rz(0.090242537) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9245042) q[0];
sx q[0];
rz(-2.6441433) q[0];
sx q[0];
rz(0.50636943) q[0];
rz(-pi) q[1];
rz(-0.87553067) q[2];
sx q[2];
rz(-1.9771431) q[2];
sx q[2];
rz(-1.6900495) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3857194) q[1];
sx q[1];
rz(-1.8579646) q[1];
sx q[1];
rz(1.9906635) q[1];
rz(0.19321038) q[3];
sx q[3];
rz(-1.8487559) q[3];
sx q[3];
rz(1.6707458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51284853) q[2];
sx q[2];
rz(-0.57400846) q[2];
sx q[2];
rz(-0.64344704) q[2];
rz(-0.63055864) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.028320463) q[0];
sx q[0];
rz(-1.7009108) q[0];
sx q[0];
rz(-1.2233618) q[0];
rz(-0.89448294) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(-1.9953413) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2558738) q[0];
sx q[0];
rz(-2.4303655) q[0];
sx q[0];
rz(-2.2125986) q[0];
rz(-pi) q[1];
rz(-0.38739631) q[2];
sx q[2];
rz(-1.9700389) q[2];
sx q[2];
rz(-2.7170167) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8889956) q[1];
sx q[1];
rz(-1.8114144) q[1];
sx q[1];
rz(0.29901286) q[1];
x q[2];
rz(2.1560539) q[3];
sx q[3];
rz(-1.7904591) q[3];
sx q[3];
rz(-1.0128164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.869183) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(1.6669081) q[2];
rz(-0.21371755) q[3];
sx q[3];
rz(-1.4054047) q[3];
sx q[3];
rz(0.86764446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84233442) q[0];
sx q[0];
rz(-0.61115757) q[0];
sx q[0];
rz(-2.4364731) q[0];
rz(0.57535386) q[1];
sx q[1];
rz(-1.4082785) q[1];
sx q[1];
rz(0.56952482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1181208) q[0];
sx q[0];
rz(-3.0295353) q[0];
sx q[0];
rz(-0.81046354) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6274983) q[2];
sx q[2];
rz(-2.5894458) q[2];
sx q[2];
rz(0.29468003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3635253) q[1];
sx q[1];
rz(-0.91202867) q[1];
sx q[1];
rz(-2.3450065) q[1];
rz(-1.2678746) q[3];
sx q[3];
rz(-1.1008796) q[3];
sx q[3];
rz(-2.5055714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9465955) q[2];
sx q[2];
rz(-0.94299287) q[2];
sx q[2];
rz(-2.1627964) q[2];
rz(-2.7106674) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(-2.0144958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8579213) q[0];
sx q[0];
rz(-0.019275276) q[0];
sx q[0];
rz(1.4625782) q[0];
rz(2.1025533) q[1];
sx q[1];
rz(-1.3835399) q[1];
sx q[1];
rz(-1.9336112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4911352) q[0];
sx q[0];
rz(-1.2692599) q[0];
sx q[0];
rz(-1.9670606) q[0];
rz(3.1080148) q[2];
sx q[2];
rz(-1.0462282) q[2];
sx q[2];
rz(-1.7941504) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0145899) q[1];
sx q[1];
rz(-0.63204884) q[1];
sx q[1];
rz(-1.8188502) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1753756) q[3];
sx q[3];
rz(-1.6334157) q[3];
sx q[3];
rz(0.6400125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98564467) q[2];
sx q[2];
rz(-2.1375103) q[2];
sx q[2];
rz(1.3204302) q[2];
rz(0.350658) q[3];
sx q[3];
rz(-0.81736332) q[3];
sx q[3];
rz(2.5613274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.25528) q[0];
sx q[0];
rz(-1.1863615) q[0];
sx q[0];
rz(-0.86082389) q[0];
rz(-1.4844004) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(-0.17694494) q[2];
sx q[2];
rz(-2.7569303) q[2];
sx q[2];
rz(-3.0349572) q[2];
rz(-1.2091533) q[3];
sx q[3];
rz(-2.4222838) q[3];
sx q[3];
rz(-2.1771976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
