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
rz(0.83528432) q[1];
sx q[1];
rz(-1.3568027) q[1];
sx q[1];
rz(-2.2396483) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.038739) q[0];
sx q[0];
rz(-1.7243288) q[0];
sx q[0];
rz(-2.5262336) q[0];
rz(-pi) q[1];
rz(1.3839533) q[2];
sx q[2];
rz(-1.9451491) q[2];
sx q[2];
rz(-2.6461305) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8613667) q[1];
sx q[1];
rz(-1.2806483) q[1];
sx q[1];
rz(2.9600083) q[1];
rz(2.4113184) q[3];
sx q[3];
rz(-2.090082) q[3];
sx q[3];
rz(-1.9744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96282643) q[2];
sx q[2];
rz(-2.1361394) q[2];
sx q[2];
rz(2.315305) q[2];
rz(-1.7360342) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(2.3981986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5590782) q[0];
sx q[0];
rz(-0.40575108) q[0];
sx q[0];
rz(-1.7412809) q[0];
rz(1.1236313) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(2.0434911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208291) q[0];
sx q[0];
rz(-1.4614824) q[0];
sx q[0];
rz(0.1569605) q[0];
rz(-1.2400572) q[2];
sx q[2];
rz(-1.8914701) q[2];
sx q[2];
rz(1.641524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2929165) q[1];
sx q[1];
rz(-1.6458938) q[1];
sx q[1];
rz(-0.1212705) q[1];
rz(-pi) q[2];
rz(2.9750175) q[3];
sx q[3];
rz(-1.1984571) q[3];
sx q[3];
rz(2.6083715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0850247) q[2];
sx q[2];
rz(-0.3173863) q[2];
sx q[2];
rz(1.8918096) q[2];
rz(-1.7791087) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(-1.484163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63610858) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(0.27534494) q[0];
rz(-2.6823726) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(-3.088248) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032565462) q[0];
sx q[0];
rz(-1.7836187) q[0];
sx q[0];
rz(1.4906917) q[0];
x q[1];
rz(-0.01845999) q[2];
sx q[2];
rz(-1.109913) q[2];
sx q[2];
rz(-2.8783952) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.054823067) q[1];
sx q[1];
rz(-2.5882698) q[1];
sx q[1];
rz(-3.064108) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9333604) q[3];
sx q[3];
rz(-1.2484115) q[3];
sx q[3];
rz(0.65402189) q[3];
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
rz(-1.5041171) q[2];
rz(-1.641168) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(2.8659081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212873) q[0];
sx q[0];
rz(-2.5581701) q[0];
sx q[0];
rz(0.48459184) q[0];
rz(1.8991607) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(2.7925083) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8139127) q[0];
sx q[0];
rz(-2.591344) q[0];
sx q[0];
rz(-1.489086) q[0];
x q[1];
rz(1.9140117) q[2];
sx q[2];
rz(-1.6417393) q[2];
sx q[2];
rz(-0.22954839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.0047093948) q[1];
sx q[1];
rz(-1.3420958) q[1];
sx q[1];
rz(1.3761702) q[1];
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
rz(-pi/2) q[1];
rz(-0.35859534) q[2];
sx q[2];
rz(-1.2025669) q[2];
sx q[2];
rz(2.6915468) q[2];
rz(-0.83886823) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(0.1990327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6435476) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(-3.041748) q[0];
rz(1.3205344) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(0.60755306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927101) q[0];
sx q[0];
rz(-1.6537871) q[0];
sx q[0];
rz(-1.4957499) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54536087) q[2];
sx q[2];
rz(-1.4358476) q[2];
sx q[2];
rz(-0.0098740059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1888426) q[1];
sx q[1];
rz(-2.5217068) q[1];
sx q[1];
rz(0.37921885) q[1];
x q[2];
rz(3.0436375) q[3];
sx q[3];
rz(-1.8281015) q[3];
sx q[3];
rz(-2.3549781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1943835) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-0.53059951) q[2];
rz(-0.44140205) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43668231) q[0];
sx q[0];
rz(-1.3446151) q[0];
sx q[0];
rz(-0.5740903) q[0];
rz(0.87999815) q[1];
sx q[1];
rz(-1.1209542) q[1];
sx q[1];
rz(-1.3425739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18595565) q[0];
sx q[0];
rz(-1.7711207) q[0];
sx q[0];
rz(2.0485417) q[0];
rz(-pi) q[1];
rz(0.75127496) q[2];
sx q[2];
rz(-2.5715368) q[2];
sx q[2];
rz(1.7715724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1426265) q[1];
sx q[1];
rz(-1.2104831) q[1];
sx q[1];
rz(2.8887038) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4785512) q[3];
sx q[3];
rz(-1.460307) q[3];
sx q[3];
rz(-2.7821531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67419702) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(2.6944842) q[2];
rz(-2.1111264) q[3];
sx q[3];
rz(-1.6298529) q[3];
sx q[3];
rz(0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717188) q[0];
sx q[0];
rz(-1.8444703) q[0];
sx q[0];
rz(2.7287927) q[0];
rz(2.0665456) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(0.80361754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5301054) q[0];
sx q[0];
rz(-0.75615935) q[0];
sx q[0];
rz(0.79571457) q[0];
x q[1];
rz(0.56411479) q[2];
sx q[2];
rz(-1.0061641) q[2];
sx q[2];
rz(-2.5052793) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6576801) q[1];
sx q[1];
rz(-1.5566623) q[1];
sx q[1];
rz(-0.78426984) q[1];
rz(1.4945696) q[3];
sx q[3];
rz(-1.3919358) q[3];
sx q[3];
rz(3.0146527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86218086) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(1.6406406) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(-1.6525432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.415446) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(0.55794445) q[0];
rz(-2.5803512) q[1];
sx q[1];
rz(-1.7702425) q[1];
sx q[1];
rz(1.4161313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.95074) q[0];
sx q[0];
rz(-1.2878364) q[0];
sx q[0];
rz(1.1038194) q[0];
x q[1];
rz(-0.92500706) q[2];
sx q[2];
rz(-0.46621284) q[2];
sx q[2];
rz(-0.031161873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76991612) q[1];
sx q[1];
rz(-1.8990108) q[1];
sx q[1];
rz(3.0382243) q[1];
rz(0.21702311) q[3];
sx q[3];
rz(-2.3742834) q[3];
sx q[3];
rz(1.2485152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8884362) q[2];
sx q[2];
rz(-0.77740589) q[2];
sx q[2];
rz(2.1412795) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.5001985) q[3];
sx q[3];
rz(2.9848671) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67824739) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(3.1372702) q[0];
rz(2.9777572) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(-1.77553) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8238433) q[0];
sx q[0];
rz(-1.5429521) q[0];
sx q[0];
rz(2.2180024) q[0];
rz(1.2605569) q[2];
sx q[2];
rz(-1.9266591) q[2];
sx q[2];
rz(2.9663393) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7228702) q[1];
sx q[1];
rz(-1.8218177) q[1];
sx q[1];
rz(1.6081564) q[1];
rz(-2.4716464) q[3];
sx q[3];
rz(-2.1996227) q[3];
sx q[3];
rz(2.3184702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0020478) q[2];
sx q[2];
rz(-2.1740972) q[2];
sx q[2];
rz(0.61919332) q[2];
rz(-2.5891417) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(0.14510554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3510975) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(-0.47875324) q[0];
rz(1.152773) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(2.1943888) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21175948) q[0];
sx q[0];
rz(-0.18778983) q[0];
sx q[0];
rz(-1.7540356) q[0];
x q[1];
rz(-0.98299111) q[2];
sx q[2];
rz(-1.776256) q[2];
sx q[2];
rz(-0.062177156) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0997206) q[1];
sx q[1];
rz(-2.2769089) q[1];
sx q[1];
rz(-0.40634917) q[1];
rz(-pi) q[2];
rz(0.35263108) q[3];
sx q[3];
rz(-1.336634) q[3];
sx q[3];
rz(-1.4591097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(-1.0673149) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(-2.1783569) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0530479) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(-1.0531986) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(1.5679892) q[2];
sx q[2];
rz(-1.3288099) q[2];
sx q[2];
rz(-2.228683) q[2];
rz(1.7942747) q[3];
sx q[3];
rz(-1.0031932) q[3];
sx q[3];
rz(-0.72465988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
