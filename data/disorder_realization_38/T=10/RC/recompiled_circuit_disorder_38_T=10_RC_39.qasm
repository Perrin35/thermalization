OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8367856) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(0.027713393) q[0];
x q[1];
rz(-2.7031541) q[2];
sx q[2];
rz(-1.0796094) q[2];
sx q[2];
rz(3.0868798) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5397415) q[1];
sx q[1];
rz(-1.9508233) q[1];
sx q[1];
rz(2.3439581) q[1];
rz(-pi) q[2];
rz(1.689765) q[3];
sx q[3];
rz(-0.61421466) q[3];
sx q[3];
rz(-3.0005232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41574079) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(-1.6289904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166116) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(-0.38831098) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5674044) q[2];
sx q[2];
rz(-1.2870103) q[2];
sx q[2];
rz(-2.0397182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7035547) q[1];
sx q[1];
rz(-0.8982656) q[1];
sx q[1];
rz(-2.9707675) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23962044) q[3];
sx q[3];
rz(-1.2926971) q[3];
sx q[3];
rz(0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(-2.3382323) q[2];
rz(2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-0.025618205) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(0.74584109) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94507664) q[0];
sx q[0];
rz(-1.3498422) q[0];
sx q[0];
rz(-3.0483079) q[0];
rz(1.7350115) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(2.8525713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5902139) q[1];
sx q[1];
rz(-0.56400245) q[1];
sx q[1];
rz(2.6022634) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8196705) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(-2.4528465) q[2];
rz(3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.21928366) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(-3.02137) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(-2.8682958) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7057987) q[0];
sx q[0];
rz(-1.6596284) q[0];
sx q[0];
rz(1.8949132) q[0];
rz(-1.0068514) q[2];
sx q[2];
rz(-1.0014921) q[2];
sx q[2];
rz(-3.0727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30444333) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(1.3346346) q[1];
x q[2];
rz(-2.9076505) q[3];
sx q[3];
rz(-0.66286874) q[3];
sx q[3];
rz(-0.58594221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79167241) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(-2.1544429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7894831) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(-2.2228918) q[0];
rz(0.012533112) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(-2.8963793) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.6189177) q[1];
sx q[1];
rz(-2.4213311) q[1];
sx q[1];
rz(0.48536761) q[1];
rz(-pi) q[2];
rz(-1.4694674) q[3];
sx q[3];
rz(-2.7469027) q[3];
sx q[3];
rz(0.24995382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(-0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.8621559) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813542) q[0];
sx q[0];
rz(-1.8525774) q[0];
sx q[0];
rz(-3.1262585) q[0];
x q[1];
rz(-2.5336669) q[2];
sx q[2];
rz(-1.330901) q[2];
sx q[2];
rz(1.8976256) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93950677) q[1];
sx q[1];
rz(-1.0491976) q[1];
sx q[1];
rz(-0.034394666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.389723) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8906158) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(0.77077579) q[2];
rz(-1.6714913) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(-0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(3.1380222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033747) q[0];
sx q[0];
rz(-0.57045454) q[0];
sx q[0];
rz(-2.4128782) q[0];
rz(1.1793169) q[2];
sx q[2];
rz(-0.85748312) q[2];
sx q[2];
rz(-1.6187402) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.057295416) q[1];
sx q[1];
rz(-0.97505403) q[1];
sx q[1];
rz(-0.36235313) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7690897) q[3];
sx q[3];
rz(-1.4962422) q[3];
sx q[3];
rz(-1.2411181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-0.41771093) q[3];
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
rz(-0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(-0.67869854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6060651) q[0];
sx q[0];
rz(-1.8341656) q[0];
sx q[0];
rz(2.8019964) q[0];
x q[1];
rz(-2.3085824) q[2];
sx q[2];
rz(-1.2150803) q[2];
sx q[2];
rz(-2.9290111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2862257) q[1];
sx q[1];
rz(-1.6337506) q[1];
sx q[1];
rz(2.0952058) q[1];
rz(-pi) q[2];
rz(2.7160866) q[3];
sx q[3];
rz(-1.4262428) q[3];
sx q[3];
rz(-1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(2.2237681) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-0.25979364) q[0];
rz(-2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(0.52694595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5848815) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(0.35993872) q[0];
x q[1];
rz(-0.91295816) q[2];
sx q[2];
rz(-1.3557938) q[2];
sx q[2];
rz(2.1195597) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.37104169) q[1];
sx q[1];
rz(-2.7275804) q[1];
sx q[1];
rz(2.4310334) q[1];
x q[2];
rz(-1.8606887) q[3];
sx q[3];
rz(-1.1416417) q[3];
sx q[3];
rz(3.1115301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(2.1561484) q[0];
rz(-0.5685637) q[1];
sx q[1];
rz(-1.111258) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11689582) q[0];
sx q[0];
rz(-1.7381867) q[0];
sx q[0];
rz(0.025339729) q[0];
rz(-2.0759517) q[2];
sx q[2];
rz(-2.897122) q[2];
sx q[2];
rz(-2.5214362) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0204748) q[1];
sx q[1];
rz(-1.7133683) q[1];
sx q[1];
rz(-2.3974182) q[1];
x q[2];
rz(0.71513128) q[3];
sx q[3];
rz(-0.22278856) q[3];
sx q[3];
rz(-2.1905394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0439904) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(0.94341755) q[2];
rz(-2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5671134) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(1.3654937) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(-0.75136649) q[3];
sx q[3];
rz(-2.1204227) q[3];
sx q[3];
rz(-1.942996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
