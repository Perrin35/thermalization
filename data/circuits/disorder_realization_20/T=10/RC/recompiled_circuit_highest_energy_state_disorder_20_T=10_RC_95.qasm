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
rz(0.96639291) q[0];
sx q[0];
rz(-1.7557431) q[0];
sx q[0];
rz(2.5461499) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(2.5505677) q[1];
sx q[1];
rz(12.443065) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714471) q[0];
sx q[0];
rz(-0.71200221) q[0];
sx q[0];
rz(2.2385236) q[0];
rz(-pi) q[1];
rz(0.79084556) q[2];
sx q[2];
rz(-2.186785) q[2];
sx q[2];
rz(-0.14033422) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5414303) q[1];
sx q[1];
rz(-1.3568391) q[1];
sx q[1];
rz(-3.0453277) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7756157) q[3];
sx q[3];
rz(-1.7133738) q[3];
sx q[3];
rz(-0.78998427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82723242) q[2];
sx q[2];
rz(-1.102697) q[2];
sx q[2];
rz(0.11191351) q[2];
rz(-1.7498451) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(-2.8351423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7258485) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(-0.14990212) q[0];
rz(2.8556178) q[1];
sx q[1];
rz(-1.3949225) q[1];
sx q[1];
rz(2.012595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0386815) q[0];
sx q[0];
rz(-1.8590218) q[0];
sx q[0];
rz(-3.1100357) q[0];
x q[1];
rz(-1.3311884) q[2];
sx q[2];
rz(-3.0802266) q[2];
sx q[2];
rz(2.1642016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4399229) q[1];
sx q[1];
rz(-1.4050802) q[1];
sx q[1];
rz(-2.7782337) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1993447) q[3];
sx q[3];
rz(-1.0744922) q[3];
sx q[3];
rz(1.7402349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67843208) q[2];
sx q[2];
rz(-0.9223991) q[2];
sx q[2];
rz(-1.9909667) q[2];
rz(1.9567418) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(-0.020603389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6612369) q[0];
sx q[0];
rz(-2.895597) q[0];
sx q[0];
rz(-0.38159698) q[0];
rz(-1.8474139) q[1];
sx q[1];
rz(-2.0353863) q[1];
sx q[1];
rz(-0.19827422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564045) q[0];
sx q[0];
rz(-1.745178) q[0];
sx q[0];
rz(-1.3402433) q[0];
rz(-pi) q[1];
rz(0.080520544) q[2];
sx q[2];
rz(-0.90290194) q[2];
sx q[2];
rz(2.5494247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5570017) q[1];
sx q[1];
rz(-1.7832568) q[1];
sx q[1];
rz(0.23834385) q[1];
rz(-pi) q[2];
rz(1.395969) q[3];
sx q[3];
rz(-1.5509737) q[3];
sx q[3];
rz(-2.1682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76564378) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(2.622733) q[2];
rz(1.2505924) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(2.5605104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4984703) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(0.44922391) q[0];
rz(1.6953702) q[1];
sx q[1];
rz(-2.4999373) q[1];
sx q[1];
rz(0.62072388) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0064471) q[0];
sx q[0];
rz(-0.9165316) q[0];
sx q[0];
rz(-1.3885137) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63346699) q[2];
sx q[2];
rz(-1.6704428) q[2];
sx q[2];
rz(-1.5758621) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1816239) q[1];
sx q[1];
rz(-1.9697497) q[1];
sx q[1];
rz(1.3514888) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4820547) q[3];
sx q[3];
rz(-2.6013298) q[3];
sx q[3];
rz(1.0080573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1064523) q[2];
sx q[2];
rz(-1.8702714) q[2];
sx q[2];
rz(0.16109666) q[2];
rz(-2.7437239) q[3];
sx q[3];
rz(-1.4784003) q[3];
sx q[3];
rz(2.9919992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905269) q[0];
sx q[0];
rz(-1.7790786) q[0];
sx q[0];
rz(-0.25704849) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-0.6404666) q[1];
sx q[1];
rz(-1.8880728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81957626) q[0];
sx q[0];
rz(-1.5347693) q[0];
sx q[0];
rz(-0.018112273) q[0];
rz(0.23373105) q[2];
sx q[2];
rz(-0.82373754) q[2];
sx q[2];
rz(-1.1540138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1977928) q[1];
sx q[1];
rz(-1.5256923) q[1];
sx q[1];
rz(-0.79774858) q[1];
rz(-pi) q[2];
rz(2.7944466) q[3];
sx q[3];
rz(-1.5724564) q[3];
sx q[3];
rz(-1.1356789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4132061) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(-1.2665292) q[2];
rz(-2.5201216) q[3];
sx q[3];
rz(-0.60756835) q[3];
sx q[3];
rz(0.5411886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40015873) q[0];
sx q[0];
rz(-2.3024237) q[0];
sx q[0];
rz(0.025064502) q[0];
rz(-1.9301682) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(2.8866344) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716041) q[0];
sx q[0];
rz(-1.5758262) q[0];
sx q[0];
rz(-0.049335376) q[0];
rz(0.58664257) q[2];
sx q[2];
rz(-1.1544747) q[2];
sx q[2];
rz(0.094791709) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3759758) q[1];
sx q[1];
rz(-1.0436397) q[1];
sx q[1];
rz(-0.65816452) q[1];
rz(-1.2473769) q[3];
sx q[3];
rz(-2.5445523) q[3];
sx q[3];
rz(3.0003594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24594626) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(-2.6825405) q[2];
rz(-1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(0.12115255) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46886214) q[0];
sx q[0];
rz(-0.48968306) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(2.22279) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(1.930621) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051832599) q[0];
sx q[0];
rz(-0.38516949) q[0];
sx q[0];
rz(2.5776107) q[0];
x q[1];
rz(-1.0031021) q[2];
sx q[2];
rz(-1.1922115) q[2];
sx q[2];
rz(0.97419435) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0414435) q[1];
sx q[1];
rz(-2.7584834) q[1];
sx q[1];
rz(-0.99975296) q[1];
x q[2];
rz(-2.7769412) q[3];
sx q[3];
rz(-1.9130008) q[3];
sx q[3];
rz(2.3221644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81898895) q[2];
sx q[2];
rz(-0.49590597) q[2];
sx q[2];
rz(-2.4244579) q[2];
rz(-2.1142193) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(0.43924847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9913919) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(-0.014658654) q[0];
rz(-0.75153366) q[1];
sx q[1];
rz(-1.2208168) q[1];
sx q[1];
rz(1.470648) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96318564) q[0];
sx q[0];
rz(-1.3826332) q[0];
sx q[0];
rz(-1.7117731) q[0];
rz(2.884955) q[2];
sx q[2];
rz(-1.6483288) q[2];
sx q[2];
rz(-0.77285779) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.60622207) q[1];
sx q[1];
rz(-2.722123) q[1];
sx q[1];
rz(1.8505881) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5603546) q[3];
sx q[3];
rz(-0.34509788) q[3];
sx q[3];
rz(-1.0853678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8810001) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(2.8680958) q[2];
rz(1.1547487) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(-0.65034136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038789373) q[0];
sx q[0];
rz(-2.0162855) q[0];
sx q[0];
rz(-1.0754841) q[0];
rz(0.12818809) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(0.34559616) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95042607) q[0];
sx q[0];
rz(-1.4109525) q[0];
sx q[0];
rz(0.9013844) q[0];
rz(-pi) q[1];
rz(2.8198411) q[2];
sx q[2];
rz(-1.0257799) q[2];
sx q[2];
rz(2.7779752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4491475) q[1];
sx q[1];
rz(-2.1112103) q[1];
sx q[1];
rz(-1.497333) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7365428) q[3];
sx q[3];
rz(-0.9791383) q[3];
sx q[3];
rz(1.4480643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0547611) q[2];
sx q[2];
rz(-1.0826449) q[2];
sx q[2];
rz(-1.6205988) q[2];
rz(-1.7704376) q[3];
sx q[3];
rz(-2.3833279) q[3];
sx q[3];
rz(-0.41485205) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3138251) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-0.23571043) q[0];
rz(2.56855) q[1];
sx q[1];
rz(-2.133281) q[1];
sx q[1];
rz(-1.7726353) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.154207) q[0];
sx q[0];
rz(-1.9999749) q[0];
sx q[0];
rz(-2.3230419) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4659285) q[2];
sx q[2];
rz(-2.0620637) q[2];
sx q[2];
rz(-1.5671687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0693501) q[1];
sx q[1];
rz(-1.4003039) q[1];
sx q[1];
rz(-2.1428698) q[1];
x q[2];
rz(-1.4769796) q[3];
sx q[3];
rz(-1.929885) q[3];
sx q[3];
rz(-2.1916568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77091757) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(3.0885922) q[2];
rz(0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5175405) q[0];
sx q[0];
rz(-2.1255827) q[0];
sx q[0];
rz(-0.95638635) q[0];
rz(-1.7908295) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(2.4443632) q[2];
sx q[2];
rz(-1.4705428) q[2];
sx q[2];
rz(-1.85208) q[2];
rz(2.3334669) q[3];
sx q[3];
rz(-1.1302409) q[3];
sx q[3];
rz(-1.7779868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
