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
rz(0.44952196) q[0];
sx q[0];
rz(4.8604453) q[0];
sx q[0];
rz(7.3331375) q[0];
rz(2.6226251) q[1];
sx q[1];
rz(1.5331886) q[1];
sx q[1];
rz(9.4756995) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.19456) q[0];
sx q[0];
rz(-1.8935391) q[0];
sx q[0];
rz(2.6850924) q[0];
x q[1];
rz(-2.9257183) q[2];
sx q[2];
rz(-1.2435438) q[2];
sx q[2];
rz(-0.33282166) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66827735) q[1];
sx q[1];
rz(-0.44679579) q[1];
sx q[1];
rz(2.3832641) q[1];
rz(-pi) q[2];
rz(2.256278) q[3];
sx q[3];
rz(-1.3847794) q[3];
sx q[3];
rz(1.3929588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8672987) q[2];
sx q[2];
rz(-2.2653502) q[2];
sx q[2];
rz(-1.3977316) q[2];
rz(-2.6869669) q[3];
sx q[3];
rz(-0.8780829) q[3];
sx q[3];
rz(1.7129869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4947263) q[0];
sx q[0];
rz(-2.2541101) q[0];
sx q[0];
rz(0.60321641) q[0];
rz(1.8869205) q[1];
sx q[1];
rz(-0.61379495) q[1];
sx q[1];
rz(-2.5616554) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7042687) q[0];
sx q[0];
rz(-2.5781729) q[0];
sx q[0];
rz(-2.826344) q[0];
rz(1.1179148) q[2];
sx q[2];
rz(-1.0476026) q[2];
sx q[2];
rz(2.4459237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1973639) q[1];
sx q[1];
rz(-1.6104638) q[1];
sx q[1];
rz(-0.96645379) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9469684) q[3];
sx q[3];
rz(-2.1056089) q[3];
sx q[3];
rz(0.68262284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6569528) q[2];
sx q[2];
rz(-0.75289774) q[2];
sx q[2];
rz(2.7638655) q[2];
rz(-1.7083302) q[3];
sx q[3];
rz(-1.6323615) q[3];
sx q[3];
rz(-1.9506955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.79320532) q[0];
sx q[0];
rz(-3.0866525) q[0];
sx q[0];
rz(-1.6435664) q[0];
rz(-1.9801697) q[1];
sx q[1];
rz(-1.8893416) q[1];
sx q[1];
rz(0.84322554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0204165) q[0];
sx q[0];
rz(-1.5112229) q[0];
sx q[0];
rz(-2.1942433) q[0];
x q[1];
rz(-2.0328838) q[2];
sx q[2];
rz(-0.46056754) q[2];
sx q[2];
rz(2.2572287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9084042) q[1];
sx q[1];
rz(-2.452525) q[1];
sx q[1];
rz(1.3702964) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2349818) q[3];
sx q[3];
rz(-0.89794176) q[3];
sx q[3];
rz(-1.3300174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8680385) q[2];
sx q[2];
rz(-2.4783583) q[2];
sx q[2];
rz(0.794945) q[2];
rz(-0.76977175) q[3];
sx q[3];
rz(-2.218518) q[3];
sx q[3];
rz(2.9845089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.58761) q[0];
sx q[0];
rz(-1.3251323) q[0];
sx q[0];
rz(2.8532568) q[0];
rz(-2.4544857) q[1];
sx q[1];
rz(-1.6485051) q[1];
sx q[1];
rz(-1.4289325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9754959) q[0];
sx q[0];
rz(-1.5642484) q[0];
sx q[0];
rz(0.01207821) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4304787) q[2];
sx q[2];
rz(-1.9694917) q[2];
sx q[2];
rz(1.2600419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9058207) q[1];
sx q[1];
rz(-1.8265939) q[1];
sx q[1];
rz(1.5033493) q[1];
x q[2];
rz(-0.63923423) q[3];
sx q[3];
rz(-1.0022638) q[3];
sx q[3];
rz(2.6871339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.571542) q[2];
sx q[2];
rz(-0.96633458) q[2];
sx q[2];
rz(-2.9759882) q[2];
rz(3.0770732) q[3];
sx q[3];
rz(-3.0118628) q[3];
sx q[3];
rz(-1.8701514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9128543) q[0];
sx q[0];
rz(-1.1930635) q[0];
sx q[0];
rz(-0.91113973) q[0];
rz(-2.1707824) q[1];
sx q[1];
rz(-1.3263005) q[1];
sx q[1];
rz(0.86311805) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1939094) q[0];
sx q[0];
rz(-1.860431) q[0];
sx q[0];
rz(-2.0783483) q[0];
x q[1];
rz(1.3023071) q[2];
sx q[2];
rz(-1.5483861) q[2];
sx q[2];
rz(-0.42699285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6533901) q[1];
sx q[1];
rz(-2.0987256) q[1];
sx q[1];
rz(1.8406244) q[1];
rz(-pi) q[2];
rz(1.7105402) q[3];
sx q[3];
rz(-2.826353) q[3];
sx q[3];
rz(1.3073352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0866278) q[2];
sx q[2];
rz(-1.2530155) q[2];
sx q[2];
rz(0.48913726) q[2];
rz(-2.1431811) q[3];
sx q[3];
rz(-2.9676134) q[3];
sx q[3];
rz(3.0424931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84422207) q[0];
sx q[0];
rz(-0.64592823) q[0];
sx q[0];
rz(-2.5883664) q[0];
rz(2.3440701) q[1];
sx q[1];
rz(-0.99634606) q[1];
sx q[1];
rz(1.9901989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40276179) q[0];
sx q[0];
rz(-0.58332764) q[0];
sx q[0];
rz(2.8999694) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8676574) q[2];
sx q[2];
rz(-0.56800743) q[2];
sx q[2];
rz(-0.97962475) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4321255) q[1];
sx q[1];
rz(-1.8399939) q[1];
sx q[1];
rz(1.3666735) q[1];
x q[2];
rz(-2.2362333) q[3];
sx q[3];
rz(-1.2631772) q[3];
sx q[3];
rz(-0.24505982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3449751) q[2];
sx q[2];
rz(-1.4130219) q[2];
sx q[2];
rz(-0.83149347) q[2];
rz(1.8031395) q[3];
sx q[3];
rz(-1.0748539) q[3];
sx q[3];
rz(-0.89053806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7688585) q[0];
sx q[0];
rz(-1.9356198) q[0];
sx q[0];
rz(-0.15705577) q[0];
rz(1.318469) q[1];
sx q[1];
rz(-0.942197) q[1];
sx q[1];
rz(2.7838321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3387332) q[0];
sx q[0];
rz(-2.0293616) q[0];
sx q[0];
rz(-2.0684956) q[0];
rz(-pi) q[1];
rz(-2.6951201) q[2];
sx q[2];
rz(-0.95281592) q[2];
sx q[2];
rz(-0.22836049) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4688229) q[1];
sx q[1];
rz(-2.2190418) q[1];
sx q[1];
rz(0.85521568) q[1];
rz(-pi) q[2];
rz(0.57447432) q[3];
sx q[3];
rz(-1.7041429) q[3];
sx q[3];
rz(0.80783081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99178189) q[2];
sx q[2];
rz(-1.8842183) q[2];
sx q[2];
rz(-3.120976) q[2];
rz(1.8990382) q[3];
sx q[3];
rz(-0.81762448) q[3];
sx q[3];
rz(0.29153618) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108516) q[0];
sx q[0];
rz(-1.956097) q[0];
sx q[0];
rz(-0.32671842) q[0];
rz(-2.2945981) q[1];
sx q[1];
rz(-1.0495443) q[1];
sx q[1];
rz(1.0583896) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21789078) q[0];
sx q[0];
rz(-2.0837221) q[0];
sx q[0];
rz(3.0615527) q[0];
rz(-pi) q[1];
rz(-2.5024274) q[2];
sx q[2];
rz(-1.4286388) q[2];
sx q[2];
rz(-2.7246812) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10524532) q[1];
sx q[1];
rz(-2.7509442) q[1];
sx q[1];
rz(1.3643144) q[1];
x q[2];
rz(-1.0048546) q[3];
sx q[3];
rz(-1.3019239) q[3];
sx q[3];
rz(-3.1389126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1191001) q[2];
sx q[2];
rz(-2.7969226) q[2];
sx q[2];
rz(-1.9692839) q[2];
rz(0.0017702866) q[3];
sx q[3];
rz(-0.5032731) q[3];
sx q[3];
rz(-2.1625904) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86972791) q[0];
sx q[0];
rz(-2.3379022) q[0];
sx q[0];
rz(-1.3386238) q[0];
rz(1.1805234) q[1];
sx q[1];
rz(-1.0931284) q[1];
sx q[1];
rz(2.6611633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6201037) q[0];
sx q[0];
rz(-2.0248981) q[0];
sx q[0];
rz(-0.2257077) q[0];
rz(-pi) q[1];
rz(-1.6132203) q[2];
sx q[2];
rz(-1.7406751) q[2];
sx q[2];
rz(-2.1720048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94225223) q[1];
sx q[1];
rz(-2.128731) q[1];
sx q[1];
rz(-0.11726484) q[1];
x q[2];
rz(-0.55175169) q[3];
sx q[3];
rz(-0.41424879) q[3];
sx q[3];
rz(-1.5442314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7052445) q[2];
sx q[2];
rz(-2.1465325) q[2];
sx q[2];
rz(1.0850517) q[2];
rz(-2.0294225) q[3];
sx q[3];
rz(-2.2472436) q[3];
sx q[3];
rz(-0.81356847) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8535264) q[0];
sx q[0];
rz(-0.94678322) q[0];
sx q[0];
rz(1.0119525) q[0];
rz(-2.2400253) q[1];
sx q[1];
rz(-2.8755867) q[1];
sx q[1];
rz(-2.576135) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8249567) q[0];
sx q[0];
rz(-2.4554283) q[0];
sx q[0];
rz(1.9030722) q[0];
rz(-pi) q[1];
rz(0.26763518) q[2];
sx q[2];
rz(-1.9861172) q[2];
sx q[2];
rz(-3.0618947) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31670891) q[1];
sx q[1];
rz(-1.4780586) q[1];
sx q[1];
rz(-2.2888661) q[1];
x q[2];
rz(1.8875445) q[3];
sx q[3];
rz(-1.4707139) q[3];
sx q[3];
rz(2.7730178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7908343) q[2];
sx q[2];
rz(-1.5459205) q[2];
sx q[2];
rz(-1.2606384) q[2];
rz(2.6993921) q[3];
sx q[3];
rz(-2.111777) q[3];
sx q[3];
rz(1.7650167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5762536) q[0];
sx q[0];
rz(-2.2591142) q[0];
sx q[0];
rz(2.2478065) q[0];
rz(1.5743938) q[1];
sx q[1];
rz(-1.4566474) q[1];
sx q[1];
rz(1.7172071) q[1];
rz(-0.42548634) q[2];
sx q[2];
rz(-1.31447) q[2];
sx q[2];
rz(-0.18438495) q[2];
rz(1.1461729) q[3];
sx q[3];
rz(-2.7594447) q[3];
sx q[3];
rz(-1.2654163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
