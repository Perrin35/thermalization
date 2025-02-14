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
rz(-0.72467726) q[0];
sx q[0];
rz(4.4530498) q[0];
sx q[0];
rz(11.761576) q[0];
rz(0.59504396) q[1];
sx q[1];
rz(-2.9437328) q[1];
sx q[1];
rz(1.2208389) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0752995) q[0];
sx q[0];
rz(-1.4508712) q[0];
sx q[0];
rz(-1.6060702) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8684565) q[2];
sx q[2];
rz(-0.47800666) q[2];
sx q[2];
rz(-0.82551685) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28452595) q[1];
sx q[1];
rz(-1.6150805) q[1];
sx q[1];
rz(-0.74378777) q[1];
x q[2];
rz(-0.63296707) q[3];
sx q[3];
rz(-1.6817012) q[3];
sx q[3];
rz(-1.1735011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7120984) q[2];
sx q[2];
rz(-2.6338989) q[2];
sx q[2];
rz(-1.3362159) q[2];
rz(1.3974238) q[3];
sx q[3];
rz(-1.6349399) q[3];
sx q[3];
rz(-1.5857182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76032388) q[0];
sx q[0];
rz(-1.6345432) q[0];
sx q[0];
rz(2.4632578) q[0];
rz(0.61596576) q[1];
sx q[1];
rz(-2.1149642) q[1];
sx q[1];
rz(-0.14332992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1122637) q[0];
sx q[0];
rz(-1.2486132) q[0];
sx q[0];
rz(-0.16561819) q[0];
rz(-pi) q[1];
rz(1.4392339) q[2];
sx q[2];
rz(-2.4189197) q[2];
sx q[2];
rz(-0.45718788) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7166087) q[1];
sx q[1];
rz(-1.0657446) q[1];
sx q[1];
rz(1.5165971) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7160104) q[3];
sx q[3];
rz(-1.2493361) q[3];
sx q[3];
rz(-1.5150573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6501288) q[2];
sx q[2];
rz(-2.3626707) q[2];
sx q[2];
rz(-1.9219363) q[2];
rz(2.8236112) q[3];
sx q[3];
rz(-2.3173083) q[3];
sx q[3];
rz(-2.4021751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3984482) q[0];
sx q[0];
rz(-0.92107451) q[0];
sx q[0];
rz(0.63968023) q[0];
rz(1.8094874) q[1];
sx q[1];
rz(-2.2001241) q[1];
sx q[1];
rz(-2.926362) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23519653) q[0];
sx q[0];
rz(-1.0825038) q[0];
sx q[0];
rz(-1.3286235) q[0];
rz(-pi) q[1];
rz(3.136342) q[2];
sx q[2];
rz(-1.8551146) q[2];
sx q[2];
rz(2.9470987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1647268) q[1];
sx q[1];
rz(-1.9661207) q[1];
sx q[1];
rz(-3.1302451) q[1];
rz(1.8120469) q[3];
sx q[3];
rz(-0.46866207) q[3];
sx q[3];
rz(-1.7023583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3053863) q[2];
sx q[2];
rz(-0.47577327) q[2];
sx q[2];
rz(-2.9212908) q[2];
rz(-2.3588755) q[3];
sx q[3];
rz(-1.3681151) q[3];
sx q[3];
rz(2.9960184) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7102605) q[0];
sx q[0];
rz(-0.98589698) q[0];
sx q[0];
rz(1.3166991) q[0];
rz(-2.6915164) q[1];
sx q[1];
rz(-1.7066259) q[1];
sx q[1];
rz(-3.0302474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836783) q[0];
sx q[0];
rz(-0.92843572) q[0];
sx q[0];
rz(1.8877554) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0527332) q[2];
sx q[2];
rz(-1.8919164) q[2];
sx q[2];
rz(3.113609) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9604019) q[1];
sx q[1];
rz(-0.63397303) q[1];
sx q[1];
rz(2.0540103) q[1];
rz(-pi) q[2];
rz(-2.4347859) q[3];
sx q[3];
rz(-2.4549559) q[3];
sx q[3];
rz(-1.5070908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4075809) q[2];
sx q[2];
rz(-1.13052) q[2];
sx q[2];
rz(0.13047516) q[2];
rz(2.981251) q[3];
sx q[3];
rz(-0.66157833) q[3];
sx q[3];
rz(2.9714238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57123667) q[0];
sx q[0];
rz(-0.2061051) q[0];
sx q[0];
rz(-0.11827949) q[0];
rz(2.2453399) q[1];
sx q[1];
rz(-1.6629442) q[1];
sx q[1];
rz(-0.55312696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0118857) q[0];
sx q[0];
rz(-1.0965523) q[0];
sx q[0];
rz(-1.9408731) q[0];
rz(-pi) q[1];
rz(-2.4055355) q[2];
sx q[2];
rz(-1.9070574) q[2];
sx q[2];
rz(-2.2081809) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0748079) q[1];
sx q[1];
rz(-2.5321567) q[1];
sx q[1];
rz(1.273073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0712264) q[3];
sx q[3];
rz(-2.4132846) q[3];
sx q[3];
rz(1.0906995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3502189) q[2];
sx q[2];
rz(-0.72177902) q[2];
sx q[2];
rz(1.0132033) q[2];
rz(2.8241217) q[3];
sx q[3];
rz(-1.443913) q[3];
sx q[3];
rz(2.1875994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1182227) q[0];
sx q[0];
rz(-0.88339266) q[0];
sx q[0];
rz(0.50514847) q[0];
rz(0.39816868) q[1];
sx q[1];
rz(-1.265637) q[1];
sx q[1];
rz(-1.3291976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564541) q[0];
sx q[0];
rz(-0.49665652) q[0];
sx q[0];
rz(-1.9209548) q[0];
rz(-pi) q[1];
rz(-2.5754314) q[2];
sx q[2];
rz(-2.5306411) q[2];
sx q[2];
rz(-2.6662155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3214282) q[1];
sx q[1];
rz(-1.0099995) q[1];
sx q[1];
rz(3.012077) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81327) q[3];
sx q[3];
rz(-2.7754042) q[3];
sx q[3];
rz(1.0367108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.087611) q[2];
sx q[2];
rz(-1.3302646) q[2];
sx q[2];
rz(-2.3114253) q[2];
rz(1.4812482) q[3];
sx q[3];
rz(-1.0426499) q[3];
sx q[3];
rz(-1.9771077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2723349) q[0];
sx q[0];
rz(-2.2255958) q[0];
sx q[0];
rz(-1.9328312) q[0];
rz(0.2441949) q[1];
sx q[1];
rz(-1.9701651) q[1];
sx q[1];
rz(-0.29921439) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7231185) q[0];
sx q[0];
rz(-1.5445853) q[0];
sx q[0];
rz(-1.7869201) q[0];
rz(1.0396378) q[2];
sx q[2];
rz(-2.1130145) q[2];
sx q[2];
rz(1.7249267) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7119042) q[1];
sx q[1];
rz(-1.5221847) q[1];
sx q[1];
rz(2.2232988) q[1];
rz(-2.6664566) q[3];
sx q[3];
rz(-2.3804602) q[3];
sx q[3];
rz(-0.5055529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.19305912) q[2];
sx q[2];
rz(-1.4024573) q[2];
sx q[2];
rz(2.194727) q[2];
rz(-1.7638505) q[3];
sx q[3];
rz(-0.54072127) q[3];
sx q[3];
rz(2.4989959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4544025) q[0];
sx q[0];
rz(-0.43455046) q[0];
sx q[0];
rz(-0.53892556) q[0];
rz(2.9680805) q[1];
sx q[1];
rz(-0.75650802) q[1];
sx q[1];
rz(0.38421806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095188524) q[0];
sx q[0];
rz(-2.6676819) q[0];
sx q[0];
rz(2.7907811) q[0];
rz(-pi) q[1];
rz(-0.60224763) q[2];
sx q[2];
rz(-1.411323) q[2];
sx q[2];
rz(-0.98724706) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3612566) q[1];
sx q[1];
rz(-1.5118264) q[1];
sx q[1];
rz(-3.1153989) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.027111091) q[3];
sx q[3];
rz(-2.1246006) q[3];
sx q[3];
rz(-1.1983271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82181278) q[2];
sx q[2];
rz(-2.3723448) q[2];
sx q[2];
rz(2.5452851) q[2];
rz(-2.1042306) q[3];
sx q[3];
rz(-2.3682902) q[3];
sx q[3];
rz(1.9572702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3768815) q[0];
sx q[0];
rz(-2.3887971) q[0];
sx q[0];
rz(-2.9994614) q[0];
rz(1.1171974) q[1];
sx q[1];
rz(-1.486472) q[1];
sx q[1];
rz(1.2275009) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9949029) q[0];
sx q[0];
rz(-1.0953383) q[0];
sx q[0];
rz(-1.882193) q[0];
rz(-1.6864482) q[2];
sx q[2];
rz(-3.0031548) q[2];
sx q[2];
rz(-1.1756681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.043316226) q[1];
sx q[1];
rz(-2.195561) q[1];
sx q[1];
rz(-0.079778133) q[1];
x q[2];
rz(-2.7015721) q[3];
sx q[3];
rz(-2.2624863) q[3];
sx q[3];
rz(-1.7753851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.050798945) q[2];
sx q[2];
rz(-2.8992081) q[2];
sx q[2];
rz(-2.3432815) q[2];
rz(-1.5797041) q[3];
sx q[3];
rz(-1.3825994) q[3];
sx q[3];
rz(0.17670512) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.056219) q[0];
sx q[0];
rz(-2.4985785) q[0];
sx q[0];
rz(-3.1363078) q[0];
rz(0.082848631) q[1];
sx q[1];
rz(-1.1033892) q[1];
sx q[1];
rz(-2.8740035) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0326313) q[0];
sx q[0];
rz(-1.8825476) q[0];
sx q[0];
rz(0.70272081) q[0];
rz(-pi) q[1];
rz(-0.29163714) q[2];
sx q[2];
rz(-1.1979629) q[2];
sx q[2];
rz(-2.458984) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53658518) q[1];
sx q[1];
rz(-1.3571902) q[1];
sx q[1];
rz(1.2475292) q[1];
rz(-pi) q[2];
rz(-1.0571207) q[3];
sx q[3];
rz(-1.9909715) q[3];
sx q[3];
rz(3.1204124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.573632) q[2];
sx q[2];
rz(-0.48144123) q[2];
sx q[2];
rz(1.9285704) q[2];
rz(-1.8968449) q[3];
sx q[3];
rz(-0.48143482) q[3];
sx q[3];
rz(-1.9511694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5382814) q[0];
sx q[0];
rz(-2.17157) q[0];
sx q[0];
rz(-3.0441913) q[0];
rz(-3.1311323) q[1];
sx q[1];
rz(-2.2757826) q[1];
sx q[1];
rz(-0.59715685) q[1];
rz(-2.9928906) q[2];
sx q[2];
rz(-0.54051334) q[2];
sx q[2];
rz(1.9220026) q[2];
rz(-0.98164557) q[3];
sx q[3];
rz(-1.4245778) q[3];
sx q[3];
rz(2.1733976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
