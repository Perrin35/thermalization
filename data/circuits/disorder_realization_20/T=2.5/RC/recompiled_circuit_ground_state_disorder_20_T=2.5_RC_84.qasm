OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29937509) q[0];
sx q[0];
rz(-2.8111281) q[0];
sx q[0];
rz(-1.0634896) q[0];
rz(-0.039634135) q[1];
sx q[1];
rz(-0.57365817) q[1];
sx q[1];
rz(-1.080245) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40908694) q[0];
sx q[0];
rz(-1.5686146) q[0];
sx q[0];
rz(1.1078784) q[0];
rz(-pi) q[1];
rz(3.1267794) q[2];
sx q[2];
rz(-1.3990566) q[2];
sx q[2];
rz(-1.1352676) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1138405) q[1];
sx q[1];
rz(-0.58958399) q[1];
sx q[1];
rz(-0.47926183) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61038252) q[3];
sx q[3];
rz(-2.3993407) q[3];
sx q[3];
rz(0.28377747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0237191) q[2];
sx q[2];
rz(-1.4154075) q[2];
sx q[2];
rz(2.7685557) q[2];
rz(0.55752623) q[3];
sx q[3];
rz(-1.5159461) q[3];
sx q[3];
rz(-2.663234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28906223) q[0];
sx q[0];
rz(-1.4339148) q[0];
sx q[0];
rz(1.0093932) q[0];
rz(2.8139662) q[1];
sx q[1];
rz(-0.3265003) q[1];
sx q[1];
rz(2.0351298) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94325667) q[0];
sx q[0];
rz(-0.13511629) q[0];
sx q[0];
rz(-0.19395239) q[0];
rz(1.5554788) q[2];
sx q[2];
rz(-0.88587447) q[2];
sx q[2];
rz(-0.1372125) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7944462) q[1];
sx q[1];
rz(-2.3443188) q[1];
sx q[1];
rz(-2.2496012) q[1];
rz(2.7170466) q[3];
sx q[3];
rz(-2.694482) q[3];
sx q[3];
rz(-3.0101484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.893225) q[2];
sx q[2];
rz(-1.9577586) q[2];
sx q[2];
rz(0.66967213) q[2];
rz(1.9855965) q[3];
sx q[3];
rz(-0.35231927) q[3];
sx q[3];
rz(-0.34935752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356693) q[0];
sx q[0];
rz(-2.7140129) q[0];
sx q[0];
rz(1.3315573) q[0];
rz(1.9412899) q[1];
sx q[1];
rz(-0.2526865) q[1];
sx q[1];
rz(-0.67218626) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.975544) q[0];
sx q[0];
rz(-0.48915809) q[0];
sx q[0];
rz(2.8821095) q[0];
x q[1];
rz(-1.9483445) q[2];
sx q[2];
rz(-0.43127764) q[2];
sx q[2];
rz(1.4403421) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7461639) q[1];
sx q[1];
rz(-2.6708467) q[1];
sx q[1];
rz(-0.37396273) q[1];
x q[2];
rz(3.0447118) q[3];
sx q[3];
rz(-1.5200429) q[3];
sx q[3];
rz(-1.0602204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4868698) q[2];
sx q[2];
rz(-1.9630311) q[2];
sx q[2];
rz(-2.2910924) q[2];
rz(0.9507829) q[3];
sx q[3];
rz(-2.6362004) q[3];
sx q[3];
rz(-0.90184414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1281857) q[0];
sx q[0];
rz(-0.81939092) q[0];
sx q[0];
rz(0.73981458) q[0];
rz(0.20135227) q[1];
sx q[1];
rz(-1.9722152) q[1];
sx q[1];
rz(-2.9154725) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27557954) q[0];
sx q[0];
rz(-3.1034307) q[0];
sx q[0];
rz(1.1665795) q[0];
rz(1.1925542) q[2];
sx q[2];
rz(-1.8762445) q[2];
sx q[2];
rz(0.39337197) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2831336) q[1];
sx q[1];
rz(-1.0609416) q[1];
sx q[1];
rz(0.63555952) q[1];
x q[2];
rz(-1.4410517) q[3];
sx q[3];
rz(-1.1888388) q[3];
sx q[3];
rz(0.48473919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40253338) q[2];
sx q[2];
rz(-0.89503521) q[2];
sx q[2];
rz(0.99739897) q[2];
rz(-0.41336695) q[3];
sx q[3];
rz(-1.9316542) q[3];
sx q[3];
rz(1.0443784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1362374) q[0];
sx q[0];
rz(-0.69664609) q[0];
sx q[0];
rz(-0.053939017) q[0];
rz(-2.9593762) q[1];
sx q[1];
rz(-0.91582623) q[1];
sx q[1];
rz(2.838476) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5526472) q[0];
sx q[0];
rz(-0.56845462) q[0];
sx q[0];
rz(-0.89737749) q[0];
x q[1];
rz(-1.8218173) q[2];
sx q[2];
rz(-0.66749882) q[2];
sx q[2];
rz(-2.4199744) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50072155) q[1];
sx q[1];
rz(-2.5815775) q[1];
sx q[1];
rz(-1.1997833) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3833586) q[3];
sx q[3];
rz(-1.544974) q[3];
sx q[3];
rz(1.7968221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4431241) q[2];
sx q[2];
rz(-1.8734525) q[2];
sx q[2];
rz(2.7847086) q[2];
rz(2.3667864) q[3];
sx q[3];
rz(-1.509343) q[3];
sx q[3];
rz(-2.7029412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2543432) q[0];
sx q[0];
rz(-1.9312504) q[0];
sx q[0];
rz(2.3526225) q[0];
rz(-1.7987159) q[1];
sx q[1];
rz(-1.7306381) q[1];
sx q[1];
rz(-2.3416187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8784284) q[0];
sx q[0];
rz(-1.5845926) q[0];
sx q[0];
rz(-2.2565802) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7943013) q[2];
sx q[2];
rz(-2.433521) q[2];
sx q[2];
rz(-2.023165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3331078) q[1];
sx q[1];
rz(-1.6792381) q[1];
sx q[1];
rz(-2.7641349) q[1];
x q[2];
rz(0.39928945) q[3];
sx q[3];
rz(-1.1762816) q[3];
sx q[3];
rz(-2.7189915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0903025) q[2];
sx q[2];
rz(-0.74959576) q[2];
sx q[2];
rz(-1.3517815) q[2];
rz(-2.6734062) q[3];
sx q[3];
rz(-1.0951833) q[3];
sx q[3];
rz(-0.93402544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4139597) q[0];
sx q[0];
rz(-0.41828823) q[0];
sx q[0];
rz(0.002451238) q[0];
rz(-0.0062423627) q[1];
sx q[1];
rz(-2.7793482) q[1];
sx q[1];
rz(-2.7856766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9166989) q[0];
sx q[0];
rz(-2.421764) q[0];
sx q[0];
rz(-1.0309321) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8628274) q[2];
sx q[2];
rz(-2.27072) q[2];
sx q[2];
rz(-2.2577159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.98470515) q[1];
sx q[1];
rz(-0.92159373) q[1];
sx q[1];
rz(-2.8125416) q[1];
rz(-pi) q[2];
rz(0.53419279) q[3];
sx q[3];
rz(-1.6592798) q[3];
sx q[3];
rz(-2.8719287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9956841) q[2];
sx q[2];
rz(-1.2331139) q[2];
sx q[2];
rz(2.0705059) q[2];
rz(-0.60761014) q[3];
sx q[3];
rz(-2.0376318) q[3];
sx q[3];
rz(1.5266017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0869658) q[0];
sx q[0];
rz(-0.93769756) q[0];
sx q[0];
rz(1.1338393) q[0];
rz(1.0519823) q[1];
sx q[1];
rz(-1.4583505) q[1];
sx q[1];
rz(-0.7410616) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4008557) q[0];
sx q[0];
rz(-2.1261524) q[0];
sx q[0];
rz(-1.8779138) q[0];
rz(-0.95693077) q[2];
sx q[2];
rz(-1.9939878) q[2];
sx q[2];
rz(2.500071) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8246671) q[1];
sx q[1];
rz(-0.52579907) q[1];
sx q[1];
rz(0.047917685) q[1];
rz(-pi) q[2];
rz(-2.2005733) q[3];
sx q[3];
rz(-1.9004603) q[3];
sx q[3];
rz(2.8235112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1715601) q[2];
sx q[2];
rz(-0.67983183) q[2];
sx q[2];
rz(-0.46530923) q[2];
rz(0.71632898) q[3];
sx q[3];
rz(-1.497044) q[3];
sx q[3];
rz(1.959645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5301836) q[0];
sx q[0];
rz(-2.1937328) q[0];
sx q[0];
rz(1.1780257) q[0];
rz(-1.0097965) q[1];
sx q[1];
rz(-1.3346846) q[1];
sx q[1];
rz(-3.03426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4383746) q[0];
sx q[0];
rz(-2.1461663) q[0];
sx q[0];
rz(0.98585016) q[0];
rz(1.4194518) q[2];
sx q[2];
rz(-0.79337315) q[2];
sx q[2];
rz(-1.3248487) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8047528) q[1];
sx q[1];
rz(-0.64267413) q[1];
sx q[1];
rz(-0.39037946) q[1];
rz(-pi) q[2];
rz(-1.1501649) q[3];
sx q[3];
rz(-2.5007043) q[3];
sx q[3];
rz(0.057614652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44059077) q[2];
sx q[2];
rz(-0.89459449) q[2];
sx q[2];
rz(-0.44006285) q[2];
rz(-0.20982404) q[3];
sx q[3];
rz(-1.194229) q[3];
sx q[3];
rz(-2.8275209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1485727) q[0];
sx q[0];
rz(-2.5555389) q[0];
sx q[0];
rz(1.3826189) q[0];
rz(-2.4644409) q[1];
sx q[1];
rz(-1.4644198) q[1];
sx q[1];
rz(-2.009353) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5234194) q[0];
sx q[0];
rz(-1.3480524) q[0];
sx q[0];
rz(-1.3370418) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8381565) q[2];
sx q[2];
rz(-0.79435452) q[2];
sx q[2];
rz(1.8092138) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59076485) q[1];
sx q[1];
rz(-0.71749291) q[1];
sx q[1];
rz(-0.036542459) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72158815) q[3];
sx q[3];
rz(-2.3240328) q[3];
sx q[3];
rz(-2.8374654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9558692) q[2];
sx q[2];
rz(-0.74803868) q[2];
sx q[2];
rz(-1.424074) q[2];
rz(-0.87797034) q[3];
sx q[3];
rz(-2.9214171) q[3];
sx q[3];
rz(-0.018208114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2672742) q[0];
sx q[0];
rz(-1.6813288) q[0];
sx q[0];
rz(1.9718476) q[0];
rz(2.7727903) q[1];
sx q[1];
rz(-1.4099051) q[1];
sx q[1];
rz(-3.1261408) q[1];
rz(2.7013576) q[2];
sx q[2];
rz(-1.8320089) q[2];
sx q[2];
rz(1.6377891) q[2];
rz(0.5795547) q[3];
sx q[3];
rz(-1.8892291) q[3];
sx q[3];
rz(-1.5665594) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
