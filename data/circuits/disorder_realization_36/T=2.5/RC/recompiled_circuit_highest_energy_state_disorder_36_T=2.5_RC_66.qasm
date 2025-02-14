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
rz(-0.041115887) q[0];
sx q[0];
rz(4.3423131) q[0];
sx q[0];
rz(11.095476) q[0];
rz(2.7566551) q[1];
sx q[1];
rz(-2.5632783) q[1];
sx q[1];
rz(-2.2338423) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9258869) q[0];
sx q[0];
rz(-2.3952847) q[0];
sx q[0];
rz(0.67912905) q[0];
x q[1];
rz(-0.29405354) q[2];
sx q[2];
rz(-1.8584677) q[2];
sx q[2];
rz(-1.2534382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1536394) q[1];
sx q[1];
rz(-1.0864294) q[1];
sx q[1];
rz(1.7820249) q[1];
rz(-3.0201245) q[3];
sx q[3];
rz(-1.5163912) q[3];
sx q[3];
rz(1.7732946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82751983) q[2];
sx q[2];
rz(-1.6562409) q[2];
sx q[2];
rz(0.32076389) q[2];
rz(-0.00094207923) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(-0.48085406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8568521) q[0];
sx q[0];
rz(-2.4212403) q[0];
sx q[0];
rz(-2.5208933) q[0];
rz(1.2868098) q[1];
sx q[1];
rz(-1.6544673) q[1];
sx q[1];
rz(2.1466045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5704203) q[0];
sx q[0];
rz(-1.7670399) q[0];
sx q[0];
rz(1.061383) q[0];
rz(-2.0109977) q[2];
sx q[2];
rz(-1.2792247) q[2];
sx q[2];
rz(-2.7879813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2590005) q[1];
sx q[1];
rz(-2.0152806) q[1];
sx q[1];
rz(-0.20436546) q[1];
rz(1.0424588) q[3];
sx q[3];
rz(-2.3064559) q[3];
sx q[3];
rz(0.56921116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2113125) q[2];
sx q[2];
rz(-1.899753) q[2];
sx q[2];
rz(2.9105183) q[2];
rz(1.6127582) q[3];
sx q[3];
rz(-1.684618) q[3];
sx q[3];
rz(0.55135623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424161) q[0];
sx q[0];
rz(-1.1671952) q[0];
sx q[0];
rz(-2.7447847) q[0];
rz(-1.9467758) q[1];
sx q[1];
rz(-1.8794941) q[1];
sx q[1];
rz(-1.4746812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75735649) q[0];
sx q[0];
rz(-3.0529733) q[0];
sx q[0];
rz(1.4411974) q[0];
x q[1];
rz(0.72309317) q[2];
sx q[2];
rz(-0.68487072) q[2];
sx q[2];
rz(1.9346646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4941229) q[1];
sx q[1];
rz(-2.2241333) q[1];
sx q[1];
rz(0.43237574) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9649404) q[3];
sx q[3];
rz(-1.5598462) q[3];
sx q[3];
rz(-1.8979069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84612238) q[2];
sx q[2];
rz(-2.2961871) q[2];
sx q[2];
rz(1.6579312) q[2];
rz(2.1814003) q[3];
sx q[3];
rz(-2.5064597) q[3];
sx q[3];
rz(2.4471762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33261499) q[0];
sx q[0];
rz(-1.806145) q[0];
sx q[0];
rz(0.72840011) q[0];
rz(-1.867846) q[1];
sx q[1];
rz(-1.1796917) q[1];
sx q[1];
rz(1.5625928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1611164) q[0];
sx q[0];
rz(-0.022810629) q[0];
sx q[0];
rz(-1.423832) q[0];
rz(-pi) q[1];
rz(-1.0228297) q[2];
sx q[2];
rz(-1.5336516) q[2];
sx q[2];
rz(-1.1556582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3660903) q[1];
sx q[1];
rz(-1.3555158) q[1];
sx q[1];
rz(-1.4332814) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2402917) q[3];
sx q[3];
rz(-1.309295) q[3];
sx q[3];
rz(-1.4217564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.031521156) q[2];
sx q[2];
rz(-1.4036125) q[2];
sx q[2];
rz(-0.20315476) q[2];
rz(1.3567694) q[3];
sx q[3];
rz(-1.2099096) q[3];
sx q[3];
rz(-1.320896) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2585231) q[0];
sx q[0];
rz(-1.8061545) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(-0.37172231) q[1];
sx q[1];
rz(-2.5860791) q[1];
sx q[1];
rz(2.1518167) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3759339) q[0];
sx q[0];
rz(-1.3816557) q[0];
sx q[0];
rz(-0.93171691) q[0];
rz(1.6838275) q[2];
sx q[2];
rz(-0.80146061) q[2];
sx q[2];
rz(-0.67455268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6367148) q[1];
sx q[1];
rz(-1.2413919) q[1];
sx q[1];
rz(0.18784951) q[1];
x q[2];
rz(-2.6017461) q[3];
sx q[3];
rz(-0.55571908) q[3];
sx q[3];
rz(-0.49969765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.524579) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(3.0830834) q[2];
rz(-1.8578364) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(-1.9730998) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002895) q[0];
sx q[0];
rz(-1.0463725) q[0];
sx q[0];
rz(-2.6562279) q[0];
rz(-0.46733388) q[1];
sx q[1];
rz(-2.554437) q[1];
sx q[1];
rz(2.4553518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9017603) q[0];
sx q[0];
rz(-2.2621381) q[0];
sx q[0];
rz(-1.9907238) q[0];
rz(1.2755308) q[2];
sx q[2];
rz(-0.64744189) q[2];
sx q[2];
rz(-1.9290627) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8778692) q[1];
sx q[1];
rz(-2.4179966) q[1];
sx q[1];
rz(2.3797026) q[1];
rz(-pi) q[2];
rz(2.3090906) q[3];
sx q[3];
rz(-1.0847391) q[3];
sx q[3];
rz(1.6532236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.301959) q[2];
sx q[2];
rz(-2.0531824) q[2];
sx q[2];
rz(-0.23183091) q[2];
rz(-2.2202282) q[3];
sx q[3];
rz(-1.5270343) q[3];
sx q[3];
rz(-3.11619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63474083) q[0];
sx q[0];
rz(-2.1379819) q[0];
sx q[0];
rz(-2.4941709) q[0];
rz(1.0415152) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(-1.1127068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96760945) q[0];
sx q[0];
rz(-2.4436473) q[0];
sx q[0];
rz(-0.83794727) q[0];
rz(1.9693807) q[2];
sx q[2];
rz(-2.9867134) q[2];
sx q[2];
rz(-1.9914371) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80363336) q[1];
sx q[1];
rz(-2.5542959) q[1];
sx q[1];
rz(2.8226167) q[1];
rz(-1.5912368) q[3];
sx q[3];
rz(-0.4651362) q[3];
sx q[3];
rz(-2.6488369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5645494) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(-0.88456336) q[2];
rz(2.1029643) q[3];
sx q[3];
rz(-2.0197155) q[3];
sx q[3];
rz(-1.60892) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3603947) q[0];
sx q[0];
rz(-2.9304657) q[0];
sx q[0];
rz(-1.1398844) q[0];
rz(0.79849157) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(0.07447695) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2452364) q[0];
sx q[0];
rz(-2.470861) q[0];
sx q[0];
rz(1.9696495) q[0];
x q[1];
rz(-0.76272599) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(0.21108) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.18596126) q[1];
sx q[1];
rz(-2.4096386) q[1];
sx q[1];
rz(0.47853985) q[1];
rz(-2.8387855) q[3];
sx q[3];
rz(-1.2330556) q[3];
sx q[3];
rz(-0.95888058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3227587) q[2];
sx q[2];
rz(-2.5311311) q[2];
sx q[2];
rz(0.16767137) q[2];
rz(-0.75336114) q[3];
sx q[3];
rz(-1.916901) q[3];
sx q[3];
rz(1.0838449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7849671) q[0];
sx q[0];
rz(-1.6870808) q[0];
sx q[0];
rz(1.0474569) q[0];
rz(0.52182237) q[1];
sx q[1];
rz(-1.9597041) q[1];
sx q[1];
rz(-2.3375915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88946901) q[0];
sx q[0];
rz(-1.7458545) q[0];
sx q[0];
rz(-1.1451707) q[0];
x q[1];
rz(1.2928771) q[2];
sx q[2];
rz(-0.57363311) q[2];
sx q[2];
rz(-0.81338993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2968353) q[1];
sx q[1];
rz(-1.8342819) q[1];
sx q[1];
rz(0.46864639) q[1];
rz(-pi) q[2];
rz(0.45553546) q[3];
sx q[3];
rz(-0.75596327) q[3];
sx q[3];
rz(2.8433133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1890586) q[2];
sx q[2];
rz(-1.6510094) q[2];
sx q[2];
rz(-0.35183364) q[2];
rz(-3.0686038) q[3];
sx q[3];
rz(-1.9609541) q[3];
sx q[3];
rz(2.4019057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9214582) q[0];
sx q[0];
rz(-2.298482) q[0];
sx q[0];
rz(-1.5716918) q[0];
rz(-0.72228557) q[1];
sx q[1];
rz(-2.0000358) q[1];
sx q[1];
rz(-1.3977745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98508747) q[0];
sx q[0];
rz(-2.6596375) q[0];
sx q[0];
rz(-2.5768649) q[0];
x q[1];
rz(0.70571122) q[2];
sx q[2];
rz(-1.7448261) q[2];
sx q[2];
rz(2.0311914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8885615) q[1];
sx q[1];
rz(-1.0200046) q[1];
sx q[1];
rz(1.4768449) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2820551) q[3];
sx q[3];
rz(-2.7061314) q[3];
sx q[3];
rz(2.9277152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7950644) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(1.8468924) q[2];
rz(2.0857701) q[3];
sx q[3];
rz(-1.0603696) q[3];
sx q[3];
rz(2.136371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.1718564) q[0];
sx q[0];
rz(-0.91005253) q[0];
sx q[0];
rz(-1.4685312) q[0];
rz(2.3079474) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(2.0019309) q[2];
sx q[2];
rz(-2.3464602) q[2];
sx q[2];
rz(-2.4577557) q[2];
rz(-2.007768) q[3];
sx q[3];
rz(-1.1935607) q[3];
sx q[3];
rz(0.85535645) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
