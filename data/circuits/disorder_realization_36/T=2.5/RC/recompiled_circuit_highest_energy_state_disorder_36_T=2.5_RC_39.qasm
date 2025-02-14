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
rz(3.1004768) q[0];
sx q[0];
rz(-1.2007204) q[0];
sx q[0];
rz(-1.6706985) q[0];
rz(-0.38493758) q[1];
sx q[1];
rz(-0.5783143) q[1];
sx q[1];
rz(2.2338423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9258869) q[0];
sx q[0];
rz(-0.74630794) q[0];
sx q[0];
rz(2.4624636) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8706246) q[2];
sx q[2];
rz(-1.8524287) q[2];
sx q[2];
rz(-0.23164888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.585184) q[1];
sx q[1];
rz(-2.6165462) q[1];
sx q[1];
rz(-2.7624112) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0201245) q[3];
sx q[3];
rz(-1.6252015) q[3];
sx q[3];
rz(-1.7732946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82751983) q[2];
sx q[2];
rz(-1.6562409) q[2];
sx q[2];
rz(-0.32076389) q[2];
rz(0.00094207923) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(0.48085406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568521) q[0];
sx q[0];
rz(-2.4212403) q[0];
sx q[0];
rz(2.5208933) q[0];
rz(1.8547828) q[1];
sx q[1];
rz(-1.6544673) q[1];
sx q[1];
rz(-2.1466045) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5704203) q[0];
sx q[0];
rz(-1.3745527) q[0];
sx q[0];
rz(-1.061383) q[0];
rz(-pi) q[1];
rz(0.32032712) q[2];
sx q[2];
rz(-1.1503845) q[2];
sx q[2];
rz(-1.0826031) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.4006808) q[1];
sx q[1];
rz(-1.7550635) q[1];
sx q[1];
rz(-2.0234985) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6334559) q[3];
sx q[3];
rz(-2.2656815) q[3];
sx q[3];
rz(0.14665237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2113125) q[2];
sx q[2];
rz(-1.899753) q[2];
sx q[2];
rz(-0.23107432) q[2];
rz(-1.6127582) q[3];
sx q[3];
rz(-1.4569747) q[3];
sx q[3];
rz(-2.5902364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.424161) q[0];
sx q[0];
rz(-1.9743974) q[0];
sx q[0];
rz(0.396808) q[0];
rz(1.1948168) q[1];
sx q[1];
rz(-1.2620986) q[1];
sx q[1];
rz(1.4746812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68434381) q[0];
sx q[0];
rz(-1.5822344) q[0];
sx q[0];
rz(-1.4829163) q[0];
x q[1];
rz(0.54946396) q[2];
sx q[2];
rz(-2.0026752) q[2];
sx q[2];
rz(0.23572505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6474698) q[1];
sx q[1];
rz(-0.9174594) q[1];
sx q[1];
rz(2.7092169) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1766522) q[3];
sx q[3];
rz(-1.5598462) q[3];
sx q[3];
rz(-1.2436858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2954703) q[2];
sx q[2];
rz(-0.84540558) q[2];
sx q[2];
rz(-1.6579312) q[2];
rz(-2.1814003) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(-0.6944164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8089777) q[0];
sx q[0];
rz(-1.806145) q[0];
sx q[0];
rz(-2.4131925) q[0];
rz(-1.2737466) q[1];
sx q[1];
rz(-1.961901) q[1];
sx q[1];
rz(-1.5789998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9804763) q[0];
sx q[0];
rz(-0.022810629) q[0];
sx q[0];
rz(-1.423832) q[0];
rz(-pi) q[1];
rz(-2.118763) q[2];
sx q[2];
rz(-1.5336516) q[2];
sx q[2];
rz(1.1556582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3660903) q[1];
sx q[1];
rz(-1.3555158) q[1];
sx q[1];
rz(-1.7083113) q[1];
rz(-2.8126841) q[3];
sx q[3];
rz(-2.2136627) q[3];
sx q[3];
rz(2.7907284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1100715) q[2];
sx q[2];
rz(-1.7379802) q[2];
sx q[2];
rz(2.9384379) q[2];
rz(1.7848232) q[3];
sx q[3];
rz(-1.2099096) q[3];
sx q[3];
rz(1.320896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88306952) q[0];
sx q[0];
rz(-1.8061545) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(0.37172231) q[1];
sx q[1];
rz(-0.55551353) q[1];
sx q[1];
rz(-0.98977596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5842642) q[0];
sx q[0];
rz(-2.4788878) q[0];
sx q[0];
rz(1.2602379) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3690574) q[2];
sx q[2];
rz(-1.651911) q[2];
sx q[2];
rz(-0.81743956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.02640662) q[1];
sx q[1];
rz(-2.7640928) q[1];
sx q[1];
rz(2.0707692) q[1];
rz(2.6017461) q[3];
sx q[3];
rz(-2.5858736) q[3];
sx q[3];
rz(-0.49969765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.524579) q[2];
sx q[2];
rz(-2.9144574) q[2];
sx q[2];
rz(3.0830834) q[2];
rz(1.2837563) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(-1.9730998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44130317) q[0];
sx q[0];
rz(-2.0952201) q[0];
sx q[0];
rz(0.48536479) q[0];
rz(-0.46733388) q[1];
sx q[1];
rz(-0.58715564) q[1];
sx q[1];
rz(0.68624085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94636665) q[0];
sx q[0];
rz(-1.2513046) q[0];
sx q[0];
rz(0.736306) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9250018) q[2];
sx q[2];
rz(-2.1859043) q[2];
sx q[2];
rz(1.5647581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8778692) q[1];
sx q[1];
rz(-0.72359607) q[1];
sx q[1];
rz(2.3797026) q[1];
rz(-pi) q[2];
rz(-2.3090906) q[3];
sx q[3];
rz(-2.0568536) q[3];
sx q[3];
rz(1.6532236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83963362) q[2];
sx q[2];
rz(-1.0884103) q[2];
sx q[2];
rz(-0.23183091) q[2];
rz(-2.2202282) q[3];
sx q[3];
rz(-1.5270343) q[3];
sx q[3];
rz(0.025402633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63474083) q[0];
sx q[0];
rz(-1.0036108) q[0];
sx q[0];
rz(-2.4941709) q[0];
rz(-2.1000775) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(2.0288859) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10202399) q[0];
sx q[0];
rz(-2.0687851) q[0];
sx q[0];
rz(2.6302393) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9693807) q[2];
sx q[2];
rz(-0.1548793) q[2];
sx q[2];
rz(-1.9914371) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3379593) q[1];
sx q[1];
rz(-0.58729672) q[1];
sx q[1];
rz(-0.31897591) q[1];
x q[2];
rz(1.5912368) q[3];
sx q[3];
rz(-2.6764565) q[3];
sx q[3];
rz(-2.6488369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.57704321) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(0.88456336) q[2];
rz(2.1029643) q[3];
sx q[3];
rz(-1.1218772) q[3];
sx q[3];
rz(-1.5326726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3603947) q[0];
sx q[0];
rz(-0.21112694) q[0];
sx q[0];
rz(-2.0017083) q[0];
rz(-2.3431011) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(-3.0671157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2452364) q[0];
sx q[0];
rz(-2.470861) q[0];
sx q[0];
rz(1.1719431) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3788667) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(0.21108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9556314) q[1];
sx q[1];
rz(-2.4096386) q[1];
sx q[1];
rz(0.47853985) q[1];
x q[2];
rz(-1.2182323) q[3];
sx q[3];
rz(-1.2855913) q[3];
sx q[3];
rz(2.6328301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8188339) q[2];
sx q[2];
rz(-2.5311311) q[2];
sx q[2];
rz(-2.9739213) q[2];
rz(-2.3882315) q[3];
sx q[3];
rz(-1.916901) q[3];
sx q[3];
rz(2.0577478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7849671) q[0];
sx q[0];
rz(-1.4545119) q[0];
sx q[0];
rz(1.0474569) q[0];
rz(-2.6197703) q[1];
sx q[1];
rz(-1.9597041) q[1];
sx q[1];
rz(0.80400115) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76011953) q[0];
sx q[0];
rz(-1.1520885) q[0];
sx q[0];
rz(0.19180723) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1754403) q[2];
sx q[2];
rz(-2.1198303) q[2];
sx q[2];
rz(-2.6556478) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2968353) q[1];
sx q[1];
rz(-1.8342819) q[1];
sx q[1];
rz(2.6729463) q[1];
rz(0.70254962) q[3];
sx q[3];
rz(-1.8773729) q[3];
sx q[3];
rz(2.2114913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1890586) q[2];
sx q[2];
rz(-1.6510094) q[2];
sx q[2];
rz(-0.35183364) q[2];
rz(-0.072988836) q[3];
sx q[3];
rz(-1.1806386) q[3];
sx q[3];
rz(2.4019057) q[3];
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
rz(-0.22013448) q[0];
sx q[0];
rz(-0.84311068) q[0];
sx q[0];
rz(1.5716918) q[0];
rz(-0.72228557) q[1];
sx q[1];
rz(-1.1415569) q[1];
sx q[1];
rz(1.3977745) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36436468) q[0];
sx q[0];
rz(-1.1684864) q[0];
sx q[0];
rz(-1.2978294) q[0];
rz(-pi) q[1];
rz(2.8768853) q[2];
sx q[2];
rz(-2.4183344) q[2];
sx q[2];
rz(-0.25991752) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.074894431) q[1];
sx q[1];
rz(-0.5579307) q[1];
sx q[1];
rz(2.9900223) q[1];
rz(-pi) q[2];
rz(0.13171036) q[3];
sx q[3];
rz(-1.9870821) q[3];
sx q[3];
rz(-0.10271969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7950644) q[2];
sx q[2];
rz(-1.8300191) q[2];
sx q[2];
rz(1.8468924) q[2];
rz(-2.0857701) q[3];
sx q[3];
rz(-1.0603696) q[3];
sx q[3];
rz(1.0052217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1718564) q[0];
sx q[0];
rz(-2.2315401) q[0];
sx q[0];
rz(1.6730614) q[0];
rz(-0.83364529) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(-0.40281725) q[2];
sx q[2];
rz(-2.2765712) q[2];
sx q[2];
rz(1.2650875) q[2];
rz(2.7293495) q[3];
sx q[3];
rz(-1.9751679) q[3];
sx q[3];
rz(-0.54505596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
