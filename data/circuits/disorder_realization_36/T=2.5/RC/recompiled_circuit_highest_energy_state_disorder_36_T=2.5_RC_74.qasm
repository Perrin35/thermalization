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
rz(-0.38493758) q[1];
sx q[1];
rz(2.5632783) q[1];
sx q[1];
rz(10.332528) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9258869) q[0];
sx q[0];
rz(-2.3952847) q[0];
sx q[0];
rz(-2.4624636) q[0];
rz(1.8706246) q[2];
sx q[2];
rz(-1.8524287) q[2];
sx q[2];
rz(0.23164888) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31763916) q[1];
sx q[1];
rz(-1.7574218) q[1];
sx q[1];
rz(2.6479048) q[1];
rz(0.42239093) q[3];
sx q[3];
rz(-0.13304079) q[3];
sx q[3];
rz(0.62158004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82751983) q[2];
sx q[2];
rz(-1.6562409) q[2];
sx q[2];
rz(0.32076389) q[2];
rz(3.1406506) q[3];
sx q[3];
rz(-0.92206803) q[3];
sx q[3];
rz(2.6607386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8568521) q[0];
sx q[0];
rz(-0.72035235) q[0];
sx q[0];
rz(2.5208933) q[0];
rz(-1.2868098) q[1];
sx q[1];
rz(-1.4871253) q[1];
sx q[1];
rz(2.1466045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10886562) q[0];
sx q[0];
rz(-1.0720729) q[0];
sx q[0];
rz(-2.9176955) q[0];
rz(-pi) q[1];
rz(-0.95717818) q[2];
sx q[2];
rz(-0.52268302) q[2];
sx q[2];
rz(1.7651287) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7409119) q[1];
sx q[1];
rz(-1.7550635) q[1];
sx q[1];
rz(-2.0234985) q[1];
rz(-pi) q[2];
rz(1.0424588) q[3];
sx q[3];
rz(-2.3064559) q[3];
sx q[3];
rz(-2.5723815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2113125) q[2];
sx q[2];
rz(-1.899753) q[2];
sx q[2];
rz(-2.9105183) q[2];
rz(-1.5288345) q[3];
sx q[3];
rz(-1.4569747) q[3];
sx q[3];
rz(2.5902364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.7174317) q[0];
sx q[0];
rz(-1.1671952) q[0];
sx q[0];
rz(0.396808) q[0];
rz(-1.1948168) q[1];
sx q[1];
rz(-1.2620986) q[1];
sx q[1];
rz(1.6669115) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75735649) q[0];
sx q[0];
rz(-0.088619329) q[0];
sx q[0];
rz(-1.4411974) q[0];
rz(-pi) q[1];
rz(-2.5921287) q[2];
sx q[2];
rz(-2.0026752) q[2];
sx q[2];
rz(0.23572505) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7914305) q[1];
sx q[1];
rz(-1.2315948) q[1];
sx q[1];
rz(2.2712525) q[1];
x q[2];
rz(3.1297333) q[3];
sx q[3];
rz(-1.1766772) q[3];
sx q[3];
rz(0.32255641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2954703) q[2];
sx q[2];
rz(-2.2961871) q[2];
sx q[2];
rz(-1.4836614) q[2];
rz(-2.1814003) q[3];
sx q[3];
rz(-0.63513297) q[3];
sx q[3];
rz(-0.6944164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33261499) q[0];
sx q[0];
rz(-1.3354477) q[0];
sx q[0];
rz(-2.4131925) q[0];
rz(-1.867846) q[1];
sx q[1];
rz(-1.1796917) q[1];
sx q[1];
rz(-1.5789998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.584986) q[0];
sx q[0];
rz(-1.5741363) q[0];
sx q[0];
rz(-1.5933611) q[0];
x q[1];
rz(-0.043508675) q[2];
sx q[2];
rz(-1.0232506) q[2];
sx q[2];
rz(0.39247733) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3660903) q[1];
sx q[1];
rz(-1.7860768) q[1];
sx q[1];
rz(-1.4332814) q[1];
x q[2];
rz(-2.2402917) q[3];
sx q[3];
rz(-1.309295) q[3];
sx q[3];
rz(-1.4217564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.031521156) q[2];
sx q[2];
rz(-1.4036125) q[2];
sx q[2];
rz(-0.20315476) q[2];
rz(-1.7848232) q[3];
sx q[3];
rz(-1.9316831) q[3];
sx q[3];
rz(1.320896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2585231) q[0];
sx q[0];
rz(-1.3354381) q[0];
sx q[0];
rz(-1.6254599) q[0];
rz(-0.37172231) q[1];
sx q[1];
rz(-2.5860791) q[1];
sx q[1];
rz(-0.98977596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55732841) q[0];
sx q[0];
rz(-2.4788878) q[0];
sx q[0];
rz(1.2602379) q[0];
rz(2.3690574) q[2];
sx q[2];
rz(-1.4896817) q[2];
sx q[2];
rz(0.81743956) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.02640662) q[1];
sx q[1];
rz(-0.37749981) q[1];
sx q[1];
rz(1.0708234) q[1];
rz(-1.8797713) q[3];
sx q[3];
rz(-2.0403962) q[3];
sx q[3];
rz(-2.0275786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.524579) q[2];
sx q[2];
rz(-0.22713529) q[2];
sx q[2];
rz(0.058509286) q[2];
rz(1.8578364) q[3];
sx q[3];
rz(-2.3144898) q[3];
sx q[3];
rz(1.9730998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7002895) q[0];
sx q[0];
rz(-2.0952201) q[0];
sx q[0];
rz(-0.48536479) q[0];
rz(-2.6742588) q[1];
sx q[1];
rz(-0.58715564) q[1];
sx q[1];
rz(2.4553518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2398324) q[0];
sx q[0];
rz(-0.87945451) q[0];
sx q[0];
rz(-1.9907238) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9250018) q[2];
sx q[2];
rz(-2.1859043) q[2];
sx q[2];
rz(-1.5647581) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8279449) q[1];
sx q[1];
rz(-1.0961431) q[1];
sx q[1];
rz(-0.56874911) q[1];
x q[2];
rz(0.83250203) q[3];
sx q[3];
rz(-1.0847391) q[3];
sx q[3];
rz(-1.6532236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83963362) q[2];
sx q[2];
rz(-1.0884103) q[2];
sx q[2];
rz(0.23183091) q[2];
rz(-0.92136446) q[3];
sx q[3];
rz(-1.6145584) q[3];
sx q[3];
rz(0.025402633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5068518) q[0];
sx q[0];
rz(-2.1379819) q[0];
sx q[0];
rz(-0.64742175) q[0];
rz(-2.1000775) q[1];
sx q[1];
rz(-1.2728649) q[1];
sx q[1];
rz(2.0288859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1739832) q[0];
sx q[0];
rz(-2.4436473) q[0];
sx q[0];
rz(0.83794727) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0810705) q[2];
sx q[2];
rz(-1.7134482) q[2];
sx q[2];
rz(1.5530653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7156421) q[1];
sx q[1];
rz(-1.0167158) q[1];
sx q[1];
rz(1.7765846) q[1];
x q[2];
rz(1.5912368) q[3];
sx q[3];
rz(-0.4651362) q[3];
sx q[3];
rz(2.6488369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5645494) q[2];
sx q[2];
rz(-2.3304522) q[2];
sx q[2];
rz(2.2570293) q[2];
rz(1.0386284) q[3];
sx q[3];
rz(-1.1218772) q[3];
sx q[3];
rz(1.5326726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3603947) q[0];
sx q[0];
rz(-2.9304657) q[0];
sx q[0];
rz(-2.0017083) q[0];
rz(-2.3431011) q[1];
sx q[1];
rz(-1.4356828) q[1];
sx q[1];
rz(-3.0671157) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89635623) q[0];
sx q[0];
rz(-2.470861) q[0];
sx q[0];
rz(-1.9696495) q[0];
rz(-0.76272599) q[2];
sx q[2];
rz(-1.8755091) q[2];
sx q[2];
rz(-2.9305127) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1250389) q[1];
sx q[1];
rz(-1.8836262) q[1];
sx q[1];
rz(-0.67322634) q[1];
rz(-pi) q[2];
rz(-1.9233604) q[3];
sx q[3];
rz(-1.2855913) q[3];
sx q[3];
rz(0.5087626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3227587) q[2];
sx q[2];
rz(-0.61046159) q[2];
sx q[2];
rz(-0.16767137) q[2];
rz(2.3882315) q[3];
sx q[3];
rz(-1.2246917) q[3];
sx q[3];
rz(-1.0838449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.3566256) q[0];
sx q[0];
rz(-1.6870808) q[0];
sx q[0];
rz(-1.0474569) q[0];
rz(-2.6197703) q[1];
sx q[1];
rz(-1.9597041) q[1];
sx q[1];
rz(-2.3375915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31452824) q[0];
sx q[0];
rz(-0.4581764) q[0];
sx q[0];
rz(1.9755152) q[0];
rz(2.9661524) q[2];
sx q[2];
rz(-1.0217624) q[2];
sx q[2];
rz(-2.6556478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9404133) q[1];
sx q[1];
rz(-0.53279725) q[1];
sx q[1];
rz(-0.53838191) q[1];
x q[2];
rz(-1.9639765) q[3];
sx q[3];
rz(-2.2344974) q[3];
sx q[3];
rz(0.8908602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1890586) q[2];
sx q[2];
rz(-1.6510094) q[2];
sx q[2];
rz(-2.789759) q[2];
rz(3.0686038) q[3];
sx q[3];
rz(-1.9609541) q[3];
sx q[3];
rz(-2.4019057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22013448) q[0];
sx q[0];
rz(-0.84311068) q[0];
sx q[0];
rz(1.5716918) q[0];
rz(-2.4193071) q[1];
sx q[1];
rz(-1.1415569) q[1];
sx q[1];
rz(1.7438181) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.04434) q[0];
sx q[0];
rz(-1.821479) q[0];
sx q[0];
rz(-0.41608019) q[0];
x q[1];
rz(2.4358814) q[2];
sx q[2];
rz(-1.3967665) q[2];
sx q[2];
rz(-1.1104012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0666982) q[1];
sx q[1];
rz(-2.583662) q[1];
sx q[1];
rz(-0.15157031) q[1];
rz(-pi) q[2];
rz(3.0098823) q[3];
sx q[3];
rz(-1.1545106) q[3];
sx q[3];
rz(-0.10271969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7950644) q[2];
sx q[2];
rz(-1.3115735) q[2];
sx q[2];
rz(1.8468924) q[2];
rz(-2.0857701) q[3];
sx q[3];
rz(-2.081223) q[3];
sx q[3];
rz(-1.0052217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1718564) q[0];
sx q[0];
rz(-0.91005253) q[0];
sx q[0];
rz(-1.4685312) q[0];
rz(-0.83364529) q[1];
sx q[1];
rz(-1.2794762) q[1];
sx q[1];
rz(1.5247482) q[1];
rz(2.7387754) q[2];
sx q[2];
rz(-2.2765712) q[2];
sx q[2];
rz(1.2650875) q[2];
rz(-2.7293495) q[3];
sx q[3];
rz(-1.1664248) q[3];
sx q[3];
rz(2.5965367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
