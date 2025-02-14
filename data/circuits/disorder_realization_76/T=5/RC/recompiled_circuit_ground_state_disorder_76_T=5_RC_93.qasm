OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4467093) q[0];
sx q[0];
rz(-2.0599685) q[0];
sx q[0];
rz(-2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(1.4256328) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.254296) q[0];
sx q[0];
rz(-1.0440517) q[0];
sx q[0];
rz(2.8077784) q[0];
rz(-1.8427467) q[2];
sx q[2];
rz(-2.0520979) q[2];
sx q[2];
rz(-2.0331733) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1358429) q[1];
sx q[1];
rz(-1.7187566) q[1];
sx q[1];
rz(-0.078385014) q[1];
x q[2];
rz(1.6207148) q[3];
sx q[3];
rz(-2.422159) q[3];
sx q[3];
rz(1.0843383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90929675) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(-0.24093957) q[2];
rz(-3.1124034) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(1.9624814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772684) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(0.48467317) q[0];
rz(-0.67972216) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(-1.1601123) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77185936) q[0];
sx q[0];
rz(-1.5892649) q[0];
sx q[0];
rz(2.8375677) q[0];
x q[1];
rz(2.4258575) q[2];
sx q[2];
rz(-1.01075) q[2];
sx q[2];
rz(2.3079688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94494941) q[1];
sx q[1];
rz(-1.3232051) q[1];
sx q[1];
rz(0.61243122) q[1];
rz(-2.6061222) q[3];
sx q[3];
rz(-0.59038416) q[3];
sx q[3];
rz(1.930742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29740563) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(-2.6837132) q[2];
rz(-1.1229905) q[3];
sx q[3];
rz(-0.99959683) q[3];
sx q[3];
rz(-0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726525) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(-0.031524468) q[0];
rz(-0.28745502) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.9256176) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8376802) q[0];
sx q[0];
rz(-1.0420351) q[0];
sx q[0];
rz(-2.3882475) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7599665) q[2];
sx q[2];
rz(-2.6256621) q[2];
sx q[2];
rz(2.6037773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.029473631) q[1];
sx q[1];
rz(-2.748317) q[1];
sx q[1];
rz(0.78177364) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9312882) q[3];
sx q[3];
rz(-2.1289325) q[3];
sx q[3];
rz(-1.0071514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1027801) q[2];
sx q[2];
rz(-1.5993885) q[2];
sx q[2];
rz(2.3056324) q[2];
rz(-1.3577667) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(-2.3939705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8589856) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(3.0614241) q[0];
rz(2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(1.3287883) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.778085) q[0];
sx q[0];
rz(-0.62113614) q[0];
sx q[0];
rz(-0.58025324) q[0];
x q[1];
rz(1.7432418) q[2];
sx q[2];
rz(-1.2524464) q[2];
sx q[2];
rz(0.47296745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7419974) q[1];
sx q[1];
rz(-2.6873702) q[1];
sx q[1];
rz(1.9059559) q[1];
rz(-0.33578797) q[3];
sx q[3];
rz(-1.3779252) q[3];
sx q[3];
rz(1.646281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(-0.67029101) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(-1.7214187) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-0.43148828) q[0];
rz(1.1314499) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-0.44050899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1898855) q[0];
sx q[0];
rz(-2.0312956) q[0];
sx q[0];
rz(3.0104464) q[0];
rz(-pi) q[1];
rz(-1.126664) q[2];
sx q[2];
rz(-0.45373771) q[2];
sx q[2];
rz(0.96378381) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7311598) q[1];
sx q[1];
rz(-1.9010547) q[1];
sx q[1];
rz(-1.6124658) q[1];
rz(-2.308485) q[3];
sx q[3];
rz(-2.3596977) q[3];
sx q[3];
rz(1.1469008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(0.41395536) q[2];
rz(1.3881989) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(-3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269161) q[0];
sx q[0];
rz(-1.7358945) q[0];
sx q[0];
rz(-2.2077014) q[0];
rz(2.4712708) q[1];
sx q[1];
rz(-1.8972081) q[1];
sx q[1];
rz(-1.3353039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0224755) q[0];
sx q[0];
rz(-0.60998218) q[0];
sx q[0];
rz(1.7458292) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9839758) q[2];
sx q[2];
rz(-0.15963563) q[2];
sx q[2];
rz(0.0052513382) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28999968) q[1];
sx q[1];
rz(-1.6834752) q[1];
sx q[1];
rz(3.1188909) q[1];
x q[2];
rz(0.63702668) q[3];
sx q[3];
rz(-1.5777794) q[3];
sx q[3];
rz(-2.1050326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2495217) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(1.2707233) q[2];
rz(-1.3696085) q[3];
sx q[3];
rz(-1.2415875) q[3];
sx q[3];
rz(0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524566) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(2.9879046) q[0];
rz(0.20248374) q[1];
sx q[1];
rz(-1.2362365) q[1];
sx q[1];
rz(-0.94857803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21009203) q[0];
sx q[0];
rz(-2.4903959) q[0];
sx q[0];
rz(0.33588846) q[0];
rz(-pi) q[1];
rz(-1.9948629) q[2];
sx q[2];
rz(-0.68409195) q[2];
sx q[2];
rz(2.7685194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2413509) q[1];
sx q[1];
rz(-0.35120041) q[1];
sx q[1];
rz(0.83937774) q[1];
rz(-1.7465215) q[3];
sx q[3];
rz(-2.4843289) q[3];
sx q[3];
rz(-2.6906309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.004868) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(-3.1401805) q[2];
rz(-3.0794365) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(0.96127659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75981265) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(2.3836721) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(-0.55799276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787631) q[0];
sx q[0];
rz(-2.7192594) q[0];
sx q[0];
rz(0.93393737) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5900389) q[2];
sx q[2];
rz(-1.4467014) q[2];
sx q[2];
rz(-0.53260224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26436603) q[1];
sx q[1];
rz(-2.5747882) q[1];
sx q[1];
rz(1.0344489) q[1];
rz(-pi) q[2];
rz(-3.1147546) q[3];
sx q[3];
rz(-1.2937814) q[3];
sx q[3];
rz(-2.7153496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5743635) q[2];
sx q[2];
rz(-1.9476674) q[2];
sx q[2];
rz(-0.26926678) q[2];
rz(-1.9783798) q[3];
sx q[3];
rz(-0.60267699) q[3];
sx q[3];
rz(-0.24063024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729386) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(-2.0682251) q[0];
rz(-1.1277554) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(-2.5849297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23281413) q[0];
sx q[0];
rz(-1.4289843) q[0];
sx q[0];
rz(-1.2491666) q[0];
x q[1];
rz(-2.5684909) q[2];
sx q[2];
rz(-2.0656043) q[2];
sx q[2];
rz(-3.0798517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5334415) q[1];
sx q[1];
rz(-1.9284298) q[1];
sx q[1];
rz(0.74700345) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4707556) q[3];
sx q[3];
rz(-1.3373858) q[3];
sx q[3];
rz(-1.4518713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90159455) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(-1.979801) q[2];
rz(2.6084172) q[3];
sx q[3];
rz(-2.6170001) q[3];
sx q[3];
rz(0.1575135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48851442) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(0.59481204) q[0];
rz(2.5626903) q[1];
sx q[1];
rz(-1.6659104) q[1];
sx q[1];
rz(-2.4822809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4730277) q[0];
sx q[0];
rz(-2.1921232) q[0];
sx q[0];
rz(-1.4061808) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5393007) q[2];
sx q[2];
rz(-2.4458439) q[2];
sx q[2];
rz(-0.93590036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.092992) q[1];
sx q[1];
rz(-1.8151746) q[1];
sx q[1];
rz(2.9952769) q[1];
rz(-1.8439383) q[3];
sx q[3];
rz(-1.8255263) q[3];
sx q[3];
rz(1.2431074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.351563) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(-3.0885546) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.2767295) q[3];
sx q[3];
rz(-3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85528436) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(0.098943624) q[1];
sx q[1];
rz(-1.3700486) q[1];
sx q[1];
rz(0.055421967) q[1];
rz(-1.5621875) q[2];
sx q[2];
rz(-1.9762403) q[2];
sx q[2];
rz(2.5019675) q[2];
rz(2.4002038) q[3];
sx q[3];
rz(-1.4546053) q[3];
sx q[3];
rz(-1.4447053) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
