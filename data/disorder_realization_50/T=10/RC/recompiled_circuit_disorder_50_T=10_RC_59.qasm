OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(3.1403132) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(4.2445634) q[1];
sx q[1];
rz(7.0581262) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2431508) q[0];
sx q[0];
rz(-1.4282707) q[0];
sx q[0];
rz(-1.6210763) q[0];
x q[1];
rz(1.0983659) q[2];
sx q[2];
rz(-2.35765) q[2];
sx q[2];
rz(-0.69477615) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1577658) q[1];
sx q[1];
rz(-2.9938934) q[1];
sx q[1];
rz(-1.8450792) q[1];
rz(-2.6585456) q[3];
sx q[3];
rz(-2.8222198) q[3];
sx q[3];
rz(-1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(-1.1958896) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130106) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(-2.9257665) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0963418) q[2];
sx q[2];
rz(-1.8088733) q[2];
sx q[2];
rz(2.0766052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3017722) q[1];
sx q[1];
rz(-2.7313247) q[1];
sx q[1];
rz(-2.6778585) q[1];
rz(-2.1249173) q[3];
sx q[3];
rz(-1.2614054) q[3];
sx q[3];
rz(-0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.2352357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9440207) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(0.59528415) q[0];
x q[1];
rz(-1.8232934) q[2];
sx q[2];
rz(-1.9324979) q[2];
sx q[2];
rz(0.30644882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5879844) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(-1.8334465) q[1];
x q[2];
rz(2.5123185) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(-0.34793684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(-1.7049449) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(0.11985699) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0940995) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(-0.043461965) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4255964) q[2];
sx q[2];
rz(-1.928066) q[2];
sx q[2];
rz(-1.4841339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.90384342) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(-3.1234427) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1060137) q[3];
sx q[3];
rz(-1.03973) q[3];
sx q[3];
rz(0.5862743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(0.37483254) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(1.5572085) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.24057062) q[2];
sx q[2];
rz(-1.3153937) q[2];
sx q[2];
rz(2.0672928) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1357437) q[1];
sx q[1];
rz(-1.6372576) q[1];
sx q[1];
rz(2.5367516) q[1];
x q[2];
rz(-0.9661071) q[3];
sx q[3];
rz(-1.1629512) q[3];
sx q[3];
rz(0.84433489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(-2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(-0.99037209) q[0];
rz(3.0888427) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(-1.6606768) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(1.7646164) q[0];
x q[1];
rz(2.0159857) q[2];
sx q[2];
rz(-1.5714958) q[2];
sx q[2];
rz(-0.47847754) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78696886) q[1];
sx q[1];
rz(-1.0458535) q[1];
sx q[1];
rz(0.24137361) q[1];
x q[2];
rz(2.2495066) q[3];
sx q[3];
rz(-1.5329554) q[3];
sx q[3];
rz(-0.91176646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0084373077) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(-2.0810614) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.3247103) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0990471) q[0];
sx q[0];
rz(-1.5934266) q[0];
sx q[0];
rz(-0.12734539) q[0];
x q[1];
rz(-2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(1.2457459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.054338) q[1];
sx q[1];
rz(-1.195765) q[1];
sx q[1];
rz(0.43941984) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38447325) q[3];
sx q[3];
rz(-0.95015804) q[3];
sx q[3];
rz(1.4409161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(2.1929072) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.8519648) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063110654) q[0];
sx q[0];
rz(-2.1266298) q[0];
sx q[0];
rz(-2.1102935) q[0];
x q[1];
rz(-0.82215674) q[2];
sx q[2];
rz(-0.90688721) q[2];
sx q[2];
rz(-0.31809959) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(0.11282632) q[1];
x q[2];
rz(0.21059489) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(1.3099561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(-3.0454214) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44336549) q[0];
sx q[0];
rz(-0.8015612) q[0];
sx q[0];
rz(-0.65936868) q[0];
rz(-pi) q[1];
rz(-0.56844175) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(2.8708411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54284755) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(-2.7302242) q[1];
rz(1.4275527) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(-2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(0.50869554) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1334575) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(-2.125678) q[0];
rz(2.290906) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(1.016664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21916418) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(1.0132755) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52313157) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(-0.79062068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(2.534261) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(-2.5088359) q[2];
sx q[2];
rz(-0.7706332) q[2];
sx q[2];
rz(2.3494233) q[2];
rz(1.3689465) q[3];
sx q[3];
rz(-1.1931843) q[3];
sx q[3];
rz(-1.9184792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
