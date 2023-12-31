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
rz(-0.0012794415) q[0];
rz(-1.6969504) q[1];
sx q[1];
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9026731) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(-2.8047049) q[0];
rz(-pi) q[1];
rz(-1.0983659) q[2];
sx q[2];
rz(-0.78394267) q[2];
sx q[2];
rz(-0.69477615) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2831813) q[1];
sx q[1];
rz(-1.6106669) q[1];
sx q[1];
rz(-1.4285425) q[1];
x q[2];
rz(2.8567021) q[3];
sx q[3];
rz(-1.7171515) q[3];
sx q[3];
rz(2.871606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(-1.7529863) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7213223) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401706) q[0];
sx q[0];
rz(-1.4835843) q[0];
sx q[0];
rz(2.7313822) q[0];
rz(-1.7550811) q[2];
sx q[2];
rz(-0.24225907) q[2];
sx q[2];
rz(-2.266303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3017722) q[1];
sx q[1];
rz(-2.7313247) q[1];
sx q[1];
rz(2.6778585) q[1];
x q[2];
rz(-1.0248915) q[3];
sx q[3];
rz(-0.62666577) q[3];
sx q[3];
rz(2.6729667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.7958612) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(2.0080163) q[3];
sx q[3];
rz(-0.19530345) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(1.906357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19757195) q[0];
sx q[0];
rz(-2.4320142) q[0];
sx q[0];
rz(0.59528415) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58358242) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(-0.32354087) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55360824) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(-1.8334465) q[1];
rz(-pi) q[2];
rz(2.0577621) q[3];
sx q[3];
rz(-2.1420797) q[3];
sx q[3];
rz(-0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-2.5333372) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.5290715) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9712517) q[0];
sx q[0];
rz(-2.5345638) q[0];
sx q[0];
rz(1.6334565) q[0];
x q[1];
rz(0.36074952) q[2];
sx q[2];
rz(-1.7067688) q[2];
sx q[2];
rz(3.0038358) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2377492) q[1];
sx q[1];
rz(-1.8029873) q[1];
sx q[1];
rz(-0.018149908) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0355789) q[3];
sx q[3];
rz(-1.03973) q[3];
sx q[3];
rz(2.5553184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(-0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(-1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(0.48686349) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(1.9015076) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70532521) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(-1.5572085) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3081595) q[2];
sx q[2];
rz(-1.8034168) q[2];
sx q[2];
rz(2.5831985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.33927321) q[1];
sx q[1];
rz(-0.60802751) q[1];
sx q[1];
rz(-0.1165216) q[1];
x q[2];
rz(2.6579882) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(0.70303482) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96419656) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(-0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(-1.4809158) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.7646164) q[0];
rz(-0.00077498282) q[2];
sx q[2];
rz(-2.0159855) q[2];
sx q[2];
rz(-2.0489401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90658014) q[1];
sx q[1];
rz(-1.3624411) q[1];
sx q[1];
rz(-2.1085897) q[1];
rz(-0.89208608) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(0.91176646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(-1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.4690171) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.3247103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0990471) q[0];
sx q[0];
rz(-1.5934266) q[0];
sx q[0];
rz(-0.12734539) q[0];
rz(-pi) q[1];
rz(-2.6323694) q[2];
sx q[2];
rz(-2.074713) q[2];
sx q[2];
rz(-3.0798806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8285671) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(-1.9810956) q[1];
x q[2];
rz(-0.38447325) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(-1.2896279) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78532082) q[0];
sx q[0];
rz(-0.7542146) q[0];
sx q[0];
rz(0.69099364) q[0];
rz(-pi) q[1];
rz(-2.3233534) q[2];
sx q[2];
rz(-1.0050251) q[2];
sx q[2];
rz(-0.73275369) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.928927) q[1];
sx q[1];
rz(-1.3839098) q[1];
sx q[1];
rz(0.11282632) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4054221) q[3];
sx q[3];
rz(-0.91455063) q[3];
sx q[3];
rz(1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(-3.0454214) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-2.8589378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3960421) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(1.0066443) q[0];
x q[1];
rz(0.56844175) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(-2.8708411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.54284755) q[1];
sx q[1];
rz(-1.3951021) q[1];
sx q[1];
rz(0.41136841) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88528021) q[3];
sx q[3];
rz(-1.6618528) q[3];
sx q[3];
rz(-2.1212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(0.7412509) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(0.68181109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1334575) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(2.125678) q[0];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-1.0444813) q[1];
sx q[1];
rz(-2.0437355) q[1];
sx q[1];
rz(2.5333235) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52313157) q[3];
sx q[3];
rz(-1.0732068) q[3];
sx q[3];
rz(0.79062068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(0.63275679) q[2];
sx q[2];
rz(-0.7706332) q[2];
sx q[2];
rz(2.3494233) q[2];
rz(-1.7726462) q[3];
sx q[3];
rz(-1.1931843) q[3];
sx q[3];
rz(-1.9184792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
