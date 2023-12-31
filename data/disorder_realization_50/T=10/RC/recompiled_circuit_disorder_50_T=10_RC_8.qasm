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
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(2.3666518) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9026731) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(2.8047049) q[0];
x q[1];
rz(1.0983659) q[2];
sx q[2];
rz(-0.78394267) q[2];
sx q[2];
rz(-2.4468165) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4349164) q[1];
sx q[1];
rz(-1.7129363) q[1];
sx q[1];
rz(-3.1013156) q[1];
rz(2.6585456) q[3];
sx q[3];
rz(-0.31937283) q[3];
sx q[3];
rz(1.3787624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9346416) q[0];
sx q[0];
rz(-1.1622381) q[0];
sx q[0];
rz(1.6658528) q[0];
rz(1.3865115) q[2];
sx q[2];
rz(-0.24225907) q[2];
sx q[2];
rz(0.87528961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83982044) q[1];
sx q[1];
rz(-2.7313247) q[1];
sx q[1];
rz(-0.46373414) q[1];
rz(-2.7820285) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89861682) q[0];
sx q[0];
rz(-1.9447864) q[0];
sx q[0];
rz(0.6181194) q[0];
x q[1];
rz(0.37249506) q[2];
sx q[2];
rz(-1.8066346) q[2];
sx q[2];
rz(-1.3553938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.55360824) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(-1.3081461) q[1];
rz(2.0577621) q[3];
sx q[3];
rz(-0.999513) q[3];
sx q[3];
rz(0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76413313) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(-1.4366478) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(-0.60825545) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-0.80379379) q[0];
rz(-0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(3.0217357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6896497) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(-2.1769051) q[0];
x q[1];
rz(-2.7718133) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(-1.0880926) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67112982) q[1];
sx q[1];
rz(-1.5531335) q[1];
sx q[1];
rz(-1.8030241) q[1];
x q[2];
rz(0.5991163) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.1981296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(0.48686349) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70532521) q[0];
sx q[0];
rz(-1.3316532) q[0];
sx q[0];
rz(1.5843841) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83112049) q[2];
sx q[2];
rz(-2.7925425) q[2];
sx q[2];
rz(2.8380053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6606635) q[1];
sx q[1];
rz(-2.1741121) q[1];
sx q[1];
rz(1.4900581) q[1];
x q[2];
rz(0.48360444) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-0.70303482) q[2];
rz(1.4098343) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(-0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-2.2149448) q[1];
sx q[1];
rz(1.6606768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4617417) q[0];
sx q[0];
rz(-1.4977507) q[0];
sx q[0];
rz(-1.1887656) q[0];
rz(-pi) q[1];
rz(3.1408177) q[2];
sx q[2];
rz(-2.0159855) q[2];
sx q[2];
rz(-2.0489401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2350125) q[1];
sx q[1];
rz(-1.3624411) q[1];
sx q[1];
rz(-2.1085897) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2495066) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(2.2298262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.19787191) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0990471) q[0];
sx q[0];
rz(-1.5934266) q[0];
sx q[0];
rz(3.0142473) q[0];
rz(-2.6323694) q[2];
sx q[2];
rz(-2.074713) q[2];
sx q[2];
rz(-3.0798806) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3130256) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(-1.9810956) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7571194) q[3];
sx q[3];
rz(-0.95015804) q[3];
sx q[3];
rz(1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.8519648) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3279009) q[0];
sx q[0];
rz(-1.1192338) q[0];
sx q[0];
rz(0.62664647) q[0];
rz(2.3233534) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(-0.73275369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3370812) q[1];
sx q[1];
rz(-1.6816499) q[1];
sx q[1];
rz(-1.3827419) q[1];
rz(1.7361705) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(-2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-2.8073231) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-2.8589378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219014) q[0];
sx q[0];
rz(-2.0265409) q[0];
sx q[0];
rz(2.4569608) q[0];
rz(-pi) q[1];
rz(0.33090584) q[2];
sx q[2];
rz(-1.7754284) q[2];
sx q[2];
rz(-0.76588878) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0375367) q[1];
sx q[1];
rz(-1.1661342) q[1];
sx q[1];
rz(-1.3794823) q[1];
rz(-pi) q[2];
rz(-0.11741365) q[3];
sx q[3];
rz(-2.252929) q[3];
sx q[3];
rz(0.62473245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21215542) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(0.50869554) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(-2.4597816) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081351) q[0];
sx q[0];
rz(-1.4032161) q[0];
sx q[0];
rz(-2.125678) q[0];
rz(2.3186458) q[2];
sx q[2];
rz(-2.1087077) q[2];
sx q[2];
rz(3.092098) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1055923) q[1];
sx q[1];
rz(-2.3899374) q[1];
sx q[1];
rz(-2.4113301) q[1];
x q[2];
rz(-0.82716771) q[3];
sx q[3];
rz(-2.4359772) q[3];
sx q[3];
rz(-0.08882113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(-1.0836481) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(-0.60733168) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(-2.4773459) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
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
