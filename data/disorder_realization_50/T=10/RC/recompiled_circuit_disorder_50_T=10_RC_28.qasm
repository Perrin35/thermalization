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
rz(2.476386) q[0];
sx q[0];
rz(-1.6205661) q[0];
sx q[0];
rz(2.9988891) q[0];
x q[1];
rz(1.0983659) q[2];
sx q[2];
rz(-2.35765) q[2];
sx q[2];
rz(-0.69477615) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98382681) q[1];
sx q[1];
rz(-0.14769927) q[1];
sx q[1];
rz(1.8450792) q[1];
x q[2];
rz(-1.418387) q[3];
sx q[3];
rz(-1.8525575) q[3];
sx q[3];
rz(-1.8834653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.9457031) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.3886064) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4202704) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.8700245) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401706) q[0];
sx q[0];
rz(-1.6580083) q[0];
sx q[0];
rz(-0.4102104) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3324845) q[2];
sx q[2];
rz(-1.5268227) q[2];
sx q[2];
rz(-2.6251052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83982044) q[1];
sx q[1];
rz(-0.41026792) q[1];
sx q[1];
rz(2.6778585) q[1];
rz(-pi) q[2];
rz(-0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-1.1121225) q[3];
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
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-0.24060732) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.2352357) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2152527) q[0];
sx q[0];
rz(-2.140576) q[0];
sx q[0];
rz(1.1220054) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3182993) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(-2.8351438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55360824) q[1];
sx q[1];
rz(-1.8318892) q[1];
sx q[1];
rz(1.8334465) q[1];
rz(-pi) q[2];
rz(-0.62927411) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(2.7936558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(2.5333372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(3.0217357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.170341) q[0];
sx q[0];
rz(-0.60702885) q[0];
sx q[0];
rz(1.5081362) q[0];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(2.0535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67112982) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(1.8030241) q[1];
rz(-pi) q[2];
rz(0.5991163) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(-2.7667601) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(0.48686349) q[0];
rz(2.411719) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(1.9015076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4362674) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(-1.5572085) q[0];
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
rz(2.6606635) q[1];
sx q[1];
rz(-2.1741121) q[1];
sx q[1];
rz(-1.4900581) q[1];
x q[2];
rz(-0.9661071) q[3];
sx q[3];
rz(-1.9786414) q[3];
sx q[3];
rz(2.2972578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(-2.4385578) q[2];
rz(1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(-0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(-1.4809158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.3769763) q[0];
rz(-pi) q[1];
rz(2.0159857) q[2];
sx q[2];
rz(-1.5700969) q[2];
sx q[2];
rz(-2.6631151) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(0.24137361) q[1];
x q[2];
rz(-2.2495066) q[3];
sx q[3];
rz(-1.5329554) q[3];
sx q[3];
rz(-2.2298262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
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
rz(-1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0425456) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(0.12734539) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0075188) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(1.2457459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99609375) q[1];
sx q[1];
rz(-0.56963527) q[1];
sx q[1];
rz(0.74665229) q[1];
rz(-pi) q[2];
rz(-2.7571194) q[3];
sx q[3];
rz(-0.95015804) q[3];
sx q[3];
rz(1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(-0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396486) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(1.2896279) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.078482) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(-2.1102935) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3233534) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(2.408839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3813078) q[1];
sx q[1];
rz(-2.9236301) q[1];
sx q[1];
rz(-1.0337619) q[1];
rz(-pi) q[2];
rz(1.7361705) q[3];
sx q[3];
rz(-0.91455063) q[3];
sx q[3];
rz(-1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9926247) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(-0.096171245) q[3];
sx q[3];
rz(-1.440719) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(2.8589378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196913) q[0];
sx q[0];
rz(-2.0265409) q[0];
sx q[0];
rz(-2.4569608) q[0];
rz(-pi) q[1];
rz(-1.7868144) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(-0.87460364) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.104056) q[1];
sx q[1];
rz(-1.9754585) q[1];
sx q[1];
rz(1.3794823) q[1];
x q[2];
rz(-2.2563124) q[3];
sx q[3];
rz(-1.4797398) q[3];
sx q[3];
rz(-1.0202927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-0.58399502) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.042645) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(-0.11518654) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-0.68181109) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9671229) q[0];
sx q[0];
rz(-0.57708626) q[0];
sx q[0];
rz(1.8814927) q[0];
rz(-pi) q[1];
rz(2.3186458) q[2];
sx q[2];
rz(-2.1087077) q[2];
sx q[2];
rz(3.092098) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9224285) q[1];
sx q[1];
rz(-1.0370967) q[1];
sx q[1];
rz(-2.1283172) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0107332) q[3];
sx q[3];
rz(-1.1162973) q[3];
sx q[3];
rz(-1.0487995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(2.0579445) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1682128) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(1.0495937) q[2];
sx q[2];
rz(-0.97432077) q[2];
sx q[2];
rz(-1.5885098) q[2];
rz(2.6736005) q[3];
sx q[3];
rz(-2.7157126) q[3];
sx q[3];
rz(-1.4117905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];