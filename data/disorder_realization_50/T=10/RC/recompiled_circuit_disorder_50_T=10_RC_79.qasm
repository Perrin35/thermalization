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
rz(-2.0386219) q[1];
sx q[1];
rz(-2.3666518) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9026731) q[0];
sx q[0];
rz(-2.9905149) q[0];
sx q[0];
rz(2.8047049) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84471976) q[2];
sx q[2];
rz(-1.8978999) q[2];
sx q[2];
rz(-2.6127882) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4349164) q[1];
sx q[1];
rz(-1.4286563) q[1];
sx q[1];
rz(3.1013156) q[1];
x q[2];
rz(-2.8567021) q[3];
sx q[3];
rz(-1.4244411) q[3];
sx q[3];
rz(2.871606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7213223) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(1.3756479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1130106) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(0.21582614) q[0];
rz(0.045250821) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(1.0649875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1609636) q[1];
sx q[1];
rz(-1.7501608) q[1];
sx q[1];
rz(2.7705926) q[1];
rz(-2.1249173) q[3];
sx q[3];
rz(-1.8801873) q[3];
sx q[3];
rz(0.644899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21330825) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(-1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(2.074923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60004822) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.2352357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2429758) q[0];
sx q[0];
rz(-1.1968062) q[0];
sx q[0];
rz(-2.5234733) q[0];
x q[1];
rz(-1.3182993) q[2];
sx q[2];
rz(-1.2090948) q[2];
sx q[2];
rz(-2.8351438) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55360824) q[1];
sx q[1];
rz(-1.3097035) q[1];
sx q[1];
rz(1.3081461) q[1];
rz(1.0838305) q[3];
sx q[3];
rz(-2.1420797) q[3];
sx q[3];
rz(-2.7146102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(-1.4012339) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(-3.0217357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9712517) q[0];
sx q[0];
rz(-2.5345638) q[0];
sx q[0];
rz(-1.6334565) q[0];
rz(2.7718133) q[2];
sx q[2];
rz(-2.7571207) q[2];
sx q[2];
rz(-1.0880926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3164697) q[1];
sx q[1];
rz(-2.9087062) q[1];
sx q[1];
rz(1.6474001) q[1];
x q[2];
rz(2.4265392) q[3];
sx q[3];
rz(-2.4063769) q[3];
sx q[3];
rz(0.27763593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(0.072337739) q[2];
rz(-0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(1.9015076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86868984) q[0];
sx q[0];
rz(-1.5839974) q[0];
sx q[0];
rz(-2.9024283) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8334332) q[2];
sx q[2];
rz(-1.3381759) q[2];
sx q[2];
rz(2.5831985) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48092914) q[1];
sx q[1];
rz(-2.1741121) q[1];
sx q[1];
rz(-1.6515345) q[1];
x q[2];
rz(0.9661071) q[3];
sx q[3];
rz(-1.1629512) q[3];
sx q[3];
rz(-0.84433489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(0.15771244) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96419656) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(-2.1512206) q[0];
rz(-3.0888427) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(-1.6606768) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67985095) q[0];
sx q[0];
rz(-1.6438419) q[0];
sx q[0];
rz(1.1887656) q[0];
x q[1];
rz(1.5724206) q[2];
sx q[2];
rz(-0.44518984) q[2];
sx q[2];
rz(-1.0908529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2350125) q[1];
sx q[1];
rz(-1.3624411) q[1];
sx q[1];
rz(2.1085897) q[1];
rz(-pi) q[2];
rz(0.89208608) q[3];
sx q[3];
rz(-1.5329554) q[3];
sx q[3];
rz(0.91176646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(0.56224242) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.3247103) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46885195) q[0];
sx q[0];
rz(-1.4434837) q[0];
sx q[0];
rz(-1.5936113) q[0];
rz(-pi) q[1];
rz(-1.0075188) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(-1.2457459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.054338) q[1];
sx q[1];
rz(-1.195765) q[1];
sx q[1];
rz(2.7021728) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7571194) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(1.7006765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7396486) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(-3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.8519648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063110654) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(-1.0312992) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3194359) q[2];
sx q[2];
rz(-2.2347054) q[2];
sx q[2];
rz(-2.8234931) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21266567) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(0.11282632) q[1];
rz(1.7361705) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(-2.0986433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.148968) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6982272) q[0];
sx q[0];
rz(-0.8015612) q[0];
sx q[0];
rz(-2.482224) q[0];
x q[1];
rz(-0.56844175) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(2.8708411) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5987451) q[1];
sx q[1];
rz(-1.3951021) q[1];
sx q[1];
rz(-2.7302242) q[1];
rz(-2.2563124) q[3];
sx q[3];
rz(-1.6618528) q[3];
sx q[3];
rz(-2.1212999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9294372) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(-2.6397928) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9671229) q[0];
sx q[0];
rz(-2.5645064) q[0];
sx q[0];
rz(-1.2601) q[0];
rz(-0.82294686) q[2];
sx q[2];
rz(-2.1087077) q[2];
sx q[2];
rz(3.092098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0444813) q[1];
sx q[1];
rz(-1.0978571) q[1];
sx q[1];
rz(0.60826917) q[1];
rz(-pi) q[2];
rz(-2.6184611) q[3];
sx q[3];
rz(-2.0683859) q[3];
sx q[3];
rz(-0.79062068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(2.8005023) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(-0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682128) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(0.66424673) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
rz(1.7726462) q[3];
sx q[3];
rz(-1.9484083) q[3];
sx q[3];
rz(1.2231135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
