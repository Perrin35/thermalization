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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(3.1290913) q[1];
sx q[1];
rz(-1.8379509) q[1];
sx q[1];
rz(-1.5703896) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16906315) q[0];
sx q[0];
rz(-2.1816945) q[0];
sx q[0];
rz(-1.1283895) q[0];
rz(1.7292132) q[2];
sx q[2];
rz(-2.7374501) q[2];
sx q[2];
rz(2.2805285) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.058152288) q[1];
sx q[1];
rz(-1.6279164) q[1];
sx q[1];
rz(-1.2099464) q[1];
rz(-pi) q[2];
rz(-1.8800354) q[3];
sx q[3];
rz(-1.6282363) q[3];
sx q[3];
rz(2.8185237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6115173) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(2.7008936) q[2];
rz(-0.079744451) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(1.124148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1752862) q[0];
sx q[0];
rz(-3.0847302) q[0];
sx q[0];
rz(2.9735907) q[0];
rz(3.1213144) q[1];
sx q[1];
rz(-0.30826491) q[1];
sx q[1];
rz(-1.6049989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.080606) q[0];
sx q[0];
rz(-1.8053375) q[0];
sx q[0];
rz(-1.3479665) q[0];
rz(-pi) q[1];
rz(-1.5606784) q[2];
sx q[2];
rz(-1.5732906) q[2];
sx q[2];
rz(-3.1163355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5913493) q[1];
sx q[1];
rz(-3.1367932) q[1];
sx q[1];
rz(-1.8018434) q[1];
rz(-pi) q[2];
x q[2];
rz(0.042293799) q[3];
sx q[3];
rz(-1.5797857) q[3];
sx q[3];
rz(0.90153722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7961879) q[2];
sx q[2];
rz(-0.91190839) q[2];
sx q[2];
rz(-1.7339285) q[2];
rz(2.0965072) q[3];
sx q[3];
rz(-0.049523517) q[3];
sx q[3];
rz(-0.27200562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8318091) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(2.580544) q[0];
rz(-2.8647515) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(1.3078088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5878752) q[0];
sx q[0];
rz(-0.29969117) q[0];
sx q[0];
rz(1.4315579) q[0];
x q[1];
rz(1.5708013) q[2];
sx q[2];
rz(-1.5781856) q[2];
sx q[2];
rz(2.0781197) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0289291) q[1];
sx q[1];
rz(-1.6293007) q[1];
sx q[1];
rz(-2.1444291) q[1];
rz(0.4588608) q[3];
sx q[3];
rz(-2.1956177) q[3];
sx q[3];
rz(0.030908728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7967367) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(0.5893839) q[2];
rz(1.899259) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(1.7813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.1347443) q[0];
sx q[0];
rz(-2.6297748) q[0];
sx q[0];
rz(1.7929329) q[0];
rz(-3.1349365) q[1];
sx q[1];
rz(-1.8204047) q[1];
sx q[1];
rz(3.1080918) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071246192) q[0];
sx q[0];
rz(-1.7435562) q[0];
sx q[0];
rz(0.91342775) q[0];
rz(1.5752931) q[2];
sx q[2];
rz(-1.6904313) q[2];
sx q[2];
rz(-0.41882354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12043145) q[1];
sx q[1];
rz(-2.8725) q[1];
sx q[1];
rz(1.5325559) q[1];
rz(-1.4064404) q[3];
sx q[3];
rz(-1.6911611) q[3];
sx q[3];
rz(0.09935483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.032430705) q[2];
sx q[2];
rz(-3.1350632) q[2];
sx q[2];
rz(-0.29864857) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(-0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3700767) q[0];
sx q[0];
rz(-1.6271485) q[0];
sx q[0];
rz(2.694743) q[0];
rz(2.9554548) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(1.4131379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22767775) q[0];
sx q[0];
rz(-0.23328885) q[0];
sx q[0];
rz(-1.4248217) q[0];
rz(-pi) q[1];
rz(-1.1158285) q[2];
sx q[2];
rz(-2.6778497) q[2];
sx q[2];
rz(-2.595903) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91292101) q[1];
sx q[1];
rz(-3.0745709) q[1];
sx q[1];
rz(2.3093501) q[1];
x q[2];
rz(1.8522315) q[3];
sx q[3];
rz(-1.3556983) q[3];
sx q[3];
rz(0.18751442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83429217) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(0.49868047) q[2];
rz(-0.57791609) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(-0.59670603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96108288) q[0];
sx q[0];
rz(-2.0208277) q[0];
sx q[0];
rz(0.34641308) q[0];
rz(0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(2.3880889) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2997871) q[0];
sx q[0];
rz(-0.25887576) q[0];
sx q[0];
rz(2.4095834) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0664938) q[2];
sx q[2];
rz(-1.4595928) q[2];
sx q[2];
rz(2.5203343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1767039) q[1];
sx q[1];
rz(-0.93550013) q[1];
sx q[1];
rz(-0.52365644) q[1];
x q[2];
rz(-1.0585045) q[3];
sx q[3];
rz(-0.18658328) q[3];
sx q[3];
rz(3.1087524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.57009131) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(1.5241148) q[2];
rz(3.0213455) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(-2.6038468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26871249) q[0];
sx q[0];
rz(-0.92395067) q[0];
sx q[0];
rz(2.9203316) q[0];
rz(-1.4606754) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(-3.0642919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54228264) q[0];
sx q[0];
rz(-3.0998383) q[0];
sx q[0];
rz(-1.6058654) q[0];
rz(-pi) q[1];
rz(3.1368106) q[2];
sx q[2];
rz(-1.5796697) q[2];
sx q[2];
rz(-2.9298669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52409808) q[1];
sx q[1];
rz(-1.3957983) q[1];
sx q[1];
rz(1.4982759) q[1];
rz(-pi) q[2];
rz(-1.037498) q[3];
sx q[3];
rz(-0.24720705) q[3];
sx q[3];
rz(-0.63106288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7920502) q[2];
sx q[2];
rz(-3.1304066) q[2];
sx q[2];
rz(-0.95996094) q[2];
rz(2.8121484) q[3];
sx q[3];
rz(-0.0080778413) q[3];
sx q[3];
rz(-2.2913057) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195381) q[0];
sx q[0];
rz(-2.52849) q[0];
sx q[0];
rz(3.0323113) q[0];
rz(2.7682313) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(1.2304617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7007992) q[0];
sx q[0];
rz(-2.1820118) q[0];
sx q[0];
rz(-0.81328765) q[0];
rz(-2.3745499) q[2];
sx q[2];
rz(-0.27077507) q[2];
sx q[2];
rz(-0.8052288) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8459863) q[1];
sx q[1];
rz(-1.5822268) q[1];
sx q[1];
rz(-1.5039526) q[1];
rz(-pi) q[2];
x q[2];
rz(1.968574) q[3];
sx q[3];
rz(-0.69044411) q[3];
sx q[3];
rz(-0.08716128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5665148) q[2];
sx q[2];
rz(-1.9058303) q[2];
sx q[2];
rz(-1.3285948) q[2];
rz(-1.7447507) q[3];
sx q[3];
rz(-3.1378523) q[3];
sx q[3];
rz(-2.1197135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.063865572) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(0.57300895) q[0];
rz(-2.833448) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(1.0073957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.751811) q[0];
sx q[0];
rz(-3.0349019) q[0];
sx q[0];
rz(1.5185028) q[0];
rz(-pi) q[1];
rz(-1.8780872) q[2];
sx q[2];
rz(-0.14423926) q[2];
sx q[2];
rz(2.7968614) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7511661) q[1];
sx q[1];
rz(-1.6114283) q[1];
sx q[1];
rz(-1.4875814) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5406403) q[3];
sx q[3];
rz(-1.584225) q[3];
sx q[3];
rz(1.3659988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8241626) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(-0.38995788) q[2];
rz(3.0725078) q[3];
sx q[3];
rz(-0.0091113541) q[3];
sx q[3];
rz(0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020141715) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(-0.87156975) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(-1.647324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0594306) q[0];
sx q[0];
rz(-0.63157394) q[0];
sx q[0];
rz(-0.53379121) q[0];
rz(-pi) q[1];
rz(0.062537161) q[2];
sx q[2];
rz(-0.61574575) q[2];
sx q[2];
rz(-3.1341022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20032665) q[1];
sx q[1];
rz(-1.8706053) q[1];
sx q[1];
rz(-2.7977562) q[1];
rz(-pi) q[2];
rz(3.028454) q[3];
sx q[3];
rz(-1.6125049) q[3];
sx q[3];
rz(1.0678408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5747052) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-3.1097143) q[2];
rz(-0.77830642) q[3];
sx q[3];
rz(-0.0068155546) q[3];
sx q[3];
rz(0.2955029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7185709) q[0];
sx q[0];
rz(-1.6091249) q[0];
sx q[0];
rz(-1.3269497) q[0];
rz(-0.12693916) q[1];
sx q[1];
rz(-2.9025684) q[1];
sx q[1];
rz(-2.92166) q[1];
rz(3.1359966) q[2];
sx q[2];
rz(-1.7104618) q[2];
sx q[2];
rz(0.24431123) q[2];
rz(1.6416141) q[3];
sx q[3];
rz(-0.67588617) q[3];
sx q[3];
rz(-2.5839154) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
