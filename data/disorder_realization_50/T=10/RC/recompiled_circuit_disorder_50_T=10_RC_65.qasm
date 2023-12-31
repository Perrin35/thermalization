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
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(2.3666518) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89844184) q[0];
sx q[0];
rz(-1.713322) q[0];
sx q[0];
rz(1.5205163) q[0];
rz(2.7156419) q[2];
sx q[2];
rz(-2.2507239) q[2];
sx q[2];
rz(-1.8217063) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2831813) q[1];
sx q[1];
rz(-1.5309257) q[1];
sx q[1];
rz(-1.7130501) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6585456) q[3];
sx q[3];
rz(-2.8222198) q[3];
sx q[3];
rz(-1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(-1.3886064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4202704) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.7659448) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20695105) q[0];
sx q[0];
rz(-1.9793545) q[0];
sx q[0];
rz(1.6658528) q[0];
rz(-1.3865115) q[2];
sx q[2];
rz(-2.8993336) q[2];
sx q[2];
rz(0.87528961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8010506) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(-1.3786475) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7820285) q[3];
sx q[3];
rz(-2.095795) q[3];
sx q[3];
rz(1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-2.6518872) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-2.074923) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-0.24060732) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-2.1668285) q[1];
sx q[1];
rz(-1.2352357) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19757195) q[0];
sx q[0];
rz(-0.70957843) q[0];
sx q[0];
rz(-2.5463085) q[0];
x q[1];
rz(0.58358242) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(-0.32354087) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2521378) q[1];
sx q[1];
rz(-0.368202) q[1];
sx q[1];
rz(-0.77106573) q[1];
x q[2];
rz(0.62883212) q[3];
sx q[3];
rz(-1.1662081) q[3];
sx q[3];
rz(-1.7189327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(-1.7049449) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.876372) q[3];
sx q[3];
rz(0.60825545) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.065598) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(3.0217357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.170341) q[0];
sx q[0];
rz(-0.60702885) q[0];
sx q[0];
rz(-1.6334565) q[0];
rz(2.7808431) q[2];
sx q[2];
rz(-1.7067688) q[2];
sx q[2];
rz(-3.0038358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67112982) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(-1.8030241) q[1];
rz(-pi) q[2];
rz(-0.5991163) q[3];
sx q[3];
rz(-2.0261507) q[3];
sx q[3];
rz(1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5084761) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(3.0692549) q[2];
rz(-0.37483254) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(-1.9434631) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(-0.72987366) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2729028) q[0];
sx q[0];
rz(-1.5575952) q[0];
sx q[0];
rz(2.9024283) q[0];
x q[1];
rz(-1.8334332) q[2];
sx q[2];
rz(-1.3381759) q[2];
sx q[2];
rz(2.5831985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6606635) q[1];
sx q[1];
rz(-0.96748057) q[1];
sx q[1];
rz(1.4900581) q[1];
rz(-pi) q[2];
rz(2.6579882) q[3];
sx q[3];
rz(-2.1198453) q[3];
sx q[3];
rz(-2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.038625) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(0.70303482) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(0.15771244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
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
rz(2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.7646164) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0159857) q[2];
sx q[2];
rz(-1.5714958) q[2];
sx q[2];
rz(-0.47847754) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78696886) q[1];
sx q[1];
rz(-2.0957392) q[1];
sx q[1];
rz(-0.24137361) q[1];
x q[2];
rz(0.89208608) q[3];
sx q[3];
rz(-1.6086372) q[3];
sx q[3];
rz(2.2298262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.4921654) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.7892388) q[0];
sx q[0];
rz(1.6725756) q[0];
rz(-1.0143657) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(1.8168824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46885195) q[0];
sx q[0];
rz(-1.6981089) q[0];
sx q[0];
rz(-1.5936113) q[0];
rz(-0.5092233) q[2];
sx q[2];
rz(-1.0668796) q[2];
sx q[2];
rz(0.06171209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99609375) q[1];
sx q[1];
rz(-2.5719574) q[1];
sx q[1];
rz(0.74665229) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2277101) q[3];
sx q[3];
rz(-1.8808639) q[3];
sx q[3];
rz(-0.1012181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(1.8519648) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8136918) q[0];
sx q[0];
rz(-2.0223589) q[0];
sx q[0];
rz(-0.62664647) q[0];
rz(-pi) q[1];
rz(-0.81823924) q[2];
sx q[2];
rz(-1.0050251) q[2];
sx q[2];
rz(-2.408839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3813078) q[1];
sx q[1];
rz(-2.9236301) q[1];
sx q[1];
rz(1.0337619) q[1];
x q[2];
rz(2.9309978) q[3];
sx q[3];
rz(-0.67376332) q[3];
sx q[3];
rz(1.8316366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(1.0827433) q[2];
rz(3.0454214) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.2307897) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(2.8073231) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-2.8589378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219014) q[0];
sx q[0];
rz(-2.0265409) q[0];
sx q[0];
rz(0.68463188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5731509) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(-2.8708411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5987451) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(2.7302242) q[1];
rz(-0.11741365) q[3];
sx q[3];
rz(-0.88866361) q[3];
sx q[3];
rz(-0.62473245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.1634469) q[2];
sx q[2];
rz(-2.4003417) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
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
x q[3];
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
rz(1.0989477) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-0.50869554) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(-0.68181109) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9671229) q[0];
sx q[0];
rz(-0.57708626) q[0];
sx q[0];
rz(1.2601) q[0];
rz(-pi) q[1];
x q[1];
rz(2.290906) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(-2.1249287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0971113) q[1];
sx q[1];
rz(-2.0437355) q[1];
sx q[1];
rz(0.60826917) q[1];
rz(-pi) q[2];
rz(-2.6184611) q[3];
sx q[3];
rz(-1.0732068) q[3];
sx q[3];
rz(0.79062068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.5073551) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9733799) q[0];
sx q[0];
rz(-1.6642878) q[0];
sx q[0];
rz(-0.99933495) q[0];
rz(2.534261) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(2.091999) q[2];
sx q[2];
rz(-2.1672719) q[2];
sx q[2];
rz(1.5530829) q[2];
rz(2.7568983) q[3];
sx q[3];
rz(-1.7582498) q[3];
sx q[3];
rz(2.869217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
