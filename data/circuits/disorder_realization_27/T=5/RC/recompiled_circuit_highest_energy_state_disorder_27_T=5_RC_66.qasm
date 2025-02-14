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
rz(-1.902154) q[0];
sx q[0];
rz(-1.3286989) q[0];
sx q[0];
rz(0.21188307) q[0];
rz(1.7243241) q[1];
sx q[1];
rz(-0.53242004) q[1];
sx q[1];
rz(2.7636757) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53055313) q[0];
sx q[0];
rz(-1.5612649) q[0];
sx q[0];
rz(3.1411489) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7013999) q[2];
sx q[2];
rz(-0.72840103) q[2];
sx q[2];
rz(-3.13934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18641414) q[1];
sx q[1];
rz(-0.70611533) q[1];
sx q[1];
rz(0.43118711) q[1];
x q[2];
rz(-2.3232949) q[3];
sx q[3];
rz(-0.94486559) q[3];
sx q[3];
rz(1.4656386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8523031) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-0.079785384) q[2];
rz(-0.96528178) q[3];
sx q[3];
rz(-1.691317) q[3];
sx q[3];
rz(2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(1.7695919) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(1.1044097) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8744132) q[0];
sx q[0];
rz(-2.4738418) q[0];
sx q[0];
rz(-0.37607583) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49141617) q[2];
sx q[2];
rz(-2.2960536) q[2];
sx q[2];
rz(0.27057901) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1762878) q[1];
sx q[1];
rz(-1.4215648) q[1];
sx q[1];
rz(-1.1231827) q[1];
x q[2];
rz(-1.4711963) q[3];
sx q[3];
rz(-1.0336813) q[3];
sx q[3];
rz(0.84505862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.264512) q[2];
sx q[2];
rz(-1.1254213) q[2];
sx q[2];
rz(1.4667) q[2];
rz(3.1125715) q[3];
sx q[3];
rz(-1.0163739) q[3];
sx q[3];
rz(-2.4513054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324209) q[0];
sx q[0];
rz(-3.0747774) q[0];
sx q[0];
rz(0.27798852) q[0];
rz(1.4311283) q[1];
sx q[1];
rz(-2.2093096) q[1];
sx q[1];
rz(-2.7412282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5516282) q[0];
sx q[0];
rz(-1.1271521) q[0];
sx q[0];
rz(-0.69728627) q[0];
rz(0.30817356) q[2];
sx q[2];
rz(-0.71012596) q[2];
sx q[2];
rz(-1.9911204) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6975721) q[1];
sx q[1];
rz(-2.6409147) q[1];
sx q[1];
rz(-0.33659192) q[1];
rz(-0.72680803) q[3];
sx q[3];
rz(-2.6977959) q[3];
sx q[3];
rz(-0.82647317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6843188) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(2.686783) q[2];
rz(1.8999892) q[3];
sx q[3];
rz(-0.95507115) q[3];
sx q[3];
rz(2.4212867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95551816) q[0];
sx q[0];
rz(-1.3294514) q[0];
sx q[0];
rz(3.1410826) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-2.3020703) q[1];
sx q[1];
rz(2.9972163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3568319) q[0];
sx q[0];
rz(-2.5481173) q[0];
sx q[0];
rz(1.1952101) q[0];
rz(1.7208485) q[2];
sx q[2];
rz(-1.320082) q[2];
sx q[2];
rz(1.5992407) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3444654) q[1];
sx q[1];
rz(-0.68736156) q[1];
sx q[1];
rz(0.3631773) q[1];
x q[2];
rz(1.4636643) q[3];
sx q[3];
rz(-2.1153573) q[3];
sx q[3];
rz(-0.0060826172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1986177) q[2];
sx q[2];
rz(-1.696442) q[2];
sx q[2];
rz(0.96251881) q[2];
rz(-1.5396384) q[3];
sx q[3];
rz(-1.4055777) q[3];
sx q[3];
rz(-2.8008154) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9050423) q[0];
sx q[0];
rz(-1.6860697) q[0];
sx q[0];
rz(-1.9566253) q[0];
rz(-2.9153337) q[1];
sx q[1];
rz(-2.2626651) q[1];
sx q[1];
rz(2.102898) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0112891) q[0];
sx q[0];
rz(-1.4892007) q[0];
sx q[0];
rz(2.5528583) q[0];
x q[1];
rz(-0.67579999) q[2];
sx q[2];
rz(-0.86188176) q[2];
sx q[2];
rz(1.0647286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8983137) q[1];
sx q[1];
rz(-2.2479575) q[1];
sx q[1];
rz(-2.2449298) q[1];
x q[2];
rz(-1.709278) q[3];
sx q[3];
rz(-0.44381986) q[3];
sx q[3];
rz(-0.5139155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(-1.4656674) q[3];
sx q[3];
rz(-1.1301872) q[3];
sx q[3];
rz(1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6370711) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(-1.1035408) q[0];
rz(2.6612813) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(2.5441817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99973122) q[0];
sx q[0];
rz(-1.7787361) q[0];
sx q[0];
rz(-0.19241649) q[0];
rz(1.3296606) q[2];
sx q[2];
rz(-0.84000194) q[2];
sx q[2];
rz(-2.3868449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.061042) q[1];
sx q[1];
rz(-1.4054728) q[1];
sx q[1];
rz(0.25456659) q[1];
rz(1.1998981) q[3];
sx q[3];
rz(-1.6565595) q[3];
sx q[3];
rz(-0.24468064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9395113) q[2];
sx q[2];
rz(-1.6263522) q[2];
sx q[2];
rz(1.0763947) q[2];
rz(1.1427897) q[3];
sx q[3];
rz(-0.79311526) q[3];
sx q[3];
rz(-0.49016652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6550605) q[0];
sx q[0];
rz(-1.9831816) q[0];
sx q[0];
rz(0.038473815) q[0];
rz(-3.0768652) q[1];
sx q[1];
rz(-1.3836626) q[1];
sx q[1];
rz(-0.23385349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.192894) q[0];
sx q[0];
rz(-1.7561551) q[0];
sx q[0];
rz(0.22484397) q[0];
x q[1];
rz(-2.269417) q[2];
sx q[2];
rz(-0.46469122) q[2];
sx q[2];
rz(1.8605491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0044627) q[1];
sx q[1];
rz(-2.4871768) q[1];
sx q[1];
rz(-0.14738247) q[1];
rz(-pi) q[2];
rz(-1.5130299) q[3];
sx q[3];
rz(-1.7004629) q[3];
sx q[3];
rz(-0.71021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48876277) q[2];
sx q[2];
rz(-1.5831999) q[2];
sx q[2];
rz(2.5433507) q[2];
rz(3.0277142) q[3];
sx q[3];
rz(-1.7555534) q[3];
sx q[3];
rz(2.2940476) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050196) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(0.2463499) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-0.49682239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2082191) q[0];
sx q[0];
rz(-0.95666203) q[0];
sx q[0];
rz(-2.0858411) q[0];
x q[1];
rz(0.70513983) q[2];
sx q[2];
rz(-1.459483) q[2];
sx q[2];
rz(1.706858) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0223479) q[1];
sx q[1];
rz(-1.3356707) q[1];
sx q[1];
rz(2.0658595) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9867191) q[3];
sx q[3];
rz(-2.057977) q[3];
sx q[3];
rz(-0.39184141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7184427) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(-1.144484) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(-2.9437039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941187) q[0];
sx q[0];
rz(-2.9070774) q[0];
sx q[0];
rz(-0.06614729) q[0];
rz(1.9937531) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(-0.60417169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3124638) q[0];
sx q[0];
rz(-1.7760217) q[0];
sx q[0];
rz(0.28136307) q[0];
rz(1.4140903) q[2];
sx q[2];
rz(-2.092166) q[2];
sx q[2];
rz(-2.7336313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0688404) q[1];
sx q[1];
rz(-1.3786331) q[1];
sx q[1];
rz(0.56984624) q[1];
rz(-pi) q[2];
rz(-0.87852134) q[3];
sx q[3];
rz(-2.1067348) q[3];
sx q[3];
rz(-1.5978447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.26312795) q[2];
sx q[2];
rz(-2.2915514) q[2];
sx q[2];
rz(2.7086332) q[2];
rz(-1.2290907) q[3];
sx q[3];
rz(-1.2058328) q[3];
sx q[3];
rz(1.3273201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.1935254) q[0];
sx q[0];
rz(-2.0663517) q[0];
sx q[0];
rz(-1.2731592) q[0];
rz(-0.46514568) q[1];
sx q[1];
rz(-1.7756614) q[1];
sx q[1];
rz(-2.926631) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5453542) q[0];
sx q[0];
rz(-0.69850498) q[0];
sx q[0];
rz(0.39909382) q[0];
rz(-pi) q[1];
rz(0.6647756) q[2];
sx q[2];
rz(-0.6419581) q[2];
sx q[2];
rz(0.73981111) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4196288) q[1];
sx q[1];
rz(-2.2441494) q[1];
sx q[1];
rz(-0.90760214) q[1];
rz(-pi) q[2];
rz(2.5796579) q[3];
sx q[3];
rz(-1.7896381) q[3];
sx q[3];
rz(0.17499017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2994069) q[2];
sx q[2];
rz(-2.3278548) q[2];
sx q[2];
rz(-2.1827533) q[2];
rz(-1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(1.6963262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5398298) q[0];
sx q[0];
rz(-1.6015263) q[0];
sx q[0];
rz(1.482561) q[0];
rz(-2.4330347) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(0.56253994) q[2];
sx q[2];
rz(-1.8481135) q[2];
sx q[2];
rz(-2.6279966) q[2];
rz(-0.72193969) q[3];
sx q[3];
rz(-1.7122713) q[3];
sx q[3];
rz(0.28817978) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
