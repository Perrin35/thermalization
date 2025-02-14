OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(4.058429) q[0];
sx q[0];
rz(9.8596758) q[0];
rz(1.5974367) q[1];
sx q[1];
rz(-2.6225852) q[1];
sx q[1];
rz(2.5595698) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25044826) q[0];
sx q[0];
rz(-1.8316557) q[0];
sx q[0];
rz(2.9820739) q[0];
rz(-0.46635038) q[2];
sx q[2];
rz(-0.98578054) q[2];
sx q[2];
rz(-2.0908751) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4705321) q[1];
sx q[1];
rz(-1.3249505) q[1];
sx q[1];
rz(-1.125294) q[1];
rz(-pi) q[2];
rz(2.1277027) q[3];
sx q[3];
rz(-2.8134973) q[3];
sx q[3];
rz(-1.8522287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5883098) q[2];
sx q[2];
rz(-1.2574235) q[2];
sx q[2];
rz(-2.8498939) q[2];
rz(-0.45082539) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(0.60744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.70158231) q[0];
sx q[0];
rz(-2.1511183) q[0];
sx q[0];
rz(2.0462346) q[0];
rz(-1.9502684) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(-2.2390168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0046834) q[0];
sx q[0];
rz(-2.1666514) q[0];
sx q[0];
rz(2.3803902) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.90229374) q[2];
sx q[2];
rz(-1.8040787) q[2];
sx q[2];
rz(-1.4177314) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2486156) q[1];
sx q[1];
rz(-1.6910403) q[1];
sx q[1];
rz(2.8477232) q[1];
x q[2];
rz(-0.5036854) q[3];
sx q[3];
rz(-1.0102444) q[3];
sx q[3];
rz(-0.47893804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9537182) q[2];
sx q[2];
rz(-0.97430054) q[2];
sx q[2];
rz(-0.11889674) q[2];
rz(2.4925354) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(-1.4547179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.651799) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(2.6258262) q[0];
rz(2.0491397) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(1.2752424) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37947734) q[0];
sx q[0];
rz(-1.5889052) q[0];
sx q[0];
rz(-3.1000231) q[0];
x q[1];
rz(0.31351201) q[2];
sx q[2];
rz(-1.3517153) q[2];
sx q[2];
rz(-0.30922019) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4750649) q[1];
sx q[1];
rz(-1.3747066) q[1];
sx q[1];
rz(-0.88351078) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6728064) q[3];
sx q[3];
rz(-2.3622741) q[3];
sx q[3];
rz(-1.23097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0755997) q[2];
sx q[2];
rz(-1.8241355) q[2];
sx q[2];
rz(-1.9817748) q[2];
rz(-2.6521111) q[3];
sx q[3];
rz(-1.5533841) q[3];
sx q[3];
rz(0.71973962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2751665) q[0];
sx q[0];
rz(-2.8186099) q[0];
sx q[0];
rz(-2.6599356) q[0];
rz(-0.9306759) q[1];
sx q[1];
rz(-1.8447256) q[1];
sx q[1];
rz(-0.01893386) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1279432) q[0];
sx q[0];
rz(-3.0873723) q[0];
sx q[0];
rz(1.9677866) q[0];
x q[1];
rz(-1.6274979) q[2];
sx q[2];
rz(-2.4783274) q[2];
sx q[2];
rz(-0.5723638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9908617) q[1];
sx q[1];
rz(-2.1580937) q[1];
sx q[1];
rz(2.2022922) q[1];
rz(-pi) q[2];
rz(-0.83796699) q[3];
sx q[3];
rz(-2.1893756) q[3];
sx q[3];
rz(-0.092242084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9229752) q[2];
sx q[2];
rz(-2.7851892) q[2];
sx q[2];
rz(-2.3962928) q[2];
rz(-2.7636012) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(-2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9208263) q[0];
sx q[0];
rz(-2.8296318) q[0];
sx q[0];
rz(1.1908603) q[0];
rz(-2.6530755) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(2.3337505) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90276001) q[0];
sx q[0];
rz(-1.897398) q[0];
sx q[0];
rz(-3.045911) q[0];
rz(-pi) q[1];
rz(1.135181) q[2];
sx q[2];
rz(-2.6382338) q[2];
sx q[2];
rz(1.8653009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1792308) q[1];
sx q[1];
rz(-2.6768502) q[1];
sx q[1];
rz(2.9443113) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1937815) q[3];
sx q[3];
rz(-1.4075507) q[3];
sx q[3];
rz(-2.3415274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1390344) q[2];
sx q[2];
rz(-0.99271861) q[2];
sx q[2];
rz(0.16415088) q[2];
rz(-2.3955591) q[3];
sx q[3];
rz(-1.371871) q[3];
sx q[3];
rz(-3.0677838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98110759) q[0];
sx q[0];
rz(-1.2657413) q[0];
sx q[0];
rz(-2.4142081) q[0];
rz(-0.47850594) q[1];
sx q[1];
rz(-1.452927) q[1];
sx q[1];
rz(-2.4047638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3872358) q[0];
sx q[0];
rz(-0.14850907) q[0];
sx q[0];
rz(1.1867619) q[0];
rz(1.1428035) q[2];
sx q[2];
rz(-0.30747947) q[2];
sx q[2];
rz(2.3124419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.068253156) q[1];
sx q[1];
rz(-2.1756449) q[1];
sx q[1];
rz(2.6920435) q[1];
rz(-0.013929587) q[3];
sx q[3];
rz(-1.1831814) q[3];
sx q[3];
rz(0.25957169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9937146) q[2];
sx q[2];
rz(-1.0241877) q[2];
sx q[2];
rz(0.68515879) q[2];
rz(2.1327175) q[3];
sx q[3];
rz(-0.69005552) q[3];
sx q[3];
rz(-2.8907997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93609) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(-1.093338) q[0];
rz(-3.1178442) q[1];
sx q[1];
rz(-2.4045585) q[1];
sx q[1];
rz(-0.49560961) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28195295) q[0];
sx q[0];
rz(-3.0744327) q[0];
sx q[0];
rz(-0.18224506) q[0];
rz(-1.7772113) q[2];
sx q[2];
rz(-1.141618) q[2];
sx q[2];
rz(0.51233722) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.520185) q[1];
sx q[1];
rz(-1.9193135) q[1];
sx q[1];
rz(1.7640616) q[1];
x q[2];
rz(-0.43066671) q[3];
sx q[3];
rz(-1.1209295) q[3];
sx q[3];
rz(-1.7783742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89639837) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(0.29655656) q[2];
rz(-0.5101997) q[3];
sx q[3];
rz(-1.5254131) q[3];
sx q[3];
rz(2.8701674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(1.9940014) q[0];
sx q[0];
rz(-0.95518249) q[0];
sx q[0];
rz(1.9872794) q[0];
rz(-1.1555903) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(1.6824228) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66508913) q[0];
sx q[0];
rz(-2.4263093) q[0];
sx q[0];
rz(-0.034937783) q[0];
rz(-pi) q[1];
rz(0.84648561) q[2];
sx q[2];
rz(-1.5716388) q[2];
sx q[2];
rz(-0.049132012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0399122) q[1];
sx q[1];
rz(-1.7500568) q[1];
sx q[1];
rz(-1.0615361) q[1];
x q[2];
rz(-0.63520875) q[3];
sx q[3];
rz(-2.016474) q[3];
sx q[3];
rz(3.095568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3079754) q[2];
sx q[2];
rz(-0.54215446) q[2];
sx q[2];
rz(1.3758434) q[2];
rz(3.0552676) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(-1.255792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1983222) q[0];
sx q[0];
rz(-0.69381303) q[0];
sx q[0];
rz(2.2721403) q[0];
rz(-2.4122639) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(2.9076911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9651523) q[0];
sx q[0];
rz(-1.6185221) q[0];
sx q[0];
rz(-0.41712572) q[0];
x q[1];
rz(2.7178784) q[2];
sx q[2];
rz(-2.9117081) q[2];
sx q[2];
rz(-2.4073441) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4471376) q[1];
sx q[1];
rz(-1.7181953) q[1];
sx q[1];
rz(0.44575341) q[1];
rz(-pi) q[2];
rz(-2.8222705) q[3];
sx q[3];
rz(-0.87933137) q[3];
sx q[3];
rz(-2.6806074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0539315) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(-1.1762478) q[2];
rz(2.4704399) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(-1.6254856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072356) q[0];
sx q[0];
rz(-2.2800627) q[0];
sx q[0];
rz(-1.0435411) q[0];
rz(0.097298233) q[1];
sx q[1];
rz(-1.0568591) q[1];
sx q[1];
rz(1.1204488) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0138884) q[0];
sx q[0];
rz(-0.59895016) q[0];
sx q[0];
rz(1.4636135) q[0];
rz(-pi) q[1];
rz(-2.3612622) q[2];
sx q[2];
rz(-0.95338168) q[2];
sx q[2];
rz(1.8898026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0617025) q[1];
sx q[1];
rz(-1.7482867) q[1];
sx q[1];
rz(-1.5212785) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2213732) q[3];
sx q[3];
rz(-2.6333678) q[3];
sx q[3];
rz(-2.1247381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4572767) q[2];
sx q[2];
rz(-0.61351675) q[2];
sx q[2];
rz(-1.786001) q[2];
rz(-1.5654303) q[3];
sx q[3];
rz(-1.0836982) q[3];
sx q[3];
rz(2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23715699) q[0];
sx q[0];
rz(-2.5813527) q[0];
sx q[0];
rz(-0.78975633) q[0];
rz(2.4850028) q[1];
sx q[1];
rz(-2.0273392) q[1];
sx q[1];
rz(-0.56710342) q[1];
rz(-2.1661027) q[2];
sx q[2];
rz(-0.75303034) q[2];
sx q[2];
rz(-1.2122214) q[2];
rz(2.8511467) q[3];
sx q[3];
rz(-1.757156) q[3];
sx q[3];
rz(0.7745756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
