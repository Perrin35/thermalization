OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(-1.2195725) q[0];
rz(3.976877) q[1];
sx q[1];
rz(4.4983954) q[1];
sx q[1];
rz(8.5228336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.565958) q[0];
sx q[0];
rz(-2.1778641) q[0];
sx q[0];
rz(1.7580887) q[0];
rz(2.7612258) q[2];
sx q[2];
rz(-1.7445626) q[2];
sx q[2];
rz(1.1443537) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28022596) q[1];
sx q[1];
rz(-1.8609443) q[1];
sx q[1];
rz(-2.9600083) q[1];
x q[2];
rz(-2.4331081) q[3];
sx q[3];
rz(-0.86743858) q[3];
sx q[3];
rz(-3.0391703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1787662) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(0.82628769) q[2];
rz(-1.7360342) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(-0.74339408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5590782) q[0];
sx q[0];
rz(-0.40575108) q[0];
sx q[0];
rz(-1.4003117) q[0];
rz(-1.1236313) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(-2.0434911) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208291) q[0];
sx q[0];
rz(-1.4614824) q[0];
sx q[0];
rz(-2.9846322) q[0];
x q[1];
rz(1.2400572) q[2];
sx q[2];
rz(-1.2501226) q[2];
sx q[2];
rz(-1.5000686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.73126331) q[1];
sx q[1];
rz(-1.6917233) q[1];
sx q[1];
rz(1.6464473) q[1];
x q[2];
rz(-0.16657515) q[3];
sx q[3];
rz(-1.9431356) q[3];
sx q[3];
rz(-2.6083715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0565679) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(-1.8918096) q[2];
rz(1.7791087) q[3];
sx q[3];
rz(-2.1407849) q[3];
sx q[3];
rz(1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63610858) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(-2.8662477) q[0];
rz(2.6823726) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(-0.053344639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5551852) q[0];
sx q[0];
rz(-1.6490899) q[0];
sx q[0];
rz(-0.21348641) q[0];
rz(-pi) q[1];
rz(1.6079536) q[2];
sx q[2];
rz(-2.6803663) q[2];
sx q[2];
rz(-0.30468582) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.054823067) q[1];
sx q[1];
rz(-0.55332282) q[1];
sx q[1];
rz(0.077484681) q[1];
rz(-pi) q[2];
rz(-1.2417913) q[3];
sx q[3];
rz(-1.3734372) q[3];
sx q[3];
rz(0.98361919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.82103819) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(1.5041171) q[2];
rz(1.641168) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(0.27568451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212873) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(2.6570008) q[0];
rz(1.2424319) q[1];
sx q[1];
rz(-0.74598765) q[1];
sx q[1];
rz(2.7925083) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.718133) q[0];
sx q[0];
rz(-1.0225931) q[0];
sx q[0];
rz(-0.050027442) q[0];
rz(-pi) q[1];
rz(0.075320638) q[2];
sx q[2];
rz(-1.2284797) q[2];
sx q[2];
rz(1.3665733) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1368833) q[1];
sx q[1];
rz(-1.7994969) q[1];
sx q[1];
rz(-1.7654224) q[1];
rz(2.146017) q[3];
sx q[3];
rz(-0.68505008) q[3];
sx q[3];
rz(0.12888651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35859534) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-2.6915468) q[2];
rz(-0.83886823) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(0.1990327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435476) q[0];
sx q[0];
rz(-1.6723375) q[0];
sx q[0];
rz(-0.099844649) q[0];
rz(1.3205344) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(0.60755306) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71220416) q[0];
sx q[0];
rz(-0.11183248) q[0];
sx q[0];
rz(-0.73356103) q[0];
x q[1];
rz(-1.728292) q[2];
sx q[2];
rz(-1.0309426) q[2];
sx q[2];
rz(1.4992204) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1888426) q[1];
sx q[1];
rz(-2.5217068) q[1];
sx q[1];
rz(2.7623738) q[1];
rz(-pi) q[2];
rz(-1.8292866) q[3];
sx q[3];
rz(-1.6655169) q[3];
sx q[3];
rz(-2.3324089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1943835) q[2];
sx q[2];
rz(-0.76346976) q[2];
sx q[2];
rz(-0.53059951) q[2];
rz(-2.7001906) q[3];
sx q[3];
rz(-0.97359052) q[3];
sx q[3];
rz(-2.1921564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7049103) q[0];
sx q[0];
rz(-1.3446151) q[0];
sx q[0];
rz(2.5675024) q[0];
rz(-2.2615945) q[1];
sx q[1];
rz(-1.1209542) q[1];
sx q[1];
rz(1.7990187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8594151) q[0];
sx q[0];
rz(-1.1033774) q[0];
sx q[0];
rz(-2.9168081) q[0];
rz(-pi) q[1];
rz(-1.1583352) q[2];
sx q[2];
rz(-1.165373) q[2];
sx q[2];
rz(2.6087922) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33730971) q[1];
sx q[1];
rz(-1.8071257) q[1];
sx q[1];
rz(-1.9419036) q[1];
rz(-pi) q[2];
x q[2];
rz(1.430949) q[3];
sx q[3];
rz(-2.2290843) q[3];
sx q[3];
rz(2.0161395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-0.28417045) q[2];
sx q[2];
rz(-2.6944842) q[2];
rz(1.0304662) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(-0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717188) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(0.41279992) q[0];
rz(2.0665456) q[1];
sx q[1];
rz(-1.1378891) q[1];
sx q[1];
rz(-0.80361754) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5980854) q[0];
sx q[0];
rz(-1.0585365) q[0];
sx q[0];
rz(0.58337636) q[0];
rz(-pi) q[1];
rz(0.86974025) q[2];
sx q[2];
rz(-0.77574965) q[2];
sx q[2];
rz(2.9091331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4839125) q[1];
sx q[1];
rz(-1.5849304) q[1];
sx q[1];
rz(-2.3573228) q[1];
rz(1.647023) q[3];
sx q[3];
rz(-1.7496568) q[3];
sx q[3];
rz(-0.12693996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.86218086) q[2];
sx q[2];
rz(-1.923424) q[2];
sx q[2];
rz(1.5009521) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-2.3707135) q[3];
sx q[3];
rz(-1.4890495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.415446) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(0.55794445) q[0];
rz(-0.56124148) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(1.4161313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6218044) q[0];
sx q[0];
rz(-1.1237696) q[0];
sx q[0];
rz(-2.8267953) q[0];
x q[1];
rz(0.92500706) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(-0.031161873) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0603283) q[1];
sx q[1];
rz(-0.34354478) q[1];
sx q[1];
rz(-1.8650123) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.775557) q[3];
sx q[3];
rz(-2.3157218) q[3];
sx q[3];
rz(-1.5958169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8884362) q[2];
sx q[2];
rz(-0.77740589) q[2];
sx q[2];
rz(1.0003132) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(0.15672556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67824739) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(3.1372702) q[0];
rz(-0.16383544) q[1];
sx q[1];
rz(-0.42525735) q[1];
sx q[1];
rz(-1.3660627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28986606) q[0];
sx q[0];
rz(-0.64771876) q[0];
sx q[0];
rz(-1.6169548) q[0];
rz(-pi) q[1];
rz(0.37213426) q[2];
sx q[2];
rz(-1.8610125) q[2];
sx q[2];
rz(-1.8572825) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1613579) q[1];
sx q[1];
rz(-1.606985) q[1];
sx q[1];
rz(0.25118942) q[1];
x q[2];
rz(-2.3187997) q[3];
sx q[3];
rz(-1.0446915) q[3];
sx q[3];
rz(1.9578809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0020478) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(-0.61919332) q[2];
rz(2.5891417) q[3];
sx q[3];
rz(-1.8360454) q[3];
sx q[3];
rz(2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3510975) q[0];
sx q[0];
rz(-0.18432291) q[0];
sx q[0];
rz(0.47875324) q[0];
rz(1.152773) q[1];
sx q[1];
rz(-0.29007998) q[1];
sx q[1];
rz(2.1943888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1789492) q[0];
sx q[0];
rz(-1.5367723) q[0];
sx q[0];
rz(1.7555139) q[0];
x q[1];
rz(2.8962063) q[2];
sx q[2];
rz(-0.99694251) q[2];
sx q[2];
rz(1.4978486) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.041872) q[1];
sx q[1];
rz(-2.2769089) q[1];
sx q[1];
rz(0.40634917) q[1];
rz(1.8197038) q[3];
sx q[3];
rz(-1.2281903) q[3];
sx q[3];
rz(-0.026503868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(-2.0742778) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08854475) q[0];
sx q[0];
rz(-1.6041258) q[0];
sx q[0];
rz(-1.5969101) q[0];
rz(-1.0531986) q[1];
sx q[1];
rz(-1.2146626) q[1];
sx q[1];
rz(0.76520898) q[1];
rz(-2.8996053) q[2];
sx q[2];
rz(-1.568071) q[2];
sx q[2];
rz(-0.65721401) q[2];
rz(-2.5625145) q[3];
sx q[3];
rz(-1.758772) q[3];
sx q[3];
rz(0.72455345) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
