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
rz(2.4636318) q[0];
sx q[0];
rz(-2.8673661) q[0];
sx q[0];
rz(-1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3363269) q[0];
sx q[0];
rz(-1.299843) q[0];
sx q[0];
rz(1.8150369) q[0];
rz(2.847509) q[2];
sx q[2];
rz(-0.80162797) q[2];
sx q[2];
rz(1.5844567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5668402) q[1];
sx q[1];
rz(-2.2623203) q[1];
sx q[1];
rz(-1.994446) q[1];
rz(-0.94120195) q[3];
sx q[3];
rz(-1.1735907) q[3];
sx q[3];
rz(-1.8091699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(2.2402666) q[2];
rz(2.6775635) q[3];
sx q[3];
rz(-1.2814859) q[3];
sx q[3];
rz(0.06981167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0483911) q[0];
sx q[0];
rz(-3.1049325) q[0];
sx q[0];
rz(-1.2777591) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-2.6779046) q[1];
sx q[1];
rz(-1.5346079) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8746895) q[0];
sx q[0];
rz(-2.1703566) q[0];
sx q[0];
rz(3.1019109) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8707451) q[2];
sx q[2];
rz(-0.19057628) q[2];
sx q[2];
rz(1.1491164) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6834641) q[1];
sx q[1];
rz(-2.5741815) q[1];
sx q[1];
rz(0.50000425) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7369698) q[3];
sx q[3];
rz(-1.1841067) q[3];
sx q[3];
rz(-0.85961941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.719912) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(2.9812532) q[2];
rz(-2.4028589) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(1.7140478) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62683231) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(-0.50977388) q[0];
rz(-1.7104507) q[1];
sx q[1];
rz(-1.1809843) q[1];
sx q[1];
rz(0.74657718) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7187481) q[0];
sx q[0];
rz(-1.5507409) q[0];
sx q[0];
rz(2.123318) q[0];
rz(-pi) q[1];
rz(-2.4116729) q[2];
sx q[2];
rz(-2.0509208) q[2];
sx q[2];
rz(1.7034704) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4164734) q[1];
sx q[1];
rz(-1.2095569) q[1];
sx q[1];
rz(0.84049417) q[1];
rz(0.88991092) q[3];
sx q[3];
rz(-1.2632252) q[3];
sx q[3];
rz(0.5336844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9868682) q[2];
sx q[2];
rz(-2.8747323) q[2];
sx q[2];
rz(-1.223684) q[2];
rz(-2.0679421) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(-2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2735485) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(-2.7622188) q[0];
rz(-2.9647478) q[1];
sx q[1];
rz(-1.481448) q[1];
sx q[1];
rz(0.79536974) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.11957) q[0];
sx q[0];
rz(-1.9384465) q[0];
sx q[0];
rz(0.099438503) q[0];
rz(-2.3248939) q[2];
sx q[2];
rz(-1.3743322) q[2];
sx q[2];
rz(2.8981371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6713555) q[1];
sx q[1];
rz(-2.0908326) q[1];
sx q[1];
rz(-2.8440471) q[1];
x q[2];
rz(-1.627911) q[3];
sx q[3];
rz(-1.4453672) q[3];
sx q[3];
rz(-2.7386087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.4312326) q[2];
sx q[2];
rz(-1.3453588) q[2];
sx q[2];
rz(1.7809407) q[2];
rz(-1.0294754) q[3];
sx q[3];
rz(-1.4808713) q[3];
sx q[3];
rz(1.8195389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4904356) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(-1.5486451) q[0];
rz(-0.40052888) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.4097479) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1574788) q[0];
sx q[0];
rz(-1.7703919) q[0];
sx q[0];
rz(2.4019813) q[0];
rz(-pi) q[1];
rz(0.35946741) q[2];
sx q[2];
rz(-1.4233982) q[2];
sx q[2];
rz(-2.0161674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6741989) q[1];
sx q[1];
rz(-2.0122897) q[1];
sx q[1];
rz(-2.2236094) q[1];
x q[2];
rz(2.1928113) q[3];
sx q[3];
rz(-1.6262486) q[3];
sx q[3];
rz(0.90622073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3132402) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(-2.7102615) q[2];
rz(-3.0865772) q[3];
sx q[3];
rz(-0.31156817) q[3];
sx q[3];
rz(0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738275) q[0];
sx q[0];
rz(-0.25280935) q[0];
sx q[0];
rz(1.025169) q[0];
rz(0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-0.90389171) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55306474) q[0];
sx q[0];
rz(-2.4143973) q[0];
sx q[0];
rz(-2.0447879) q[0];
x q[1];
rz(-2.2270062) q[2];
sx q[2];
rz(-2.0981972) q[2];
sx q[2];
rz(1.5114776) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76662582) q[1];
sx q[1];
rz(-1.300525) q[1];
sx q[1];
rz(-0.75281669) q[1];
x q[2];
rz(0.17136161) q[3];
sx q[3];
rz(-1.4079861) q[3];
sx q[3];
rz(-2.5622615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8338833) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(1.5007796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-2.0348771) q[0];
sx q[0];
rz(-2.4427781) q[0];
rz(-1.1837333) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(-2.2898477) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7296556) q[0];
sx q[0];
rz(-2.1067527) q[0];
sx q[0];
rz(0.044117731) q[0];
x q[1];
rz(1.5450655) q[2];
sx q[2];
rz(-2.1479969) q[2];
sx q[2];
rz(-2.4943697) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70515436) q[1];
sx q[1];
rz(-2.2346304) q[1];
sx q[1];
rz(-2.0363807) q[1];
rz(-pi) q[2];
rz(0.61148879) q[3];
sx q[3];
rz(-0.73966714) q[3];
sx q[3];
rz(-2.8304493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.304004) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(2.6336929) q[2];
rz(2.1128283) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(2.540551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4429338) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(3.1319295) q[0];
rz(-1.6161605) q[1];
sx q[1];
rz(-1.7639672) q[1];
sx q[1];
rz(-2.756871) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5670861) q[0];
sx q[0];
rz(-0.28038803) q[0];
sx q[0];
rz(-1.0100421) q[0];
rz(1.184395) q[2];
sx q[2];
rz(-1.7586054) q[2];
sx q[2];
rz(1.9822497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8249121) q[1];
sx q[1];
rz(-2.7928565) q[1];
sx q[1];
rz(2.2168) q[1];
rz(-pi) q[2];
rz(-2.3703897) q[3];
sx q[3];
rz(-1.1386385) q[3];
sx q[3];
rz(1.766253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.559451) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(-0.01586308) q[2];
rz(-1.2265685) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(0.62172833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7568307) q[0];
sx q[0];
rz(-0.66487304) q[0];
sx q[0];
rz(-1.0913947) q[0];
rz(0.81820828) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(-2.4726726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10169928) q[0];
sx q[0];
rz(-0.98838663) q[0];
sx q[0];
rz(2.3132313) q[0];
rz(0.42538504) q[2];
sx q[2];
rz(-0.90406075) q[2];
sx q[2];
rz(-2.3053856) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.131625) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(1.4373006) q[1];
x q[2];
rz(2.0505597) q[3];
sx q[3];
rz(-2.5026813) q[3];
sx q[3];
rz(1.0425488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5074671) q[2];
sx q[2];
rz(-2.6125245) q[2];
sx q[2];
rz(0.96366209) q[2];
rz(1.4108747) q[3];
sx q[3];
rz(-1.7129292) q[3];
sx q[3];
rz(1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016634781) q[0];
sx q[0];
rz(-0.5683012) q[0];
sx q[0];
rz(-1.6868663) q[0];
rz(2.0330873) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(-2.813521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1301918) q[0];
sx q[0];
rz(-1.5131803) q[0];
sx q[0];
rz(1.685623) q[0];
rz(-pi) q[1];
rz(2.5727082) q[2];
sx q[2];
rz(-2.9279104) q[2];
sx q[2];
rz(1.3485707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0340424) q[1];
sx q[1];
rz(-0.75354105) q[1];
sx q[1];
rz(2.5667002) q[1];
rz(-pi) q[2];
rz(-0.107268) q[3];
sx q[3];
rz(-1.7431431) q[3];
sx q[3];
rz(1.337826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6803117) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(2.7247562) q[2];
rz(-2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1998491) q[0];
sx q[0];
rz(-1.2171634) q[0];
sx q[0];
rz(-1.4456277) q[0];
rz(-2.4327714) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(2.9403953) q[2];
sx q[2];
rz(-0.94379776) q[2];
sx q[2];
rz(-1.1371053) q[2];
rz(-2.3117456) q[3];
sx q[3];
rz(-2.0086858) q[3];
sx q[3];
rz(0.84926844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
