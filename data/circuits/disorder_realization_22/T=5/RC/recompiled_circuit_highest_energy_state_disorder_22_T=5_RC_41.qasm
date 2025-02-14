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
rz(-0.67996112) q[0];
sx q[0];
rz(-0.90587076) q[0];
sx q[0];
rz(-1.0933956) q[0];
rz(-1.6361341) q[1];
sx q[1];
rz(-1.1727762) q[1];
sx q[1];
rz(2.7170031) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68442518) q[0];
sx q[0];
rz(-2.9363605) q[0];
sx q[0];
rz(1.7590909) q[0];
rz(-pi) q[1];
rz(1.7248254) q[2];
sx q[2];
rz(-1.5206778) q[2];
sx q[2];
rz(-0.94238867) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5988785) q[1];
sx q[1];
rz(-1.492036) q[1];
sx q[1];
rz(-1.3244737) q[1];
x q[2];
rz(0.13386676) q[3];
sx q[3];
rz(-2.2765719) q[3];
sx q[3];
rz(5.2701252e-05) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7666185) q[2];
sx q[2];
rz(-1.5593636) q[2];
sx q[2];
rz(-1.1589104) q[2];
rz(0.22289395) q[3];
sx q[3];
rz(-1.4054207) q[3];
sx q[3];
rz(-2.8515653) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952154) q[0];
sx q[0];
rz(-1.615849) q[0];
sx q[0];
rz(-1.079153) q[0];
rz(1.7201299) q[1];
sx q[1];
rz(-2.454897) q[1];
sx q[1];
rz(3.0770643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33693211) q[0];
sx q[0];
rz(-1.48359) q[0];
sx q[0];
rz(0.052636458) q[0];
x q[1];
rz(-1.5987492) q[2];
sx q[2];
rz(-1.0947026) q[2];
sx q[2];
rz(-2.3700489) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77212438) q[1];
sx q[1];
rz(-1.8126948) q[1];
sx q[1];
rz(-2.3568704) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3202479) q[3];
sx q[3];
rz(-1.4331766) q[3];
sx q[3];
rz(2.2087165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1530389) q[2];
sx q[2];
rz(-1.7861853) q[2];
sx q[2];
rz(2.0241731) q[2];
rz(2.2828263) q[3];
sx q[3];
rz(-2.358181) q[3];
sx q[3];
rz(2.7045238) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3808463) q[0];
sx q[0];
rz(-0.28457156) q[0];
sx q[0];
rz(0.17307702) q[0];
rz(-0.95678798) q[1];
sx q[1];
rz(-2.1134977) q[1];
sx q[1];
rz(-2.079336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040471023) q[0];
sx q[0];
rz(-0.72325828) q[0];
sx q[0];
rz(1.5803807) q[0];
rz(-0.71346475) q[2];
sx q[2];
rz(-2.3629945) q[2];
sx q[2];
rz(0.32003357) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.838321) q[1];
sx q[1];
rz(-2.1049989) q[1];
sx q[1];
rz(0.81373416) q[1];
rz(-pi) q[2];
rz(0.62543243) q[3];
sx q[3];
rz(-2.1533826) q[3];
sx q[3];
rz(0.39804493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.7303077) q[2];
sx q[2];
rz(-1.2016808) q[2];
sx q[2];
rz(-2.3968706) q[2];
rz(-2.5005285) q[3];
sx q[3];
rz(-1.791626) q[3];
sx q[3];
rz(0.49066576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.6915801) q[0];
sx q[0];
rz(-0.55679655) q[0];
sx q[0];
rz(-2.2963754) q[0];
rz(-2.7382964) q[1];
sx q[1];
rz(-1.6019628) q[1];
sx q[1];
rz(-2.8135615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75719365) q[0];
sx q[0];
rz(-1.4310657) q[0];
sx q[0];
rz(-1.3881042) q[0];
x q[1];
rz(-2.9035527) q[2];
sx q[2];
rz(-1.5063707) q[2];
sx q[2];
rz(3.079133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33095903) q[1];
sx q[1];
rz(-1.0077268) q[1];
sx q[1];
rz(-1.710239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27356903) q[3];
sx q[3];
rz(-2.0988587) q[3];
sx q[3];
rz(2.7298965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24568096) q[2];
sx q[2];
rz(-2.3978105) q[2];
sx q[2];
rz(-2.9823859) q[2];
rz(3.1253452) q[3];
sx q[3];
rz(-1.0280321) q[3];
sx q[3];
rz(2.4102559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943587) q[0];
sx q[0];
rz(-1.9488229) q[0];
sx q[0];
rz(-2.0514945) q[0];
rz(-3.0116426) q[1];
sx q[1];
rz(-2.3268685) q[1];
sx q[1];
rz(-1.5257588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9620044) q[0];
sx q[0];
rz(-1.9854913) q[0];
sx q[0];
rz(2.5432822) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1431057) q[2];
sx q[2];
rz(-2.9614355) q[2];
sx q[2];
rz(-3.0017972) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33541778) q[1];
sx q[1];
rz(-1.5433784) q[1];
sx q[1];
rz(-0.63372483) q[1];
rz(-pi) q[2];
x q[2];
rz(0.035338621) q[3];
sx q[3];
rz(-1.3900847) q[3];
sx q[3];
rz(2.904195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7502363) q[2];
sx q[2];
rz(-0.82166925) q[2];
sx q[2];
rz(0.4772805) q[2];
rz(-0.20398772) q[3];
sx q[3];
rz(-2.9450649) q[3];
sx q[3];
rz(-2.5175214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9452962) q[0];
sx q[0];
rz(-0.81552234) q[0];
sx q[0];
rz(2.3639009) q[0];
rz(-0.6174736) q[1];
sx q[1];
rz(-1.6505046) q[1];
sx q[1];
rz(-2.2106574) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8645267) q[0];
sx q[0];
rz(-2.0210938) q[0];
sx q[0];
rz(-0.45876512) q[0];
x q[1];
rz(-1.0920877) q[2];
sx q[2];
rz(-1.1943814) q[2];
sx q[2];
rz(2.7680754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8691784) q[1];
sx q[1];
rz(-0.34124103) q[1];
sx q[1];
rz(-1.7029087) q[1];
rz(2.191675) q[3];
sx q[3];
rz(-1.2395879) q[3];
sx q[3];
rz(-0.036605926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65333873) q[2];
sx q[2];
rz(-1.3335756) q[2];
sx q[2];
rz(-0.2674357) q[2];
rz(-0.93287647) q[3];
sx q[3];
rz(-1.8578015) q[3];
sx q[3];
rz(-0.4944087) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6701508) q[0];
sx q[0];
rz(-1.143456) q[0];
sx q[0];
rz(-0.66584051) q[0];
rz(2.937607) q[1];
sx q[1];
rz(-0.98122707) q[1];
sx q[1];
rz(1.9542255) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68722938) q[0];
sx q[0];
rz(-0.47894127) q[0];
sx q[0];
rz(-2.9402551) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1631931) q[2];
sx q[2];
rz(-2.3401988) q[2];
sx q[2];
rz(1.8155542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5449808) q[1];
sx q[1];
rz(-0.89681936) q[1];
sx q[1];
rz(-1.9121714) q[1];
rz(-1.7660308) q[3];
sx q[3];
rz(-0.25849202) q[3];
sx q[3];
rz(0.067276567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5953956) q[2];
sx q[2];
rz(-0.51837817) q[2];
sx q[2];
rz(-2.4714244) q[2];
rz(-0.3012805) q[3];
sx q[3];
rz(-1.8471085) q[3];
sx q[3];
rz(0.41845751) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30348521) q[0];
sx q[0];
rz(-2.1391588) q[0];
sx q[0];
rz(1.0656892) q[0];
rz(2.4844555) q[1];
sx q[1];
rz(-2.4333351) q[1];
sx q[1];
rz(1.4297952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098786548) q[0];
sx q[0];
rz(-2.1273861) q[0];
sx q[0];
rz(-0.3996398) q[0];
rz(-pi) q[1];
rz(-2.6515342) q[2];
sx q[2];
rz(-2.769751) q[2];
sx q[2];
rz(-1.8143559) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89420477) q[1];
sx q[1];
rz(-0.46890837) q[1];
sx q[1];
rz(-0.50007485) q[1];
rz(-2.4938857) q[3];
sx q[3];
rz(-0.2362902) q[3];
sx q[3];
rz(-1.7648197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26057217) q[2];
sx q[2];
rz(-1.2994956) q[2];
sx q[2];
rz(2.1232429) q[2];
rz(-1.5923422) q[3];
sx q[3];
rz(-1.5038303) q[3];
sx q[3];
rz(-0.75756592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8757979) q[0];
sx q[0];
rz(-0.51545155) q[0];
sx q[0];
rz(-2.1739668) q[0];
rz(0.34010092) q[1];
sx q[1];
rz(-2.3022771) q[1];
sx q[1];
rz(-2.1077154) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144185) q[0];
sx q[0];
rz(-2.5670739) q[0];
sx q[0];
rz(1.3142725) q[0];
rz(1.1393093) q[2];
sx q[2];
rz(-2.1528006) q[2];
sx q[2];
rz(-1.7503357) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6996697) q[1];
sx q[1];
rz(-1.1811387) q[1];
sx q[1];
rz(0.81894919) q[1];
x q[2];
rz(1.7179862) q[3];
sx q[3];
rz(-1.7924252) q[3];
sx q[3];
rz(2.2972884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10244441) q[2];
sx q[2];
rz(-1.4463964) q[2];
sx q[2];
rz(0.038912494) q[2];
rz(3.0012722) q[3];
sx q[3];
rz(-0.084429927) q[3];
sx q[3];
rz(1.3154359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4843531) q[0];
sx q[0];
rz(-1.0056647) q[0];
sx q[0];
rz(-1.2580309) q[0];
rz(0.21615061) q[1];
sx q[1];
rz(-0.70544568) q[1];
sx q[1];
rz(-0.36453882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1572239) q[0];
sx q[0];
rz(-1.8175392) q[0];
sx q[0];
rz(-1.1086199) q[0];
rz(2.6993518) q[2];
sx q[2];
rz(-0.85390515) q[2];
sx q[2];
rz(-0.88767641) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8382524) q[1];
sx q[1];
rz(-1.5152182) q[1];
sx q[1];
rz(1.5838632) q[1];
rz(0.057897827) q[3];
sx q[3];
rz(-1.6030884) q[3];
sx q[3];
rz(-2.3253289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2898499) q[2];
sx q[2];
rz(-1.936329) q[2];
sx q[2];
rz(0.64129788) q[2];
rz(-1.4801721) q[3];
sx q[3];
rz(-1.8931754) q[3];
sx q[3];
rz(-0.67354584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4079473) q[0];
sx q[0];
rz(-2.6483364) q[0];
sx q[0];
rz(-0.37089621) q[0];
rz(-2.1570878) q[1];
sx q[1];
rz(-1.6592818) q[1];
sx q[1];
rz(2.0936113) q[1];
rz(1.0417837) q[2];
sx q[2];
rz(-2.9137901) q[2];
sx q[2];
rz(0.13005039) q[2];
rz(1.0897286) q[3];
sx q[3];
rz(-1.0952598) q[3];
sx q[3];
rz(2.918603) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
