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
rz(-2.5290639) q[0];
sx q[0];
rz(-0.6331442) q[0];
sx q[0];
rz(2.7781515) q[0];
rz(2.5545622) q[1];
sx q[1];
rz(-0.5293923) q[1];
sx q[1];
rz(1.6338978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8315805) q[0];
sx q[0];
rz(-0.10619199) q[0];
sx q[0];
rz(-2.0640762) q[0];
rz(-pi) q[1];
rz(2.426067) q[2];
sx q[2];
rz(-1.160485) q[2];
sx q[2];
rz(0.35718756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22544293) q[1];
sx q[1];
rz(-2.5366211) q[1];
sx q[1];
rz(-1.1892645) q[1];
rz(-pi) q[2];
rz(-0.81256688) q[3];
sx q[3];
rz(-1.9877533) q[3];
sx q[3];
rz(-1.6171232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2734566) q[2];
sx q[2];
rz(-0.09082219) q[2];
sx q[2];
rz(2.691213) q[2];
rz(-0.069933794) q[3];
sx q[3];
rz(-1.115333) q[3];
sx q[3];
rz(-0.44225881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0792585) q[0];
sx q[0];
rz(-2.313995) q[0];
sx q[0];
rz(0.29912478) q[0];
rz(-2.4807855) q[1];
sx q[1];
rz(-1.1102763) q[1];
sx q[1];
rz(2.2181236) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7851334) q[0];
sx q[0];
rz(-1.3970427) q[0];
sx q[0];
rz(1.2453402) q[0];
x q[1];
rz(-2.9060676) q[2];
sx q[2];
rz(-1.356429) q[2];
sx q[2];
rz(1.216824) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3287433) q[1];
sx q[1];
rz(-1.5069811) q[1];
sx q[1];
rz(-2.7743474) q[1];
rz(-pi) q[2];
rz(-2.1577182) q[3];
sx q[3];
rz(-1.3345342) q[3];
sx q[3];
rz(2.7191702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6580761) q[2];
sx q[2];
rz(-2.7407756) q[2];
sx q[2];
rz(0.083902396) q[2];
rz(0.99533254) q[3];
sx q[3];
rz(-1.9845125) q[3];
sx q[3];
rz(-1.6307433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8453269) q[0];
sx q[0];
rz(-1.9060598) q[0];
sx q[0];
rz(-0.96257019) q[0];
rz(2.8341017) q[1];
sx q[1];
rz(-0.88491416) q[1];
sx q[1];
rz(-2.9226551) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1277351) q[0];
sx q[0];
rz(-1.225219) q[0];
sx q[0];
rz(1.6446949) q[0];
x q[1];
rz(-1.8586765) q[2];
sx q[2];
rz(-2.4426443) q[2];
sx q[2];
rz(-1.3290249) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.98597279) q[1];
sx q[1];
rz(-0.39955214) q[1];
sx q[1];
rz(-1.5635207) q[1];
x q[2];
rz(2.5265273) q[3];
sx q[3];
rz(-0.99216539) q[3];
sx q[3];
rz(2.8556153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4585877) q[2];
sx q[2];
rz(-1.747749) q[2];
sx q[2];
rz(1.638442) q[2];
rz(-0.25034869) q[3];
sx q[3];
rz(-0.72385794) q[3];
sx q[3];
rz(2.943128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64971626) q[0];
sx q[0];
rz(-2.8174077) q[0];
sx q[0];
rz(2.8247996) q[0];
rz(-3.1253539) q[1];
sx q[1];
rz(-0.79252487) q[1];
sx q[1];
rz(-2.8932755) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93706276) q[0];
sx q[0];
rz(-2.2872675) q[0];
sx q[0];
rz(-0.69964827) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9944579) q[2];
sx q[2];
rz(-2.6700041) q[2];
sx q[2];
rz(-1.4855282) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9435472) q[1];
sx q[1];
rz(-1.24453) q[1];
sx q[1];
rz(-1.7034143) q[1];
rz(-pi) q[2];
rz(-1.3518356) q[3];
sx q[3];
rz(-1.7585424) q[3];
sx q[3];
rz(-2.2098372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.877964) q[2];
sx q[2];
rz(-0.10136494) q[2];
sx q[2];
rz(-1.144009) q[2];
rz(0.7621724) q[3];
sx q[3];
rz(-2.6165369) q[3];
sx q[3];
rz(-2.150056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85631973) q[0];
sx q[0];
rz(-3.0506436) q[0];
sx q[0];
rz(-1.0524622) q[0];
rz(0.4902803) q[1];
sx q[1];
rz(-1.0888381) q[1];
sx q[1];
rz(3.0163684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76716048) q[0];
sx q[0];
rz(-1.1781791) q[0];
sx q[0];
rz(0.36628337) q[0];
x q[1];
rz(-2.5936761) q[2];
sx q[2];
rz(-1.3306018) q[2];
sx q[2];
rz(-1.5870767) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0908806) q[1];
sx q[1];
rz(-1.6428324) q[1];
sx q[1];
rz(1.8211831) q[1];
rz(-pi) q[2];
rz(0.85876561) q[3];
sx q[3];
rz(-1.7384343) q[3];
sx q[3];
rz(-1.110838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6584955) q[2];
sx q[2];
rz(-0.53762895) q[2];
sx q[2];
rz(-1.5888692) q[2];
rz(3.0850278) q[3];
sx q[3];
rz(-1.6274692) q[3];
sx q[3];
rz(-0.27730832) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5796984) q[0];
sx q[0];
rz(-2.6772006) q[0];
sx q[0];
rz(-0.42771801) q[0];
rz(-2.600421) q[1];
sx q[1];
rz(-2.7310557) q[1];
sx q[1];
rz(1.7778832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83374524) q[0];
sx q[0];
rz(-2.5302093) q[0];
sx q[0];
rz(1.7627384) q[0];
x q[1];
rz(-0.42714675) q[2];
sx q[2];
rz(-1.5183798) q[2];
sx q[2];
rz(0.77709955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1178255) q[1];
sx q[1];
rz(-2.3741873) q[1];
sx q[1];
rz(1.1249705) q[1];
rz(-pi) q[2];
x q[2];
rz(1.545752) q[3];
sx q[3];
rz(-0.4129172) q[3];
sx q[3];
rz(-0.027137952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0987739) q[2];
sx q[2];
rz(-2.6218178) q[2];
sx q[2];
rz(0.40936142) q[2];
rz(2.9218819) q[3];
sx q[3];
rz(-1.5922056) q[3];
sx q[3];
rz(1.4896013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4396502) q[0];
sx q[0];
rz(-2.7203429) q[0];
sx q[0];
rz(-0.14081328) q[0];
rz(-0.5367865) q[1];
sx q[1];
rz(-1.355492) q[1];
sx q[1];
rz(0.37667325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0842782) q[0];
sx q[0];
rz(-1.8207826) q[0];
sx q[0];
rz(-2.9765073) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67492731) q[2];
sx q[2];
rz(-1.2528193) q[2];
sx q[2];
rz(-0.392158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.46493775) q[1];
sx q[1];
rz(-2.3466517) q[1];
sx q[1];
rz(-1.388768) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34318681) q[3];
sx q[3];
rz(-2.7036441) q[3];
sx q[3];
rz(-1.393895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.10304) q[2];
sx q[2];
rz(-1.1662177) q[2];
sx q[2];
rz(2.340359) q[2];
rz(-1.0771105) q[3];
sx q[3];
rz(-0.25682768) q[3];
sx q[3];
rz(0.1424772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.24256663) q[0];
sx q[0];
rz(-2.9809256) q[0];
sx q[0];
rz(-2.2005079) q[0];
rz(-2.2451027) q[1];
sx q[1];
rz(-1.4184378) q[1];
sx q[1];
rz(-0.81472188) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3692432) q[0];
sx q[0];
rz(-0.54074484) q[0];
sx q[0];
rz(2.3750099) q[0];
rz(0.28622921) q[2];
sx q[2];
rz(-0.17795086) q[2];
sx q[2];
rz(-2.9432757) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9174609) q[1];
sx q[1];
rz(-2.1814709) q[1];
sx q[1];
rz(0.4179722) q[1];
rz(-0.37121935) q[3];
sx q[3];
rz(-2.1069848) q[3];
sx q[3];
rz(-0.49193383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61646378) q[2];
sx q[2];
rz(-2.5069671) q[2];
sx q[2];
rz(-2.5647707) q[2];
rz(2.6680464) q[3];
sx q[3];
rz(-1.1487995) q[3];
sx q[3];
rz(-0.043341652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43805495) q[0];
sx q[0];
rz(-0.75551581) q[0];
sx q[0];
rz(-0.087906539) q[0];
rz(-1.5880623) q[1];
sx q[1];
rz(-0.5843662) q[1];
sx q[1];
rz(-0.25762525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15730298) q[0];
sx q[0];
rz(-2.7098795) q[0];
sx q[0];
rz(-1.2832416) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1981774) q[2];
sx q[2];
rz(-1.6281272) q[2];
sx q[2];
rz(-1.0939404) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0637197) q[1];
sx q[1];
rz(-1.8053836) q[1];
sx q[1];
rz(-0.44198702) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78225953) q[3];
sx q[3];
rz(-0.74499797) q[3];
sx q[3];
rz(1.9375325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0660925) q[2];
sx q[2];
rz(-1.6801715) q[2];
sx q[2];
rz(-0.077787682) q[2];
rz(-1.0545688) q[3];
sx q[3];
rz(-2.4013897) q[3];
sx q[3];
rz(-0.65092248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.2442378) q[0];
sx q[0];
rz(-0.033441823) q[0];
sx q[0];
rz(0.43575409) q[0];
rz(-0.032595366) q[1];
sx q[1];
rz(-2.4041912) q[1];
sx q[1];
rz(-1.3479412) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7028026) q[0];
sx q[0];
rz(-2.3910671) q[0];
sx q[0];
rz(-3.1311959) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6145337) q[2];
sx q[2];
rz(-1.3599476) q[2];
sx q[2];
rz(1.2817037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7256335) q[1];
sx q[1];
rz(-1.6068649) q[1];
sx q[1];
rz(-3.1276324) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3147821) q[3];
sx q[3];
rz(-1.3482932) q[3];
sx q[3];
rz(0.088929847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2160448) q[2];
sx q[2];
rz(-1.0341158) q[2];
sx q[2];
rz(2.3852868) q[2];
rz(-3.0594399) q[3];
sx q[3];
rz(-2.7492456) q[3];
sx q[3];
rz(-2.4503585) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77770752) q[0];
sx q[0];
rz(-0.73962279) q[0];
sx q[0];
rz(2.2864322) q[0];
rz(-1.1126407) q[1];
sx q[1];
rz(-2.009195) q[1];
sx q[1];
rz(2.7809273) q[1];
rz(1.3203717) q[2];
sx q[2];
rz(-1.741521) q[2];
sx q[2];
rz(2.3157816) q[2];
rz(2.6769911) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(-0.28649022) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
