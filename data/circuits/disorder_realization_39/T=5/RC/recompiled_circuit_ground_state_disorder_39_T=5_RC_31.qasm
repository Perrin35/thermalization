OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.47867632) q[0];
sx q[0];
rz(7.8362099) q[0];
sx q[0];
rz(10.683164) q[0];
rz(-5.0212669) q[1];
sx q[1];
rz(6.8016383) q[1];
sx q[1];
rz(6.3432884) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3451884) q[0];
sx q[0];
rz(-0.095392659) q[0];
sx q[0];
rz(-0.32890396) q[0];
rz(-2.6491935) q[2];
sx q[2];
rz(-1.4352893) q[2];
sx q[2];
rz(2.8092334) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2138252) q[1];
sx q[1];
rz(-0.37437427) q[1];
sx q[1];
rz(-1.8273749) q[1];
rz(0.93688776) q[3];
sx q[3];
rz(-2.4408615) q[3];
sx q[3];
rz(2.7752293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1897159) q[2];
sx q[2];
rz(-2.3592301) q[2];
sx q[2];
rz(-2.961109) q[2];
rz(-0.30432025) q[3];
sx q[3];
rz(-2.2173827) q[3];
sx q[3];
rz(-2.9144104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3825398) q[0];
sx q[0];
rz(-2.5681684) q[0];
sx q[0];
rz(-3.0376814) q[0];
rz(-0.78481627) q[1];
sx q[1];
rz(-1.3094614) q[1];
sx q[1];
rz(-0.58194247) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73954813) q[0];
sx q[0];
rz(-0.53240896) q[0];
sx q[0];
rz(0.043921434) q[0];
rz(-pi) q[1];
x q[1];
rz(1.21851) q[2];
sx q[2];
rz(-2.0535882) q[2];
sx q[2];
rz(2.5329075) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4247893) q[1];
sx q[1];
rz(-1.5967073) q[1];
sx q[1];
rz(2.0477363) q[1];
rz(-1.9690912) q[3];
sx q[3];
rz(-1.952926) q[3];
sx q[3];
rz(2.7492461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4265784) q[2];
sx q[2];
rz(-2.7535186) q[2];
sx q[2];
rz(1.0283872) q[2];
rz(0.55108023) q[3];
sx q[3];
rz(-1.9469399) q[3];
sx q[3];
rz(-0.50857956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55117637) q[0];
sx q[0];
rz(-1.506378) q[0];
sx q[0];
rz(-1.2580385) q[0];
rz(-0.60351562) q[1];
sx q[1];
rz(-1.3110833) q[1];
sx q[1];
rz(2.9579732) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1574137) q[0];
sx q[0];
rz(-1.3065728) q[0];
sx q[0];
rz(-1.1452894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25213653) q[2];
sx q[2];
rz(-0.49415961) q[2];
sx q[2];
rz(1.3141989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0634707) q[1];
sx q[1];
rz(-2.2666711) q[1];
sx q[1];
rz(0.9982001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19422005) q[3];
sx q[3];
rz(-2.9306539) q[3];
sx q[3];
rz(-0.4926745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6151578) q[2];
sx q[2];
rz(-1.3765455) q[2];
sx q[2];
rz(-2.5395425) q[2];
rz(0.43271068) q[3];
sx q[3];
rz(-1.5475464) q[3];
sx q[3];
rz(1.6656779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.3848307) q[0];
sx q[0];
rz(-2.6669406) q[0];
sx q[0];
rz(2.5153644) q[0];
rz(-1.8734044) q[1];
sx q[1];
rz(-1.7363997) q[1];
sx q[1];
rz(-0.38571206) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2572266) q[0];
sx q[0];
rz(-0.8943087) q[0];
sx q[0];
rz(2.1348597) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5587444) q[2];
sx q[2];
rz(-1.1839317) q[2];
sx q[2];
rz(-0.033535784) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9358237) q[1];
sx q[1];
rz(-0.21142928) q[1];
sx q[1];
rz(0.32306674) q[1];
rz(-pi) q[2];
rz(-1.3703501) q[3];
sx q[3];
rz(-0.55209898) q[3];
sx q[3];
rz(-0.31262661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6993774) q[2];
sx q[2];
rz(-2.0065887) q[2];
sx q[2];
rz(-0.29083148) q[2];
rz(1.95131) q[3];
sx q[3];
rz(-2.5255327) q[3];
sx q[3];
rz(1.0440089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6749343) q[0];
sx q[0];
rz(-3.0789154) q[0];
sx q[0];
rz(1.4516996) q[0];
rz(2.2845204) q[1];
sx q[1];
rz(-1.8430201) q[1];
sx q[1];
rz(-0.42795408) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33855406) q[0];
sx q[0];
rz(-1.1341652) q[0];
sx q[0];
rz(-3.0019041) q[0];
x q[1];
rz(-2.0796135) q[2];
sx q[2];
rz(-2.2428838) q[2];
sx q[2];
rz(0.88269688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.9792069) q[1];
sx q[1];
rz(-2.2507994) q[1];
sx q[1];
rz(1.1080145) q[1];
x q[2];
rz(-2.1367504) q[3];
sx q[3];
rz(-2.1814846) q[3];
sx q[3];
rz(-0.37328675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2032623) q[2];
sx q[2];
rz(-1.0469971) q[2];
sx q[2];
rz(-0.14673512) q[2];
rz(-2.9444368) q[3];
sx q[3];
rz(-1.6517755) q[3];
sx q[3];
rz(-0.76876918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7305304) q[0];
sx q[0];
rz(-1.0163607) q[0];
sx q[0];
rz(2.9732669) q[0];
rz(2.7492145) q[1];
sx q[1];
rz(-0.43586755) q[1];
sx q[1];
rz(-1.6900774) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0453059) q[0];
sx q[0];
rz(-1.1340965) q[0];
sx q[0];
rz(-0.023604579) q[0];
rz(-pi) q[1];
rz(2.2915527) q[2];
sx q[2];
rz(-1.3728567) q[2];
sx q[2];
rz(0.70591656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5785006) q[1];
sx q[1];
rz(-0.76116409) q[1];
sx q[1];
rz(-0.20462402) q[1];
rz(-pi) q[2];
rz(-1.1216954) q[3];
sx q[3];
rz(-2.3074352) q[3];
sx q[3];
rz(-2.4569805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3671941) q[2];
sx q[2];
rz(-0.44782475) q[2];
sx q[2];
rz(-1.1773342) q[2];
rz(0.92159739) q[3];
sx q[3];
rz(-0.84601837) q[3];
sx q[3];
rz(0.53068501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0717936) q[0];
sx q[0];
rz(-1.8853747) q[0];
sx q[0];
rz(1.0569093) q[0];
rz(0.22843703) q[1];
sx q[1];
rz(-0.7799131) q[1];
sx q[1];
rz(2.8905919) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042648249) q[0];
sx q[0];
rz(-1.3146719) q[0];
sx q[0];
rz(0.29447181) q[0];
rz(-pi) q[1];
rz(2.2728738) q[2];
sx q[2];
rz(-1.0660604) q[2];
sx q[2];
rz(2.2622893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.31902203) q[1];
sx q[1];
rz(-2.9290556) q[1];
sx q[1];
rz(0.6937278) q[1];
rz(-1.9422705) q[3];
sx q[3];
rz(-1.087093) q[3];
sx q[3];
rz(-2.7886645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.935427) q[2];
sx q[2];
rz(-0.82038227) q[2];
sx q[2];
rz(2.7247735) q[2];
rz(-2.5721512) q[3];
sx q[3];
rz(-2.0408401) q[3];
sx q[3];
rz(2.4429564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5702629) q[0];
sx q[0];
rz(-1.2705734) q[0];
sx q[0];
rz(-0.50462333) q[0];
rz(2.8847671) q[1];
sx q[1];
rz(-0.80373126) q[1];
sx q[1];
rz(1.9907192) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68783142) q[0];
sx q[0];
rz(-2.6021299) q[0];
sx q[0];
rz(-1.1629172) q[0];
rz(-pi) q[1];
rz(0.93575259) q[2];
sx q[2];
rz(-1.2798736) q[2];
sx q[2];
rz(-2.3036602) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2765892) q[1];
sx q[1];
rz(-2.114608) q[1];
sx q[1];
rz(1.0629884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4126115) q[3];
sx q[3];
rz(-1.9831267) q[3];
sx q[3];
rz(1.3875426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5095832) q[2];
sx q[2];
rz(-0.42319599) q[2];
sx q[2];
rz(-2.7022341) q[2];
rz(1.1257233) q[3];
sx q[3];
rz(-1.6042387) q[3];
sx q[3];
rz(1.841898) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706962) q[0];
sx q[0];
rz(-2.3211711) q[0];
sx q[0];
rz(-2.6743555) q[0];
rz(2.1029419) q[1];
sx q[1];
rz(-0.50913441) q[1];
sx q[1];
rz(-1.4899303) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56541967) q[0];
sx q[0];
rz(-1.3955981) q[0];
sx q[0];
rz(-2.6433667) q[0];
x q[1];
rz(-1.4760255) q[2];
sx q[2];
rz(-1.7168003) q[2];
sx q[2];
rz(-2.3340747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6789602) q[1];
sx q[1];
rz(-0.81607807) q[1];
sx q[1];
rz(3.1051293) q[1];
x q[2];
rz(0.64405264) q[3];
sx q[3];
rz(-0.67094147) q[3];
sx q[3];
rz(1.3356127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6963639) q[2];
sx q[2];
rz(-2.7859521) q[2];
sx q[2];
rz(-1.5544372) q[2];
rz(0.98443952) q[3];
sx q[3];
rz(-0.90960228) q[3];
sx q[3];
rz(-1.4136774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5928891) q[0];
sx q[0];
rz(-0.84243542) q[0];
sx q[0];
rz(0.62826759) q[0];
rz(-2.9382622) q[1];
sx q[1];
rz(-2.0275828) q[1];
sx q[1];
rz(0.59439739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0258163) q[0];
sx q[0];
rz(-0.62827841) q[0];
sx q[0];
rz(1.216287) q[0];
rz(-pi) q[1];
rz(1.3514694) q[2];
sx q[2];
rz(-2.0989053) q[2];
sx q[2];
rz(-2.2161432) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0600784) q[1];
sx q[1];
rz(-1.3003254) q[1];
sx q[1];
rz(-2.8287776) q[1];
rz(0.61305586) q[3];
sx q[3];
rz(-1.5125015) q[3];
sx q[3];
rz(-1.6976274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9336885) q[2];
sx q[2];
rz(-0.7242569) q[2];
sx q[2];
rz(2.9653463) q[2];
rz(1.8267953) q[3];
sx q[3];
rz(-2.3254471) q[3];
sx q[3];
rz(-2.3740681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15277302) q[0];
sx q[0];
rz(-1.7012699) q[0];
sx q[0];
rz(-1.758601) q[0];
rz(-0.11463166) q[1];
sx q[1];
rz(-0.75420598) q[1];
sx q[1];
rz(2.4892714) q[1];
rz(2.9822275) q[2];
sx q[2];
rz(-2.9014316) q[2];
sx q[2];
rz(0.17221774) q[2];
rz(-0.38269855) q[3];
sx q[3];
rz(-2.0321587) q[3];
sx q[3];
rz(1.7176499) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
