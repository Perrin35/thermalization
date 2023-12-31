OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55460632) q[0];
sx q[0];
rz(2.245683) q[0];
sx q[0];
rz(10.829344) q[0];
rz(3.8430619) q[1];
sx q[1];
rz(3.7624533) q[1];
sx q[1];
rz(7.9384595) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0151057) q[0];
sx q[0];
rz(-1.4529714) q[0];
sx q[0];
rz(-0.18527041) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15847023) q[2];
sx q[2];
rz(-2.7140693) q[2];
sx q[2];
rz(-2.1536364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9857786) q[1];
sx q[1];
rz(-1.6667546) q[1];
sx q[1];
rz(2.7491261) q[1];
rz(-pi) q[2];
rz(-0.078943723) q[3];
sx q[3];
rz(-1.5642089) q[3];
sx q[3];
rz(1.5353257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0861337) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(0.80111516) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.79008094) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(0.077985667) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(2.7819113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65950899) q[0];
sx q[0];
rz(-2.2523899) q[0];
sx q[0];
rz(0.73007749) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7799136) q[2];
sx q[2];
rz(-0.86161999) q[2];
sx q[2];
rz(1.4271971) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3188394) q[1];
sx q[1];
rz(-2.0944632) q[1];
sx q[1];
rz(-0.19284064) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5248604) q[3];
sx q[3];
rz(-2.4270298) q[3];
sx q[3];
rz(3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(-1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(2.673972) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(2.8288249) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395961) q[0];
sx q[0];
rz(-1.6507286) q[0];
sx q[0];
rz(-1.3586112) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8094693) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(2.3534564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.410027) q[1];
sx q[1];
rz(-1.7419852) q[1];
sx q[1];
rz(2.9982135) q[1];
rz(-pi) q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(-1.0696166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(0.30511937) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6453648) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(-0.67273295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3161635) q[0];
sx q[0];
rz(-2.5220693) q[0];
sx q[0];
rz(-2.6114458) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8276617) q[2];
sx q[2];
rz(-1.0360403) q[2];
sx q[2];
rz(-1.9321439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2697301) q[1];
sx q[1];
rz(-1.4566263) q[1];
sx q[1];
rz(-1.8037379) q[1];
x q[2];
rz(2.1898515) q[3];
sx q[3];
rz(-1.4674867) q[3];
sx q[3];
rz(-0.88224525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(-0.76081863) q[2];
rz(-0.86756724) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(2.0779579) q[0];
rz(-1.4798374) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(0.038287727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69890755) q[0];
sx q[0];
rz(-1.2783056) q[0];
sx q[0];
rz(-0.12683503) q[0];
rz(-0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(-0.40508168) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8748862) q[1];
sx q[1];
rz(-1.6883856) q[1];
sx q[1];
rz(1.7456732) q[1];
rz(-pi) q[2];
rz(-1.0682265) q[3];
sx q[3];
rz(-0.52191075) q[3];
sx q[3];
rz(-3.0038602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(2.0495474) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(-0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096481) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(-0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(0.19097701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.068114) q[0];
sx q[0];
rz(-1.3604135) q[0];
sx q[0];
rz(2.99899) q[0];
x q[1];
rz(-0.35031788) q[2];
sx q[2];
rz(-0.91674524) q[2];
sx q[2];
rz(-1.0612812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8543507) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(-2.1139305) q[1];
rz(-2.0885002) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(0.5047441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(-2.7065275) q[2];
rz(-2.2502031) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1973535) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(-2.7213668) q[0];
rz(0.48121437) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(1.0302672) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7468443) q[0];
sx q[0];
rz(-2.0472102) q[0];
sx q[0];
rz(-1.9917411) q[0];
rz(3.0646965) q[2];
sx q[2];
rz(-1.352407) q[2];
sx q[2];
rz(1.0181392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3863341) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(-0.30524409) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7610735) q[3];
sx q[3];
rz(-1.7576302) q[3];
sx q[3];
rz(0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.9497005) q[0];
rz(0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(1.0620767) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5545643) q[0];
sx q[0];
rz(-1.4438859) q[0];
sx q[0];
rz(2.1112604) q[0];
rz(-pi) q[1];
rz(-2.383963) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(2.4766395) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8909104) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(-2.2834884) q[1];
rz(-0.90537269) q[3];
sx q[3];
rz(-1.5983799) q[3];
sx q[3];
rz(-1.8679004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(0.051076802) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(2.0843845) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99701571) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-2.2924246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6309755) q[0];
sx q[0];
rz(-1.346934) q[0];
sx q[0];
rz(-1.1918681) q[0];
x q[1];
rz(-1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(0.6347444) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2850212) q[1];
sx q[1];
rz(-1.3339431) q[1];
sx q[1];
rz(-2.9076438) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9407223) q[3];
sx q[3];
rz(-1.4088653) q[3];
sx q[3];
rz(-2.014132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(-2.3000439) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(2.9291005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8520078) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(-1.1803124) q[0];
rz(-pi) q[1];
x q[1];
rz(0.041065644) q[2];
sx q[2];
rz(-1.8827056) q[2];
sx q[2];
rz(1.2817771) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6007235) q[1];
sx q[1];
rz(-1.0717372) q[1];
sx q[1];
rz(-0.869511) q[1];
rz(-pi) q[2];
rz(-2.859697) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(-2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.4546222) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385948) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(-0.075642792) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(1.5008925) q[2];
sx q[2];
rz(-1.5887661) q[2];
sx q[2];
rz(-0.80815732) q[2];
rz(2.4235519) q[3];
sx q[3];
rz(-0.64648872) q[3];
sx q[3];
rz(1.7512376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
