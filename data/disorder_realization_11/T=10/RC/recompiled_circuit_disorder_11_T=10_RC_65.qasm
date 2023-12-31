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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0044999997) q[0];
sx q[0];
rz(-2.922393) q[0];
sx q[0];
rz(-2.5704434) q[0];
rz(1.4990184) q[2];
sx q[2];
rz(-1.1489747) q[2];
sx q[2];
rz(-2.3274802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3753429) q[1];
sx q[1];
rz(-1.9613593) q[1];
sx q[1];
rz(1.4669963) q[1];
rz(-pi) q[2];
x q[2];
rz(0.078943723) q[3];
sx q[3];
rz(-1.5773838) q[3];
sx q[3];
rz(-1.6062669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(0.80111516) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-3.063607) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(2.7819113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8436444) q[0];
sx q[0];
rz(-0.95373646) q[0];
sx q[0];
rz(2.258837) q[0];
rz(0.72007911) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(3.135323) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45047346) q[1];
sx q[1];
rz(-2.5866661) q[1];
sx q[1];
rz(1.2503442) q[1];
rz(-pi) q[2];
rz(3.1017786) q[3];
sx q[3];
rz(-2.2844444) q[3];
sx q[3];
rz(-3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(2.673972) q[0];
rz(-0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62374672) q[0];
sx q[0];
rz(-0.22652921) q[0];
sx q[0];
rz(-1.2073327) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8094693) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(-0.78813625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2777729) q[1];
sx q[1];
rz(-1.4295271) q[1];
sx q[1];
rz(1.743725) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7160277) q[3];
sx q[3];
rz(-1.8947453) q[3];
sx q[3];
rz(2.5026929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9096845) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(-2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(0.028045068) q[0];
rz(-1.6758945) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(2.4688597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3161635) q[0];
sx q[0];
rz(-0.61952335) q[0];
sx q[0];
rz(-2.6114458) q[0];
rz(-pi) q[1];
rz(0.54949923) q[2];
sx q[2];
rz(-1.7911583) q[2];
sx q[2];
rz(0.22827497) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67191254) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(3.0242821) q[1];
x q[2];
rz(-1.7476051) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(0.76081863) q[2];
rz(2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(-0.76604617) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(1.0636348) q[0];
rz(1.6617552) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-0.038287727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1153039) q[0];
sx q[0];
rz(-0.3180819) q[0];
sx q[0];
rz(1.9684857) q[0];
rz(-pi) q[1];
rz(1.080004) q[2];
sx q[2];
rz(-2.12825) q[2];
sx q[2];
rz(1.4412396) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89025154) q[1];
sx q[1];
rz(-2.9311935) q[1];
sx q[1];
rz(0.9743147) q[1];
rz(-pi) q[2];
rz(2.0733662) q[3];
sx q[3];
rz(-0.52191075) q[3];
sx q[3];
rz(0.13773242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15497196) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(2.0495474) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(2.3585414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.18096481) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(-3.0788105) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(-0.19097701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67577772) q[0];
sx q[0];
rz(-2.8880279) q[0];
sx q[0];
rz(0.98357865) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.149876) q[2];
sx q[2];
rz(-0.72962609) q[2];
sx q[2];
rz(2.6211477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8543507) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(-1.0276621) q[1];
x q[2];
rz(-2.7950068) q[3];
sx q[3];
rz(-2.1086018) q[3];
sx q[3];
rz(2.0227405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-0.43506518) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(-0.42022589) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-1.0302672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39474836) q[0];
sx q[0];
rz(-2.0472102) q[0];
sx q[0];
rz(1.1498515) q[0];
rz(-1.7898125) q[2];
sx q[2];
rz(-1.6458626) q[2];
sx q[2];
rz(-0.53596562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.269304) q[1];
sx q[1];
rz(-1.270712) q[1];
sx q[1];
rz(-1.7606723) q[1];
rz(-pi) q[2];
rz(-0.19018634) q[3];
sx q[3];
rz(-1.38387) q[3];
sx q[3];
rz(-0.72894064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(1.9497005) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-2.079516) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5870283) q[0];
sx q[0];
rz(-1.6977068) q[0];
sx q[0];
rz(-1.0303322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7161712) q[2];
sx q[2];
rz(-1.9429978) q[2];
sx q[2];
rz(-2.9120955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2506822) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(2.2834884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5261371) q[3];
sx q[3];
rz(-2.4756845) q[3];
sx q[3];
rz(0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.186211) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(0.051076802) q[2];
rz(0.71930277) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99701571) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-2.4966519) q[0];
rz(1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97201559) q[0];
sx q[0];
rz(-1.2017842) q[0];
sx q[0];
rz(-0.24032648) q[0];
rz(-pi) q[1];
rz(0.37372132) q[2];
sx q[2];
rz(-1.4288651) q[2];
sx q[2];
rz(2.2609557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2850212) q[1];
sx q[1];
rz(-1.8076496) q[1];
sx q[1];
rz(-0.23394886) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4056021) q[3];
sx q[3];
rz(-1.372589) q[3];
sx q[3];
rz(-2.7310731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(-0.84154877) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(0.21249214) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1157163) q[0];
sx q[0];
rz(-1.8159144) q[0];
sx q[0];
rz(0.91761597) q[0];
rz(3.100527) q[2];
sx q[2];
rz(-1.2588871) q[2];
sx q[2];
rz(1.2817771) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4957628) q[1];
sx q[1];
rz(-2.1730575) q[1];
sx q[1];
rz(-2.5219003) q[1];
rz(0.28189567) q[3];
sx q[3];
rz(-2.3882139) q[3];
sx q[3];
rz(-0.58338141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96182573) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(-1.6869705) q[2];
rz(1.9466594) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-0.71804071) q[3];
sx q[3];
rz(-0.64648872) q[3];
sx q[3];
rz(1.7512376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
