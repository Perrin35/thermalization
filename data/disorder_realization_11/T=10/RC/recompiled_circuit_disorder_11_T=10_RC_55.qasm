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
rz(-0.89590961) q[0];
sx q[0];
rz(1.7370261) q[0];
rz(3.8430619) q[1];
sx q[1];
rz(3.7624533) q[1];
sx q[1];
rz(7.9384595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5777187) q[0];
sx q[0];
rz(-1.3868252) q[0];
sx q[0];
rz(1.4509393) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42278554) q[2];
sx q[2];
rz(-1.6362731) q[2];
sx q[2];
rz(-0.72725429) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4991551) q[1];
sx q[1];
rz(-0.40343522) q[1];
sx q[1];
rz(0.24654504) q[1];
rz(-pi) q[2];
rz(0.08333929) q[3];
sx q[3];
rz(-3.0623751) q[3];
sx q[3];
rz(3.0230429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.7479744) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(0.077985667) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(-0.35968131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4820837) q[0];
sx q[0];
rz(-2.2523899) q[0];
sx q[0];
rz(0.73007749) q[0];
x q[1];
rz(0.72007911) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(-0.0062696487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82275326) q[1];
sx q[1];
rz(-2.0944632) q[1];
sx q[1];
rz(-0.19284064) q[1];
rz(1.6167323) q[3];
sx q[3];
rz(-0.71456281) q[3];
sx q[3];
rz(3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(-0.94334156) q[2];
rz(1.2603166) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(-0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.48370353) q[1];
sx q[1];
rz(2.8288249) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(1.3586112) q[0];
rz(-pi) q[1];
rz(0.33212338) q[2];
sx q[2];
rz(-0.87288522) q[2];
sx q[2];
rz(-2.3534564) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2777729) q[1];
sx q[1];
rz(-1.4295271) q[1];
sx q[1];
rz(-1.743725) q[1];
rz(-2.7346482) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(2.0719761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(-1.336608) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(-0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(2.4688597) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8254291) q[0];
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
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87186253) q[1];
sx q[1];
rz(-1.6849663) q[1];
sx q[1];
rz(1.8037379) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1898515) q[3];
sx q[3];
rz(-1.674106) q[3];
sx q[3];
rz(0.88224525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(-2.380774) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(0.76604617) q[3];
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
rz(2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-2.0779579) q[0];
rz(1.4798374) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(3.1033049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3064561) q[0];
sx q[0];
rz(-1.4493754) q[0];
sx q[0];
rz(1.8655213) q[0];
rz(-pi) q[1];
rz(-2.4942336) q[2];
sx q[2];
rz(-0.72509662) q[2];
sx q[2];
rz(-0.65078306) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8748862) q[1];
sx q[1];
rz(-1.453207) q[1];
sx q[1];
rz(-1.7456732) q[1];
rz(-pi) q[2];
rz(-0.27023817) q[3];
sx q[3];
rz(-2.022937) q[3];
sx q[3];
rz(-2.7142392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(1.0920452) q[2];
rz(-1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096481) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(1.3076179) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(2.9506156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6742498) q[0];
sx q[0];
rz(-1.7102339) q[0];
sx q[0];
rz(1.3583202) q[0];
rz(-pi) q[1];
rz(1.9917166) q[2];
sx q[2];
rz(-0.72962609) q[2];
sx q[2];
rz(-0.52044496) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2872419) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(1.0276621) q[1];
rz(-pi) q[2];
rz(-2.1359547) q[3];
sx q[3];
rz(-1.2747545) q[3];
sx q[3];
rz(2.5067096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-0.43506518) q[2];
rz(-2.2502031) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(2.7213668) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(1.0302672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(-2.6269873) q[0];
rz(-3.0646965) q[2];
sx q[2];
rz(-1.7891857) q[2];
sx q[2];
rz(1.0181392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.269304) q[1];
sx q[1];
rz(-1.270712) q[1];
sx q[1];
rz(1.7606723) q[1];
x q[2];
rz(2.9514063) q[3];
sx q[3];
rz(-1.38387) q[3];
sx q[3];
rz(-0.72894064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(1.1918921) q[0];
rz(0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(1.0620767) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22411352) q[0];
sx q[0];
rz(-0.55372059) q[0];
sx q[0];
rz(-1.327716) q[0];
rz(-1.1658792) q[2];
sx q[2];
rz(-1.9654044) q[2];
sx q[2];
rz(1.6369866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3036365) q[1];
sx q[1];
rz(-2.084823) q[1];
sx q[1];
rz(2.282826) q[1];
x q[2];
rz(-3.1065337) q[3];
sx q[3];
rz(-2.2359214) q[3];
sx q[3];
rz(-2.866131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.186211) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(0.051076802) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99701571) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061714) q[0];
sx q[0];
rz(-1.346934) q[0];
sx q[0];
rz(-1.9497245) q[0];
x q[1];
rz(-0.37306771) q[2];
sx q[2];
rz(-0.39857769) q[2];
sx q[2];
rz(2.1052436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2850212) q[1];
sx q[1];
rz(-1.3339431) q[1];
sx q[1];
rz(-2.9076438) q[1];
x q[2];
rz(2.455515) q[3];
sx q[3];
rz(-0.25732532) q[3];
sx q[3];
rz(-1.1130594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7944305) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(-2.9629422) q[2];
rz(0.84154877) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(2.2264746) q[0];
rz(1.2471584) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(-2.9291005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8520078) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(1.1803124) q[0];
rz(1.4441522) q[2];
sx q[2];
rz(-2.8270792) q[2];
sx q[2];
rz(-1.7267137) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.64582981) q[1];
sx q[1];
rz(-0.96853515) q[1];
sx q[1];
rz(2.5219003) q[1];
rz(-pi) q[2];
rz(0.73331613) q[3];
sx q[3];
rz(-1.3793257) q[3];
sx q[3];
rz(1.9460033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1797669) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(1.6869705) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4029978) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(0.075642792) q[1];
sx q[1];
rz(-1.5374001) q[1];
sx q[1];
rz(-1.7061445) q[1];
rz(-3.1235789) q[2];
sx q[2];
rz(-1.6406888) q[2];
sx q[2];
rz(-2.3802118) q[2];
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
