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
rz(-2.4401234) q[1];
sx q[1];
rz(-2.520732) q[1];
sx q[1];
rz(-1.4863185) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0151057) q[0];
sx q[0];
rz(-1.4529714) q[0];
sx q[0];
rz(-0.18527041) q[0];
rz(-pi) q[1];
rz(1.4990184) q[2];
sx q[2];
rz(-1.992618) q[2];
sx q[2];
rz(2.3274802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15581407) q[1];
sx q[1];
rz(-1.6667546) q[1];
sx q[1];
rz(2.7491261) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0626489) q[3];
sx q[3];
rz(-1.5642089) q[3];
sx q[3];
rz(-1.5353257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79008094) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-3.063607) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(0.35968131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247503) q[0];
sx q[0];
rz(-2.1152088) q[0];
sx q[0];
rz(2.398688) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9042269) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(1.3993625) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2962131) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(2.1026033) q[1];
rz(-pi) q[2];
rz(3.1017786) q[3];
sx q[3];
rz(-0.85714825) q[3];
sx q[3];
rz(-0.1114705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(-1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-0.31276774) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8395961) q[0];
sx q[0];
rz(-1.6507286) q[0];
sx q[0];
rz(1.3586112) q[0];
rz(2.2964988) q[2];
sx q[2];
rz(-1.8232864) q[2];
sx q[2];
rz(-1.0007728) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.410027) q[1];
sx q[1];
rz(-1.7419852) q[1];
sx q[1];
rz(2.9982135) q[1];
x q[2];
rz(-1.425565) q[3];
sx q[3];
rz(-1.2468474) q[3];
sx q[3];
rz(2.5026929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(-1.336608) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(0.46521503) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(1.4656981) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(-2.4688597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69220316) q[0];
sx q[0];
rz(-1.0461079) q[0];
sx q[0];
rz(-1.9169109) q[0];
rz(0.54949923) q[2];
sx q[2];
rz(-1.7911583) q[2];
sx q[2];
rz(0.22827497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9949918) q[1];
sx q[1];
rz(-0.25895893) q[1];
sx q[1];
rz(2.0318356) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3939875) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(-0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(-0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(-1.4798374) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(3.1033049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4426851) q[0];
sx q[0];
rz(-1.2783056) q[0];
sx q[0];
rz(0.12683503) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64735909) q[2];
sx q[2];
rz(-0.72509662) q[2];
sx q[2];
rz(2.4908096) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28336477) q[1];
sx q[1];
rz(-1.3971395) q[1];
sx q[1];
rz(-0.11939343) q[1];
rz(-0.27023817) q[3];
sx q[3];
rz(-1.1186557) q[3];
sx q[3];
rz(2.7142392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(-1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(2.3585414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606278) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(-0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(0.19097701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.068114) q[0];
sx q[0];
rz(-1.3604135) q[0];
sx q[0];
rz(-0.14260261) q[0];
x q[1];
rz(0.35031788) q[2];
sx q[2];
rz(-0.91674524) q[2];
sx q[2];
rz(1.0612812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8543507) q[1];
sx q[1];
rz(-1.5503503) q[1];
sx q[1];
rz(-2.1139305) q[1];
rz(0.34658587) q[3];
sx q[3];
rz(-1.0329909) q[3];
sx q[3];
rz(1.1188521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.442231) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(-0.43506518) q[2];
rz(-2.2502031) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.94423914) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(2.7213668) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(2.1113254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1680465) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(-0.51460534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3517802) q[2];
sx q[2];
rz(-1.6458626) q[2];
sx q[2];
rz(0.53596562) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.269304) q[1];
sx q[1];
rz(-1.8708806) q[1];
sx q[1];
rz(1.7606723) q[1];
rz(-pi) q[2];
rz(-0.19018634) q[3];
sx q[3];
rz(-1.7577226) q[3];
sx q[3];
rz(0.72894064) q[3];
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
rz(-2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-2.079516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5545643) q[0];
sx q[0];
rz(-1.4438859) q[0];
sx q[0];
rz(1.0303322) q[0];
x q[1];
rz(-1.1658792) q[2];
sx q[2];
rz(-1.9654044) q[2];
sx q[2];
rz(1.6369866) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2506822) q[1];
sx q[1];
rz(-2.2905596) q[1];
sx q[1];
rz(0.85810424) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1065337) q[3];
sx q[3];
rz(-0.9056712) q[3];
sx q[3];
rz(2.866131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1445769) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-2.4966519) q[0];
rz(1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5727974) q[0];
sx q[0];
rz(-0.43734567) q[0];
sx q[0];
rz(-1.019078) q[0];
rz(-pi) q[1];
rz(2.7685249) q[2];
sx q[2];
rz(-0.39857769) q[2];
sx q[2];
rz(-1.0363491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6497242) q[1];
sx q[1];
rz(-0.33136156) q[1];
sx q[1];
rz(2.3359873) q[1];
x q[2];
rz(-2.9407223) q[3];
sx q[3];
rz(-1.7327274) q[3];
sx q[3];
rz(-2.014132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(0.84154877) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(-0.21249214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2895848) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(1.9612802) q[0];
x q[1];
rz(-1.882952) q[2];
sx q[2];
rz(-1.6098795) q[2];
sx q[2];
rz(2.8651819) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64582981) q[1];
sx q[1];
rz(-0.96853515) q[1];
sx q[1];
rz(-0.61969238) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3155762) q[3];
sx q[3];
rz(-2.2877684) q[3];
sx q[3];
rz(2.9361801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(-1.6869705) q[2];
rz(1.9466594) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(-0.52836829) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(1.6407001) q[2];
sx q[2];
rz(-1.5528266) q[2];
sx q[2];
rz(2.3334353) q[2];
rz(2.6247737) q[3];
sx q[3];
rz(-1.1632945) q[3];
sx q[3];
rz(-2.3522285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];