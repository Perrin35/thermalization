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
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(1.4863185) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0151057) q[0];
sx q[0];
rz(-1.6886212) q[0];
sx q[0];
rz(-0.18527041) q[0];
x q[1];
rz(-1.6425743) q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-2.9857786) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(2.7491261) q[1];
rz(0.078943723) q[3];
sx q[3];
rz(-1.5773838) q[3];
sx q[3];
rz(-1.6062669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0554589) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(-2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(-1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-3.063607) q[0];
rz(-0.69411913) q[1];
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
x q[1];
rz(-2.4215135) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(-0.0062696487) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.82275326) q[1];
sx q[1];
rz(-1.0471294) q[1];
sx q[1];
rz(2.948752) q[1];
rz(-pi) q[2];
rz(0.85675591) q[3];
sx q[3];
rz(-1.6008915) q[3];
sx q[3];
rz(1.4853958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(-0.94334156) q[2];
rz(-1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-2.4356306) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.869732) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(-2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-0.31276774) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(1.7829814) q[0];
x q[1];
rz(-2.8094693) q[2];
sx q[2];
rz(-0.87288522) q[2];
sx q[2];
rz(-2.3534564) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.410027) q[1];
sx q[1];
rz(-1.7419852) q[1];
sx q[1];
rz(-2.9982135) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32715601) q[3];
sx q[3];
rz(-1.7084242) q[3];
sx q[3];
rz(-0.88537346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9096845) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.6453648) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(2.4688597) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4417277) q[0];
sx q[0];
rz(-1.2727951) q[0];
sx q[0];
rz(2.5900048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54949923) q[2];
sx q[2];
rz(-1.3504343) q[2];
sx q[2];
rz(-0.22827497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9949918) q[1];
sx q[1];
rz(-0.25895893) q[1];
sx q[1];
rz(-2.0318356) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1898515) q[3];
sx q[3];
rz(-1.674106) q[3];
sx q[3];
rz(2.2593474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68040401) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-0.86756724) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(2.0779579) q[0];
rz(-1.6617552) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(-0.038287727) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4426851) q[0];
sx q[0];
rz(-1.2783056) q[0];
sx q[0];
rz(0.12683503) q[0];
rz(-2.5262881) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(0.40508168) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8582279) q[1];
sx q[1];
rz(-1.7444532) q[1];
sx q[1];
rz(3.0221992) q[1];
rz(-pi) q[2];
rz(-2.0733662) q[3];
sx q[3];
rz(-2.6196819) q[3];
sx q[3];
rz(0.13773242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
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
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(-2.9506156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.068114) q[0];
sx q[0];
rz(-1.7811791) q[0];
sx q[0];
rz(0.14260261) q[0];
rz(-pi) q[1];
rz(1.9917166) q[2];
sx q[2];
rz(-2.4119666) q[2];
sx q[2];
rz(-2.6211477) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8543507) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(-2.1139305) q[1];
rz(-pi) q[2];
rz(2.0885002) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(-0.5047441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6993616) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(2.2502031) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(-2.7213668) q[0];
rz(-0.48121437) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-2.1113254) q[1];
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
rz(-pi) q[1];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.7891857) q[2];
sx q[2];
rz(-1.0181392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.269304) q[1];
sx q[1];
rz(-1.8708806) q[1];
sx q[1];
rz(1.7606723) q[1];
rz(-pi) q[2];
rz(1.3805192) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(-0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(-2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(-1.0620767) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05957219) q[0];
sx q[0];
rz(-2.1064415) q[0];
sx q[0];
rz(0.14772149) q[0];
rz(-pi) q[1];
rz(-0.42542142) q[2];
sx q[2];
rz(-1.9429978) q[2];
sx q[2];
rz(0.22949716) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2506822) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(-0.85810424) q[1];
x q[2];
rz(-3.1065337) q[3];
sx q[3];
rz(-0.9056712) q[3];
sx q[3];
rz(-0.2754617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(3.0905159) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1445769) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-2.2924246) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97201559) q[0];
sx q[0];
rz(-1.2017842) q[0];
sx q[0];
rz(-0.24032648) q[0];
rz(-pi) q[1];
rz(-0.37306771) q[2];
sx q[2];
rz(-2.743015) q[2];
sx q[2];
rz(1.0363491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2850212) q[1];
sx q[1];
rz(-1.3339431) q[1];
sx q[1];
rz(-0.23394886) q[1];
rz(1.7359906) q[3];
sx q[3];
rz(-1.372589) q[3];
sx q[3];
rz(2.7310731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7944305) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(-2.9629422) q[2];
rz(-2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-0.49079045) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(-1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1157163) q[0];
sx q[0];
rz(-1.3256782) q[0];
sx q[0];
rz(-2.2239767) q[0];
rz(-pi) q[1];
rz(-1.4441522) q[2];
sx q[2];
rz(-2.8270792) q[2];
sx q[2];
rz(-1.414879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5960658) q[1];
sx q[1];
rz(-0.83546987) q[1];
sx q[1];
rz(-2.2722785) q[1];
rz(0.28189567) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(0.58338141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1797669) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.6869705) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385948) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(0.075642792) q[1];
sx q[1];
rz(-1.5374001) q[1];
sx q[1];
rz(-1.7061445) q[1];
rz(0.018013714) q[2];
sx q[2];
rz(-1.6406888) q[2];
sx q[2];
rz(-2.3802118) q[2];
rz(0.71804071) q[3];
sx q[3];
rz(-2.4951039) q[3];
sx q[3];
rz(-1.390355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
