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
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(1.4863185) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1264869) q[0];
sx q[0];
rz(-1.4529714) q[0];
sx q[0];
rz(-0.18527041) q[0];
x q[1];
rz(0.15847023) q[2];
sx q[2];
rz(-0.42752334) q[2];
sx q[2];
rz(-2.1536364) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6424375) q[1];
sx q[1];
rz(-0.40343522) q[1];
sx q[1];
rz(-0.24654504) q[1];
rz(0.08333929) q[3];
sx q[3];
rz(-0.079217521) q[3];
sx q[3];
rz(-3.0230429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0861337) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(-1.3936183) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(-1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79008094) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-0.077985667) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(-2.7819113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4820837) q[0];
sx q[0];
rz(-2.2523899) q[0];
sx q[0];
rz(0.73007749) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4215135) q[2];
sx q[2];
rz(-1.7290001) q[2];
sx q[2];
rz(3.135323) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2962131) q[1];
sx q[1];
rz(-1.4040596) q[1];
sx q[1];
rz(-1.0389894) q[1];
rz(-pi) q[2];
x q[2];
rz(0.039814063) q[3];
sx q[3];
rz(-0.85714825) q[3];
sx q[3];
rz(-3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89547196) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(-1.2603166) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(-2.673972) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(-1.7829814) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84509387) q[2];
sx q[2];
rz(-1.3183062) q[2];
sx q[2];
rz(-2.1408199) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1130502) q[1];
sx q[1];
rz(-2.9187435) q[1];
sx q[1];
rz(0.88009665) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7160277) q[3];
sx q[3];
rz(-1.8947453) q[3];
sx q[3];
rz(-0.63889972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(-0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-0.67273295) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4417277) q[0];
sx q[0];
rz(-1.2727951) q[0];
sx q[0];
rz(-2.5900048) q[0];
x q[1];
rz(1.3139309) q[2];
sx q[2];
rz(-2.1055524) q[2];
sx q[2];
rz(1.9321439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9949918) q[1];
sx q[1];
rz(-2.8826337) q[1];
sx q[1];
rz(-2.0318356) q[1];
rz(1.7476051) q[3];
sx q[3];
rz(-0.62649957) q[3];
sx q[3];
rz(-2.3092983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.68040401) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(-2.380774) q[2];
rz(-0.86756724) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(-0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5701533) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(2.0779579) q[0];
rz(-1.4798374) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-0.038287727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69890755) q[0];
sx q[0];
rz(-1.8632871) q[0];
sx q[0];
rz(-0.12683503) q[0];
rz(-pi) q[1];
x q[1];
rz(1.080004) q[2];
sx q[2];
rz(-1.0133427) q[2];
sx q[2];
rz(-1.4412396) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28336477) q[1];
sx q[1];
rz(-1.7444532) q[1];
sx q[1];
rz(-3.0221992) q[1];
rz(-1.1039626) q[3];
sx q[3];
rz(-1.3282913) q[3];
sx q[3];
rz(1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(1.0920452) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(2.9506156) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.068114) q[0];
sx q[0];
rz(-1.7811791) q[0];
sx q[0];
rz(-2.99899) q[0];
rz(-pi) q[1];
rz(-2.2553308) q[2];
sx q[2];
rz(-1.8466511) q[2];
sx q[2];
rz(0.72826284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.891891) q[1];
sx q[1];
rz(-0.54348031) q[1];
sx q[1];
rz(-1.5312503) q[1];
rz(2.1359547) q[3];
sx q[3];
rz(-1.8668381) q[3];
sx q[3];
rz(-0.63488301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6993616) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(-2.2502031) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94423914) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(2.7213668) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(1.0302672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7468443) q[0];
sx q[0];
rz(-2.0472102) q[0];
sx q[0];
rz(1.9917411) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7898125) q[2];
sx q[2];
rz(-1.4957301) q[2];
sx q[2];
rz(-2.605627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3863341) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(2.8363486) q[1];
x q[2];
rz(-2.9514063) q[3];
sx q[3];
rz(-1.38387) q[3];
sx q[3];
rz(-2.412652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(-2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(2.172327) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2675562) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-2.2299956) q[1];
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
rz(0.05957219) q[0];
sx q[0];
rz(-2.1064415) q[0];
sx q[0];
rz(2.9938712) q[0];
rz(-pi) q[1];
rz(0.42542142) q[2];
sx q[2];
rz(-1.9429978) q[2];
sx q[2];
rz(2.9120955) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.810079) q[1];
sx q[1];
rz(-2.1760097) q[1];
sx q[1];
rz(0.64085754) q[1];
rz(-pi) q[2];
rz(-3.1065337) q[3];
sx q[3];
rz(-0.9056712) q[3];
sx q[3];
rz(2.866131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(3.0905159) q[2];
rz(-0.71930277) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(2.0843845) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99701571) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(0.64494079) q[0];
rz(2.0058517) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-2.2924246) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5727974) q[0];
sx q[0];
rz(-2.704247) q[0];
sx q[0];
rz(2.1225147) q[0];
rz(-pi) q[1];
rz(2.7678713) q[2];
sx q[2];
rz(-1.7127275) q[2];
sx q[2];
rz(-0.88063699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6497242) q[1];
sx q[1];
rz(-0.33136156) q[1];
sx q[1];
rz(-0.80560537) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20087033) q[3];
sx q[3];
rz(-1.4088653) q[3];
sx q[3];
rz(2.014132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(-1.2471584) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(-0.21249214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025876306) q[0];
sx q[0];
rz(-1.3256782) q[0];
sx q[0];
rz(-0.91761597) q[0];
rz(1.6974405) q[2];
sx q[2];
rz(-2.8270792) q[2];
sx q[2];
rz(1.7267137) q[2];
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
rz(-2.0698554) q[1];
sx q[1];
rz(0.869511) q[1];
rz(-1.8260164) q[3];
sx q[3];
rz(-0.85382429) q[3];
sx q[3];
rz(0.20541257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1797669) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(-1.6869705) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(-2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4029978) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(3.0659499) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(-3.1235789) q[2];
sx q[2];
rz(-1.6406888) q[2];
sx q[2];
rz(-2.3802118) q[2];
rz(-2.4235519) q[3];
sx q[3];
rz(-2.4951039) q[3];
sx q[3];
rz(-1.390355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
