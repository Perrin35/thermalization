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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0151057) q[0];
sx q[0];
rz(-1.6886212) q[0];
sx q[0];
rz(-0.18527041) q[0];
x q[1];
rz(-2.7188071) q[2];
sx q[2];
rz(-1.6362731) q[2];
sx q[2];
rz(0.72725429) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15581407) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(-2.7491261) q[1];
rz(-pi) q[2];
rz(-3.0582534) q[3];
sx q[3];
rz(-3.0623751) q[3];
sx q[3];
rz(3.0230429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(0.80111516) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79008094) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-0.077985667) q[0];
rz(-0.69411913) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(-0.35968131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8436444) q[0];
sx q[0];
rz(-2.1878562) q[0];
sx q[0];
rz(-0.88275568) q[0];
rz(0.72007911) q[2];
sx q[2];
rz(-1.7290001) q[2];
sx q[2];
rz(0.0062696487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8453796) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(-1.0389894) q[1];
rz(-2.2848367) q[3];
sx q[3];
rz(-1.5407011) q[3];
sx q[3];
rz(1.6561968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27186069) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-0.31276774) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3019965) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(1.7829814) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8094693) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(-0.78813625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.410027) q[1];
sx q[1];
rz(-1.3996074) q[1];
sx q[1];
rz(0.14337916) q[1];
rz(-2.7346482) q[3];
sx q[3];
rz(-0.35396468) q[3];
sx q[3];
rz(-2.0719761) q[3];
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
rz(-1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(-0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(2.4688597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69986491) q[0];
sx q[0];
rz(-1.8687975) q[0];
sx q[0];
rz(0.55158789) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8276617) q[2];
sx q[2];
rz(-1.0360403) q[2];
sx q[2];
rz(-1.2094487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87186253) q[1];
sx q[1];
rz(-1.4566263) q[1];
sx q[1];
rz(-1.3378548) q[1];
rz(-3.014971) q[3];
sx q[3];
rz(-0.95553482) q[3];
sx q[3];
rz(-0.61520731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4611886) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(-0.76081863) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(-0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(1.4798374) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(3.1033049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69890755) q[0];
sx q[0];
rz(-1.8632871) q[0];
sx q[0];
rz(3.0147576) q[0];
rz(-pi) q[1];
rz(-0.61530453) q[2];
sx q[2];
rz(-1.9822789) q[2];
sx q[2];
rz(-2.736511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2667065) q[1];
sx q[1];
rz(-1.6883856) q[1];
sx q[1];
rz(1.7456732) q[1];
x q[2];
rz(-2.0376301) q[3];
sx q[3];
rz(-1.3282913) q[3];
sx q[3];
rz(1.8777101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096481) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(1.8339748) q[0];
rz(-3.0788105) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(-0.19097701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4658149) q[0];
sx q[0];
rz(-2.8880279) q[0];
sx q[0];
rz(2.158014) q[0];
rz(2.2553308) q[2];
sx q[2];
rz(-1.2949416) q[2];
sx q[2];
rz(0.72826284) q[2];
rz(-pi) q[3];
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
x q[2];
rz(-2.7950068) q[3];
sx q[3];
rz(-2.1086018) q[3];
sx q[3];
rz(2.0227405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6993616) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-0.43506518) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(-0.88164202) q[1];
sx q[1];
rz(1.0302672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1681322) q[0];
sx q[0];
rz(-0.62481835) q[0];
sx q[0];
rz(-2.4718651) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7898125) q[2];
sx q[2];
rz(-1.4957301) q[2];
sx q[2];
rz(2.605627) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.75525857) q[1];
sx q[1];
rz(-1.7520906) q[1];
sx q[1];
rz(-0.30524409) q[1];
rz(-0.78564268) q[3];
sx q[3];
rz(-0.26587405) q[3];
sx q[3];
rz(0.074093821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(-1.0366084) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(-0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.8740365) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.9497005) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(-2.079516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05957219) q[0];
sx q[0];
rz(-2.1064415) q[0];
sx q[0];
rz(-2.9938712) q[0];
rz(2.7161712) q[2];
sx q[2];
rz(-1.1985949) q[2];
sx q[2];
rz(2.9120955) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2506822) q[1];
sx q[1];
rz(-2.2905596) q[1];
sx q[1];
rz(0.85810424) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5261371) q[3];
sx q[3];
rz(-0.66590819) q[3];
sx q[3];
rz(-0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(-1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(0.64494079) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(-2.2924246) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97201559) q[0];
sx q[0];
rz(-1.9398085) q[0];
sx q[0];
rz(-0.24032648) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(-0.6347444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2850212) q[1];
sx q[1];
rz(-1.8076496) q[1];
sx q[1];
rz(-2.9076438) q[1];
rz(-pi) q[2];
rz(-1.4056021) q[3];
sx q[3];
rz(-1.7690036) q[3];
sx q[3];
rz(0.41051958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(-0.84154877) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(2.2264746) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(-0.21249214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025876306) q[0];
sx q[0];
rz(-1.8159144) q[0];
sx q[0];
rz(-2.2239767) q[0];
rz(-pi) q[1];
rz(-1.882952) q[2];
sx q[2];
rz(-1.5317132) q[2];
sx q[2];
rz(-2.8651819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4957628) q[1];
sx q[1];
rz(-0.96853515) q[1];
sx q[1];
rz(-0.61969238) q[1];
x q[2];
rz(-1.8260164) q[3];
sx q[3];
rz(-0.85382429) q[3];
sx q[3];
rz(-2.9361801) q[3];
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
rz(-1.6869705) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385948) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
rz(3.0659499) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(3.1235789) q[2];
sx q[2];
rz(-1.5009038) q[2];
sx q[2];
rz(0.76138087) q[2];
rz(2.0316488) q[3];
sx q[3];
rz(-1.0999332) q[3];
sx q[3];
rz(-0.55988452) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];