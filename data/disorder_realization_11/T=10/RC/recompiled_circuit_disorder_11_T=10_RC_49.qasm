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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5777187) q[0];
sx q[0];
rz(-1.7547675) q[0];
sx q[0];
rz(1.4509393) q[0];
rz(1.6425743) q[2];
sx q[2];
rz(-1.1489747) q[2];
sx q[2];
rz(2.3274802) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6424375) q[1];
sx q[1];
rz(-0.40343522) q[1];
sx q[1];
rz(-0.24654504) q[1];
rz(-3.0582534) q[3];
sx q[3];
rz(-0.079217521) q[3];
sx q[3];
rz(0.11854974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.7479744) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(2.7819113) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65950899) q[0];
sx q[0];
rz(-2.2523899) q[0];
sx q[0];
rz(-0.73007749) q[0];
rz(-2.9042269) q[2];
sx q[2];
rz(-2.4073905) q[2];
sx q[2];
rz(-1.7422301) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2962131) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(-1.0389894) q[1];
x q[2];
rz(1.5248604) q[3];
sx q[3];
rz(-0.71456281) q[3];
sx q[3];
rz(-3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.48370353) q[1];
sx q[1];
rz(-2.8288249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5178459) q[0];
sx q[0];
rz(-2.9150634) q[0];
sx q[0];
rz(-1.2073327) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84509387) q[2];
sx q[2];
rz(-1.8232864) q[2];
sx q[2];
rz(-1.0007728) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.410027) q[1];
sx q[1];
rz(-1.3996074) q[1];
sx q[1];
rz(0.14337916) q[1];
rz(-pi) q[2];
rz(-1.7160277) q[3];
sx q[3];
rz(-1.2468474) q[3];
sx q[3];
rz(0.63889972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(1.336608) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(0.46521503) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(0.67273295) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3161635) q[0];
sx q[0];
rz(-0.61952335) q[0];
sx q[0];
rz(0.53014689) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40517278) q[2];
sx q[2];
rz(-0.58779683) q[2];
sx q[2];
rz(1.6853465) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4696801) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(-3.0242821) q[1];
rz(-pi) q[2];
rz(-1.3939875) q[3];
sx q[3];
rz(-0.62649957) q[3];
sx q[3];
rz(-2.3092983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.68040401) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(-2.380774) q[2];
rz(2.2740254) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(-2.3755465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5701533) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-2.0779579) q[0];
rz(-1.6617552) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-3.1033049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69890755) q[0];
sx q[0];
rz(-1.8632871) q[0];
sx q[0];
rz(3.0147576) q[0];
rz(-pi) q[1];
rz(-0.64735909) q[2];
sx q[2];
rz(-2.416496) q[2];
sx q[2];
rz(-0.65078306) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.28336477) q[1];
sx q[1];
rz(-1.7444532) q[1];
sx q[1];
rz(-3.0221992) q[1];
rz(1.0682265) q[3];
sx q[3];
rz(-2.6196819) q[3];
sx q[3];
rz(-3.0038602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(2.0495474) q[2];
rz(-1.6070222) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.9654013) q[1];
sx q[1];
rz(0.19097701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.068114) q[0];
sx q[0];
rz(-1.3604135) q[0];
sx q[0];
rz(-2.99899) q[0];
rz(-pi) q[1];
rz(-2.2553308) q[2];
sx q[2];
rz(-1.8466511) q[2];
sx q[2];
rz(0.72826284) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29589614) q[1];
sx q[1];
rz(-1.0277883) q[1];
sx q[1];
rz(0.023881749) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0530924) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(2.6368486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6993616) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(-2.7065275) q[2];
rz(2.2502031) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1973535) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(-2.7213668) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(2.1113254) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(-2.6269873) q[0];
rz(0.076896197) q[2];
sx q[2];
rz(-1.352407) q[2];
sx q[2];
rz(2.1234535) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8458) q[1];
sx q[1];
rz(-0.35357057) q[1];
sx q[1];
rz(2.5938864) q[1];
x q[2];
rz(-1.7610735) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(-0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(1.0366084) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2675562) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-1.0620767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5545643) q[0];
sx q[0];
rz(-1.6977068) q[0];
sx q[0];
rz(1.0303322) q[0];
rz(-pi) q[1];
rz(-0.42542142) q[2];
sx q[2];
rz(-1.9429978) q[2];
sx q[2];
rz(0.22949716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3036365) q[1];
sx q[1];
rz(-2.084823) q[1];
sx q[1];
rz(-2.282826) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90537269) q[3];
sx q[3];
rz(-1.5432127) q[3];
sx q[3];
rz(-1.2736922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.186211) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(3.0905159) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1445769) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(2.2924246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061714) q[0];
sx q[0];
rz(-1.7946587) q[0];
sx q[0];
rz(-1.1918681) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37372132) q[2];
sx q[2];
rz(-1.4288651) q[2];
sx q[2];
rz(-2.2609557) q[2];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.455515) q[3];
sx q[3];
rz(-2.8842673) q[3];
sx q[3];
rz(1.1130594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(2.9629422) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-0.21249214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025876306) q[0];
sx q[0];
rz(-1.8159144) q[0];
sx q[0];
rz(-0.91761597) q[0];
rz(-pi) q[1];
rz(3.100527) q[2];
sx q[2];
rz(-1.8827056) q[2];
sx q[2];
rz(1.8598156) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6007235) q[1];
sx q[1];
rz(-2.0698554) q[1];
sx q[1];
rz(0.869511) q[1];
x q[2];
rz(0.28189567) q[3];
sx q[3];
rz(-2.3882139) q[3];
sx q[3];
rz(2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.96182573) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(1.4546222) q[2];
rz(1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.5008925) q[2];
sx q[2];
rz(-1.5528266) q[2];
sx q[2];
rz(2.3334353) q[2];
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