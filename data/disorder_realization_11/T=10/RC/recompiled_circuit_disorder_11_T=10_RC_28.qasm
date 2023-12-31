OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5869863) q[0];
sx q[0];
rz(-2.245683) q[0];
sx q[0];
rz(-1.7370261) q[0];
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(1.4863185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1264869) q[0];
sx q[0];
rz(-1.6886212) q[0];
sx q[0];
rz(-0.18527041) q[0];
x q[1];
rz(2.7188071) q[2];
sx q[2];
rz(-1.6362731) q[2];
sx q[2];
rz(-0.72725429) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3753429) q[1];
sx q[1];
rz(-1.9613593) q[1];
sx q[1];
rz(1.4669963) q[1];
rz(1.5774044) q[3];
sx q[3];
rz(-1.4918543) q[3];
sx q[3];
rz(-0.034949485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0861337) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(3.063607) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(2.7819113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247503) q[0];
sx q[0];
rz(-2.1152088) q[0];
sx q[0];
rz(-0.74290468) q[0];
rz(-0.23736575) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(-1.7422301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45047346) q[1];
sx q[1];
rz(-0.55492655) q[1];
sx q[1];
rz(1.2503442) q[1];
rz(-0.039814063) q[3];
sx q[3];
rz(-0.85714825) q[3];
sx q[3];
rz(-0.1114705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(-1.881276) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(-0.46762064) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(-1.3586112) q[0];
rz(2.2964988) q[2];
sx q[2];
rz(-1.3183062) q[2];
sx q[2];
rz(1.0007728) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1130502) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(0.88009665) q[1];
rz(-1.425565) q[3];
sx q[3];
rz(-1.8947453) q[3];
sx q[3];
rz(0.63889972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(0.30511937) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(-1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(-2.4688597) q[1];
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
rz(0.55158789) q[0];
x q[1];
rz(-1.3139309) q[2];
sx q[2];
rz(-1.0360403) q[2];
sx q[2];
rz(-1.2094487) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.67191254) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(-0.11731053) q[1];
rz(-pi) q[2];
rz(-0.9517412) q[3];
sx q[3];
rz(-1.4674867) q[3];
sx q[3];
rz(2.2593474) q[3];
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
rz(2.2740254) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5701533) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(-2.0779579) q[0];
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
rz(2.3064561) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(1.2760713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61530453) q[2];
sx q[2];
rz(-1.9822789) q[2];
sx q[2];
rz(2.736511) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8582279) q[1];
sx q[1];
rz(-1.3971395) q[1];
sx q[1];
rz(0.11939343) q[1];
x q[2];
rz(0.27023817) q[3];
sx q[3];
rz(-2.022937) q[3];
sx q[3];
rz(2.7142392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(-2.0495474) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-2.3585414) q[3];
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
rz(pi/2) q[0];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67577772) q[0];
sx q[0];
rz(-0.2535648) q[0];
sx q[0];
rz(2.158014) q[0];
x q[1];
rz(-1.9917166) q[2];
sx q[2];
rz(-2.4119666) q[2];
sx q[2];
rz(2.6211477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29589614) q[1];
sx q[1];
rz(-1.0277883) q[1];
sx q[1];
rz(3.1177109) q[1];
rz(0.34658587) q[3];
sx q[3];
rz(-2.1086018) q[3];
sx q[3];
rz(-1.1188521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6993616) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-2.7065275) q[2];
rz(2.2502031) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(-1.58266) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(2.7213668) q[0];
rz(-0.48121437) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(0.51460534) q[0];
x q[1];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.7891857) q[2];
sx q[2];
rz(2.1234535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29579269) q[1];
sx q[1];
rz(-0.35357057) q[1];
sx q[1];
rz(2.5938864) q[1];
rz(-pi) q[2];
rz(-1.3805192) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8740365) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.9497005) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-2.079516) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05957219) q[0];
sx q[0];
rz(-2.1064415) q[0];
sx q[0];
rz(0.14772149) q[0];
x q[1];
rz(2.383963) q[2];
sx q[2];
rz(-2.5839351) q[2];
sx q[2];
rz(2.4766395) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8909104) q[1];
sx q[1];
rz(-2.2905596) q[1];
sx q[1];
rz(0.85810424) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5261371) q[3];
sx q[3];
rz(-2.4756845) q[3];
sx q[3];
rz(-2.8093616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(0.71930277) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(0.64494079) q[0];
rz(2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(-0.84916806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51061714) q[0];
sx q[0];
rz(-1.7946587) q[0];
sx q[0];
rz(-1.9497245) q[0];
x q[1];
rz(0.37372132) q[2];
sx q[2];
rz(-1.7127275) q[2];
sx q[2];
rz(0.88063699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2850212) q[1];
sx q[1];
rz(-1.8076496) q[1];
sx q[1];
rz(-0.23394886) q[1];
rz(-pi) q[2];
rz(0.68607761) q[3];
sx q[3];
rz(-2.8842673) q[3];
sx q[3];
rz(2.0285332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-2.9629422) q[2];
rz(-2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8520078) q[0];
sx q[0];
rz(-2.4502909) q[0];
sx q[0];
rz(1.1803124) q[0];
rz(-pi) q[1];
rz(1.2586406) q[2];
sx q[2];
rz(-1.5317132) q[2];
sx q[2];
rz(0.27641073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5455268) q[1];
sx q[1];
rz(-0.83546987) q[1];
sx q[1];
rz(-2.2722785) q[1];
rz(-0.73331613) q[3];
sx q[3];
rz(-1.7622669) q[3];
sx q[3];
rz(1.9460033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1797669) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4029978) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
rz(-0.075642792) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(3.1235789) q[2];
sx q[2];
rz(-1.5009038) q[2];
sx q[2];
rz(0.76138087) q[2];
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
