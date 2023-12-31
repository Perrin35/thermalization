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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5777187) q[0];
sx q[0];
rz(-1.3868252) q[0];
sx q[0];
rz(-1.4509393) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7188071) q[2];
sx q[2];
rz(-1.5053195) q[2];
sx q[2];
rz(0.72725429) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6424375) q[1];
sx q[1];
rz(-0.40343522) q[1];
sx q[1];
rz(-0.24654504) q[1];
x q[2];
rz(0.078943723) q[3];
sx q[3];
rz(-1.5773838) q[3];
sx q[3];
rz(-1.6062669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.7479744) q[2];
rz(0.80111516) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(0.077985667) q[0];
rz(-2.4474735) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(2.7819113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.2799727) q[2];
sx q[2];
rz(1.7143955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8453796) q[1];
sx q[1];
rz(-1.4040596) q[1];
sx q[1];
rz(-1.0389894) q[1];
rz(-0.85675591) q[3];
sx q[3];
rz(-1.6008915) q[3];
sx q[3];
rz(1.6561968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(-2.673972) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(2.8288249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(-1.3586112) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9415641) q[2];
sx q[2];
rz(-2.3808378) q[2];
sx q[2];
rz(0.2955546) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.410027) q[1];
sx q[1];
rz(-1.3996074) q[1];
sx q[1];
rz(2.9982135) q[1];
rz(0.4069445) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(-1.0696166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(-1.336608) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(-2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(0.67273295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3161635) q[0];
sx q[0];
rz(-0.61952335) q[0];
sx q[0];
rz(-0.53014689) q[0];
x q[1];
rz(-2.5920934) q[2];
sx q[2];
rz(-1.3504343) q[2];
sx q[2];
rz(2.9133177) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87186253) q[1];
sx q[1];
rz(-1.6849663) q[1];
sx q[1];
rz(-1.8037379) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3939875) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(2.3092983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(-2.380774) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(-0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69890755) q[0];
sx q[0];
rz(-1.8632871) q[0];
sx q[0];
rz(3.0147576) q[0];
rz(-0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(2.736511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8582279) q[1];
sx q[1];
rz(-1.7444532) q[1];
sx q[1];
rz(0.11939343) q[1];
rz(-pi) q[2];
rz(-2.8713545) q[3];
sx q[3];
rz(-1.1186557) q[3];
sx q[3];
rz(0.42735344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(2.0495474) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(-2.3585414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(1.8339748) q[0];
rz(3.0788105) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(2.9506156) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4673429) q[0];
sx q[0];
rz(-1.4313587) q[0];
sx q[0];
rz(-1.7832725) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9917166) q[2];
sx q[2];
rz(-2.4119666) q[2];
sx q[2];
rz(-2.6211477) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2872419) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(1.0276621) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0530924) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(-2.6368486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(0.89138952) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(-0.42022589) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(1.0302672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9734605) q[0];
sx q[0];
rz(-2.5167743) q[0];
sx q[0];
rz(-0.66972759) q[0];
x q[1];
rz(-1.3517802) q[2];
sx q[2];
rz(-1.6458626) q[2];
sx q[2];
rz(0.53596562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29579269) q[1];
sx q[1];
rz(-0.35357057) q[1];
sx q[1];
rz(-0.54770623) q[1];
rz(2.35595) q[3];
sx q[3];
rz(-0.26587405) q[3];
sx q[3];
rz(-3.0674988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(-0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.079516) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9174791) q[0];
sx q[0];
rz(-2.5878721) q[0];
sx q[0];
rz(1.327716) q[0];
rz(-2.7161712) q[2];
sx q[2];
rz(-1.9429978) q[2];
sx q[2];
rz(-0.22949716) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.810079) q[1];
sx q[1];
rz(-0.96558297) q[1];
sx q[1];
rz(-0.64085754) q[1];
x q[2];
rz(-1.6154556) q[3];
sx q[3];
rz(-2.4756845) q[3];
sx q[3];
rz(-2.8093616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(0.71930277) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(-1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(0.84916806) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51061714) q[0];
sx q[0];
rz(-1.346934) q[0];
sx q[0];
rz(-1.1918681) q[0];
rz(-pi) q[1];
rz(-1.7230942) q[2];
sx q[2];
rz(-1.9405775) q[2];
sx q[2];
rz(-0.6347444) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.9076438) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20087033) q[3];
sx q[3];
rz(-1.7327274) q[3];
sx q[3];
rz(1.1274606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(2.9629422) q[2];
rz(0.84154877) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(2.2264746) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(2.9291005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3613113) q[0];
sx q[0];
rz(-0.94029501) q[0];
sx q[0];
rz(-0.30514858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.882952) q[2];
sx q[2];
rz(-1.6098795) q[2];
sx q[2];
rz(-0.27641073) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4957628) q[1];
sx q[1];
rz(-0.96853515) q[1];
sx q[1];
rz(-0.61969238) q[1];
rz(-pi) q[2];
rz(1.3155762) q[3];
sx q[3];
rz(-0.85382429) q[3];
sx q[3];
rz(-2.9361801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-0.075642792) q[1];
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
