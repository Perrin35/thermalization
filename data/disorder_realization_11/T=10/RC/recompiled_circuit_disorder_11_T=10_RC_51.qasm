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
rz(3.8430619) q[1];
sx q[1];
rz(3.7624533) q[1];
sx q[1];
rz(7.9384595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.563874) q[0];
sx q[0];
rz(-1.3868252) q[0];
sx q[0];
rz(-1.4509393) q[0];
x q[1];
rz(-2.7188071) q[2];
sx q[2];
rz(-1.6362731) q[2];
sx q[2];
rz(-2.4143384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6424375) q[1];
sx q[1];
rz(-2.7381574) q[1];
sx q[1];
rz(-2.8950476) q[1];
x q[2];
rz(-1.5641883) q[3];
sx q[3];
rz(-1.6497383) q[3];
sx q[3];
rz(0.034949485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(-1.3936183) q[2];
rz(-0.80111516) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(1.0629517) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(0.077985667) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(-0.35968131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168424) q[0];
sx q[0];
rz(-1.0263838) q[0];
sx q[0];
rz(2.398688) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23736575) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(1.7422301) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8453796) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(2.1026033) q[1];
rz(-pi) q[2];
rz(0.039814063) q[3];
sx q[3];
rz(-2.2844444) q[3];
sx q[3];
rz(-0.1114705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(1.2603166) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(-0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(-0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8899925) q[0];
sx q[0];
rz(-1.359299) q[0];
sx q[0];
rz(0.081758008) q[0];
rz(-1.2000285) q[2];
sx q[2];
rz(-0.76075483) q[2];
sx q[2];
rz(2.8460381) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86381972) q[1];
sx q[1];
rz(-1.4295271) q[1];
sx q[1];
rz(-1.743725) q[1];
rz(-pi) q[2];
rz(-0.4069445) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(-2.0719761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(1.336608) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(-1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(-2.4688597) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3161635) q[0];
sx q[0];
rz(-0.61952335) q[0];
sx q[0];
rz(-0.53014689) q[0];
x q[1];
rz(-1.8276617) q[2];
sx q[2];
rz(-2.1055524) q[2];
sx q[2];
rz(-1.2094487) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4696801) q[1];
sx q[1];
rz(-1.3393991) q[1];
sx q[1];
rz(-3.0242821) q[1];
x q[2];
rz(-2.1898515) q[3];
sx q[3];
rz(-1.4674867) q[3];
sx q[3];
rz(0.88224525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(0.86756724) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(2.3755465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5701533) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(1.6617552) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(3.1033049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83513658) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(1.8655213) q[0];
rz(-pi) q[1];
x q[1];
rz(1.080004) q[2];
sx q[2];
rz(-2.12825) q[2];
sx q[2];
rz(-1.7003531) q[2];
rz(-pi) q[3];
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
rz(-0.11939343) q[1];
rz(-pi) q[2];
rz(-0.27023817) q[3];
sx q[3];
rz(-1.1186557) q[3];
sx q[3];
rz(2.7142392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-2.0495474) q[2];
rz(-1.6070222) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-2.3585414) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(-0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(0.19097701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6742498) q[0];
sx q[0];
rz(-1.4313587) q[0];
sx q[0];
rz(-1.3583202) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2553308) q[2];
sx q[2];
rz(-1.2949416) q[2];
sx q[2];
rz(-2.4133298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8543507) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(1.0276621) q[1];
x q[2];
rz(2.7950068) q[3];
sx q[3];
rz(-1.0329909) q[3];
sx q[3];
rz(-1.1188521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-0.43506518) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.2599506) q[1];
sx q[1];
rz(2.1113254) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681322) q[0];
sx q[0];
rz(-2.5167743) q[0];
sx q[0];
rz(0.66972759) q[0];
x q[1];
rz(-1.7898125) q[2];
sx q[2];
rz(-1.4957301) q[2];
sx q[2];
rz(-2.605627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87228862) q[1];
sx q[1];
rz(-1.8708806) q[1];
sx q[1];
rz(1.7606723) q[1];
rz(-1.7610735) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(2.2639757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.39945012) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(1.9497005) q[0];
rz(2.2299956) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-1.0620767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0820205) q[0];
sx q[0];
rz(-1.0351512) q[0];
sx q[0];
rz(2.9938712) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9757134) q[2];
sx q[2];
rz(-1.9654044) q[2];
sx q[2];
rz(-1.6369866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.810079) q[1];
sx q[1];
rz(-2.1760097) q[1];
sx q[1];
rz(-2.5007351) q[1];
rz(-1.6154556) q[3];
sx q[3];
rz(-0.66590819) q[3];
sx q[3];
rz(-0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(-0.71930277) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(1.0572082) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99701571) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-2.4966519) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(0.84916806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6309755) q[0];
sx q[0];
rz(-1.7946587) q[0];
sx q[0];
rz(1.9497245) q[0];
rz(2.7678713) q[2];
sx q[2];
rz(-1.7127275) q[2];
sx q[2];
rz(-0.88063699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6497242) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(0.80560537) q[1];
x q[2];
rz(-0.20087033) q[3];
sx q[3];
rz(-1.4088653) q[3];
sx q[3];
rz(1.1274606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7944305) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(0.84154877) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(0.21249214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2895848) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(1.9612802) q[0];
x q[1];
rz(0.041065644) q[2];
sx q[2];
rz(-1.8827056) q[2];
sx q[2];
rz(1.2817771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5960658) q[1];
sx q[1];
rz(-2.3061228) q[1];
sx q[1];
rz(2.2722785) q[1];
rz(-pi) q[2];
rz(-0.73331613) q[3];
sx q[3];
rz(-1.7622669) q[3];
sx q[3];
rz(-1.1955893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.4546222) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
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
rz(1.8226345) q[2];
sx q[2];
rz(-0.072172879) q[2];
sx q[2];
rz(-2.1277438) q[2];
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
