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
rz(-1.6552742) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5777187) q[0];
sx q[0];
rz(-1.7547675) q[0];
sx q[0];
rz(1.6906533) q[0];
rz(2.7188071) q[2];
sx q[2];
rz(-1.6362731) q[2];
sx q[2];
rz(-0.72725429) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6424375) q[1];
sx q[1];
rz(-2.7381574) q[1];
sx q[1];
rz(-0.24654504) q[1];
rz(1.5641883) q[3];
sx q[3];
rz(-1.6497383) q[3];
sx q[3];
rz(-0.034949485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0554589) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(0.80111516) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-0.077985667) q[0];
rz(-0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(-2.7819113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29794824) q[0];
sx q[0];
rz(-2.1878562) q[0];
sx q[0];
rz(-2.258837) q[0];
rz(0.72007911) q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6911192) q[1];
sx q[1];
rz(-0.55492655) q[1];
sx q[1];
rz(1.2503442) q[1];
x q[2];
rz(-0.039814063) q[3];
sx q[3];
rz(-0.85714825) q[3];
sx q[3];
rz(3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(-2.1982511) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(-0.46762064) q[0];
rz(-2.6768661) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(0.31276774) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3019965) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(1.7829814) q[0];
x q[1];
rz(-0.84509387) q[2];
sx q[2];
rz(-1.8232864) q[2];
sx q[2];
rz(-1.0007728) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1130502) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(0.88009665) q[1];
rz(-2.8144366) q[3];
sx q[3];
rz(-1.4331685) q[3];
sx q[3];
rz(-0.88537346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(1.336608) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(-0.46521503) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49622789) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(-1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(0.67273295) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3161635) q[0];
sx q[0];
rz(-2.5220693) q[0];
sx q[0];
rz(-0.53014689) q[0];
rz(-0.54949923) q[2];
sx q[2];
rz(-1.3504343) q[2];
sx q[2];
rz(0.22827497) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.67191254) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(0.11731053) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.014971) q[3];
sx q[3];
rz(-0.95553482) q[3];
sx q[3];
rz(2.5263853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(-0.76081863) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(0.76604617) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-2.0779579) q[0];
rz(-1.4798374) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-0.038287727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83513658) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(1.2760713) q[0];
rz(-pi) q[1];
rz(1.080004) q[2];
sx q[2];
rz(-1.0133427) q[2];
sx q[2];
rz(1.7003531) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89025154) q[1];
sx q[1];
rz(-2.9311935) q[1];
sx q[1];
rz(2.167278) q[1];
rz(-pi) q[2];
rz(1.1039626) q[3];
sx q[3];
rz(-1.3282913) q[3];
sx q[3];
rz(-1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.15497196) q[2];
sx q[2];
rz(-1.203323) q[2];
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
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096481) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(-0.062782137) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(2.9506156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6742498) q[0];
sx q[0];
rz(-1.7102339) q[0];
sx q[0];
rz(-1.3583202) q[0];
x q[1];
rz(-1.149876) q[2];
sx q[2];
rz(-0.72962609) q[2];
sx q[2];
rz(-0.52044496) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2497017) q[1];
sx q[1];
rz(-0.54348031) q[1];
sx q[1];
rz(-1.6103423) q[1];
rz(1.0530924) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(-2.6368486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6993616) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(-1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681322) q[0];
sx q[0];
rz(-0.62481835) q[0];
sx q[0];
rz(-0.66972759) q[0];
rz(1.3517802) q[2];
sx q[2];
rz(-1.4957301) q[2];
sx q[2];
rz(0.53596562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29579269) q[1];
sx q[1];
rz(-2.7880221) q[1];
sx q[1];
rz(-2.5938864) q[1];
rz(-pi) q[2];
rz(-2.35595) q[3];
sx q[3];
rz(-2.8757186) q[3];
sx q[3];
rz(-3.0674988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(-0.31526652) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(-0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[1];
rz(2.7161712) q[2];
sx q[2];
rz(-1.1985949) q[2];
sx q[2];
rz(-0.22949716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.33151367) q[1];
sx q[1];
rz(-0.96558297) q[1];
sx q[1];
rz(2.5007351) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5261371) q[3];
sx q[3];
rz(-0.66590819) q[3];
sx q[3];
rz(0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(3.0905159) q[2];
rz(0.71930277) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.99701571) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(-0.84916806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061714) q[0];
sx q[0];
rz(-1.346934) q[0];
sx q[0];
rz(-1.9497245) q[0];
rz(2.7678713) q[2];
sx q[2];
rz(-1.7127275) q[2];
sx q[2];
rz(2.2609557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85657142) q[1];
sx q[1];
rz(-1.3339431) q[1];
sx q[1];
rz(-2.9076438) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9407223) q[3];
sx q[3];
rz(-1.4088653) q[3];
sx q[3];
rz(2.014132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(-0.84154877) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0586044) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3613113) q[0];
sx q[0];
rz(-2.2012976) q[0];
sx q[0];
rz(0.30514858) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.041065644) q[2];
sx q[2];
rz(-1.8827056) q[2];
sx q[2];
rz(1.8598156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6007235) q[1];
sx q[1];
rz(-2.0698554) q[1];
sx q[1];
rz(0.869511) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.859697) q[3];
sx q[3];
rz(-2.3882139) q[3];
sx q[3];
rz(2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.6869705) q[2];
rz(-1.1949332) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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