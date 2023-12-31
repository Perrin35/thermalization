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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.563874) q[0];
sx q[0];
rz(-1.7547675) q[0];
sx q[0];
rz(-1.4509393) q[0];
x q[1];
rz(-0.15847023) q[2];
sx q[2];
rz(-0.42752334) q[2];
sx q[2];
rz(2.1536364) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9857786) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(-0.39246651) q[1];
x q[2];
rz(3.0582534) q[3];
sx q[3];
rz(-3.0623751) q[3];
sx q[3];
rz(-3.0230429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0861337) q[2];
sx q[2];
rz(-1.7666585) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(3.063607) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(-2.7819113) q[1];
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
rz(2.9042269) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(1.3993625) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2962131) q[1];
sx q[1];
rz(-1.737533) q[1];
sx q[1];
rz(-1.0389894) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5248604) q[3];
sx q[3];
rz(-2.4270298) q[3];
sx q[3];
rz(-3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(-2.1982511) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-2.8288249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62374672) q[0];
sx q[0];
rz(-2.9150634) q[0];
sx q[0];
rz(-1.2073327) q[0];
rz(-1.9415641) q[2];
sx q[2];
rz(-2.3808378) q[2];
sx q[2];
rz(-0.2955546) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.028542472) q[1];
sx q[1];
rz(-2.9187435) q[1];
sx q[1];
rz(-0.88009665) q[1];
rz(-pi) q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(2.0719761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(-1.336608) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49622789) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(-1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(0.67273295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493895) q[0];
sx q[0];
rz(-2.0954847) q[0];
sx q[0];
rz(-1.9169109) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54949923) q[2];
sx q[2];
rz(-1.7911583) q[2];
sx q[2];
rz(-0.22827497) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67191254) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(0.11731053) q[1];
rz(1.3939875) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68040401) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(-0.76081863) q[2];
rz(2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(-0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.5714394) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(1.0636348) q[0];
rz(-1.6617552) q[1];
sx q[1];
rz(-2.1789443) q[1];
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
rz(-1.8655213) q[0];
rz(-pi) q[1];
rz(0.64735909) q[2];
sx q[2];
rz(-0.72509662) q[2];
sx q[2];
rz(2.4908096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2513411) q[1];
sx q[1];
rz(-2.9311935) q[1];
sx q[1];
rz(-0.9743147) q[1];
rz(-pi) q[2];
rz(1.1039626) q[3];
sx q[3];
rz(-1.3282913) q[3];
sx q[3];
rz(1.8777101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-2.0495474) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096481) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(-0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(0.19097701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4673429) q[0];
sx q[0];
rz(-1.7102339) q[0];
sx q[0];
rz(-1.3583202) q[0];
rz(-pi) q[1];
rz(2.7912748) q[2];
sx q[2];
rz(-2.2248474) q[2];
sx q[2];
rz(-2.0803114) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2872419) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(-2.1139305) q[1];
rz(-pi) q[2];
rz(-2.0885002) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(-2.6368486) q[3];
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
rz(0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.58266) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94423914) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(0.42022589) q[0];
rz(2.6603783) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(-1.0302672) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(-2.6269873) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9040362) q[2];
sx q[2];
rz(-2.9102647) q[2];
sx q[2];
rz(-1.3598134) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3863341) q[1];
sx q[1];
rz(-1.7520906) q[1];
sx q[1];
rz(2.8363486) q[1];
x q[2];
rz(-1.3805192) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(-1.0366084) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(0.96926564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(1.9497005) q[0];
rz(0.9115971) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(2.079516) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5545643) q[0];
sx q[0];
rz(-1.4438859) q[0];
sx q[0];
rz(-1.0303322) q[0];
rz(-pi) q[1];
rz(2.383963) q[2];
sx q[2];
rz(-2.5839351) q[2];
sx q[2];
rz(-0.66495313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.33151367) q[1];
sx q[1];
rz(-0.96558297) q[1];
sx q[1];
rz(-0.64085754) q[1];
x q[2];
rz(-1.6154556) q[3];
sx q[3];
rz(-0.66590819) q[3];
sx q[3];
rz(-0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.186211) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(0.051076802) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99701571) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(2.2924246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51061714) q[0];
sx q[0];
rz(-1.346934) q[0];
sx q[0];
rz(1.1918681) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
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
rz(1.4918684) q[1];
sx q[1];
rz(-0.33136156) q[1];
sx q[1];
rz(-2.3359873) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.455515) q[3];
sx q[3];
rz(-2.8842673) q[3];
sx q[3];
rz(-1.1130594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(0.21249214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2895848) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(-1.9612802) q[0];
x q[1];
rz(-1.6974405) q[2];
sx q[2];
rz(-0.31451348) q[2];
sx q[2];
rz(-1.414879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64582981) q[1];
sx q[1];
rz(-0.96853515) q[1];
sx q[1];
rz(2.5219003) q[1];
rz(2.859697) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(2.5582112) q[3];
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
rz(1.9466594) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(-0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.6407001) q[2];
sx q[2];
rz(-1.5887661) q[2];
sx q[2];
rz(-0.80815732) q[2];
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
