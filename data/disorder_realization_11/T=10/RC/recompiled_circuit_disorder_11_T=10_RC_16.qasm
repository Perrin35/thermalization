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
rz(2.1264869) q[0];
sx q[0];
rz(-1.4529714) q[0];
sx q[0];
rz(-0.18527041) q[0];
rz(-pi) q[1];
rz(-1.4990184) q[2];
sx q[2];
rz(-1.992618) q[2];
sx q[2];
rz(0.81411241) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6424375) q[1];
sx q[1];
rz(-0.40343522) q[1];
sx q[1];
rz(0.24654504) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5774044) q[3];
sx q[3];
rz(-1.6497383) q[3];
sx q[3];
rz(3.1066432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(-2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(-1.0629517) q[3];
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
rz(-pi) q[3];
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
rz(0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-0.077985667) q[0];
rz(-0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(0.35968131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4247503) q[0];
sx q[0];
rz(-2.1152088) q[0];
sx q[0];
rz(0.74290468) q[0];
x q[1];
rz(-2.4215135) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(3.135323) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3188394) q[1];
sx q[1];
rz(-1.0471294) q[1];
sx q[1];
rz(2.948752) q[1];
rz(-pi) q[2];
rz(-1.6167323) q[3];
sx q[3];
rz(-0.71456281) q[3];
sx q[3];
rz(0.050690953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(2.673972) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8899925) q[0];
sx q[0];
rz(-1.359299) q[0];
sx q[0];
rz(3.0598346) q[0];
x q[1];
rz(-0.33212338) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(0.78813625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.028542472) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(0.88009665) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(-1.0696166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(-2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(-0.67273295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493895) q[0];
sx q[0];
rz(-2.0954847) q[0];
sx q[0];
rz(-1.9169109) q[0];
rz(0.40517278) q[2];
sx q[2];
rz(-2.5537958) q[2];
sx q[2];
rz(-1.6853465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9949918) q[1];
sx q[1];
rz(-0.25895893) q[1];
sx q[1];
rz(-2.0318356) q[1];
rz(-0.12662162) q[3];
sx q[3];
rz(-2.1860578) q[3];
sx q[3];
rz(2.5263853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4611886) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(2.3755465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0262888) q[0];
sx q[0];
rz(-2.8235108) q[0];
sx q[0];
rz(1.9684857) q[0];
rz(-pi) q[1];
rz(-0.64735909) q[2];
sx q[2];
rz(-0.72509662) q[2];
sx q[2];
rz(0.65078306) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2513411) q[1];
sx q[1];
rz(-2.9311935) q[1];
sx q[1];
rz(0.9743147) q[1];
rz(-pi) q[2];
rz(-1.1039626) q[3];
sx q[3];
rz(-1.8133014) q[3];
sx q[3];
rz(1.8777101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(-1.5345705) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(-3.0788105) q[1];
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
rz(-1.3604135) q[0];
sx q[0];
rz(-0.14260261) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7912748) q[2];
sx q[2];
rz(-2.2248474) q[2];
sx q[2];
rz(-2.0803114) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.891891) q[1];
sx q[1];
rz(-2.5981123) q[1];
sx q[1];
rz(1.5312503) q[1];
rz(1.005638) q[3];
sx q[3];
rz(-1.2747545) q[3];
sx q[3];
rz(-0.63488301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.442231) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(-2.7065275) q[2];
rz(0.89138952) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.94423914) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(-0.42022589) q[0];
rz(-2.6603783) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680465) q[0];
sx q[0];
rz(-1.9424115) q[0];
sx q[0];
rz(2.6269873) q[0];
x q[1];
rz(1.9040362) q[2];
sx q[2];
rz(-0.2313279) q[2];
sx q[2];
rz(-1.7817792) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.269304) q[1];
sx q[1];
rz(-1.270712) q[1];
sx q[1];
rz(1.3809204) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9514063) q[3];
sx q[3];
rz(-1.7577226) q[3];
sx q[3];
rz(0.72894064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8740365) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.9497005) q[0];
rz(2.2299956) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-1.0620767) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5545643) q[0];
sx q[0];
rz(-1.6977068) q[0];
sx q[0];
rz(2.1112604) q[0];
rz(-pi) q[1];
rz(-0.75762962) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(0.66495313) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.810079) q[1];
sx q[1];
rz(-0.96558297) q[1];
sx q[1];
rz(-2.5007351) q[1];
rz(-pi) q[2];
rz(-0.90537269) q[3];
sx q[3];
rz(-1.5983799) q[3];
sx q[3];
rz(1.2736922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-0.28327981) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.4184985) q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-1.6497242) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(-0.80560537) q[1];
x q[2];
rz(-2.9407223) q[3];
sx q[3];
rz(-1.4088653) q[3];
sx q[3];
rz(-1.1274606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(2.3000439) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(-0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08298824) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(2.9291005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2895848) q[0];
sx q[0];
rz(-0.6913018) q[0];
sx q[0];
rz(-1.9612802) q[0];
rz(-pi) q[1];
rz(3.100527) q[2];
sx q[2];
rz(-1.8827056) q[2];
sx q[2];
rz(1.8598156) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5960658) q[1];
sx q[1];
rz(-0.83546987) q[1];
sx q[1];
rz(2.2722785) q[1];
rz(-pi) q[2];
rz(2.4082765) q[3];
sx q[3];
rz(-1.7622669) q[3];
sx q[3];
rz(1.9460033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1797669) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.6869705) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385948) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(-0.075642792) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(1.5008925) q[2];
sx q[2];
rz(-1.5887661) q[2];
sx q[2];
rz(-0.80815732) q[2];
rz(1.1099439) q[3];
sx q[3];
rz(-2.0416595) q[3];
sx q[3];
rz(2.5817081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
