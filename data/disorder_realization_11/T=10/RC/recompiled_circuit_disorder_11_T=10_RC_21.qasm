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
rz(-1.4045665) q[0];
rz(3.8430619) q[1];
sx q[1];
rz(3.7624533) q[1];
sx q[1];
rz(7.9384595) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0044999997) q[0];
sx q[0];
rz(-0.21919964) q[0];
sx q[0];
rz(-2.5704434) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9831224) q[2];
sx q[2];
rz(-0.42752334) q[2];
sx q[2];
rz(-0.98795623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15581407) q[1];
sx q[1];
rz(-1.6667546) q[1];
sx q[1];
rz(-2.7491261) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0626489) q[3];
sx q[3];
rz(-1.5642089) q[3];
sx q[3];
rz(-1.6062669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(-1.3936183) q[2];
rz(-2.3404775) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(0.077985667) q[0];
rz(-0.69411913) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(2.7819113) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168424) q[0];
sx q[0];
rz(-2.1152088) q[0];
sx q[0];
rz(0.74290468) q[0];
rz(-pi) q[1];
rz(-0.72007911) q[2];
sx q[2];
rz(-1.7290001) q[2];
sx q[2];
rz(3.135323) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6911192) q[1];
sx q[1];
rz(-0.55492655) q[1];
sx q[1];
rz(-1.2503442) q[1];
rz(1.6167323) q[3];
sx q[3];
rz(-2.4270298) q[3];
sx q[3];
rz(-3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(2.673972) q[0];
rz(-0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(0.31276774) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5178459) q[0];
sx q[0];
rz(-2.9150634) q[0];
sx q[0];
rz(1.2073327) q[0];
x q[1];
rz(-0.84509387) q[2];
sx q[2];
rz(-1.3183062) q[2];
sx q[2];
rz(1.0007728) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1130502) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(2.261496) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9096845) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(-1.336608) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(2.4688597) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69220316) q[0];
sx q[0];
rz(-1.0461079) q[0];
sx q[0];
rz(-1.9169109) q[0];
rz(-pi) q[1];
rz(2.7364199) q[2];
sx q[2];
rz(-2.5537958) q[2];
sx q[2];
rz(1.6853465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1466009) q[1];
sx q[1];
rz(-2.8826337) q[1];
sx q[1];
rz(1.109757) q[1];
rz(-1.7476051) q[3];
sx q[3];
rz(-0.62649957) q[3];
sx q[3];
rz(-0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.5701533) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(1.6617552) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(0.038287727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83513658) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(1.2760713) q[0];
rz(-pi) q[1];
rz(0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(-2.736511) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89025154) q[1];
sx q[1];
rz(-0.21039911) q[1];
sx q[1];
rz(-0.9743147) q[1];
rz(-1.0682265) q[3];
sx q[3];
rz(-0.52191075) q[3];
sx q[3];
rz(-3.0038602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.15497196) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(-2.0495474) q[2];
rz(-1.6070222) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096481) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(2.9506156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4673429) q[0];
sx q[0];
rz(-1.4313587) q[0];
sx q[0];
rz(1.3583202) q[0];
rz(-1.149876) q[2];
sx q[2];
rz(-0.72962609) q[2];
sx q[2];
rz(2.6211477) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29589614) q[1];
sx q[1];
rz(-2.1138043) q[1];
sx q[1];
rz(0.023881749) q[1];
x q[2];
rz(-2.0885002) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(0.5047441) q[3];
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
rz(-1.9458709) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.2599506) q[1];
sx q[1];
rz(-1.0302672) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734605) q[0];
sx q[0];
rz(-0.62481835) q[0];
sx q[0];
rz(-0.66972759) q[0];
rz(-pi) q[1];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.352407) q[2];
sx q[2];
rz(-2.1234535) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.269304) q[1];
sx q[1];
rz(-1.270712) q[1];
sx q[1];
rz(1.3809204) q[1];
rz(1.7610735) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(-1.0366084) q[3];
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
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(1.1918921) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(-2.079516) q[1];
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
rz(0.75762962) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(2.4766395) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8909104) q[1];
sx q[1];
rz(-2.2905596) q[1];
sx q[1];
rz(-0.85810424) q[1];
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
rz(-pi/2) q[1];
rz(1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.1445769) q[0];
sx q[0];
rz(-2.7269195) q[0];
sx q[0];
rz(-0.64494079) q[0];
rz(2.0058517) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(0.84916806) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51061714) q[0];
sx q[0];
rz(-1.7946587) q[0];
sx q[0];
rz(-1.9497245) q[0];
rz(-pi) q[1];
rz(0.37306771) q[2];
sx q[2];
rz(-2.743015) q[2];
sx q[2];
rz(-1.0363491) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4918684) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(2.3359873) q[1];
rz(-1.4056021) q[3];
sx q[3];
rz(-1.372589) q[3];
sx q[3];
rz(2.7310731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-2.9629422) q[2];
rz(0.84154877) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(1.8944342) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025876306) q[0];
sx q[0];
rz(-1.8159144) q[0];
sx q[0];
rz(0.91761597) q[0];
rz(-0.041065644) q[2];
sx q[2];
rz(-1.2588871) q[2];
sx q[2];
rz(1.2817771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6007235) q[1];
sx q[1];
rz(-2.0698554) q[1];
sx q[1];
rz(-2.2720816) q[1];
rz(-2.859697) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(-2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96182573) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(1.1949332) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(0.52836829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(0.075642792) q[1];
sx q[1];
rz(-1.5374001) q[1];
sx q[1];
rz(-1.7061445) q[1];
rz(-0.018013714) q[2];
sx q[2];
rz(-1.5009038) q[2];
sx q[2];
rz(0.76138087) q[2];
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