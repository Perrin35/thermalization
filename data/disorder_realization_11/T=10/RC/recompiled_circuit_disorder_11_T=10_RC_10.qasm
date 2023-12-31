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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0044999997) q[0];
sx q[0];
rz(-2.922393) q[0];
sx q[0];
rz(-0.57114925) q[0];
rz(-2.9831224) q[2];
sx q[2];
rz(-2.7140693) q[2];
sx q[2];
rz(-0.98795623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6424375) q[1];
sx q[1];
rz(-0.40343522) q[1];
sx q[1];
rz(-2.8950476) q[1];
rz(-pi) q[2];
rz(3.0626489) q[3];
sx q[3];
rz(-1.5773838) q[3];
sx q[3];
rz(1.6062669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.7479744) q[2];
rz(-0.80111516) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-0.077985667) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(2.7819113) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4820837) q[0];
sx q[0];
rz(-2.2523899) q[0];
sx q[0];
rz(2.4115152) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23736575) q[2];
sx q[2];
rz(-0.73420213) q[2];
sx q[2];
rz(1.3993625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.82275326) q[1];
sx q[1];
rz(-2.0944632) q[1];
sx q[1];
rz(-2.948752) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1017786) q[3];
sx q[3];
rz(-2.2844444) q[3];
sx q[3];
rz(-3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2461207) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(-0.94334156) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(-0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
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
rz(0.62374672) q[0];
sx q[0];
rz(-2.9150634) q[0];
sx q[0];
rz(-1.9342599) q[0];
x q[1];
rz(2.2964988) q[2];
sx q[2];
rz(-1.8232864) q[2];
sx q[2];
rz(2.1408199) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1130502) q[1];
sx q[1];
rz(-0.22284914) q[1];
sx q[1];
rz(2.261496) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8144366) q[3];
sx q[3];
rz(-1.4331685) q[3];
sx q[3];
rz(0.88537346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(-2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(-0.67273295) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493895) q[0];
sx q[0];
rz(-2.0954847) q[0];
sx q[0];
rz(-1.9169109) q[0];
rz(-pi) q[1];
rz(0.40517278) q[2];
sx q[2];
rz(-2.5537958) q[2];
sx q[2];
rz(-1.6853465) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67191254) q[1];
sx q[1];
rz(-1.3393991) q[1];
sx q[1];
rz(-3.0242821) q[1];
rz(-1.3939875) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(-0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4611886) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(-0.76081863) q[2];
rz(0.86756724) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(0.76604617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(-1.4798374) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(-0.038287727) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69890755) q[0];
sx q[0];
rz(-1.2783056) q[0];
sx q[0];
rz(0.12683503) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5262881) q[2];
sx q[2];
rz(-1.9822789) q[2];
sx q[2];
rz(2.736511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28336477) q[1];
sx q[1];
rz(-1.7444532) q[1];
sx q[1];
rz(3.0221992) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1039626) q[3];
sx q[3];
rz(-1.8133014) q[3];
sx q[3];
rz(-1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-2.0495474) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18096481) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(3.0788105) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(2.9506156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.068114) q[0];
sx q[0];
rz(-1.3604135) q[0];
sx q[0];
rz(0.14260261) q[0];
rz(-pi) q[1];
rz(-1.149876) q[2];
sx q[2];
rz(-0.72962609) q[2];
sx q[2];
rz(-0.52044496) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.891891) q[1];
sx q[1];
rz(-2.5981123) q[1];
sx q[1];
rz(1.5312503) q[1];
rz(-pi) q[2];
rz(2.7950068) q[3];
sx q[3];
rz(-1.0329909) q[3];
sx q[3];
rz(-1.1188521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6993616) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1973535) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(0.42022589) q[0];
rz(0.48121437) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97354613) q[0];
sx q[0];
rz(-1.1991812) q[0];
sx q[0];
rz(2.6269873) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.7891857) q[2];
sx q[2];
rz(-1.0181392) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75525857) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(2.8363486) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19018634) q[3];
sx q[3];
rz(-1.7577226) q[3];
sx q[3];
rz(2.412652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(-1.0366084) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(1.1918921) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(2.079516) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5870283) q[0];
sx q[0];
rz(-1.4438859) q[0];
sx q[0];
rz(1.0303322) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75762962) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(-0.66495313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8909104) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(-2.2834884) q[1];
x q[2];
rz(3.1065337) q[3];
sx q[3];
rz(-0.9056712) q[3];
sx q[3];
rz(0.2754617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.186211) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(0.051076802) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.6254814) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.93808162) q[1];
sx q[1];
rz(0.84916806) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5687953) q[0];
sx q[0];
rz(-2.704247) q[0];
sx q[0];
rz(1.019078) q[0];
rz(-pi) q[1];
rz(2.7678713) q[2];
sx q[2];
rz(-1.7127275) q[2];
sx q[2];
rz(-0.88063699) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4918684) q[1];
sx q[1];
rz(-0.33136156) q[1];
sx q[1];
rz(2.3359873) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4056021) q[3];
sx q[3];
rz(-1.372589) q[3];
sx q[3];
rz(-2.7310731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7944305) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(-2.9629422) q[2];
rz(2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(-1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2895848) q[0];
sx q[0];
rz(-2.4502909) q[0];
sx q[0];
rz(-1.9612802) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6974405) q[2];
sx q[2];
rz(-2.8270792) q[2];
sx q[2];
rz(-1.7267137) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5455268) q[1];
sx q[1];
rz(-2.3061228) q[1];
sx q[1];
rz(-0.86931418) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28189567) q[3];
sx q[3];
rz(-2.3882139) q[3];
sx q[3];
rz(0.58338141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.6407001) q[2];
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
