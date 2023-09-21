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
rz(-2.4401234) q[1];
sx q[1];
rz(-2.520732) q[1];
sx q[1];
rz(-1.4863185) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0044999997) q[0];
sx q[0];
rz(-2.922393) q[0];
sx q[0];
rz(-0.57114925) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15847023) q[2];
sx q[2];
rz(-2.7140693) q[2];
sx q[2];
rz(0.98795623) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.15581407) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(-2.7491261) q[1];
rz(-3.0582534) q[3];
sx q[3];
rz(-0.079217521) q[3];
sx q[3];
rz(-3.0230429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0861337) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(-0.80111516) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(-1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(3.063607) q[0];
rz(2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
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
x q[1];
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
rz(-0.82275326) q[1];
sx q[1];
rz(-2.0944632) q[1];
sx q[1];
rz(0.19284064) q[1];
rz(-pi) q[2];
rz(-1.5248604) q[3];
sx q[3];
rz(-0.71456281) q[3];
sx q[3];
rz(3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89547196) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(1.2603166) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27186069) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(-2.673972) q[0];
rz(-0.46472654) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(2.8288249) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8395961) q[0];
sx q[0];
rz(-1.490864) q[0];
sx q[0];
rz(1.7829814) q[0];
rz(1.2000285) q[2];
sx q[2];
rz(-2.3808378) q[2];
sx q[2];
rz(2.8460381) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.028542472) q[1];
sx q[1];
rz(-2.9187435) q[1];
sx q[1];
rz(-2.261496) q[1];
x q[2];
rz(-2.7346482) q[3];
sx q[3];
rz(-2.787628) q[3];
sx q[3];
rz(2.0719761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(-2.8364733) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(2.4688597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69986491) q[0];
sx q[0];
rz(-1.8687975) q[0];
sx q[0];
rz(-2.5900048) q[0];
rz(-2.5920934) q[2];
sx q[2];
rz(-1.7911583) q[2];
sx q[2];
rz(-2.9133177) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87186253) q[1];
sx q[1];
rz(-1.6849663) q[1];
sx q[1];
rz(1.8037379) q[1];
rz(2.1898515) q[3];
sx q[3];
rz(-1.4674867) q[3];
sx q[3];
rz(2.2593474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68040401) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(2.380774) q[2];
rz(2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(2.3755465) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3064561) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(1.8655213) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(-2.736511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.89025154) q[1];
sx q[1];
rz(-2.9311935) q[1];
sx q[1];
rz(2.167278) q[1];
rz(2.0376301) q[3];
sx q[3];
rz(-1.3282913) q[3];
sx q[3];
rz(1.2638826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(2.0495474) q[2];
rz(-1.6070222) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-2.3585414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096481) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(-1.8339748) q[0];
rz(-3.0788105) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(2.9506156) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073478621) q[0];
sx q[0];
rz(-1.3604135) q[0];
sx q[0];
rz(2.99899) q[0];
x q[1];
rz(2.7912748) q[2];
sx q[2];
rz(-0.91674524) q[2];
sx q[2];
rz(-1.0612812) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8543507) q[1];
sx q[1];
rz(-1.5912424) q[1];
sx q[1];
rz(-1.0276621) q[1];
x q[2];
rz(1.0530924) q[3];
sx q[3];
rz(-0.63044237) q[3];
sx q[3];
rz(-2.6368486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.442231) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(2.7065275) q[2];
rz(-2.2502031) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.58266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(2.7213668) q[0];
rz(0.48121437) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9734605) q[0];
sx q[0];
rz(-0.62481835) q[0];
sx q[0];
rz(-2.4718651) q[0];
rz(1.7898125) q[2];
sx q[2];
rz(-1.6458626) q[2];
sx q[2];
rz(0.53596562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75525857) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(-2.8363486) q[1];
rz(-pi) q[2];
rz(1.3805192) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(2.2639757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(2.8263261) q[2];
rz(-2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(-2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.9497005) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(1.0620767) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0820205) q[0];
sx q[0];
rz(-2.1064415) q[0];
sx q[0];
rz(0.14772149) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42542142) q[2];
sx q[2];
rz(-1.1985949) q[2];
sx q[2];
rz(-0.22949716) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8909104) q[1];
sx q[1];
rz(-0.85103304) q[1];
sx q[1];
rz(0.85810424) q[1];
x q[2];
rz(0.90537269) q[3];
sx q[3];
rz(-1.5983799) q[3];
sx q[3];
rz(1.8679004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(3.0905159) q[2];
rz(-2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(-1.0572082) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5687953) q[0];
sx q[0];
rz(-0.43734567) q[0];
sx q[0];
rz(-2.1225147) q[0];
rz(-0.37306771) q[2];
sx q[2];
rz(-0.39857769) q[2];
sx q[2];
rz(2.1052436) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.65836425) q[1];
sx q[1];
rz(-1.798097) q[1];
sx q[1];
rz(1.3275654) q[1];
x q[2];
rz(2.9407223) q[3];
sx q[3];
rz(-1.7327274) q[3];
sx q[3];
rz(-1.1274606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34716216) q[2];
sx q[2];
rz(-0.55134761) q[2];
sx q[2];
rz(2.9629422) q[2];
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
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(-1.2471584) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(-0.21249214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025876306) q[0];
sx q[0];
rz(-1.3256782) q[0];
sx q[0];
rz(-2.2239767) q[0];
rz(-pi) q[1];
rz(1.4441522) q[2];
sx q[2];
rz(-0.31451348) q[2];
sx q[2];
rz(-1.414879) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.64582981) q[1];
sx q[1];
rz(-2.1730575) q[1];
sx q[1];
rz(-2.5219003) q[1];
rz(-pi) q[2];
rz(0.73331613) q[3];
sx q[3];
rz(-1.7622669) q[3];
sx q[3];
rz(-1.9460033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1797669) q[2];
sx q[2];
rz(-3.0111854) q[2];
sx q[2];
rz(-1.4546222) q[2];
rz(-1.1949332) q[3];
sx q[3];
rz(-1.3062198) q[3];
sx q[3];
rz(2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(0.51681896) q[3];
sx q[3];
rz(-1.9782981) q[3];
sx q[3];
rz(0.78936418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
