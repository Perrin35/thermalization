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
rz(-2.4401234) q[1];
sx q[1];
rz(-2.520732) q[1];
sx q[1];
rz(-1.4863185) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1264869) q[0];
sx q[0];
rz(-1.4529714) q[0];
sx q[0];
rz(2.9563222) q[0];
x q[1];
rz(2.9831224) q[2];
sx q[2];
rz(-2.7140693) q[2];
sx q[2];
rz(0.98795623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9857786) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(-2.7491261) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.078943723) q[3];
sx q[3];
rz(-1.5773838) q[3];
sx q[3];
rz(-1.5353257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0861337) q[2];
sx q[2];
rz(-1.7666585) q[2];
sx q[2];
rz(-1.7479744) q[2];
rz(2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-3.063607) q[0];
rz(-2.4474735) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(0.35968131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4820837) q[0];
sx q[0];
rz(-2.2523899) q[0];
sx q[0];
rz(-0.73007749) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7799136) q[2];
sx q[2];
rz(-2.2799727) q[2];
sx q[2];
rz(1.4271971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82275326) q[1];
sx q[1];
rz(-2.0944632) q[1];
sx q[1];
rz(-2.948752) q[1];
rz(-pi) q[2];
rz(1.5248604) q[3];
sx q[3];
rz(-0.71456281) q[3];
sx q[3];
rz(-3.0909017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2461207) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(-0.94334156) q[2];
rz(-1.2603166) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(-2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.27186069) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(2.673972) q[0];
rz(0.46472654) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-0.31276774) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5178459) q[0];
sx q[0];
rz(-2.9150634) q[0];
sx q[0];
rz(1.9342599) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33212338) q[2];
sx q[2];
rz(-0.87288522) q[2];
sx q[2];
rz(-0.78813625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.028542472) q[1];
sx q[1];
rz(-2.9187435) q[1];
sx q[1];
rz(2.261496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32715601) q[3];
sx q[3];
rz(-1.4331685) q[3];
sx q[3];
rz(-0.88537346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(2.6763776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6453648) q[0];
sx q[0];
rz(-0.35571686) q[0];
sx q[0];
rz(3.1135476) q[0];
rz(1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(2.4688597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493895) q[0];
sx q[0];
rz(-1.0461079) q[0];
sx q[0];
rz(1.9169109) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3139309) q[2];
sx q[2];
rz(-2.1055524) q[2];
sx q[2];
rz(-1.2094487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.67191254) q[1];
sx q[1];
rz(-1.8021936) q[1];
sx q[1];
rz(0.11731053) q[1];
rz(-pi) q[2];
rz(-1.3939875) q[3];
sx q[3];
rz(-0.62649957) q[3];
sx q[3];
rz(0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.68040401) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(-2.380774) q[2];
rz(2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(2.3755465) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(2.0779579) q[0];
rz(1.4798374) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(-0.038287727) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4426851) q[0];
sx q[0];
rz(-1.2783056) q[0];
sx q[0];
rz(0.12683503) q[0];
x q[1];
rz(-0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(-0.40508168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2667065) q[1];
sx q[1];
rz(-1.6883856) q[1];
sx q[1];
rz(1.7456732) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27023817) q[3];
sx q[3];
rz(-1.1186557) q[3];
sx q[3];
rz(-0.42735344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9866207) q[2];
sx q[2];
rz(-1.203323) q[2];
sx q[2];
rz(1.0920452) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096481) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(-3.0788105) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(-2.9506156) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073478621) q[0];
sx q[0];
rz(-1.7811791) q[0];
sx q[0];
rz(-2.99899) q[0];
rz(-pi) q[1];
rz(-1.9917166) q[2];
sx q[2];
rz(-2.4119666) q[2];
sx q[2];
rz(2.6211477) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2497017) q[1];
sx q[1];
rz(-2.5981123) q[1];
sx q[1];
rz(-1.5312503) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1359547) q[3];
sx q[3];
rz(-1.8668381) q[3];
sx q[3];
rz(-0.63488301) q[3];
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
rz(0.89138952) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(1.5589327) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-0.42022589) q[0];
rz(2.6603783) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39474836) q[0];
sx q[0];
rz(-1.0943824) q[0];
sx q[0];
rz(1.1498515) q[0];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.352407) q[2];
sx q[2];
rz(1.0181392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75525857) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(2.8363486) q[1];
rz(-1.7610735) q[3];
sx q[3];
rz(-1.7576302) q[3];
sx q[3];
rz(0.877617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(-2.8263261) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.8866115) q[3];
sx q[3];
rz(0.96926564) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(2.2299956) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(2.079516) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22411352) q[0];
sx q[0];
rz(-0.55372059) q[0];
sx q[0];
rz(-1.327716) q[0];
x q[1];
rz(-0.75762962) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(0.66495313) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3036365) q[1];
sx q[1];
rz(-1.0567697) q[1];
sx q[1];
rz(2.282826) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90537269) q[3];
sx q[3];
rz(-1.5432127) q[3];
sx q[3];
rz(-1.8679004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9553817) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99701571) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(-2.4966519) q[0];
rz(-1.1357409) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-2.2924246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5727974) q[0];
sx q[0];
rz(-0.43734567) q[0];
sx q[0];
rz(1.019078) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7230942) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(2.5068482) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6497242) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(-2.3359873) q[1];
rz(-pi) q[2];
rz(-0.20087033) q[3];
sx q[3];
rz(-1.7327274) q[3];
sx q[3];
rz(2.014132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(2.3000439) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(-2.2264746) q[0];
rz(1.2471584) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(0.21249214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.4441522) q[2];
sx q[2];
rz(-0.31451348) q[2];
sx q[2];
rz(1.414879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5960658) q[1];
sx q[1];
rz(-2.3061228) q[1];
sx q[1];
rz(2.2722785) q[1];
rz(0.28189567) q[3];
sx q[3];
rz(-0.75337871) q[3];
sx q[3];
rz(-2.5582112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96182573) q[2];
sx q[2];
rz(-3.0111854) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-3.1235789) q[2];
sx q[2];
rz(-1.6406888) q[2];
sx q[2];
rz(-2.3802118) q[2];
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
