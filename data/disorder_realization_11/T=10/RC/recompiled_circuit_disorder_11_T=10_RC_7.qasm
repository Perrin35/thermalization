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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0044999997) q[0];
sx q[0];
rz(-2.922393) q[0];
sx q[0];
rz(2.5704434) q[0];
x q[1];
rz(1.6425743) q[2];
sx q[2];
rz(-1.1489747) q[2];
sx q[2];
rz(2.3274802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6424375) q[1];
sx q[1];
rz(-2.7381574) q[1];
sx q[1];
rz(-2.8950476) q[1];
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
rz(-pi/2) q[1];
rz(1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(-1.3936183) q[2];
rz(-2.3404775) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(1.0629517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-3.063607) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(0.35968131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29794824) q[0];
sx q[0];
rz(-2.1878562) q[0];
sx q[0];
rz(-0.88275568) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4215135) q[2];
sx q[2];
rz(-1.7290001) q[2];
sx q[2];
rz(0.0062696487) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82275326) q[1];
sx q[1];
rz(-1.0471294) q[1];
sx q[1];
rz(0.19284064) q[1];
rz(-pi) q[2];
rz(3.1017786) q[3];
sx q[3];
rz(-0.85714825) q[3];
sx q[3];
rz(3.0301222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2461207) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(2.1982511) q[2];
rz(1.2603166) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(0.70596203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(-2.673972) q[0];
rz(-2.6768661) q[1];
sx q[1];
rz(-2.6578891) q[1];
sx q[1];
rz(-0.31276774) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8899925) q[0];
sx q[0];
rz(-1.7822937) q[0];
sx q[0];
rz(0.081758008) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33212338) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(0.78813625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1130502) q[1];
sx q[1];
rz(-2.9187435) q[1];
sx q[1];
rz(0.88009665) q[1];
x q[2];
rz(2.8144366) q[3];
sx q[3];
rz(-1.4331685) q[3];
sx q[3];
rz(-2.2562192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9096845) q[2];
sx q[2];
rz(-1.5732876) q[2];
sx q[2];
rz(-0.30511937) q[2];
rz(1.8049847) q[3];
sx q[3];
rz(-2.2620585) q[3];
sx q[3];
rz(0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49622789) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(0.028045068) q[0];
rz(-1.4656981) q[1];
sx q[1];
rz(-0.60246712) q[1];
sx q[1];
rz(-2.4688597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3161635) q[0];
sx q[0];
rz(-2.5220693) q[0];
sx q[0];
rz(-2.6114458) q[0];
x q[1];
rz(-2.7364199) q[2];
sx q[2];
rz(-0.58779683) q[2];
sx q[2];
rz(1.6853465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67191254) q[1];
sx q[1];
rz(-1.3393991) q[1];
sx q[1];
rz(3.0242821) q[1];
x q[2];
rz(-1.3939875) q[3];
sx q[3];
rz(-2.5150931) q[3];
sx q[3];
rz(-0.83229438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(0.76081863) q[2];
rz(-2.2740254) q[3];
sx q[3];
rz(-1.3191728) q[3];
sx q[3];
rz(-2.3755465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-2.0779579) q[0];
rz(1.4798374) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(3.1033049) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83513658) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(-1.8655213) q[0];
rz(2.0615887) q[2];
sx q[2];
rz(-1.0133427) q[2];
sx q[2];
rz(1.4412396) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28336477) q[1];
sx q[1];
rz(-1.3971395) q[1];
sx q[1];
rz(-0.11939343) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0733662) q[3];
sx q[3];
rz(-2.6196819) q[3];
sx q[3];
rz(3.0038602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(1.0920452) q[2];
rz(1.6070222) q[3];
sx q[3];
rz(-0.97073308) q[3];
sx q[3];
rz(0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9606278) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(1.3076179) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(0.19097701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67577772) q[0];
sx q[0];
rz(-0.2535648) q[0];
sx q[0];
rz(0.98357865) q[0];
rz(-1.149876) q[2];
sx q[2];
rz(-0.72962609) q[2];
sx q[2];
rz(-0.52044496) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.891891) q[1];
sx q[1];
rz(-2.5981123) q[1];
sx q[1];
rz(1.6103423) q[1];
rz(-pi) q[2];
rz(-1.005638) q[3];
sx q[3];
rz(-1.2747545) q[3];
sx q[3];
rz(0.63488301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6993616) q[2];
sx q[2];
rz(-1.3012393) q[2];
sx q[2];
rz(-0.43506518) q[2];
rz(2.2502031) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94423914) q[0];
sx q[0];
rz(-0.22560142) q[0];
sx q[0];
rz(-0.42022589) q[0];
rz(-0.48121437) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-2.1113254) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39474836) q[0];
sx q[0];
rz(-2.0472102) q[0];
sx q[0];
rz(1.9917411) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7898125) q[2];
sx q[2];
rz(-1.6458626) q[2];
sx q[2];
rz(-2.605627) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75525857) q[1];
sx q[1];
rz(-1.3895021) q[1];
sx q[1];
rz(-0.30524409) q[1];
x q[2];
rz(-1.7610735) q[3];
sx q[3];
rz(-1.3839625) q[3];
sx q[3];
rz(2.2639757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39945012) q[2];
sx q[2];
rz(-0.43562499) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(1.0366084) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(-2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(1.1918921) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(2.079516) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9174791) q[0];
sx q[0];
rz(-0.55372059) q[0];
sx q[0];
rz(1.327716) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.383963) q[2];
sx q[2];
rz(-2.5839351) q[2];
sx q[2];
rz(0.66495313) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2506822) q[1];
sx q[1];
rz(-2.2905596) q[1];
sx q[1];
rz(0.85810424) q[1];
rz(3.1065337) q[3];
sx q[3];
rz(-2.2359214) q[3];
sx q[3];
rz(2.866131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.186211) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-0.051076802) q[2];
rz(2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(1.0572082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5687953) q[0];
sx q[0];
rz(-0.43734567) q[0];
sx q[0];
rz(1.019078) q[0];
rz(2.7678713) q[2];
sx q[2];
rz(-1.7127275) q[2];
sx q[2];
rz(2.2609557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6497242) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(-2.3359873) q[1];
x q[2];
rz(-0.68607761) q[3];
sx q[3];
rz(-2.8842673) q[3];
sx q[3];
rz(-2.0285332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(-1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-2.9291005) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3613113) q[0];
sx q[0];
rz(-2.2012976) q[0];
sx q[0];
rz(-0.30514858) q[0];
rz(1.6974405) q[2];
sx q[2];
rz(-2.8270792) q[2];
sx q[2];
rz(1.7267137) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5960658) q[1];
sx q[1];
rz(-0.83546987) q[1];
sx q[1];
rz(-2.2722785) q[1];
rz(-1.3155762) q[3];
sx q[3];
rz(-0.85382429) q[3];
sx q[3];
rz(2.9361801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1797669) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.4546222) q[2];
rz(-1.1949332) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-2.6132244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(1.7385948) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
rz(0.075642792) q[1];
sx q[1];
rz(-1.5374001) q[1];
sx q[1];
rz(-1.7061445) q[1];
rz(1.3189581) q[2];
sx q[2];
rz(-3.0694198) q[2];
sx q[2];
rz(1.0138489) q[2];
rz(-1.1099439) q[3];
sx q[3];
rz(-1.0999332) q[3];
sx q[3];
rz(-0.55988452) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
