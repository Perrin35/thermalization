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
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(1.4863185) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1370927) q[0];
sx q[0];
rz(-0.21919964) q[0];
sx q[0];
rz(2.5704434) q[0];
rz(-pi) q[1];
rz(1.6425743) q[2];
sx q[2];
rz(-1.992618) q[2];
sx q[2];
rz(0.81411241) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6424375) q[1];
sx q[1];
rz(-2.7381574) q[1];
sx q[1];
rz(-0.24654504) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.08333929) q[3];
sx q[3];
rz(-0.079217521) q[3];
sx q[3];
rz(-0.11854974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(-2.3404775) q[3];
sx q[3];
rz(-1.5812185) q[3];
sx q[3];
rz(2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79008094) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-0.077985667) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-2.0072939) q[1];
sx q[1];
rz(-2.7819113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168424) q[0];
sx q[0];
rz(-1.0263838) q[0];
sx q[0];
rz(2.398688) q[0];
x q[1];
rz(-0.72007911) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(0.0062696487) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6911192) q[1];
sx q[1];
rz(-0.55492655) q[1];
sx q[1];
rz(-1.8912485) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(2.2461207) q[2];
sx q[2];
rz(-1.8254231) q[2];
sx q[2];
rz(-2.1982511) q[2];
rz(1.881276) q[3];
sx q[3];
rz(-0.80788079) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.869732) q[0];
sx q[0];
rz(-1.4881217) q[0];
sx q[0];
rz(-0.46762064) q[0];
rz(2.6768661) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(2.8288249) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395961) q[0];
sx q[0];
rz(-1.6507286) q[0];
sx q[0];
rz(1.7829814) q[0];
rz(1.9415641) q[2];
sx q[2];
rz(-2.3808378) q[2];
sx q[2];
rz(-2.8460381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.73156563) q[1];
sx q[1];
rz(-1.7419852) q[1];
sx q[1];
rz(-0.14337916) q[1];
x q[2];
rz(-0.32715601) q[3];
sx q[3];
rz(-1.4331685) q[3];
sx q[3];
rz(-2.2562192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23190817) q[2];
sx q[2];
rz(-1.5732876) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-3.1135476) q[0];
rz(1.4656981) q[1];
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
rz(2.6114458) q[0];
rz(-pi) q[1];
rz(1.8276617) q[2];
sx q[2];
rz(-1.0360403) q[2];
sx q[2];
rz(1.9321439) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87186253) q[1];
sx q[1];
rz(-1.4566263) q[1];
sx q[1];
rz(-1.3378548) q[1];
rz(3.014971) q[3];
sx q[3];
rz(-2.1860578) q[3];
sx q[3];
rz(2.5263853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4611886) q[2];
sx q[2];
rz(-0.79128069) q[2];
sx q[2];
rz(-2.380774) q[2];
rz(0.86756724) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(2.3755465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5714394) q[0];
sx q[0];
rz(-2.0338991) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(1.6617552) q[1];
sx q[1];
rz(-2.1789443) q[1];
sx q[1];
rz(-3.1033049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3064561) q[0];
sx q[0];
rz(-1.6922173) q[0];
sx q[0];
rz(1.2760713) q[0];
rz(0.61530453) q[2];
sx q[2];
rz(-1.1593137) q[2];
sx q[2];
rz(-2.736511) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2667065) q[1];
sx q[1];
rz(-1.6883856) q[1];
sx q[1];
rz(1.3959195) q[1];
x q[2];
rz(-1.0682265) q[3];
sx q[3];
rz(-2.6196819) q[3];
sx q[3];
rz(-0.13773242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(2.0495474) q[2];
rz(-1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(-0.78305125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9606278) q[0];
sx q[0];
rz(-1.5009078) q[0];
sx q[0];
rz(1.3076179) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.1761913) q[1];
sx q[1];
rz(-0.19097701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4673429) q[0];
sx q[0];
rz(-1.4313587) q[0];
sx q[0];
rz(1.7832725) q[0];
x q[1];
rz(0.35031788) q[2];
sx q[2];
rz(-0.91674524) q[2];
sx q[2];
rz(1.0612812) q[2];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7950068) q[3];
sx q[3];
rz(-1.0329909) q[3];
sx q[3];
rz(-1.1188521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.442231) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(2.7065275) q[2];
rz(-0.89138952) q[3];
sx q[3];
rz(-1.1957217) q[3];
sx q[3];
rz(1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(0.42022589) q[0];
rz(0.48121437) q[1];
sx q[1];
rz(-0.88164202) q[1];
sx q[1];
rz(-2.1113254) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7468443) q[0];
sx q[0];
rz(-2.0472102) q[0];
sx q[0];
rz(-1.9917411) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.076896197) q[2];
sx q[2];
rz(-1.352407) q[2];
sx q[2];
rz(1.0181392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.269304) q[1];
sx q[1];
rz(-1.8708806) q[1];
sx q[1];
rz(-1.3809204) q[1];
rz(2.9514063) q[3];
sx q[3];
rz(-1.38387) q[3];
sx q[3];
rz(2.412652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7421425) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8740365) q[0];
sx q[0];
rz(-2.3836305) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-0.9115971) q[1];
sx q[1];
rz(-0.67271581) q[1];
sx q[1];
rz(-1.0620767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5870283) q[0];
sx q[0];
rz(-1.4438859) q[0];
sx q[0];
rz(-1.0303322) q[0];
rz(-2.383963) q[2];
sx q[2];
rz(-0.55765753) q[2];
sx q[2];
rz(-0.66495313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3036365) q[1];
sx q[1];
rz(-1.0567697) q[1];
sx q[1];
rz(-0.85876667) q[1];
rz(-3.1065337) q[3];
sx q[3];
rz(-2.2359214) q[3];
sx q[3];
rz(-2.866131) q[3];
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
rz(-2.4222899) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1445769) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(0.64494079) q[0];
rz(2.0058517) q[1];
sx q[1];
rz(-2.203511) q[1];
sx q[1];
rz(-0.84916806) q[1];
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
rz(0.37372132) q[2];
sx q[2];
rz(-1.4288651) q[2];
sx q[2];
rz(-0.88063699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4918684) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(-0.80560537) q[1];
rz(-2.9407223) q[3];
sx q[3];
rz(-1.7327274) q[3];
sx q[3];
rz(1.1274606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.34716216) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(-0.17865044) q[2];
rz(-2.3000439) q[3];
sx q[3];
rz(-1.0431362) q[3];
sx q[3];
rz(0.49079045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.08298824) q[0];
sx q[0];
rz(-1.2036136) q[0];
sx q[0];
rz(0.9151181) q[0];
rz(1.2471584) q[1];
sx q[1];
rz(-1.1727389) q[1];
sx q[1];
rz(-0.21249214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8520078) q[0];
sx q[0];
rz(-2.4502909) q[0];
sx q[0];
rz(-1.1803124) q[0];
rz(-pi) q[1];
rz(-1.6974405) q[2];
sx q[2];
rz(-0.31451348) q[2];
sx q[2];
rz(-1.414879) q[2];
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
rz(2.5219003) q[1];
x q[2];
rz(-2.4082765) q[3];
sx q[3];
rz(-1.7622669) q[3];
sx q[3];
rz(-1.9460033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385948) q[0];
sx q[0];
rz(-1.7398555) q[0];
sx q[0];
rz(1.6237988) q[0];
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
rz(-2.6247737) q[3];
sx q[3];
rz(-1.9782981) q[3];
sx q[3];
rz(0.78936418) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
