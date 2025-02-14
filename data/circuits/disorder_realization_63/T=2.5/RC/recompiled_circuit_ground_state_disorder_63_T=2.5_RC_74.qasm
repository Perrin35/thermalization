OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0286921) q[0];
sx q[0];
rz(-2.8960462) q[0];
sx q[0];
rz(-3.088933) q[0];
rz(1.117299) q[1];
sx q[1];
rz(-2.8007562) q[1];
sx q[1];
rz(-2.053082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3153305) q[0];
sx q[0];
rz(-1.966094) q[0];
sx q[0];
rz(3.0733103) q[0];
rz(-1.4289189) q[2];
sx q[2];
rz(-1.4617702) q[2];
sx q[2];
rz(3.1224868) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0170235) q[1];
sx q[1];
rz(-0.93734159) q[1];
sx q[1];
rz(1.3658466) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1375211) q[3];
sx q[3];
rz(-1.1511251) q[3];
sx q[3];
rz(1.3908902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6344305) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(0.6081028) q[2];
rz(1.1242584) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(2.2500136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41481498) q[0];
sx q[0];
rz(-1.5730653) q[0];
sx q[0];
rz(-2.538105) q[0];
rz(-0.51015774) q[1];
sx q[1];
rz(-2.3699103) q[1];
sx q[1];
rz(-2.1147494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48906836) q[0];
sx q[0];
rz(-1.4155672) q[0];
sx q[0];
rz(-0.17431577) q[0];
rz(-0.6925236) q[2];
sx q[2];
rz(-2.1998458) q[2];
sx q[2];
rz(-0.72450996) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4812117) q[1];
sx q[1];
rz(-2.2218895) q[1];
sx q[1];
rz(-2.2269985) q[1];
rz(-pi) q[2];
rz(-0.92225109) q[3];
sx q[3];
rz(-0.77510683) q[3];
sx q[3];
rz(0.67476455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5828731) q[2];
sx q[2];
rz(-1.1483973) q[2];
sx q[2];
rz(-0.64644512) q[2];
rz(1.2785814) q[3];
sx q[3];
rz(-0.94978142) q[3];
sx q[3];
rz(1.4518552) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.932514) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(-2.7945844) q[1];
sx q[1];
rz(-1.4748814) q[1];
sx q[1];
rz(-0.86403799) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61053172) q[0];
sx q[0];
rz(-2.5342434) q[0];
sx q[0];
rz(2.7694398) q[0];
x q[1];
rz(1.3190325) q[2];
sx q[2];
rz(-1.1472817) q[2];
sx q[2];
rz(-0.26811312) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.62926312) q[1];
sx q[1];
rz(-2.1844668) q[1];
sx q[1];
rz(-0.83028173) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0357843) q[3];
sx q[3];
rz(-2.0943301) q[3];
sx q[3];
rz(0.86712718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.049909264) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(1.4570215) q[2];
rz(1.0157887) q[3];
sx q[3];
rz(-0.53272811) q[3];
sx q[3];
rz(0.31639019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24093534) q[0];
sx q[0];
rz(-0.37739402) q[0];
sx q[0];
rz(-1.578791) q[0];
rz(2.0511625) q[1];
sx q[1];
rz(-2.5999887) q[1];
sx q[1];
rz(1.9900367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59127677) q[0];
sx q[0];
rz(-1.341686) q[0];
sx q[0];
rz(0.59058779) q[0];
rz(2.3811023) q[2];
sx q[2];
rz(-1.2142688) q[2];
sx q[2];
rz(-2.6127151) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87428625) q[1];
sx q[1];
rz(-1.5360254) q[1];
sx q[1];
rz(1.0465996) q[1];
rz(-pi) q[2];
rz(2.1776803) q[3];
sx q[3];
rz(-1.5935894) q[3];
sx q[3];
rz(-0.43302872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2942723) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(-1.4997743) q[2];
rz(-3.0070987) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(2.3597778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0624369) q[0];
sx q[0];
rz(-2.2298614) q[0];
sx q[0];
rz(-2.675918) q[0];
rz(1.2767208) q[1];
sx q[1];
rz(-2.1379505) q[1];
sx q[1];
rz(-1.0328971) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532294) q[0];
sx q[0];
rz(-1.7288701) q[0];
sx q[0];
rz(3.1309978) q[0];
x q[1];
rz(1.6711556) q[2];
sx q[2];
rz(-1.7102229) q[2];
sx q[2];
rz(3.0318236) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.031449854) q[1];
sx q[1];
rz(-1.8217161) q[1];
sx q[1];
rz(-0.67278426) q[1];
rz(-2.1139268) q[3];
sx q[3];
rz(-1.2333721) q[3];
sx q[3];
rz(-1.8130482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59139645) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-2.2097394) q[2];
rz(-3.0469117) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(-1.3100821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69022995) q[0];
sx q[0];
rz(-2.9424423) q[0];
sx q[0];
rz(3.0354101) q[0];
rz(2.5871318) q[1];
sx q[1];
rz(-2.3531871) q[1];
sx q[1];
rz(-1.6000481) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0743177) q[0];
sx q[0];
rz(-1.5467001) q[0];
sx q[0];
rz(-0.39053183) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0279858) q[2];
sx q[2];
rz(-1.3234252) q[2];
sx q[2];
rz(-0.24277011) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2682802) q[1];
sx q[1];
rz(-2.2619216) q[1];
sx q[1];
rz(1.4119488) q[1];
x q[2];
rz(-2.8951449) q[3];
sx q[3];
rz(-1.6250027) q[3];
sx q[3];
rz(-1.1185631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23325486) q[2];
sx q[2];
rz(-2.4007863) q[2];
sx q[2];
rz(-1.6439269) q[2];
rz(2.6805367) q[3];
sx q[3];
rz(-2.4888829) q[3];
sx q[3];
rz(1.3155931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(-1.0787971) q[0];
rz(-0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(2.914391) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0086454) q[0];
sx q[0];
rz(-1.7194178) q[0];
sx q[0];
rz(0.11283837) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14102139) q[2];
sx q[2];
rz(-1.900809) q[2];
sx q[2];
rz(-0.38244707) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0409483) q[1];
sx q[1];
rz(-0.64926636) q[1];
sx q[1];
rz(-2.098564) q[1];
rz(-pi) q[2];
rz(0.10031767) q[3];
sx q[3];
rz(-1.9411818) q[3];
sx q[3];
rz(-0.60587347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5220773) q[2];
sx q[2];
rz(-0.41831133) q[2];
sx q[2];
rz(-0.72959161) q[2];
rz(-1.6884165) q[3];
sx q[3];
rz(-1.4540693) q[3];
sx q[3];
rz(-2.5280473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59298092) q[0];
sx q[0];
rz(-0.91658968) q[0];
sx q[0];
rz(-2.7899637) q[0];
rz(-2.8278606) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(1.7650013) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2257654) q[0];
sx q[0];
rz(-1.6179913) q[0];
sx q[0];
rz(2.7412358) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96504546) q[2];
sx q[2];
rz(-0.20119431) q[2];
sx q[2];
rz(0.66058841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2699956) q[1];
sx q[1];
rz(-1.0078127) q[1];
sx q[1];
rz(1.3857122) q[1];
rz(1.3198765) q[3];
sx q[3];
rz(-0.38057113) q[3];
sx q[3];
rz(1.6264834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54479638) q[2];
sx q[2];
rz(-2.4185541) q[2];
sx q[2];
rz(-1.6205622) q[2];
rz(0.14791402) q[3];
sx q[3];
rz(-1.7971797) q[3];
sx q[3];
rz(1.3676876) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77383298) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.2441147) q[0];
rz(1.4338088) q[1];
sx q[1];
rz(-1.560874) q[1];
sx q[1];
rz(0.22020766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4996195) q[0];
sx q[0];
rz(-2.4110378) q[0];
sx q[0];
rz(-0.99910183) q[0];
rz(-pi) q[1];
rz(0.16248361) q[2];
sx q[2];
rz(-1.1097317) q[2];
sx q[2];
rz(-1.7562255) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9088194) q[1];
sx q[1];
rz(-2.4202883) q[1];
sx q[1];
rz(-2.2386293) q[1];
rz(1.4768019) q[3];
sx q[3];
rz(-1.6686182) q[3];
sx q[3];
rz(-2.0879012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9571017) q[2];
sx q[2];
rz(-1.9944921) q[2];
sx q[2];
rz(-1.6415049) q[2];
rz(0.084065048) q[3];
sx q[3];
rz(-1.8819239) q[3];
sx q[3];
rz(-3.0715004) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641007) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(-0.41859928) q[0];
rz(-2.1981926) q[1];
sx q[1];
rz(-1.489233) q[1];
sx q[1];
rz(-2.8387866) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26696268) q[0];
sx q[0];
rz(-1.7374674) q[0];
sx q[0];
rz(-3.002949) q[0];
rz(-pi) q[1];
rz(1.3933106) q[2];
sx q[2];
rz(-1.7587593) q[2];
sx q[2];
rz(1.6924015) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9040519) q[1];
sx q[1];
rz(-1.6142323) q[1];
sx q[1];
rz(3.0588845) q[1];
x q[2];
rz(-0.11649152) q[3];
sx q[3];
rz(-1.9132275) q[3];
sx q[3];
rz(-0.064700944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4960949) q[2];
sx q[2];
rz(-2.2690513) q[2];
sx q[2];
rz(2.3910451) q[2];
rz(-1.4613072) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(0.51699483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1207598) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(1.9137406) q[1];
sx q[1];
rz(-0.49929437) q[1];
sx q[1];
rz(1.1566537) q[1];
rz(-2.5748809) q[2];
sx q[2];
rz(-1.7253582) q[2];
sx q[2];
rz(2.2127989) q[2];
rz(1.0629366) q[3];
sx q[3];
rz(-2.776317) q[3];
sx q[3];
rz(-1.1767514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
