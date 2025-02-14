OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71894574) q[0];
sx q[0];
rz(6.8721539) q[0];
sx q[0];
rz(6.4018259) q[0];
rz(-0.9912107) q[1];
sx q[1];
rz(-2.0999496) q[1];
sx q[1];
rz(-0.51377327) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.936718) q[0];
sx q[0];
rz(-2.4002909) q[0];
sx q[0];
rz(0.78420774) q[0];
rz(-0.1206042) q[2];
sx q[2];
rz(-2.6636811) q[2];
sx q[2];
rz(1.8574886) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2456677) q[1];
sx q[1];
rz(-1.8918719) q[1];
sx q[1];
rz(1.7968057) q[1];
x q[2];
rz(3.0904055) q[3];
sx q[3];
rz(-0.59268206) q[3];
sx q[3];
rz(0.41646233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.674268) q[2];
sx q[2];
rz(-1.6215723) q[2];
sx q[2];
rz(-0.33217126) q[2];
rz(-0.25003555) q[3];
sx q[3];
rz(-1.2269521) q[3];
sx q[3];
rz(-1.2927607) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2501204) q[0];
sx q[0];
rz(-1.253506) q[0];
sx q[0];
rz(0.44723311) q[0];
rz(2.1121292) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(-0.10118016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3295022) q[0];
sx q[0];
rz(-2.0116099) q[0];
sx q[0];
rz(2.0157218) q[0];
rz(-pi) q[1];
rz(2.9415628) q[2];
sx q[2];
rz(-2.1406271) q[2];
sx q[2];
rz(-2.9404158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78522462) q[1];
sx q[1];
rz(-1.1572946) q[1];
sx q[1];
rz(1.861839) q[1];
rz(-pi) q[2];
rz(1.4199791) q[3];
sx q[3];
rz(-0.64435092) q[3];
sx q[3];
rz(2.2805813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0826565) q[2];
sx q[2];
rz(-2.3857748) q[2];
sx q[2];
rz(1.5400881) q[2];
rz(-0.61795175) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(-1.9790953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59442941) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(0.51602236) q[0];
rz(2.0891321) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(1.1722391) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2502278) q[0];
sx q[0];
rz(-0.42420441) q[0];
sx q[0];
rz(0.026861743) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6555696) q[2];
sx q[2];
rz(-2.4293278) q[2];
sx q[2];
rz(-2.149947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7683923) q[1];
sx q[1];
rz(-0.52881634) q[1];
sx q[1];
rz(1.1952728) q[1];
rz(0.013962176) q[3];
sx q[3];
rz(-2.3413071) q[3];
sx q[3];
rz(1.518702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27511328) q[2];
sx q[2];
rz(-1.8786414) q[2];
sx q[2];
rz(2.4647554) q[2];
rz(-0.46946851) q[3];
sx q[3];
rz(-2.2468086) q[3];
sx q[3];
rz(1.9278056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6770099) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(-1.1997892) q[0];
rz(-0.59263539) q[1];
sx q[1];
rz(-2.0694144) q[1];
sx q[1];
rz(-1.4592272) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8358187) q[0];
sx q[0];
rz(-0.61821067) q[0];
sx q[0];
rz(2.799116) q[0];
rz(-pi) q[1];
rz(2.8921669) q[2];
sx q[2];
rz(-2.327033) q[2];
sx q[2];
rz(3.0944173) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0778592) q[1];
sx q[1];
rz(-2.1850536) q[1];
sx q[1];
rz(2.8251404) q[1];
x q[2];
rz(1.9497037) q[3];
sx q[3];
rz(-2.2188596) q[3];
sx q[3];
rz(-0.60175446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6478641) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(-2.4596821) q[2];
rz(2.72825) q[3];
sx q[3];
rz(-1.6091434) q[3];
sx q[3];
rz(-1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84363627) q[0];
sx q[0];
rz(-1.4763259) q[0];
sx q[0];
rz(-2.8635039) q[0];
rz(-0.9043215) q[1];
sx q[1];
rz(-1.4537289) q[1];
sx q[1];
rz(-2.1542737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073719115) q[0];
sx q[0];
rz(-2.4885396) q[0];
sx q[0];
rz(-1.9685754) q[0];
rz(-pi) q[1];
rz(2.9402307) q[2];
sx q[2];
rz(-0.1397396) q[2];
sx q[2];
rz(-2.1167378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6892859) q[1];
sx q[1];
rz(-2.747162) q[1];
sx q[1];
rz(-1.4107262) q[1];
rz(-2.868949) q[3];
sx q[3];
rz(-1.050625) q[3];
sx q[3];
rz(-0.70462728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.014293369) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(2.3822752) q[2];
rz(-1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21869126) q[0];
sx q[0];
rz(-2.3221115) q[0];
sx q[0];
rz(-2.4500093) q[0];
rz(-0.85451952) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(0.42133322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3499419) q[0];
sx q[0];
rz(-1.6393108) q[0];
sx q[0];
rz(1.1380026) q[0];
rz(0.26701228) q[2];
sx q[2];
rz(-1.2894372) q[2];
sx q[2];
rz(2.6002778) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1189733) q[1];
sx q[1];
rz(-1.2105064) q[1];
sx q[1];
rz(1.6605571) q[1];
x q[2];
rz(0.38479067) q[3];
sx q[3];
rz(-2.5005385) q[3];
sx q[3];
rz(1.3271063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2104346) q[2];
sx q[2];
rz(-1.6633818) q[2];
sx q[2];
rz(-0.20273905) q[2];
rz(1.5431131) q[3];
sx q[3];
rz(-2.8212382) q[3];
sx q[3];
rz(-1.4690442) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572674) q[0];
sx q[0];
rz(-1.4284644) q[0];
sx q[0];
rz(-2.7427234) q[0];
rz(1.9786037) q[1];
sx q[1];
rz(-0.82610026) q[1];
sx q[1];
rz(0.2805447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5359176) q[0];
sx q[0];
rz(-1.1725764) q[0];
sx q[0];
rz(-1.1541973) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7107096) q[2];
sx q[2];
rz(-2.2019349) q[2];
sx q[2];
rz(2.7106896) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.043634) q[1];
sx q[1];
rz(-2.6065738) q[1];
sx q[1];
rz(0.94162264) q[1];
x q[2];
rz(-0.97784247) q[3];
sx q[3];
rz(-2.2740318) q[3];
sx q[3];
rz(-0.11444005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82105381) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(-1.4873571) q[2];
rz(-2.1982543) q[3];
sx q[3];
rz(-0.92883674) q[3];
sx q[3];
rz(1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6134375) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(0.50233895) q[0];
rz(-2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(2.6796403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0620856) q[0];
sx q[0];
rz(-1.2794198) q[0];
sx q[0];
rz(-1.0568922) q[0];
rz(-pi) q[1];
rz(1.5080323) q[2];
sx q[2];
rz(-0.6522921) q[2];
sx q[2];
rz(1.4029274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5344055) q[1];
sx q[1];
rz(-2.7173373) q[1];
sx q[1];
rz(1.3035745) q[1];
x q[2];
rz(-1.5835207) q[3];
sx q[3];
rz(-1.0403578) q[3];
sx q[3];
rz(-0.96543559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.168557) q[2];
sx q[2];
rz(-0.23739561) q[2];
sx q[2];
rz(0.17042223) q[2];
rz(0.45525822) q[3];
sx q[3];
rz(-1.596343) q[3];
sx q[3];
rz(2.4236603) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97487226) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(-0.18044743) q[0];
rz(-2.6249053) q[1];
sx q[1];
rz(-1.5770715) q[1];
sx q[1];
rz(-2.7076941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2045) q[0];
sx q[0];
rz(-2.4334416) q[0];
sx q[0];
rz(0.88237472) q[0];
rz(-2.7942068) q[2];
sx q[2];
rz(-1.1466031) q[2];
sx q[2];
rz(-0.004384282) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5499159) q[1];
sx q[1];
rz(-0.94144758) q[1];
sx q[1];
rz(0.1484599) q[1];
rz(-2.850344) q[3];
sx q[3];
rz(-2.2403803) q[3];
sx q[3];
rz(2.254829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58069289) q[2];
sx q[2];
rz(-0.43792024) q[2];
sx q[2];
rz(1.3194579) q[2];
rz(2.8699919) q[3];
sx q[3];
rz(-1.3180132) q[3];
sx q[3];
rz(2.0297089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44529799) q[0];
sx q[0];
rz(-0.91067186) q[0];
sx q[0];
rz(1.7355504) q[0];
rz(1.7156853) q[1];
sx q[1];
rz(-1.1049263) q[1];
sx q[1];
rz(1.2474733) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919148) q[0];
sx q[0];
rz(-2.5818338) q[0];
sx q[0];
rz(2.0790529) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1389473) q[2];
sx q[2];
rz(-3.0122535) q[2];
sx q[2];
rz(2.6292588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0746324) q[1];
sx q[1];
rz(-1.6116643) q[1];
sx q[1];
rz(0.34217477) q[1];
rz(-pi) q[2];
rz(-2.5188451) q[3];
sx q[3];
rz(-2.1295665) q[3];
sx q[3];
rz(1.2947242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4282816) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(-0.10078079) q[2];
rz(0.0095602592) q[3];
sx q[3];
rz(-1.5567501) q[3];
sx q[3];
rz(1.6077707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1495001) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-2.5485582) q[1];
sx q[1];
rz(-1.7094163) q[1];
sx q[1];
rz(-0.97547668) q[1];
rz(-3.0708457) q[2];
sx q[2];
rz(-1.9683217) q[2];
sx q[2];
rz(-0.34839658) q[2];
rz(0.22105263) q[3];
sx q[3];
rz(-1.5879769) q[3];
sx q[3];
rz(-1.069608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
