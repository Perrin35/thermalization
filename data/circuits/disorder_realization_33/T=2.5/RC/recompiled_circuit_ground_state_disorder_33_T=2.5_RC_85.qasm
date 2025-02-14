OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4274974) q[0];
sx q[0];
rz(-0.56199718) q[0];
sx q[0];
rz(-2.9105817) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(2.6872771) q[1];
sx q[1];
rz(8.1375577) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92316931) q[0];
sx q[0];
rz(-1.8418663) q[0];
sx q[0];
rz(2.0725033) q[0];
x q[1];
rz(-1.8051992) q[2];
sx q[2];
rz(-1.8097098) q[2];
sx q[2];
rz(2.6221681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8227777) q[1];
sx q[1];
rz(-2.9479369) q[1];
sx q[1];
rz(-2.6050287) q[1];
x q[2];
rz(-3.0662905) q[3];
sx q[3];
rz(-1.7455532) q[3];
sx q[3];
rz(2.99461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7370558) q[2];
sx q[2];
rz(-0.70713592) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(2.0170085) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(-1.2688961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8812113) q[0];
sx q[0];
rz(-2.9660048) q[0];
sx q[0];
rz(-0.65688175) q[0];
rz(-0.33292133) q[1];
sx q[1];
rz(-1.0924783) q[1];
sx q[1];
rz(0.93516707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2774076) q[0];
sx q[0];
rz(-1.5404697) q[0];
sx q[0];
rz(0.10459374) q[0];
x q[1];
rz(1.4788394) q[2];
sx q[2];
rz(-1.3146245) q[2];
sx q[2];
rz(1.0630449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8063025) q[1];
sx q[1];
rz(-2.9893901) q[1];
sx q[1];
rz(0.2633414) q[1];
rz(-1.7788497) q[3];
sx q[3];
rz(-1.2223772) q[3];
sx q[3];
rz(-1.6121685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3071345) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(-1.9251941) q[2];
rz(-2.6136716) q[3];
sx q[3];
rz(-2.1025751) q[3];
sx q[3];
rz(-1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6127748) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(0.58498996) q[0];
rz(-3.0168369) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(-1.6927208) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3901236) q[0];
sx q[0];
rz(-1.5737783) q[0];
sx q[0];
rz(3.1389159) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0694088) q[2];
sx q[2];
rz(-1.584314) q[2];
sx q[2];
rz(-1.9428321) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3347149) q[1];
sx q[1];
rz(-2.2141465) q[1];
sx q[1];
rz(-2.592464) q[1];
rz(2.8040941) q[3];
sx q[3];
rz(-1.8194767) q[3];
sx q[3];
rz(-1.6832222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2584194) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(1.6955356) q[2];
rz(2.3840733) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(0.83938804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7739173) q[0];
sx q[0];
rz(-0.92905074) q[0];
sx q[0];
rz(-2.0539334) q[0];
rz(0.37519535) q[1];
sx q[1];
rz(-1.1253858) q[1];
sx q[1];
rz(0.96022022) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1183557) q[0];
sx q[0];
rz(-2.873422) q[0];
sx q[0];
rz(-1.7212746) q[0];
rz(-pi) q[1];
rz(1.2795936) q[2];
sx q[2];
rz(-0.10379496) q[2];
sx q[2];
rz(-1.7218931) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31611495) q[1];
sx q[1];
rz(-2.0160455) q[1];
sx q[1];
rz(0.47834217) q[1];
x q[2];
rz(1.6665015) q[3];
sx q[3];
rz(-1.8391575) q[3];
sx q[3];
rz(-1.0017504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20597657) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(1.4975632) q[2];
rz(-2.2802672) q[3];
sx q[3];
rz(-0.70278168) q[3];
sx q[3];
rz(-2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(2.5373051) q[0];
rz(1.1445649) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(0.42246517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9086994) q[0];
sx q[0];
rz(-2.9541203) q[0];
sx q[0];
rz(-2.0220246) q[0];
rz(-pi) q[1];
rz(2.3461012) q[2];
sx q[2];
rz(-2.4078712) q[2];
sx q[2];
rz(-0.41662859) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3389272) q[1];
sx q[1];
rz(-0.69310729) q[1];
sx q[1];
rz(0.1130123) q[1];
x q[2];
rz(1.2286387) q[3];
sx q[3];
rz(-1.2099724) q[3];
sx q[3];
rz(2.7865041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8233238) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(2.4853415) q[2];
rz(1.6107791) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.4055279) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(0.62498012) q[0];
rz(2.3896353) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(1.5350852) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4706375) q[0];
sx q[0];
rz(-2.8482901) q[0];
sx q[0];
rz(-0.17099149) q[0];
x q[1];
rz(-0.70930945) q[2];
sx q[2];
rz(-2.0481718) q[2];
sx q[2];
rz(-2.2717182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8934763) q[1];
sx q[1];
rz(-1.3424338) q[1];
sx q[1];
rz(-0.62947692) q[1];
x q[2];
rz(-0.97134892) q[3];
sx q[3];
rz(-1.1068212) q[3];
sx q[3];
rz(0.52677192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.357888) q[2];
sx q[2];
rz(-2.2124186) q[2];
rz(0.82516986) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5850942) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(-1.7671385) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-2.3511032) q[1];
sx q[1];
rz(-2.5183831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0334367) q[0];
sx q[0];
rz(-1.9029473) q[0];
sx q[0];
rz(-2.9812814) q[0];
x q[1];
rz(2.8002732) q[2];
sx q[2];
rz(-0.53487294) q[2];
sx q[2];
rz(-0.90082263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53155073) q[1];
sx q[1];
rz(-1.4248214) q[1];
sx q[1];
rz(-3.0290305) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74759746) q[3];
sx q[3];
rz(-1.4211402) q[3];
sx q[3];
rz(0.29779321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0906543) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(-1.9742924) q[2];
rz(3.1189611) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(-2.0126655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24564329) q[0];
sx q[0];
rz(-1.1351981) q[0];
sx q[0];
rz(-2.1242712) q[0];
rz(-1.3941049) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(-0.62754935) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16495569) q[0];
sx q[0];
rz(-2.4058488) q[0];
sx q[0];
rz(-0.95860211) q[0];
x q[1];
rz(-1.6553788) q[2];
sx q[2];
rz(-2.248317) q[2];
sx q[2];
rz(1.6135474) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7146475) q[1];
sx q[1];
rz(-1.8083739) q[1];
sx q[1];
rz(1.7725043) q[1];
rz(-pi) q[2];
rz(-0.47431176) q[3];
sx q[3];
rz(-1.2657058) q[3];
sx q[3];
rz(0.69661372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5158186) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(-0.41075692) q[2];
rz(-0.050203236) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589767) q[0];
sx q[0];
rz(-2.2447383) q[0];
sx q[0];
rz(-0.26891747) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-2.0545484) q[1];
sx q[1];
rz(-3.0115829) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63142473) q[0];
sx q[0];
rz(-1.893259) q[0];
sx q[0];
rz(-2.9084951) q[0];
rz(-pi) q[1];
x q[1];
rz(1.214005) q[2];
sx q[2];
rz(-1.901682) q[2];
sx q[2];
rz(-2.7499104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0669032) q[1];
sx q[1];
rz(-0.58614158) q[1];
sx q[1];
rz(-2.8476004) q[1];
rz(2.2488908) q[3];
sx q[3];
rz(-1.9719651) q[3];
sx q[3];
rz(2.490171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93398634) q[2];
sx q[2];
rz(-2.3748368) q[2];
sx q[2];
rz(-3.0450191) q[2];
rz(-1.5038331) q[3];
sx q[3];
rz(-1.3373172) q[3];
sx q[3];
rz(-0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12829517) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(0.53303322) q[0];
rz(3.10532) q[1];
sx q[1];
rz(-0.78250042) q[1];
sx q[1];
rz(-1.4520377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9617758) q[0];
sx q[0];
rz(-0.70588934) q[0];
sx q[0];
rz(0.96700661) q[0];
rz(-1.6797941) q[2];
sx q[2];
rz(-1.2533292) q[2];
sx q[2];
rz(-0.017680971) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2382792) q[1];
sx q[1];
rz(-2.7809445) q[1];
sx q[1];
rz(-3.0566932) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6598654) q[3];
sx q[3];
rz(-1.1067821) q[3];
sx q[3];
rz(2.9651412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71296802) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(3.1089605) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7964771) q[0];
sx q[0];
rz(-0.65237541) q[0];
sx q[0];
rz(2.8271578) q[0];
rz(-2.3445917) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-0.47427872) q[2];
sx q[2];
rz(-2.5534292) q[2];
sx q[2];
rz(-3.0234887) q[2];
rz(2.3305994) q[3];
sx q[3];
rz(-1.4832433) q[3];
sx q[3];
rz(1.087838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
