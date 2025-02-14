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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79345771) q[0];
sx q[0];
rz(-2.0525816) q[0];
sx q[0];
rz(-0.30695494) q[0];
x q[1];
rz(-1.3363935) q[2];
sx q[2];
rz(-1.3318828) q[2];
sx q[2];
rz(2.6221681) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8227777) q[1];
sx q[1];
rz(-2.9479369) q[1];
sx q[1];
rz(-2.6050287) q[1];
rz(1.9735967) q[3];
sx q[3];
rz(-2.9514545) q[3];
sx q[3];
rz(-0.26241747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7370558) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-0.36049584) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-2.093061) q[3];
sx q[3];
rz(1.8726965) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8812113) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(0.65688175) q[0];
rz(-2.8086713) q[1];
sx q[1];
rz(-1.0924783) q[1];
sx q[1];
rz(2.2064256) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2774076) q[0];
sx q[0];
rz(-1.5404697) q[0];
sx q[0];
rz(-3.0369989) q[0];
rz(-pi) q[1];
rz(-0.33719535) q[2];
sx q[2];
rz(-2.8697578) q[2];
sx q[2];
rz(-1.7295009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1166653) q[1];
sx q[1];
rz(-1.6102734) q[1];
sx q[1];
rz(-2.9945606) q[1];
rz(-pi) q[2];
rz(2.6245429) q[3];
sx q[3];
rz(-2.73797) q[3];
sx q[3];
rz(2.0832186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8344581) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(1.2163986) q[2];
rz(0.52792102) q[3];
sx q[3];
rz(-1.0390176) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5288178) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(-0.58498996) q[0];
rz(-0.12475573) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(1.6927208) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0199658) q[0];
sx q[0];
rz(-0.0040071132) q[0];
sx q[0];
rz(0.83929707) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5572433) q[2];
sx q[2];
rz(-1.4986191) q[2];
sx q[2];
rz(-0.37301317) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53939143) q[1];
sx q[1];
rz(-2.3219206) q[1];
sx q[1];
rz(-0.96266268) q[1];
rz(1.307906) q[3];
sx q[3];
rz(-1.2440727) q[3];
sx q[3];
rz(-0.19858352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8831732) q[2];
sx q[2];
rz(-1.0018188) q[2];
sx q[2];
rz(1.4460571) q[2];
rz(2.3840733) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(-2.3022046) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3676753) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(1.0876592) q[0];
rz(2.7663973) q[1];
sx q[1];
rz(-1.1253858) q[1];
sx q[1];
rz(2.1813724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2743234) q[0];
sx q[0];
rz(-1.835863) q[0];
sx q[0];
rz(-0.041170711) q[0];
x q[1];
rz(-1.4713418) q[2];
sx q[2];
rz(-1.6005472) q[2];
sx q[2];
rz(-3.0029675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.580048) q[1];
sx q[1];
rz(-0.64148884) q[1];
sx q[1];
rz(-0.80342355) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4750911) q[3];
sx q[3];
rz(-1.8391575) q[3];
sx q[3];
rz(2.1398422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20597657) q[2];
sx q[2];
rz(-1.2674067) q[2];
sx q[2];
rz(1.4975632) q[2];
rz(0.86132541) q[3];
sx q[3];
rz(-0.70278168) q[3];
sx q[3];
rz(-2.6845045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1645666) q[0];
sx q[0];
rz(-2.470546) q[0];
sx q[0];
rz(-0.60428756) q[0];
rz(1.9970278) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(2.7191275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10642636) q[0];
sx q[0];
rz(-1.4894333) q[0];
sx q[0];
rz(-1.7398796) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1429135) q[2];
sx q[2];
rz(-1.0829751) q[2];
sx q[2];
rz(-1.7833379) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14478165) q[1];
sx q[1];
rz(-1.4986804) q[1];
sx q[1];
rz(2.4516289) q[1];
x q[2];
rz(-2.4147968) q[3];
sx q[3];
rz(-2.6495547) q[3];
sx q[3];
rz(2.7068822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8233238) q[2];
sx q[2];
rz(-2.4906929) q[2];
sx q[2];
rz(-0.65625119) q[2];
rz(1.5308135) q[3];
sx q[3];
rz(-1.8717513) q[3];
sx q[3];
rz(2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(2.5166125) q[0];
rz(2.3896353) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(1.5350852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4055734) q[0];
sx q[0];
rz(-1.6200119) q[0];
sx q[0];
rz(0.28926138) q[0];
x q[1];
rz(-0.70930945) q[2];
sx q[2];
rz(-1.0934208) q[2];
sx q[2];
rz(-0.8698744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5175486) q[1];
sx q[1];
rz(-2.4772948) q[1];
sx q[1];
rz(2.7655927) q[1];
rz(-pi) q[2];
rz(0.97134892) q[3];
sx q[3];
rz(-2.0347715) q[3];
sx q[3];
rz(-2.6148207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.357888) q[2];
sx q[2];
rz(0.9291741) q[2];
rz(-2.3164228) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(-2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55649844) q[0];
sx q[0];
rz(-0.80699054) q[0];
sx q[0];
rz(-1.3744542) q[0];
rz(-2.6311686) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(2.5183831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685112) q[0];
sx q[0];
rz(-0.36751908) q[0];
sx q[0];
rz(-1.1373684) q[0];
rz(-1.7665777) q[2];
sx q[2];
rz(-2.0718241) q[2];
sx q[2];
rz(2.632338) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0556866) q[1];
sx q[1];
rz(-1.6821563) q[1];
sx q[1];
rz(1.7176877) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3939952) q[3];
sx q[3];
rz(-1.7204525) q[3];
sx q[3];
rz(-0.29779321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0509384) q[2];
sx q[2];
rz(-0.81523681) q[2];
sx q[2];
rz(1.9742924) q[2];
rz(-3.1189611) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(1.1289271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8959494) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(-1.0173215) q[0];
rz(1.7474878) q[1];
sx q[1];
rz(-0.21109763) q[1];
sx q[1];
rz(0.62754935) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59333364) q[0];
sx q[0];
rz(-2.152266) q[0];
sx q[0];
rz(2.6618746) q[0];
x q[1];
rz(1.4862138) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(-1.6135474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4269451) q[1];
sx q[1];
rz(-1.8083739) q[1];
sx q[1];
rz(-1.3690884) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47431176) q[3];
sx q[3];
rz(-1.2657058) q[3];
sx q[3];
rz(2.4449789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5158186) q[2];
sx q[2];
rz(-2.5067582) q[2];
sx q[2];
rz(0.41075692) q[2];
rz(-3.0913894) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5589767) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(0.26891747) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(3.0115829) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0144671) q[0];
sx q[0];
rz(-1.7916745) q[0];
sx q[0];
rz(1.2400024) q[0];
x q[1];
rz(2.7902044) q[2];
sx q[2];
rz(-1.2341502) q[2];
sx q[2];
rz(-1.2996246) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0669032) q[1];
sx q[1];
rz(-0.58614158) q[1];
sx q[1];
rz(2.8476004) q[1];
x q[2];
rz(-0.49874108) q[3];
sx q[3];
rz(-0.95512701) q[3];
sx q[3];
rz(-1.2241253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93398634) q[2];
sx q[2];
rz(-0.76675582) q[2];
sx q[2];
rz(-0.096573528) q[2];
rz(1.5038331) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12829517) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(0.53303322) q[0];
rz(-0.036272613) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(1.4520377) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55646) q[0];
sx q[0];
rz(-2.1341354) q[0];
sx q[0];
rz(0.4507395) q[0];
rz(-pi) q[1];
rz(0.31970892) q[2];
sx q[2];
rz(-2.8065348) q[2];
sx q[2];
rz(-0.31949319) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.328985) q[1];
sx q[1];
rz(-1.9300864) q[1];
sx q[1];
rz(-1.6027662) q[1];
x q[2];
rz(0.17590268) q[3];
sx q[3];
rz(-2.6697192) q[3];
sx q[3];
rz(-0.37341082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71296802) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(0.0326322) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(-1.0802065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7964771) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(0.79700094) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(0.5354698) q[2];
sx q[2];
rz(-1.3146123) q[2];
sx q[2];
rz(-1.8563369) q[2];
rz(-0.81099323) q[3];
sx q[3];
rz(-1.4832433) q[3];
sx q[3];
rz(1.087838) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
