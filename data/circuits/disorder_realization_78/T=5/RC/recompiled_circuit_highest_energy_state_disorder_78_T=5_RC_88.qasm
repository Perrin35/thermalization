OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5155091) q[0];
sx q[0];
rz(-2.8647487) q[0];
sx q[0];
rz(-0.11685637) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(-0.35419551) q[1];
sx q[1];
rz(0.49122223) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15165937) q[0];
sx q[0];
rz(-2.058968) q[0];
sx q[0];
rz(3.0753646) q[0];
x q[1];
rz(1.7826177) q[2];
sx q[2];
rz(-1.0928177) q[2];
sx q[2];
rz(-2.4732531) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3892606) q[1];
sx q[1];
rz(-1.6698784) q[1];
sx q[1];
rz(2.5072751) q[1];
rz(-pi) q[2];
rz(-1.9141044) q[3];
sx q[3];
rz(-1.5470389) q[3];
sx q[3];
rz(0.46059617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.7626071) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(2.911705) q[2];
rz(-0.057279438) q[3];
sx q[3];
rz(-0.66904896) q[3];
sx q[3];
rz(-2.2288286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(-0.03288658) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-2.4864973) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74467662) q[0];
sx q[0];
rz(-1.370805) q[0];
sx q[0];
rz(-2.295664) q[0];
x q[1];
rz(0.35327618) q[2];
sx q[2];
rz(-0.52918079) q[2];
sx q[2];
rz(-2.8904817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.82258525) q[1];
sx q[1];
rz(-2.2862541) q[1];
sx q[1];
rz(-2.4501178) q[1];
x q[2];
rz(0.14343393) q[3];
sx q[3];
rz(-0.8415045) q[3];
sx q[3];
rz(-2.7021331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.81755) q[2];
sx q[2];
rz(-2.5174759) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(-1.921418) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(-2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2636181) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(0.61400145) q[0];
rz(-1.7738495) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(-0.49555379) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14392631) q[0];
sx q[0];
rz(-0.94519061) q[0];
sx q[0];
rz(2.1493874) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88286384) q[2];
sx q[2];
rz(-1.9830215) q[2];
sx q[2];
rz(-2.8799873) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5835201) q[1];
sx q[1];
rz(-2.8360543) q[1];
sx q[1];
rz(-2.2347441) q[1];
x q[2];
rz(2.6747781) q[3];
sx q[3];
rz(-1.8895849) q[3];
sx q[3];
rz(-2.7960645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76501781) q[2];
sx q[2];
rz(-1.0708662) q[2];
sx q[2];
rz(-0.59256727) q[2];
rz(-2.9060034) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(3.1357646) q[0];
rz(-0.5799154) q[1];
sx q[1];
rz(-2.7858851) q[1];
sx q[1];
rz(-0.52111202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91189188) q[0];
sx q[0];
rz(-1.5661514) q[0];
sx q[0];
rz(-1.5636496) q[0];
rz(-pi) q[1];
rz(-0.02946202) q[2];
sx q[2];
rz(-2.1715487) q[2];
sx q[2];
rz(-0.0044057152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5558943) q[1];
sx q[1];
rz(-2.0221077) q[1];
sx q[1];
rz(-2.0183619) q[1];
x q[2];
rz(1.6507876) q[3];
sx q[3];
rz(-1.8820401) q[3];
sx q[3];
rz(2.1224006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8899322) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(0.31822515) q[2];
rz(-2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.4819063) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707959) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(-0.036291432) q[0];
rz(2.6151784) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(-0.2821736) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2370034) q[0];
sx q[0];
rz(-2.1891199) q[0];
sx q[0];
rz(2.8008921) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5897609) q[2];
sx q[2];
rz(-0.95317344) q[2];
sx q[2];
rz(-3.0611401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1212738) q[1];
sx q[1];
rz(-1.500638) q[1];
sx q[1];
rz(-2.4362045) q[1];
rz(2.3318841) q[3];
sx q[3];
rz(-0.96862176) q[3];
sx q[3];
rz(-0.47145999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37461773) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(3.0401163) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(-0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58979708) q[0];
sx q[0];
rz(-0.17687251) q[0];
sx q[0];
rz(-1.9675323) q[0];
rz(0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(0.039965872) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5863881) q[0];
sx q[0];
rz(-2.7551536) q[0];
sx q[0];
rz(-0.9903592) q[0];
rz(2.2117679) q[2];
sx q[2];
rz(-1.8397959) q[2];
sx q[2];
rz(-0.19942927) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3921622) q[1];
sx q[1];
rz(-0.71041965) q[1];
sx q[1];
rz(-2.602052) q[1];
rz(-pi) q[2];
rz(-1.8043936) q[3];
sx q[3];
rz(-2.3832678) q[3];
sx q[3];
rz(2.9616464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2018955) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(0.8197909) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(-1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49067295) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(-1.3167205) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(0.43693158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9428944) q[0];
sx q[0];
rz(-1.6075396) q[0];
sx q[0];
rz(1.3803287) q[0];
rz(-0.11796911) q[2];
sx q[2];
rz(-1.461339) q[2];
sx q[2];
rz(2.3326959) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.97368846) q[1];
sx q[1];
rz(-1.2120795) q[1];
sx q[1];
rz(0.78679797) q[1];
rz(-pi) q[2];
rz(-0.42610069) q[3];
sx q[3];
rz(-1.0747194) q[3];
sx q[3];
rz(-2.8471198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(0.13460049) q[2];
rz(-1.2234737) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(1.9756165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(-0.18184161) q[0];
rz(0.34135154) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(-2.9827319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0358064) q[0];
sx q[0];
rz(-1.5022359) q[0];
sx q[0];
rz(2.6935656) q[0];
x q[1];
rz(0.57639052) q[2];
sx q[2];
rz(-0.49242556) q[2];
sx q[2];
rz(-0.8053636) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6312069) q[1];
sx q[1];
rz(-0.34086168) q[1];
sx q[1];
rz(2.4426961) q[1];
rz(-pi) q[2];
rz(-2.587593) q[3];
sx q[3];
rz(-1.7007728) q[3];
sx q[3];
rz(0.25773559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86026704) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(2.2567828) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-2.6917246) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29174969) q[0];
sx q[0];
rz(-1.2717286) q[0];
sx q[0];
rz(-2.1480609) q[0];
rz(-1.7811071) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(2.5885168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9270614) q[0];
sx q[0];
rz(-1.1916516) q[0];
sx q[0];
rz(1.8168455) q[0];
rz(-pi) q[1];
rz(-3.0680066) q[2];
sx q[2];
rz(-2.4017757) q[2];
sx q[2];
rz(2.432636) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.304405) q[1];
sx q[1];
rz(-0.68301979) q[1];
sx q[1];
rz(0.4153312) q[1];
rz(2.870375) q[3];
sx q[3];
rz(-2.1925049) q[3];
sx q[3];
rz(0.78628899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20455827) q[2];
sx q[2];
rz(-0.60595804) q[2];
sx q[2];
rz(-2.9045612) q[2];
rz(1.6009181) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(-1.1698394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071526214) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(-2.7023244) q[1];
sx q[1];
rz(-1.2989137) q[1];
sx q[1];
rz(-3.0781504) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251415) q[0];
sx q[0];
rz(-2.2277895) q[0];
sx q[0];
rz(-2.444066) q[0];
rz(0.62427743) q[2];
sx q[2];
rz(-1.4194402) q[2];
sx q[2];
rz(-2.4454012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52229228) q[1];
sx q[1];
rz(-1.5623153) q[1];
sx q[1];
rz(-1.897476) q[1];
rz(-0.88313734) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(-2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0110317) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(1.4981184) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-1.1199718) q[3];
sx q[3];
rz(1.2137265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7508004) q[0];
sx q[0];
rz(-1.1815429) q[0];
sx q[0];
rz(-1.0829157) q[0];
rz(-2.2687268) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(-0.93449705) q[2];
sx q[2];
rz(-0.9899944) q[2];
sx q[2];
rz(-1.8805671) q[2];
rz(-1.6428236) q[3];
sx q[3];
rz(-0.99456064) q[3];
sx q[3];
rz(-1.4069788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
