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
rz(-1.6260835) q[0];
sx q[0];
rz(-0.27684394) q[0];
sx q[0];
rz(-3.0247363) q[0];
rz(2.6729743) q[1];
sx q[1];
rz(-0.35419551) q[1];
sx q[1];
rz(-2.6503704) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1304162) q[0];
sx q[0];
rz(-2.6493086) q[0];
sx q[0];
rz(1.4468132) q[0];
rz(-0.38552706) q[2];
sx q[2];
rz(-0.51947278) q[2];
sx q[2];
rz(-2.9105718) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1894816) q[1];
sx q[1];
rz(-0.64095381) q[1];
sx q[1];
rz(-0.16619311) q[1];
rz(0.025229021) q[3];
sx q[3];
rz(-1.9140035) q[3];
sx q[3];
rz(-2.0398839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3789856) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(0.22988764) q[2];
rz(0.057279438) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(-2.2288286) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(2.4864973) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0948308) q[0];
sx q[0];
rz(-0.74709409) q[0];
sx q[0];
rz(1.8674891) q[0];
rz(-pi) q[1];
rz(2.6397471) q[2];
sx q[2];
rz(-1.3952394) q[2];
sx q[2];
rz(1.6278536) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0634708) q[1];
sx q[1];
rz(-0.95032108) q[1];
sx q[1];
rz(-2.2038151) q[1];
rz(0.14343393) q[3];
sx q[3];
rz(-2.3000882) q[3];
sx q[3];
rz(-0.4394596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32404262) q[2];
sx q[2];
rz(-2.5174759) q[2];
sx q[2];
rz(-1.2980655) q[2];
rz(-1.2201747) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(-2.4841681) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2636181) q[0];
sx q[0];
rz(-0.89732301) q[0];
sx q[0];
rz(0.61400145) q[0];
rz(-1.3677431) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(-2.6460389) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4456291) q[0];
sx q[0];
rz(-0.82484748) q[0];
sx q[0];
rz(0.64795171) q[0];
rz(-pi) q[1];
rz(-0.96769963) q[2];
sx q[2];
rz(-2.357238) q[2];
sx q[2];
rz(1.3788144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6537955) q[1];
sx q[1];
rz(-1.7572409) q[1];
sx q[1];
rz(-1.814278) q[1];
rz(-pi) q[2];
rz(0.63276799) q[3];
sx q[3];
rz(-2.583021) q[3];
sx q[3];
rz(-2.4726923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76501781) q[2];
sx q[2];
rz(-1.0708662) q[2];
sx q[2];
rz(-0.59256727) q[2];
rz(-0.23558922) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(-2.6872915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(-0.0058280514) q[0];
rz(2.5616772) q[1];
sx q[1];
rz(-2.7858851) q[1];
sx q[1];
rz(-0.52111202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91189188) q[0];
sx q[0];
rz(-1.5754412) q[0];
sx q[0];
rz(1.5636496) q[0];
rz(-pi) q[1];
rz(-3.1121306) q[2];
sx q[2];
rz(-2.1715487) q[2];
sx q[2];
rz(-3.1371869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5558943) q[1];
sx q[1];
rz(-2.0221077) q[1];
sx q[1];
rz(-1.1232308) q[1];
rz(-pi) q[2];
rz(-1.6507876) q[3];
sx q[3];
rz(-1.8820401) q[3];
sx q[3];
rz(1.019192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(2.8233675) q[2];
rz(-0.71277726) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37079674) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(3.1053012) q[0];
rz(-2.6151784) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(0.2821736) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68809915) q[0];
sx q[0];
rz(-0.69506139) q[0];
sx q[0];
rz(-1.1316677) q[0];
rz(-pi) q[1];
rz(-0.93499489) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(0.73551169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5314922) q[1];
sx q[1];
rz(-2.2740915) q[1];
sx q[1];
rz(1.478757) q[1];
x q[2];
rz(0.7871233) q[3];
sx q[3];
rz(-0.93138444) q[3];
sx q[3];
rz(-1.6358262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7669749) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(-0.10147632) q[2];
rz(2.1756419) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(3.0275893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58979708) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(-1.9675323) q[0];
rz(-2.7575098) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(-0.039965872) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9703107) q[0];
sx q[0];
rz(-1.2501646) q[0];
sx q[0];
rz(2.9220504) q[0];
rz(2.8103175) q[2];
sx q[2];
rz(-0.9563947) q[2];
sx q[2];
rz(-1.5671052) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0607073) q[1];
sx q[1];
rz(-0.97700143) q[1];
sx q[1];
rz(-1.9869366) q[1];
x q[2];
rz(-0.82620718) q[3];
sx q[3];
rz(-1.7306657) q[3];
sx q[3];
rz(1.9217971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93969718) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(2.3218018) q[2];
rz(2.36006) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509197) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(-1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(2.7046611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36501554) q[0];
sx q[0];
rz(-1.3804589) q[0];
sx q[0];
rz(0.037419407) q[0];
rz(-pi) q[1];
rz(-0.11796911) q[2];
sx q[2];
rz(-1.6802537) q[2];
sx q[2];
rz(-2.3326959) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26035151) q[1];
sx q[1];
rz(-0.8484183) q[1];
sx q[1];
rz(0.48697014) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1070621) q[3];
sx q[3];
rz(-1.1987743) q[3];
sx q[3];
rz(-2.0780502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.74066073) q[2];
sx q[2];
rz(-2.0621982) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(-1.9756165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(0.18184161) q[0];
rz(2.8002411) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-2.9827319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358064) q[0];
sx q[0];
rz(-1.5022359) q[0];
sx q[0];
rz(-2.6935656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57639052) q[2];
sx q[2];
rz(-2.6491671) q[2];
sx q[2];
rz(0.8053636) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.78217172) q[1];
sx q[1];
rz(-1.3119932) q[1];
sx q[1];
rz(1.3464298) q[1];
rz(-pi) q[2];
rz(0.24352945) q[3];
sx q[3];
rz(-2.5741044) q[3];
sx q[3];
rz(1.6220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2813256) q[2];
sx q[2];
rz(-0.93026668) q[2];
sx q[2];
rz(-0.88480985) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.849843) q[0];
sx q[0];
rz(-1.2717286) q[0];
sx q[0];
rz(-0.9935317) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(-2.5885168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9270614) q[0];
sx q[0];
rz(-1.1916516) q[0];
sx q[0];
rz(1.8168455) q[0];
x q[1];
rz(-3.0680066) q[2];
sx q[2];
rz(-0.739817) q[2];
sx q[2];
rz(-2.432636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.304405) q[1];
sx q[1];
rz(-0.68301979) q[1];
sx q[1];
rz(-0.4153312) q[1];
x q[2];
rz(2.210238) q[3];
sx q[3];
rz(-1.7903312) q[3];
sx q[3];
rz(2.1965248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-0.23703144) q[2];
rz(-1.5406746) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(1.1698394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.071526214) q[0];
sx q[0];
rz(-1.5109753) q[0];
sx q[0];
rz(0.1189098) q[0];
rz(2.7023244) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-3.0781504) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251415) q[0];
sx q[0];
rz(-2.2277895) q[0];
sx q[0];
rz(-2.444066) q[0];
rz(-pi) q[1];
rz(-0.25524885) q[2];
sx q[2];
rz(-2.5016157) q[2];
sx q[2];
rz(2.4733123) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1181142) q[1];
sx q[1];
rz(-0.32678582) q[1];
sx q[1];
rz(-1.5972196) q[1];
x q[2];
rz(-2.2584553) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.1305609) q[2];
sx q[2];
rz(-1.0728027) q[2];
sx q[2];
rz(-1.4981184) q[2];
rz(0.9324075) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(-1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(-0.87286585) q[1];
sx q[1];
rz(-2.0582336) q[1];
sx q[1];
rz(2.2179926) q[1];
rz(-2.4058061) q[2];
sx q[2];
rz(-2.308261) q[2];
sx q[2];
rz(-2.8125662) q[2];
rz(1.4987691) q[3];
sx q[3];
rz(-0.99456064) q[3];
sx q[3];
rz(-1.4069788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
