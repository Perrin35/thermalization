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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1304162) q[0];
sx q[0];
rz(-2.6493086) q[0];
sx q[0];
rz(-1.4468132) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7826177) q[2];
sx q[2];
rz(-2.048775) q[2];
sx q[2];
rz(-0.66833959) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3892606) q[1];
sx q[1];
rz(-1.6698784) q[1];
sx q[1];
rz(-2.5072751) q[1];
rz(-pi) q[2];
rz(1.5003203) q[3];
sx q[3];
rz(-0.34409663) q[3];
sx q[3];
rz(-1.1765574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3789856) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(0.22988764) q[2];
rz(3.0843132) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(2.2288286) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(3.1087061) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-2.4864973) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0467619) q[0];
sx q[0];
rz(-0.74709409) q[0];
sx q[0];
rz(-1.8674891) q[0];
rz(-pi) q[1];
rz(0.50184558) q[2];
sx q[2];
rz(-1.7463533) q[2];
sx q[2];
rz(1.6278536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8959143) q[1];
sx q[1];
rz(-2.0729471) q[1];
sx q[1];
rz(-2.4163626) q[1];
rz(1.412185) q[3];
sx q[3];
rz(-2.4008823) q[3];
sx q[3];
rz(-0.22601688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(1.921418) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(0.65742457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87797457) q[0];
sx q[0];
rz(-0.89732301) q[0];
sx q[0];
rz(-0.61400145) q[0];
rz(-1.3677431) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(-0.49555379) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6959636) q[0];
sx q[0];
rz(-0.82484748) q[0];
sx q[0];
rz(0.64795171) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88286384) q[2];
sx q[2];
rz(-1.9830215) q[2];
sx q[2];
rz(2.8799873) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0125791) q[1];
sx q[1];
rz(-1.3316175) q[1];
sx q[1];
rz(0.19197477) q[1];
rz(-pi) q[2];
rz(0.63276799) q[3];
sx q[3];
rz(-2.583021) q[3];
sx q[3];
rz(-2.4726923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3765748) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(2.5490254) q[2];
rz(2.9060034) q[3];
sx q[3];
rz(-0.75747907) q[3];
sx q[3];
rz(0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(0.0058280514) q[0];
rz(-2.5616772) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(2.6204806) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2297008) q[0];
sx q[0];
rz(-1.5661514) q[0];
sx q[0];
rz(-1.5779431) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1121306) q[2];
sx q[2];
rz(-2.1715487) q[2];
sx q[2];
rz(0.0044057152) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8940034) q[1];
sx q[1];
rz(-0.62452468) q[1];
sx q[1];
rz(2.4127059) q[1];
x q[2];
rz(0.24346015) q[3];
sx q[3];
rz(-2.820558) q[3];
sx q[3];
rz(-2.3784172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8899322) q[2];
sx q[2];
rz(-2.6106788) q[2];
sx q[2];
rz(-0.31822515) q[2];
rz(2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.6596863) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707959) q[0];
sx q[0];
rz(-2.9279121) q[0];
sx q[0];
rz(3.1053012) q[0];
rz(0.52641422) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(0.2821736) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68809915) q[0];
sx q[0];
rz(-0.69506139) q[0];
sx q[0];
rz(2.009925) q[0];
rz(-pi) q[1];
rz(0.93499489) q[2];
sx q[2];
rz(-2.3381669) q[2];
sx q[2];
rz(0.73551169) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6101004) q[1];
sx q[1];
rz(-0.86750114) q[1];
sx q[1];
rz(1.6628357) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7871233) q[3];
sx q[3];
rz(-0.93138444) q[3];
sx q[3];
rz(1.5057664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37461773) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(-0.96595079) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(-0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.9675323) q[0];
rz(0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(0.039965872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.171282) q[0];
sx q[0];
rz(-1.2501646) q[0];
sx q[0];
rz(-0.21954222) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8103175) q[2];
sx q[2];
rz(-2.185198) q[2];
sx q[2];
rz(1.5744874) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.74943042) q[1];
sx q[1];
rz(-0.71041965) q[1];
sx q[1];
rz(-2.602052) q[1];
rz(-pi) q[2];
rz(2.3153855) q[3];
sx q[3];
rz(-1.4109269) q[3];
sx q[3];
rz(-1.9217971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2018955) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(2.3218018) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6509197) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.9098234) q[1];
sx q[1];
rz(0.43693158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9428944) q[0];
sx q[0];
rz(-1.6075396) q[0];
sx q[0];
rz(-1.7612639) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11796911) q[2];
sx q[2];
rz(-1.461339) q[2];
sx q[2];
rz(-2.3326959) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1679042) q[1];
sx q[1];
rz(-1.2120795) q[1];
sx q[1];
rz(0.78679797) q[1];
rz(-pi) q[2];
rz(-0.42610069) q[3];
sx q[3];
rz(-2.0668732) q[3];
sx q[3];
rz(2.8471198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4009319) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(-1.918119) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80158919) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(-0.18184161) q[0];
rz(0.34135154) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-0.15886074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0358064) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(-0.44802702) q[0];
x q[1];
rz(1.2863289) q[2];
sx q[2];
rz(-1.9783696) q[2];
sx q[2];
rz(1.4409232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73028683) q[1];
sx q[1];
rz(-1.3540233) q[1];
sx q[1];
rz(-0.2651528) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4182866) q[3];
sx q[3];
rz(-2.1195863) q[3];
sx q[3];
rz(-1.2330518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.86026704) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(-0.88480985) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-2.116674) q[3];
sx q[3];
rz(2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29174969) q[0];
sx q[0];
rz(-1.2717286) q[0];
sx q[0];
rz(2.1480609) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(-2.5885168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21453126) q[0];
sx q[0];
rz(-1.1916516) q[0];
sx q[0];
rz(-1.3247471) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.07358609) q[2];
sx q[2];
rz(-2.4017757) q[2];
sx q[2];
rz(0.70895665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.40402624) q[1];
sx q[1];
rz(-1.3133004) q[1];
sx q[1];
rz(2.5016258) q[1];
rz(-pi) q[2];
rz(-0.93135466) q[3];
sx q[3];
rz(-1.7903312) q[3];
sx q[3];
rz(2.1965248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-2.9045612) q[2];
rz(1.6009181) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(-1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071526214) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(0.1189098) q[0];
rz(2.7023244) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-3.0781504) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3812696) q[0];
sx q[0];
rz(-1.0372235) q[0];
sx q[0];
rz(-0.78223451) q[0];
rz(-pi) q[1];
rz(1.3849867) q[2];
sx q[2];
rz(-2.1868621) q[2];
sx q[2];
rz(-2.1587929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0234785) q[1];
sx q[1];
rz(-2.8148068) q[1];
sx q[1];
rz(-1.5972196) q[1];
x q[2];
rz(2.2584553) q[3];
sx q[3];
rz(-0.68810191) q[3];
sx q[3];
rz(-2.7795252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1305609) q[2];
sx q[2];
rz(-2.06879) q[2];
sx q[2];
rz(-1.6434742) q[2];
rz(-0.9324075) q[3];
sx q[3];
rz(-1.1199718) q[3];
sx q[3];
rz(-1.9278661) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-0.73578655) q[2];
sx q[2];
rz(-0.83333165) q[2];
sx q[2];
rz(0.32902645) q[2];
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
