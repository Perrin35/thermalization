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
rz(0.49122223) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15165937) q[0];
sx q[0];
rz(-2.058968) q[0];
sx q[0];
rz(-3.0753646) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38552706) q[2];
sx q[2];
rz(-0.51947278) q[2];
sx q[2];
rz(-0.2310209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95211103) q[1];
sx q[1];
rz(-2.5006388) q[1];
sx q[1];
rz(-2.9753995) q[1];
x q[2];
rz(1.6412723) q[3];
sx q[3];
rz(-0.34409663) q[3];
sx q[3];
rz(1.1765574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3789856) q[2];
sx q[2];
rz(-0.28756046) q[2];
sx q[2];
rz(-2.911705) q[2];
rz(-3.0843132) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(-2.2288286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3463992) q[0];
sx q[0];
rz(-2.754358) q[0];
sx q[0];
rz(-2.9246395) q[0];
rz(-3.1087061) q[1];
sx q[1];
rz(-0.93371987) q[1];
sx q[1];
rz(-0.65509534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74467662) q[0];
sx q[0];
rz(-1.370805) q[0];
sx q[0];
rz(2.295664) q[0];
rz(1.7704324) q[2];
sx q[2];
rz(-2.0642274) q[2];
sx q[2];
rz(-0.1525998) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2456783) q[1];
sx q[1];
rz(-1.0686456) q[1];
sx q[1];
rz(2.4163626) q[1];
x q[2];
rz(0.14343393) q[3];
sx q[3];
rz(-2.3000882) q[3];
sx q[3];
rz(2.7021331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(1.2980655) q[2];
rz(1.921418) q[3];
sx q[3];
rz(-1.4379359) q[3];
sx q[3];
rz(-2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2636181) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(-2.5275912) q[0];
rz(1.3677431) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(-2.6460389) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9976663) q[0];
sx q[0];
rz(-2.196402) q[0];
sx q[0];
rz(2.1493874) q[0];
rz(-pi) q[1];
rz(0.51505417) q[2];
sx q[2];
rz(-2.1916789) q[2];
sx q[2];
rz(-2.1504998) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0125791) q[1];
sx q[1];
rz(-1.8099752) q[1];
sx q[1];
rz(-0.19197477) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2167778) q[3];
sx q[3];
rz(-1.1292158) q[3];
sx q[3];
rz(1.3819249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76501781) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(0.59256727) q[2];
rz(2.9060034) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(2.6872915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4950824) q[0];
sx q[0];
rz(-0.80146924) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(-0.5799154) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(0.52111202) q[1];
rz(pi/2) q[2];
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
rz(0.96984152) q[2];
sx q[2];
rz(-1.5950987) q[2];
sx q[2];
rz(-1.5497335) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.191491) q[1];
sx q[1];
rz(-1.9708212) q[1];
sx q[1];
rz(-2.6482976) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31217869) q[3];
sx q[3];
rz(-1.6469367) q[3];
sx q[3];
rz(-2.5654441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-0.53091383) q[2];
sx q[2];
rz(-0.31822515) q[2];
rz(-2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.4819063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7707959) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(3.1053012) q[0];
rz(-0.52641422) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(2.859419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53647868) q[0];
sx q[0];
rz(-1.8465586) q[0];
sx q[0];
rz(-0.92425297) q[0];
rz(-0.93499489) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(0.73551169) q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
rz(2.3544694) q[3];
sx q[3];
rz(-2.2102082) q[3];
sx q[3];
rz(1.5057664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37461773) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(3.0401163) q[2];
rz(2.1756419) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46972457) q[0];
sx q[0];
rz(-1.7789808) q[0];
sx q[0];
rz(-1.8987659) q[0];
rz(1.1388117) q[2];
sx q[2];
rz(-0.68772763) q[2];
sx q[2];
rz(1.0291531) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2474976) q[1];
sx q[1];
rz(-1.9124417) q[1];
sx q[1];
rz(-2.5057807) q[1];
rz(-pi) q[2];
rz(-2.9257366) q[3];
sx q[3];
rz(-2.3036968) q[3];
sx q[3];
rz(2.6449316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93969718) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(-2.3218018) q[2];
rz(-2.36006) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6509197) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(-1.3680178) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(2.7046611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56037134) q[0];
sx q[0];
rz(-0.19393714) q[0];
sx q[0];
rz(1.3790129) q[0];
rz(-pi) q[1];
rz(-0.75114366) q[2];
sx q[2];
rz(-0.1607543) q[2];
sx q[2];
rz(-3.1243665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2059784) q[1];
sx q[1];
rz(-2.2956478) q[1];
sx q[1];
rz(-2.0589252) q[1];
rz(1.0345305) q[3];
sx q[3];
rz(-1.1987743) q[3];
sx q[3];
rz(1.0635424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4009319) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(-3.0069922) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(1.9756165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(-2.959751) q[0];
rz(2.8002411) q[1];
sx q[1];
rz(-2.0149442) q[1];
sx q[1];
rz(2.9827319) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5349992) q[0];
sx q[0];
rz(-0.45289055) q[0];
sx q[0];
rz(2.9843828) q[0];
rz(-2.5652021) q[2];
sx q[2];
rz(-2.6491671) q[2];
sx q[2];
rz(0.8053636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3594209) q[1];
sx q[1];
rz(-1.8295994) q[1];
sx q[1];
rz(-1.3464298) q[1];
rz(0.24352945) q[3];
sx q[3];
rz(-2.5741044) q[3];
sx q[3];
rz(-1.5195888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2813256) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(-2.2567828) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(0.44986808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29174969) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(2.1480609) q[0];
rz(-1.7811071) q[1];
sx q[1];
rz(-0.88881701) q[1];
sx q[1];
rz(-0.55307585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21453126) q[0];
sx q[0];
rz(-1.1916516) q[0];
sx q[0];
rz(1.3247471) q[0];
rz(-pi) q[1];
rz(1.5037914) q[2];
sx q[2];
rz(-2.308146) q[2];
sx q[2];
rz(-0.60947567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83718761) q[1];
sx q[1];
rz(-0.68301979) q[1];
sx q[1];
rz(0.4153312) q[1];
rz(-0.93135466) q[3];
sx q[3];
rz(-1.3512615) q[3];
sx q[3];
rz(-2.1965248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20455827) q[2];
sx q[2];
rz(-0.60595804) q[2];
sx q[2];
rz(2.9045612) q[2];
rz(1.5406746) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(1.1698394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.071526214) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-3.0226829) q[0];
rz(0.43926829) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(3.0781504) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.760323) q[0];
sx q[0];
rz(-2.1043692) q[0];
sx q[0];
rz(2.3593581) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8863438) q[2];
sx q[2];
rz(-2.5016157) q[2];
sx q[2];
rz(-0.6682804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52229228) q[1];
sx q[1];
rz(-1.5792773) q[1];
sx q[1];
rz(1.897476) q[1];
rz(-pi) q[2];
rz(-2.2584553) q[3];
sx q[3];
rz(-2.4534907) q[3];
sx q[3];
rz(0.36206743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0110317) q[2];
sx q[2];
rz(-1.0728027) q[2];
sx q[2];
rz(1.6434742) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.3907923) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(-2.2687268) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(-2.4571669) q[2];
sx q[2];
rz(-2.0906856) q[2];
sx q[2];
rz(-0.69490614) q[2];
rz(1.6428236) q[3];
sx q[3];
rz(-2.147032) q[3];
sx q[3];
rz(1.7346139) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
