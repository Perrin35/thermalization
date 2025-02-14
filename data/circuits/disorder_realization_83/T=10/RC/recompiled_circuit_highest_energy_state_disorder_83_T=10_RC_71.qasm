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
rz(-0.42521617) q[0];
sx q[0];
rz(4.4758237) q[0];
sx q[0];
rz(9.0444179) q[0];
rz(1.7827787) q[1];
sx q[1];
rz(2.9549197) q[1];
sx q[1];
rz(8.5001707) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3596852) q[0];
sx q[0];
rz(-1.230002) q[0];
sx q[0];
rz(-0.73990344) q[0];
rz(-pi) q[1];
rz(1.0546646) q[2];
sx q[2];
rz(-1.4251815) q[2];
sx q[2];
rz(0.21805412) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2362931) q[1];
sx q[1];
rz(-2.2989706) q[1];
sx q[1];
rz(-3.0889126) q[1];
rz(-pi) q[2];
rz(-3.0713586) q[3];
sx q[3];
rz(-3.0133504) q[3];
sx q[3];
rz(-1.8907036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1251462) q[2];
sx q[2];
rz(-0.95637286) q[2];
sx q[2];
rz(-2.1535786) q[2];
rz(-2.9347349) q[3];
sx q[3];
rz(-1.4068973) q[3];
sx q[3];
rz(0.44001165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4619231) q[0];
sx q[0];
rz(-1.4689057) q[0];
sx q[0];
rz(-0.2521421) q[0];
rz(-3.120046) q[1];
sx q[1];
rz(-0.69308678) q[1];
sx q[1];
rz(-0.85539877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6538174) q[0];
sx q[0];
rz(-0.027055351) q[0];
sx q[0];
rz(-1.5144801) q[0];
rz(-0.62521817) q[2];
sx q[2];
rz(-2.4897235) q[2];
sx q[2];
rz(2.9232581) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.80551863) q[1];
sx q[1];
rz(-1.0401498) q[1];
sx q[1];
rz(-2.8855399) q[1];
rz(-0.40327252) q[3];
sx q[3];
rz(-0.72411116) q[3];
sx q[3];
rz(-1.6722752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97027332) q[2];
sx q[2];
rz(-3.1372941) q[2];
sx q[2];
rz(0.26101905) q[2];
rz(3.1018992) q[3];
sx q[3];
rz(-1.3986162) q[3];
sx q[3];
rz(0.54350129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71565851) q[0];
sx q[0];
rz(-1.4733227) q[0];
sx q[0];
rz(0.69450992) q[0];
rz(2.014324) q[1];
sx q[1];
rz(-0.85113168) q[1];
sx q[1];
rz(2.1515501) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3304402) q[0];
sx q[0];
rz(-1.2227204) q[0];
sx q[0];
rz(1.6578107) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98857356) q[2];
sx q[2];
rz(-1.7372088) q[2];
sx q[2];
rz(2.5022282) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5591919) q[1];
sx q[1];
rz(-1.9291996) q[1];
sx q[1];
rz(2.6646975) q[1];
rz(-pi) q[2];
rz(2.2740836) q[3];
sx q[3];
rz(-1.1425619) q[3];
sx q[3];
rz(3.0813062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29828829) q[2];
sx q[2];
rz(-0.97323209) q[2];
sx q[2];
rz(0.49015552) q[2];
rz(2.7847024) q[3];
sx q[3];
rz(-0.39339742) q[3];
sx q[3];
rz(1.2662158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0855584) q[0];
sx q[0];
rz(-1.9558676) q[0];
sx q[0];
rz(-2.2991142) q[0];
rz(0.8017686) q[1];
sx q[1];
rz(-2.8623878) q[1];
sx q[1];
rz(0.78559771) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426823) q[0];
sx q[0];
rz(-1.8225095) q[0];
sx q[0];
rz(2.4221735) q[0];
x q[1];
rz(2.9688492) q[2];
sx q[2];
rz(-0.72451353) q[2];
sx q[2];
rz(-2.8959664) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3093058) q[1];
sx q[1];
rz(-0.63696276) q[1];
sx q[1];
rz(0.5354521) q[1];
x q[2];
rz(0.81964747) q[3];
sx q[3];
rz(-1.7213806) q[3];
sx q[3];
rz(2.8898847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0164612) q[2];
sx q[2];
rz(-2.3931914) q[2];
sx q[2];
rz(-2.1330736) q[2];
rz(2.290001) q[3];
sx q[3];
rz(-2.3840756) q[3];
sx q[3];
rz(-1.0803224) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1763024) q[0];
sx q[0];
rz(-0.98133123) q[0];
sx q[0];
rz(-1.0705795) q[0];
rz(2.3640682) q[1];
sx q[1];
rz(-2.2309525) q[1];
sx q[1];
rz(-1.0106962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010653822) q[0];
sx q[0];
rz(-1.9294327) q[0];
sx q[0];
rz(0.50827311) q[0];
x q[1];
rz(-1.9254743) q[2];
sx q[2];
rz(-0.18624072) q[2];
sx q[2];
rz(2.8493488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0194369) q[1];
sx q[1];
rz(-2.2194462) q[1];
sx q[1];
rz(-2.1557356) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6701974) q[3];
sx q[3];
rz(-1.0109006) q[3];
sx q[3];
rz(-0.82370629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8396478) q[2];
sx q[2];
rz(-3.0366615) q[2];
sx q[2];
rz(-1.5860175) q[2];
rz(-1.0157061) q[3];
sx q[3];
rz(-1.1414707) q[3];
sx q[3];
rz(1.7293845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79065901) q[0];
sx q[0];
rz(-2.9582294) q[0];
sx q[0];
rz(-2.3405128) q[0];
rz(-3.0084897) q[1];
sx q[1];
rz(-1.3214654) q[1];
sx q[1];
rz(0.4020234) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19045705) q[0];
sx q[0];
rz(-1.2371089) q[0];
sx q[0];
rz(-2.3526089) q[0];
x q[1];
rz(1.5779675) q[2];
sx q[2];
rz(-2.8507887) q[2];
sx q[2];
rz(-0.70722843) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0084486246) q[1];
sx q[1];
rz(-2.5767972) q[1];
sx q[1];
rz(2.1022335) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16120637) q[3];
sx q[3];
rz(-0.66449249) q[3];
sx q[3];
rz(-2.4624069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5037527) q[2];
sx q[2];
rz(-1.736234) q[2];
sx q[2];
rz(2.3657738) q[2];
rz(-1.4555629) q[3];
sx q[3];
rz(-1.9676696) q[3];
sx q[3];
rz(-2.1337401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044947226) q[0];
sx q[0];
rz(-0.89770397) q[0];
sx q[0];
rz(0.67895472) q[0];
rz(-1.8567765) q[1];
sx q[1];
rz(-0.7130475) q[1];
sx q[1];
rz(1.3791893) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5011936) q[0];
sx q[0];
rz(-1.3974579) q[0];
sx q[0];
rz(0.79357432) q[0];
x q[1];
rz(-1.8544312) q[2];
sx q[2];
rz(-2.8511308) q[2];
sx q[2];
rz(-2.1112372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0022565) q[1];
sx q[1];
rz(-2.4508307) q[1];
sx q[1];
rz(0.73378566) q[1];
x q[2];
rz(1.3706743) q[3];
sx q[3];
rz(-1.6255696) q[3];
sx q[3];
rz(1.2204264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8046367) q[2];
sx q[2];
rz(-2.3865484) q[2];
sx q[2];
rz(1.9026559) q[2];
rz(-1.8094481) q[3];
sx q[3];
rz(-2.381031) q[3];
sx q[3];
rz(-0.24470394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.4181344) q[0];
sx q[0];
rz(-1.9754388) q[0];
sx q[0];
rz(0.13846692) q[0];
rz(0.18374099) q[1];
sx q[1];
rz(-0.67444363) q[1];
sx q[1];
rz(-0.26652452) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12611575) q[0];
sx q[0];
rz(-1.4086282) q[0];
sx q[0];
rz(1.9006985) q[0];
rz(2.8710352) q[2];
sx q[2];
rz(-1.0204401) q[2];
sx q[2];
rz(1.5824883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5268685) q[1];
sx q[1];
rz(-2.7182332) q[1];
sx q[1];
rz(1.7080659) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4146904) q[3];
sx q[3];
rz(-1.57734) q[3];
sx q[3];
rz(2.9605799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0900241) q[2];
sx q[2];
rz(-1.2568306) q[2];
sx q[2];
rz(2.9070692) q[2];
rz(-1.6656434) q[3];
sx q[3];
rz(-1.6022976) q[3];
sx q[3];
rz(0.75638151) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.804857) q[0];
sx q[0];
rz(-2.2754301) q[0];
sx q[0];
rz(2.3418703) q[0];
rz(2.5526478) q[1];
sx q[1];
rz(-2.2942693) q[1];
sx q[1];
rz(2.934093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.554858) q[0];
sx q[0];
rz(-1.7014781) q[0];
sx q[0];
rz(-2.7452637) q[0];
x q[1];
rz(2.2830711) q[2];
sx q[2];
rz(-1.8078765) q[2];
sx q[2];
rz(1.750994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3424399) q[1];
sx q[1];
rz(-0.9121597) q[1];
sx q[1];
rz(0.51167804) q[1];
x q[2];
rz(2.0405681) q[3];
sx q[3];
rz(-1.5964526) q[3];
sx q[3];
rz(-2.2750843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3346682) q[2];
sx q[2];
rz(-0.88403264) q[2];
sx q[2];
rz(-0.36726382) q[2];
rz(1.4372829) q[3];
sx q[3];
rz(-1.9673012) q[3];
sx q[3];
rz(1.1737163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2631898) q[0];
sx q[0];
rz(-2.1387565) q[0];
sx q[0];
rz(2.3314085) q[0];
rz(1.1148249) q[1];
sx q[1];
rz(-2.7572542) q[1];
sx q[1];
rz(-0.55327639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0598349) q[0];
sx q[0];
rz(-2.5011269) q[0];
sx q[0];
rz(1.6637131) q[0];
rz(-pi) q[1];
rz(0.79375879) q[2];
sx q[2];
rz(-1.1053893) q[2];
sx q[2];
rz(-2.2528354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9844891) q[1];
sx q[1];
rz(-1.2890639) q[1];
sx q[1];
rz(0.19428044) q[1];
x q[2];
rz(1.2180042) q[3];
sx q[3];
rz(-0.81656352) q[3];
sx q[3];
rz(1.0366131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1222003) q[2];
sx q[2];
rz(-2.4330752) q[2];
sx q[2];
rz(0.23966399) q[2];
rz(-1.8137118) q[3];
sx q[3];
rz(-2.0769104) q[3];
sx q[3];
rz(-2.6740668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620621) q[0];
sx q[0];
rz(-1.7876328) q[0];
sx q[0];
rz(-2.7476516) q[0];
rz(-2.486034) q[1];
sx q[1];
rz(-1.1759023) q[1];
sx q[1];
rz(-2.1942153) q[1];
rz(1.9026359) q[2];
sx q[2];
rz(-1.9658026) q[2];
sx q[2];
rz(2.6323242) q[2];
rz(1.3169133) q[3];
sx q[3];
rz(-1.0051654) q[3];
sx q[3];
rz(-1.2011423) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
