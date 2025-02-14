OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.42216766) q[0];
sx q[0];
rz(-2.4165805) q[0];
sx q[0];
rz(2.3056735) q[0];
rz(0.82756502) q[1];
sx q[1];
rz(-2.8711072) q[1];
sx q[1];
rz(-0.40721133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2094858) q[0];
sx q[0];
rz(-1.7833685) q[0];
sx q[0];
rz(1.2680156) q[0];
x q[1];
rz(2.902342) q[2];
sx q[2];
rz(-1.2206589) q[2];
sx q[2];
rz(0.020933271) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8460257) q[1];
sx q[1];
rz(-2.7252619) q[1];
sx q[1];
rz(2.8041072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3916733) q[3];
sx q[3];
rz(-1.1259606) q[3];
sx q[3];
rz(0.47575853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.779458) q[2];
sx q[2];
rz(-2.802765) q[2];
sx q[2];
rz(-0.96620488) q[2];
rz(2.2400014) q[3];
sx q[3];
rz(-0.85351557) q[3];
sx q[3];
rz(1.0454267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7411984) q[0];
sx q[0];
rz(-2.5303831) q[0];
sx q[0];
rz(-3.0389767) q[0];
rz(1.1688894) q[1];
sx q[1];
rz(-2.1714307) q[1];
sx q[1];
rz(0.51365596) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4054905) q[0];
sx q[0];
rz(-1.0879192) q[0];
sx q[0];
rz(3.0391284) q[0];
rz(2.9048237) q[2];
sx q[2];
rz(-0.54950031) q[2];
sx q[2];
rz(2.7236746) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.8674492) q[1];
sx q[1];
rz(-1.871487) q[1];
sx q[1];
rz(-2.1656028) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97764303) q[3];
sx q[3];
rz(-2.6321075) q[3];
sx q[3];
rz(2.2397985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1932842) q[2];
sx q[2];
rz(-0.19364348) q[2];
sx q[2];
rz(2.9330758) q[2];
rz(2.4559313) q[3];
sx q[3];
rz(-2.0833569) q[3];
sx q[3];
rz(-2.9764777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044428069) q[0];
sx q[0];
rz(-1.2760289) q[0];
sx q[0];
rz(-1.2318508) q[0];
rz(-2.7992898) q[1];
sx q[1];
rz(-1.1016223) q[1];
sx q[1];
rz(1.5493578) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3521455) q[0];
sx q[0];
rz(-1.88228) q[0];
sx q[0];
rz(-1.7622663) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9804055) q[2];
sx q[2];
rz(-1.4518132) q[2];
sx q[2];
rz(-2.7143402) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0682837) q[1];
sx q[1];
rz(-1.7080293) q[1];
sx q[1];
rz(-0.84899606) q[1];
rz(-pi) q[2];
rz(0.83264787) q[3];
sx q[3];
rz(-2.4111028) q[3];
sx q[3];
rz(0.60527847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0937097) q[2];
sx q[2];
rz(-1.1992531) q[2];
sx q[2];
rz(-2.1236146) q[2];
rz(1.4114981) q[3];
sx q[3];
rz(-0.74677765) q[3];
sx q[3];
rz(1.5311034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133076) q[0];
sx q[0];
rz(-1.3628553) q[0];
sx q[0];
rz(-3.0174729) q[0];
rz(2.2165551) q[1];
sx q[1];
rz(-1.3003474) q[1];
sx q[1];
rz(0.73840028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8014073) q[0];
sx q[0];
rz(-2.3871195) q[0];
sx q[0];
rz(-1.6615926) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52780789) q[2];
sx q[2];
rz(-2.4619812) q[2];
sx q[2];
rz(-2.2371593) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.55220375) q[1];
sx q[1];
rz(-0.60221803) q[1];
sx q[1];
rz(-0.37005421) q[1];
x q[2];
rz(0.57598007) q[3];
sx q[3];
rz(-1.4864731) q[3];
sx q[3];
rz(2.4654441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.063252123) q[2];
sx q[2];
rz(-1.1880778) q[2];
sx q[2];
rz(-2.038302) q[2];
rz(2.8830146) q[3];
sx q[3];
rz(-1.0753814) q[3];
sx q[3];
rz(-0.28178373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89541188) q[0];
sx q[0];
rz(-2.2639332) q[0];
sx q[0];
rz(1.3377162) q[0];
rz(0.61996639) q[1];
sx q[1];
rz(-1.679531) q[1];
sx q[1];
rz(1.5043129) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1150944) q[0];
sx q[0];
rz(-2.8015602) q[0];
sx q[0];
rz(-0.97397255) q[0];
rz(-pi) q[1];
rz(-0.049311801) q[2];
sx q[2];
rz(-2.4386536) q[2];
sx q[2];
rz(0.21156034) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6474939) q[1];
sx q[1];
rz(-1.955619) q[1];
sx q[1];
rz(-2.5267151) q[1];
rz(-1.4535459) q[3];
sx q[3];
rz(-2.4053889) q[3];
sx q[3];
rz(2.8506926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6164246) q[2];
sx q[2];
rz(-1.1526356) q[2];
sx q[2];
rz(2.0174513) q[2];
rz(0.21640402) q[3];
sx q[3];
rz(-0.83306044) q[3];
sx q[3];
rz(2.0719349) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44806099) q[0];
sx q[0];
rz(-0.8140642) q[0];
sx q[0];
rz(2.6779209) q[0];
rz(3.0060153) q[1];
sx q[1];
rz(-0.99628535) q[1];
sx q[1];
rz(-2.7574976) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6679055) q[0];
sx q[0];
rz(-2.0376922) q[0];
sx q[0];
rz(1.6993194) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21125085) q[2];
sx q[2];
rz(-2.3659035) q[2];
sx q[2];
rz(2.3084909) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64318618) q[1];
sx q[1];
rz(-1.2075042) q[1];
sx q[1];
rz(-0.17885991) q[1];
x q[2];
rz(1.0225419) q[3];
sx q[3];
rz(-1.5060079) q[3];
sx q[3];
rz(2.8509051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7163081) q[2];
sx q[2];
rz(-0.3766489) q[2];
sx q[2];
rz(2.6728805) q[2];
rz(-2.5675755) q[3];
sx q[3];
rz(-1.470397) q[3];
sx q[3];
rz(0.65042692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9531517) q[0];
sx q[0];
rz(-1.3583536) q[0];
sx q[0];
rz(-0.69001946) q[0];
rz(0.72235876) q[1];
sx q[1];
rz(-2.0004309) q[1];
sx q[1];
rz(2.3648327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0714949) q[0];
sx q[0];
rz(-2.8970708) q[0];
sx q[0];
rz(-2.8273409) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58906196) q[2];
sx q[2];
rz(-1.2576332) q[2];
sx q[2];
rz(-2.2472266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1550786) q[1];
sx q[1];
rz(-2.6342794) q[1];
sx q[1];
rz(0.55965565) q[1];
x q[2];
rz(-1.5942105) q[3];
sx q[3];
rz(-0.6802313) q[3];
sx q[3];
rz(-1.7643203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40778273) q[2];
sx q[2];
rz(-1.5668543) q[2];
sx q[2];
rz(-0.053827914) q[2];
rz(2.5194061) q[3];
sx q[3];
rz(-2.6572808) q[3];
sx q[3];
rz(3.0622862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6513885) q[0];
sx q[0];
rz(-1.9909998) q[0];
sx q[0];
rz(1.401061) q[0];
rz(0.98848629) q[1];
sx q[1];
rz(-1.8732312) q[1];
sx q[1];
rz(-2.7338457) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4238472) q[0];
sx q[0];
rz(-1.167594) q[0];
sx q[0];
rz(1.3333596) q[0];
rz(-pi) q[1];
rz(0.8959391) q[2];
sx q[2];
rz(-1.6701588) q[2];
sx q[2];
rz(-2.4473913) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8922006) q[1];
sx q[1];
rz(-2.2731051) q[1];
sx q[1];
rz(2.5125395) q[1];
x q[2];
rz(-2.678773) q[3];
sx q[3];
rz(-0.30202391) q[3];
sx q[3];
rz(-2.8800396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4150348) q[2];
sx q[2];
rz(-1.1316391) q[2];
sx q[2];
rz(-0.69990194) q[2];
rz(1.0837726) q[3];
sx q[3];
rz(-0.86728573) q[3];
sx q[3];
rz(0.23610246) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.938195) q[0];
sx q[0];
rz(-1.6827826) q[0];
sx q[0];
rz(2.7517125) q[0];
rz(-1.1979206) q[1];
sx q[1];
rz(-3.0471264) q[1];
sx q[1];
rz(-0.87497154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0194454) q[0];
sx q[0];
rz(-2.5179389) q[0];
sx q[0];
rz(1.6092997) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88458832) q[2];
sx q[2];
rz(-1.8914127) q[2];
sx q[2];
rz(-0.55701643) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81351133) q[1];
sx q[1];
rz(-1.9167148) q[1];
sx q[1];
rz(-0.4608058) q[1];
rz(2.7194831) q[3];
sx q[3];
rz(-2.096746) q[3];
sx q[3];
rz(1.6933481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8209057) q[2];
sx q[2];
rz(-2.4497538) q[2];
sx q[2];
rz(-2.4917277) q[2];
rz(-2.0248905) q[3];
sx q[3];
rz(-1.2973123) q[3];
sx q[3];
rz(-2.9445649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4138625) q[0];
sx q[0];
rz(-3.08857) q[0];
sx q[0];
rz(-2.9470288) q[0];
rz(2.3692756) q[1];
sx q[1];
rz(-1.1338502) q[1];
sx q[1];
rz(2.9639941) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1188162) q[0];
sx q[0];
rz(-1.5701862) q[0];
sx q[0];
rz(1.5691891) q[0];
rz(2.3760711) q[2];
sx q[2];
rz(-0.81350785) q[2];
sx q[2];
rz(-0.080527079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0655443) q[1];
sx q[1];
rz(-2.0418344) q[1];
sx q[1];
rz(2.440743) q[1];
rz(-pi) q[2];
x q[2];
rz(2.394049) q[3];
sx q[3];
rz(-1.0870544) q[3];
sx q[3];
rz(1.4418937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.78581587) q[2];
sx q[2];
rz(-2.5444701) q[2];
sx q[2];
rz(-0.43006483) q[2];
rz(-2.0307342) q[3];
sx q[3];
rz(-1.4579803) q[3];
sx q[3];
rz(0.41657579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149902) q[0];
sx q[0];
rz(-1.5530598) q[0];
sx q[0];
rz(-0.1834827) q[0];
rz(-1.1732187) q[1];
sx q[1];
rz(-1.219974) q[1];
sx q[1];
rz(-3.0164607) q[1];
rz(-2.8539544) q[2];
sx q[2];
rz(-2.0783744) q[2];
sx q[2];
rz(-1.7219096) q[2];
rz(-0.85632242) q[3];
sx q[3];
rz(-1.9998989) q[3];
sx q[3];
rz(-2.1378703) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
