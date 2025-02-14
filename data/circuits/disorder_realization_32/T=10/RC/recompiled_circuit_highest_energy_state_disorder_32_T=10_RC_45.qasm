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
rz(1.2122756) q[0];
sx q[0];
rz(4.2733856) q[0];
sx q[0];
rz(10.797664) q[0];
rz(2.0139439) q[1];
sx q[1];
rz(-1.7665266) q[1];
sx q[1];
rz(2.3553203) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92106956) q[0];
sx q[0];
rz(-1.2036909) q[0];
sx q[0];
rz(-2.9216032) q[0];
x q[1];
rz(0.55297466) q[2];
sx q[2];
rz(-0.99054256) q[2];
sx q[2];
rz(2.0748823) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9278501) q[1];
sx q[1];
rz(-1.8600271) q[1];
sx q[1];
rz(-0.031164483) q[1];
rz(-2.1350098) q[3];
sx q[3];
rz(-0.67970961) q[3];
sx q[3];
rz(-2.2077479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11975153) q[2];
sx q[2];
rz(-1.4234875) q[2];
sx q[2];
rz(-1.8915668) q[2];
rz(2.405808) q[3];
sx q[3];
rz(-0.17446987) q[3];
sx q[3];
rz(-1.638394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4982872) q[0];
sx q[0];
rz(-1.1730288) q[0];
sx q[0];
rz(-0.071320891) q[0];
rz(1.1530863) q[1];
sx q[1];
rz(-0.76140296) q[1];
sx q[1];
rz(-0.52871314) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3246177) q[0];
sx q[0];
rz(-1.2551089) q[0];
sx q[0];
rz(2.0633008) q[0];
x q[1];
rz(0.18888338) q[2];
sx q[2];
rz(-1.3391277) q[2];
sx q[2];
rz(-0.29555368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9829927) q[1];
sx q[1];
rz(-1.0550749) q[1];
sx q[1];
rz(0.92293585) q[1];
rz(-pi) q[2];
rz(1.0924322) q[3];
sx q[3];
rz(-0.4606263) q[3];
sx q[3];
rz(-0.13308976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66136393) q[2];
sx q[2];
rz(-2.7429548) q[2];
sx q[2];
rz(-0.028623494) q[2];
rz(0.35401595) q[3];
sx q[3];
rz(-0.75494868) q[3];
sx q[3];
rz(2.5825175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70327586) q[0];
sx q[0];
rz(-1.0030614) q[0];
sx q[0];
rz(-2.2502374) q[0];
rz(-0.50093961) q[1];
sx q[1];
rz(-1.5762065) q[1];
sx q[1];
rz(-0.058813728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0552154) q[0];
sx q[0];
rz(-0.098345938) q[0];
sx q[0];
rz(2.6340061) q[0];
rz(-pi) q[1];
rz(-1.0436613) q[2];
sx q[2];
rz(-1.322515) q[2];
sx q[2];
rz(-0.52094007) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0173024) q[1];
sx q[1];
rz(-0.92355403) q[1];
sx q[1];
rz(-0.94396293) q[1];
rz(-pi) q[2];
rz(-2.9611582) q[3];
sx q[3];
rz(-1.5431607) q[3];
sx q[3];
rz(0.35865006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7848876) q[2];
sx q[2];
rz(-2.5829743) q[2];
sx q[2];
rz(-0.2365665) q[2];
rz(2.2208354) q[3];
sx q[3];
rz(-2.0906788) q[3];
sx q[3];
rz(-1.4111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.22223602) q[0];
sx q[0];
rz(-1.284143) q[0];
sx q[0];
rz(-0.87523571) q[0];
rz(1.596176) q[1];
sx q[1];
rz(-0.86219209) q[1];
sx q[1];
rz(-0.045305591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047642853) q[0];
sx q[0];
rz(-1.4685306) q[0];
sx q[0];
rz(-0.065667466) q[0];
rz(-2.6338121) q[2];
sx q[2];
rz(-0.96755257) q[2];
sx q[2];
rz(2.809462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.873535) q[1];
sx q[1];
rz(-0.90374871) q[1];
sx q[1];
rz(-0.11063834) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1940345) q[3];
sx q[3];
rz(-2.1450419) q[3];
sx q[3];
rz(-1.6358835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.28216013) q[2];
sx q[2];
rz(-2.1789447) q[2];
sx q[2];
rz(-1.1939987) q[2];
rz(-0.89030877) q[3];
sx q[3];
rz(-1.2229536) q[3];
sx q[3];
rz(-0.94849006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0328261) q[0];
sx q[0];
rz(-2.325401) q[0];
sx q[0];
rz(-0.68242514) q[0];
rz(1.2884864) q[1];
sx q[1];
rz(-1.5792184) q[1];
sx q[1];
rz(1.3312181) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2829943) q[0];
sx q[0];
rz(-2.4147075) q[0];
sx q[0];
rz(1.5461393) q[0];
rz(-pi) q[1];
rz(-1.17679) q[2];
sx q[2];
rz(-1.3298243) q[2];
sx q[2];
rz(-0.35277982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8537123) q[1];
sx q[1];
rz(-1.9287709) q[1];
sx q[1];
rz(1.720721) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87792895) q[3];
sx q[3];
rz(-2.4932541) q[3];
sx q[3];
rz(-0.058365783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1149301) q[2];
sx q[2];
rz(-2.5422091) q[2];
sx q[2];
rz(0.66693532) q[2];
rz(-0.55495787) q[3];
sx q[3];
rz(-2.4542377) q[3];
sx q[3];
rz(-2.8426898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32894593) q[0];
sx q[0];
rz(-1.2586559) q[0];
sx q[0];
rz(2.4378648) q[0];
rz(0.16920371) q[1];
sx q[1];
rz(-1.6177982) q[1];
sx q[1];
rz(0.15883787) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7089091) q[0];
sx q[0];
rz(-1.7129717) q[0];
sx q[0];
rz(0.12169038) q[0];
x q[1];
rz(-1.6936982) q[2];
sx q[2];
rz(-1.7235061) q[2];
sx q[2];
rz(2.0982775) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.006268) q[1];
sx q[1];
rz(-0.49990955) q[1];
sx q[1];
rz(-0.27854021) q[1];
rz(-pi) q[2];
rz(-2.5037836) q[3];
sx q[3];
rz(-1.7553698) q[3];
sx q[3];
rz(-1.0705494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77330294) q[2];
sx q[2];
rz(-0.94179073) q[2];
sx q[2];
rz(-0.85757315) q[2];
rz(-1.1693906) q[3];
sx q[3];
rz(-1.8604859) q[3];
sx q[3];
rz(-2.9331971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9661949) q[0];
sx q[0];
rz(-0.76681391) q[0];
sx q[0];
rz(0.71863693) q[0];
rz(-0.28193685) q[1];
sx q[1];
rz(-1.7666631) q[1];
sx q[1];
rz(2.4729572) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40067264) q[0];
sx q[0];
rz(-2.1869279) q[0];
sx q[0];
rz(0.4192062) q[0];
x q[1];
rz(-0.35923194) q[2];
sx q[2];
rz(-0.97373913) q[2];
sx q[2];
rz(3.0199241) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9033176) q[1];
sx q[1];
rz(-1.2438477) q[1];
sx q[1];
rz(-1.9592778) q[1];
rz(1.5626146) q[3];
sx q[3];
rz(-0.93957179) q[3];
sx q[3];
rz(-0.62731987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.64688993) q[2];
sx q[2];
rz(-1.9560445) q[2];
sx q[2];
rz(1.0364214) q[2];
rz(-0.63452619) q[3];
sx q[3];
rz(-2.1536638) q[3];
sx q[3];
rz(-0.79276597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.46102872) q[0];
sx q[0];
rz(-0.11756086) q[0];
sx q[0];
rz(-0.72599894) q[0];
rz(1.9793824) q[1];
sx q[1];
rz(-1.9644968) q[1];
sx q[1];
rz(-1.1964218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40563075) q[0];
sx q[0];
rz(-2.9655632) q[0];
sx q[0];
rz(1.733863) q[0];
rz(-pi) q[1];
rz(3.093946) q[2];
sx q[2];
rz(-1.8777913) q[2];
sx q[2];
rz(0.36587151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4359522) q[1];
sx q[1];
rz(-2.0355823) q[1];
sx q[1];
rz(2.57354) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2231908) q[3];
sx q[3];
rz(-2.2112527) q[3];
sx q[3];
rz(0.15061041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2339345) q[2];
sx q[2];
rz(-1.5936759) q[2];
sx q[2];
rz(2.2080803) q[2];
rz(-2.524611) q[3];
sx q[3];
rz(-2.139293) q[3];
sx q[3];
rz(-2.2945819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9762978) q[0];
sx q[0];
rz(-3.0369861) q[0];
sx q[0];
rz(-0.053255178) q[0];
rz(2.8540197) q[1];
sx q[1];
rz(-0.92929274) q[1];
sx q[1];
rz(2.8823749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1142562) q[0];
sx q[0];
rz(-1.5948609) q[0];
sx q[0];
rz(1.5177478) q[0];
rz(-1.0561814) q[2];
sx q[2];
rz(-0.26517235) q[2];
sx q[2];
rz(2.3725703) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1500435) q[1];
sx q[1];
rz(-0.75356149) q[1];
sx q[1];
rz(-1.4808473) q[1];
rz(0.15786981) q[3];
sx q[3];
rz(-2.3693858) q[3];
sx q[3];
rz(-1.7534353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58664924) q[2];
sx q[2];
rz(-1.2536851) q[2];
sx q[2];
rz(-2.9602236) q[2];
rz(0.3950611) q[3];
sx q[3];
rz(-2.919988) q[3];
sx q[3];
rz(2.1612397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3996537) q[0];
sx q[0];
rz(-1.548865) q[0];
sx q[0];
rz(0.3821061) q[0];
rz(2.8451879) q[1];
sx q[1];
rz(-2.1370685) q[1];
sx q[1];
rz(-1.2688676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4110376) q[0];
sx q[0];
rz(-2.0417622) q[0];
sx q[0];
rz(2.6593326) q[0];
rz(-0.89752723) q[2];
sx q[2];
rz(-1.6746876) q[2];
sx q[2];
rz(1.8477224) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3281607) q[1];
sx q[1];
rz(-0.34036553) q[1];
sx q[1];
rz(-1.634738) q[1];
x q[2];
rz(-2.6878854) q[3];
sx q[3];
rz(-1.6804916) q[3];
sx q[3];
rz(0.060893313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2962013) q[2];
sx q[2];
rz(-0.56075823) q[2];
sx q[2];
rz(2.7517547) q[2];
rz(1.2975533) q[3];
sx q[3];
rz(-1.1417979) q[3];
sx q[3];
rz(-2.2977184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525477) q[0];
sx q[0];
rz(-1.4023517) q[0];
sx q[0];
rz(-1.5040816) q[0];
rz(1.3564431) q[1];
sx q[1];
rz(-2.103613) q[1];
sx q[1];
rz(1.6115859) q[1];
rz(0.46089725) q[2];
sx q[2];
rz(-2.1620038) q[2];
sx q[2];
rz(2.9481859) q[2];
rz(2.8036694) q[3];
sx q[3];
rz(-0.81339627) q[3];
sx q[3];
rz(0.42311121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
