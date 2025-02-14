OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.77914971) q[0];
sx q[0];
rz(2.7326512) q[0];
sx q[0];
rz(10.41806) q[0];
rz(-0.14708695) q[1];
sx q[1];
rz(-0.99647254) q[1];
sx q[1];
rz(-1.7239404) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0298457) q[0];
sx q[0];
rz(-0.8518712) q[0];
sx q[0];
rz(-1.0875888) q[0];
x q[1];
rz(2.4383847) q[2];
sx q[2];
rz(-2.0173912) q[2];
sx q[2];
rz(-2.869702) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0252391) q[1];
sx q[1];
rz(-2.1903746) q[1];
sx q[1];
rz(-1.2337633) q[1];
rz(-0.032790498) q[3];
sx q[3];
rz(-1.7049978) q[3];
sx q[3];
rz(-0.57574948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.950497) q[2];
sx q[2];
rz(-1.3206626) q[2];
sx q[2];
rz(-2.0868059) q[2];
rz(-0.47109207) q[3];
sx q[3];
rz(-1.6867009) q[3];
sx q[3];
rz(-2.2733222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3926113) q[0];
sx q[0];
rz(-1.457021) q[0];
sx q[0];
rz(0.23496041) q[0];
rz(1.2745534) q[1];
sx q[1];
rz(-1.3095368) q[1];
sx q[1];
rz(2.4539006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10768453) q[0];
sx q[0];
rz(-1.8982197) q[0];
sx q[0];
rz(-1.4998687) q[0];
x q[1];
rz(2.657452) q[2];
sx q[2];
rz(-2.2233367) q[2];
sx q[2];
rz(0.86838858) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.60581698) q[1];
sx q[1];
rz(-1.5652085) q[1];
sx q[1];
rz(1.5582915) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95716547) q[3];
sx q[3];
rz(-1.6464697) q[3];
sx q[3];
rz(-1.5694973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4253652) q[2];
sx q[2];
rz(-1.8391515) q[2];
sx q[2];
rz(2.9322374) q[2];
rz(-2.1685205) q[3];
sx q[3];
rz(-0.89997411) q[3];
sx q[3];
rz(-0.59276855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8415602) q[0];
sx q[0];
rz(-2.7356) q[0];
sx q[0];
rz(-5/(11*pi)) q[0];
rz(1.9035829) q[1];
sx q[1];
rz(-2.369945) q[1];
sx q[1];
rz(-3.0498116) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9973213) q[0];
sx q[0];
rz(-1.034212) q[0];
sx q[0];
rz(-2.7377605) q[0];
rz(-pi) q[1];
rz(-2.1917159) q[2];
sx q[2];
rz(-0.11813049) q[2];
sx q[2];
rz(1.0850414) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5621693) q[1];
sx q[1];
rz(-1.3797288) q[1];
sx q[1];
rz(-1.3538989) q[1];
rz(-pi) q[2];
rz(-1.1385659) q[3];
sx q[3];
rz(-1.5097268) q[3];
sx q[3];
rz(1.40738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7564275) q[2];
sx q[2];
rz(-1.6331208) q[2];
sx q[2];
rz(2.7623994) q[2];
rz(-2.3181629) q[3];
sx q[3];
rz(-0.40612602) q[3];
sx q[3];
rz(-2.8520975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8150811) q[0];
sx q[0];
rz(-0.87711763) q[0];
sx q[0];
rz(2.1429578) q[0];
rz(-0.48209349) q[1];
sx q[1];
rz(-0.95078743) q[1];
sx q[1];
rz(-1.8236209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4559193) q[0];
sx q[0];
rz(-0.64252526) q[0];
sx q[0];
rz(2.5701017) q[0];
rz(-pi) q[1];
rz(1.3236207) q[2];
sx q[2];
rz(-0.99946076) q[2];
sx q[2];
rz(0.94117576) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.4596583) q[1];
sx q[1];
rz(-1.0409969) q[1];
sx q[1];
rz(1.3270686) q[1];
rz(-1.0467806) q[3];
sx q[3];
rz(-2.2923427) q[3];
sx q[3];
rz(1.0224354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2202997) q[2];
sx q[2];
rz(-2.2396294) q[2];
sx q[2];
rz(2.656142) q[2];
rz(-2.7576533) q[3];
sx q[3];
rz(-1.5555614) q[3];
sx q[3];
rz(2.3575822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081472814) q[0];
sx q[0];
rz(-0.10426846) q[0];
sx q[0];
rz(0.01165788) q[0];
rz(-0.13721379) q[1];
sx q[1];
rz(-1.5533181) q[1];
sx q[1];
rz(1.3020017) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0903895) q[0];
sx q[0];
rz(-0.21349354) q[0];
sx q[0];
rz(2.6043686) q[0];
rz(0.88413357) q[2];
sx q[2];
rz(-2.4120286) q[2];
sx q[2];
rz(-1.8232249) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.72921371) q[1];
sx q[1];
rz(-1.679031) q[1];
sx q[1];
rz(-1.257105) q[1];
x q[2];
rz(-1.2106895) q[3];
sx q[3];
rz(-1.2146287) q[3];
sx q[3];
rz(2.3115013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7234708) q[2];
sx q[2];
rz(-1.301441) q[2];
sx q[2];
rz(2.8483025) q[2];
rz(-0.063118525) q[3];
sx q[3];
rz(-1.9189574) q[3];
sx q[3];
rz(-2.2466834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1191331) q[0];
sx q[0];
rz(-0.28894153) q[0];
sx q[0];
rz(3.1030848) q[0];
rz(-1.0844082) q[1];
sx q[1];
rz(-2.7190828) q[1];
sx q[1];
rz(-0.96673059) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9597067) q[0];
sx q[0];
rz(-0.27508914) q[0];
sx q[0];
rz(2.3692714) q[0];
x q[1];
rz(-3.128563) q[2];
sx q[2];
rz(-1.3037062) q[2];
sx q[2];
rz(2.992127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0930007) q[1];
sx q[1];
rz(-1.1735667) q[1];
sx q[1];
rz(0.36919682) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8916675) q[3];
sx q[3];
rz(-2.2062613) q[3];
sx q[3];
rz(0.15033406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66095573) q[2];
sx q[2];
rz(-0.2350685) q[2];
sx q[2];
rz(2.1601775) q[2];
rz(1.6437982) q[3];
sx q[3];
rz(-1.8794329) q[3];
sx q[3];
rz(1.5463382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3129231) q[0];
sx q[0];
rz(-0.68083119) q[0];
sx q[0];
rz(2.8269826) q[0];
rz(0.18868748) q[1];
sx q[1];
rz(-0.43470946) q[1];
sx q[1];
rz(-2.6729118) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9148736) q[0];
sx q[0];
rz(-1.813199) q[0];
sx q[0];
rz(-1.7049432) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1397083) q[2];
sx q[2];
rz(-1.4030289) q[2];
sx q[2];
rz(1.1275856) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0951765) q[1];
sx q[1];
rz(-0.79907387) q[1];
sx q[1];
rz(-2.0002736) q[1];
x q[2];
rz(-0.46814274) q[3];
sx q[3];
rz(-1.8624412) q[3];
sx q[3];
rz(0.95671747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9357052) q[2];
sx q[2];
rz(-2.2998655) q[2];
sx q[2];
rz(0.97071281) q[2];
rz(-2.6054221) q[3];
sx q[3];
rz(-1.9325117) q[3];
sx q[3];
rz(-2.4255588) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69990528) q[0];
sx q[0];
rz(-0.55460414) q[0];
sx q[0];
rz(-1.6957138) q[0];
rz(1.0384167) q[1];
sx q[1];
rz(-2.0380135) q[1];
sx q[1];
rz(-0.081092484) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86085194) q[0];
sx q[0];
rz(-0.74374712) q[0];
sx q[0];
rz(2.70983) q[0];
rz(-pi) q[1];
rz(2.7620188) q[2];
sx q[2];
rz(-2.4502769) q[2];
sx q[2];
rz(2.8065681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3861919) q[1];
sx q[1];
rz(-0.31495783) q[1];
sx q[1];
rz(1.1909199) q[1];
rz(-pi) q[2];
rz(-2.0399953) q[3];
sx q[3];
rz(-1.0476255) q[3];
sx q[3];
rz(-1.2468456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9424092) q[2];
sx q[2];
rz(-0.77745357) q[2];
sx q[2];
rz(-0.24277631) q[2];
rz(-1.5593922) q[3];
sx q[3];
rz(-0.42583164) q[3];
sx q[3];
rz(-1.8886214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.2389857) q[0];
sx q[0];
rz(-1.0317529) q[0];
sx q[0];
rz(1.8865939) q[0];
rz(-1.7651419) q[1];
sx q[1];
rz(-2.1170728) q[1];
sx q[1];
rz(2.4741516) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28813002) q[0];
sx q[0];
rz(-1.9024182) q[0];
sx q[0];
rz(1.0378855) q[0];
rz(-pi) q[1];
rz(-1.5282643) q[2];
sx q[2];
rz(-2.0669524) q[2];
sx q[2];
rz(-0.68788487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63782843) q[1];
sx q[1];
rz(-2.1385647) q[1];
sx q[1];
rz(0.34620398) q[1];
x q[2];
rz(-1.8967751) q[3];
sx q[3];
rz(-1.2080844) q[3];
sx q[3];
rz(2.6397702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61975512) q[2];
sx q[2];
rz(-2.4027368) q[2];
sx q[2];
rz(-2.6832306) q[2];
rz(1.5357664) q[3];
sx q[3];
rz(-2.063844) q[3];
sx q[3];
rz(2.2569236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.2430724) q[0];
sx q[0];
rz(-0.5721108) q[0];
sx q[0];
rz(-2.7761053) q[0];
rz(-0.99994031) q[1];
sx q[1];
rz(-1.702407) q[1];
sx q[1];
rz(1.684729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92577584) q[0];
sx q[0];
rz(-2.121637) q[0];
sx q[0];
rz(-0.23683817) q[0];
x q[1];
rz(-1.959895) q[2];
sx q[2];
rz(-0.28841296) q[2];
sx q[2];
rz(2.4450977) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7813606) q[1];
sx q[1];
rz(-1.7443027) q[1];
sx q[1];
rz(0.75299112) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0356009) q[3];
sx q[3];
rz(-0.92886954) q[3];
sx q[3];
rz(-1.2479051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0008056) q[2];
sx q[2];
rz(-1.8203338) q[2];
sx q[2];
rz(2.136039) q[2];
rz(2.4908861) q[3];
sx q[3];
rz(-1.4395827) q[3];
sx q[3];
rz(-0.85056359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.0477796) q[0];
sx q[0];
rz(-2.35738) q[0];
sx q[0];
rz(2.1398075) q[0];
rz(1.4108989) q[1];
sx q[1];
rz(-1.4374562) q[1];
sx q[1];
rz(1.1387574) q[1];
rz(0.11713709) q[2];
sx q[2];
rz(-1.8979372) q[2];
sx q[2];
rz(-2.7330782) q[2];
rz(0.10436124) q[3];
sx q[3];
rz(-1.6634533) q[3];
sx q[3];
rz(-2.8692393) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
