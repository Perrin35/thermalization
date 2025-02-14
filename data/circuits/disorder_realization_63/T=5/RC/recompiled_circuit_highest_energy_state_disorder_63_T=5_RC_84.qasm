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
rz(-1.0492078) q[0];
sx q[0];
rz(-2.3823491) q[0];
sx q[0];
rz(0.37641755) q[0];
rz(-1.5837826) q[1];
sx q[1];
rz(-2.0720785) q[1];
sx q[1];
rz(2.1853316) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96891785) q[0];
sx q[0];
rz(-2.4233074) q[0];
sx q[0];
rz(2.0613614) q[0];
rz(-pi) q[1];
rz(2.7189288) q[2];
sx q[2];
rz(-1.3939121) q[2];
sx q[2];
rz(-1.6852578) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2884644) q[1];
sx q[1];
rz(-0.83544105) q[1];
sx q[1];
rz(2.9771027) q[1];
x q[2];
rz(-0.27122916) q[3];
sx q[3];
rz(-1.800897) q[3];
sx q[3];
rz(-0.35121894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3236397) q[2];
sx q[2];
rz(-0.91922593) q[2];
sx q[2];
rz(-1.9383355) q[2];
rz(2.1008927) q[3];
sx q[3];
rz(-2.2668656) q[3];
sx q[3];
rz(-3.021595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1131209) q[0];
sx q[0];
rz(-0.57336837) q[0];
sx q[0];
rz(-1.1125125) q[0];
rz(1.6587229) q[1];
sx q[1];
rz(-0.590938) q[1];
sx q[1];
rz(-2.9842751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7027037) q[0];
sx q[0];
rz(-1.7458785) q[0];
sx q[0];
rz(-0.69467993) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4837166) q[2];
sx q[2];
rz(-1.9425937) q[2];
sx q[2];
rz(1.4390989) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8437905) q[1];
sx q[1];
rz(-2.4239967) q[1];
sx q[1];
rz(-1.9466496) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6460765) q[3];
sx q[3];
rz(-1.7930248) q[3];
sx q[3];
rz(-0.27929515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4979672) q[2];
sx q[2];
rz(-0.47933856) q[2];
sx q[2];
rz(-0.41110006) q[2];
rz(2.5655365) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(-2.1435553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75329798) q[0];
sx q[0];
rz(-1.0030712) q[0];
sx q[0];
rz(0.74111795) q[0];
rz(-1.6916212) q[1];
sx q[1];
rz(-1.8284109) q[1];
sx q[1];
rz(-0.57201874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30878371) q[0];
sx q[0];
rz(-2.7833412) q[0];
sx q[0];
rz(0.97106309) q[0];
rz(-1.3640312) q[2];
sx q[2];
rz(-1.8305147) q[2];
sx q[2];
rz(2.5842359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0894172) q[1];
sx q[1];
rz(-1.1686109) q[1];
sx q[1];
rz(-2.0341464) q[1];
rz(-pi) q[2];
rz(-2.9111262) q[3];
sx q[3];
rz(-1.8982045) q[3];
sx q[3];
rz(2.5203506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.584562) q[2];
sx q[2];
rz(-1.9494373) q[2];
sx q[2];
rz(-2.3438047) q[2];
rz(-1.3095464) q[3];
sx q[3];
rz(-1.8833501) q[3];
sx q[3];
rz(1.4057188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.13084594) q[0];
sx q[0];
rz(-2.2444785) q[0];
sx q[0];
rz(-2.77453) q[0];
rz(-1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(-0.22148111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85773112) q[0];
sx q[0];
rz(-0.22131187) q[0];
sx q[0];
rz(2.8279634) q[0];
rz(-pi) q[1];
rz(-0.367737) q[2];
sx q[2];
rz(-1.6984816) q[2];
sx q[2];
rz(2.5177296) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74879157) q[1];
sx q[1];
rz(-1.9828771) q[1];
sx q[1];
rz(-1.2296299) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52677299) q[3];
sx q[3];
rz(-1.793141) q[3];
sx q[3];
rz(0.46565817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14747846) q[2];
sx q[2];
rz(-1.9693815) q[2];
sx q[2];
rz(-1.3383024) q[2];
rz(-0.5021247) q[3];
sx q[3];
rz(-3.0410671) q[3];
sx q[3];
rz(-2.8978735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62622708) q[0];
sx q[0];
rz(-2.5318662) q[0];
sx q[0];
rz(1.6060265) q[0];
rz(-3.0961127) q[1];
sx q[1];
rz(-1.3810424) q[1];
sx q[1];
rz(-2.6754726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.202318) q[0];
sx q[0];
rz(-0.28616787) q[0];
sx q[0];
rz(-2.0757458) q[0];
rz(-pi) q[1];
rz(-1.1423436) q[2];
sx q[2];
rz(-2.4523297) q[2];
sx q[2];
rz(1.9826629) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.032253232) q[1];
sx q[1];
rz(-1.6620599) q[1];
sx q[1];
rz(-1.5534459) q[1];
rz(-pi) q[2];
rz(0.97753559) q[3];
sx q[3];
rz(-1.007742) q[3];
sx q[3];
rz(-2.7138591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8972299) q[2];
sx q[2];
rz(-1.4655317) q[2];
sx q[2];
rz(2.3587295) q[2];
rz(0.016294567) q[3];
sx q[3];
rz(-1.3085082) q[3];
sx q[3];
rz(-2.079594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29430729) q[0];
sx q[0];
rz(-0.48536244) q[0];
sx q[0];
rz(-2.7622727) q[0];
rz(2.8324221) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(-0.22383037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9681518) q[0];
sx q[0];
rz(-1.8814527) q[0];
sx q[0];
rz(-2.2947427) q[0];
rz(1.7621222) q[2];
sx q[2];
rz(-1.680964) q[2];
sx q[2];
rz(0.7502816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6675615) q[1];
sx q[1];
rz(-1.6769451) q[1];
sx q[1];
rz(-0.1072704) q[1];
x q[2];
rz(0.70586127) q[3];
sx q[3];
rz(-2.0076723) q[3];
sx q[3];
rz(-1.0918416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0763187) q[2];
sx q[2];
rz(-1.9741917) q[2];
sx q[2];
rz(-0.25263146) q[2];
rz(0.97918716) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(-0.34669909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52794367) q[0];
sx q[0];
rz(-2.376463) q[0];
sx q[0];
rz(-1.8494404) q[0];
rz(-0.53391236) q[1];
sx q[1];
rz(-1.4521867) q[1];
sx q[1];
rz(-0.36010489) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094771) q[0];
sx q[0];
rz(-1.5015389) q[0];
sx q[0];
rz(-1.2098625) q[0];
x q[1];
rz(-1.2525642) q[2];
sx q[2];
rz(-2.0352073) q[2];
sx q[2];
rz(2.0385252) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77006631) q[1];
sx q[1];
rz(-2.1202104) q[1];
sx q[1];
rz(-1.2868164) q[1];
rz(-1.8819767) q[3];
sx q[3];
rz(-0.10792416) q[3];
sx q[3];
rz(-0.51801658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1503633) q[2];
sx q[2];
rz(-1.5183307) q[2];
sx q[2];
rz(-0.72597996) q[2];
rz(2.7109072) q[3];
sx q[3];
rz(-0.78023282) q[3];
sx q[3];
rz(0.45346692) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78366572) q[0];
sx q[0];
rz(-1.4877321) q[0];
sx q[0];
rz(0.7861535) q[0];
rz(-0.17732009) q[1];
sx q[1];
rz(-1.0847849) q[1];
sx q[1];
rz(-2.1254983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8965137) q[0];
sx q[0];
rz(-1.4657146) q[0];
sx q[0];
rz(-1.8536293) q[0];
rz(-0.54722007) q[2];
sx q[2];
rz(-2.8946218) q[2];
sx q[2];
rz(2.9647567) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45276957) q[1];
sx q[1];
rz(-2.4773543) q[1];
sx q[1];
rz(3.1139873) q[1];
rz(1.949259) q[3];
sx q[3];
rz(-2.2871823) q[3];
sx q[3];
rz(2.5797957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5440386) q[2];
sx q[2];
rz(-0.95557135) q[2];
sx q[2];
rz(2.3717144) q[2];
rz(0.17635135) q[3];
sx q[3];
rz(-1.6917112) q[3];
sx q[3];
rz(-2.8281853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4835994) q[0];
sx q[0];
rz(-0.083881065) q[0];
sx q[0];
rz(1.0461079) q[0];
rz(-1.5931891) q[1];
sx q[1];
rz(-1.4175339) q[1];
sx q[1];
rz(1.2215325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51325996) q[0];
sx q[0];
rz(-1.9258537) q[0];
sx q[0];
rz(0.1116139) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1384954) q[2];
sx q[2];
rz(-2.7932248) q[2];
sx q[2];
rz(1.5846202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20925775) q[1];
sx q[1];
rz(-2.845721) q[1];
sx q[1];
rz(-0.81411718) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043552355) q[3];
sx q[3];
rz(-1.5756902) q[3];
sx q[3];
rz(-2.2088449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.017642411) q[2];
sx q[2];
rz(-1.4175748) q[2];
sx q[2];
rz(-1.9564015) q[2];
rz(2.8017398) q[3];
sx q[3];
rz(-0.9797107) q[3];
sx q[3];
rz(-3.0237696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3621984) q[0];
sx q[0];
rz(-1.8147991) q[0];
sx q[0];
rz(0.69778824) q[0];
rz(-2.4367874) q[1];
sx q[1];
rz(-2.0185399) q[1];
sx q[1];
rz(2.8885081) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20257178) q[0];
sx q[0];
rz(-0.24469412) q[0];
sx q[0];
rz(-1.8807202) q[0];
rz(-0.65480755) q[2];
sx q[2];
rz(-1.9960446) q[2];
sx q[2];
rz(-0.08928334) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4078057) q[1];
sx q[1];
rz(-1.1152667) q[1];
sx q[1];
rz(-0.84750073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8277824) q[3];
sx q[3];
rz(-0.8567613) q[3];
sx q[3];
rz(2.11588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7207328) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(0.68812686) q[2];
rz(-0.94528919) q[3];
sx q[3];
rz(-1.9663845) q[3];
sx q[3];
rz(2.2139886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59104334) q[0];
sx q[0];
rz(-1.782438) q[0];
sx q[0];
rz(1.8989643) q[0];
rz(-0.91724829) q[1];
sx q[1];
rz(-1.6684253) q[1];
sx q[1];
rz(-2.565276) q[1];
rz(-1.2082214) q[2];
sx q[2];
rz(-1.0248263) q[2];
sx q[2];
rz(-0.43760763) q[2];
rz(-1.9971725) q[3];
sx q[3];
rz(-1.6582499) q[3];
sx q[3];
rz(1.94699) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
