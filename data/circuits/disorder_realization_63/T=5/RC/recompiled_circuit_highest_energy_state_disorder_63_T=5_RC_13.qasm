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
rz(2.0923848) q[0];
sx q[0];
rz(-0.75924358) q[0];
sx q[0];
rz(2.7651751) q[0];
rz(-1.5837826) q[1];
sx q[1];
rz(4.2111068) q[1];
sx q[1];
rz(11.61011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1726748) q[0];
sx q[0];
rz(-0.71828523) q[0];
sx q[0];
rz(2.0613614) q[0];
rz(-pi) q[1];
rz(-2.7306284) q[2];
sx q[2];
rz(-2.685475) q[2];
sx q[2];
rz(-0.48738313) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8531283) q[1];
sx q[1];
rz(-2.3061516) q[1];
sx q[1];
rz(-2.9771027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.718456) q[3];
sx q[3];
rz(-2.7877533) q[3];
sx q[3];
rz(-1.906499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3236397) q[2];
sx q[2];
rz(-0.91922593) q[2];
sx q[2];
rz(-1.2032571) q[2];
rz(1.0407) q[3];
sx q[3];
rz(-2.2668656) q[3];
sx q[3];
rz(3.021595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0284718) q[0];
sx q[0];
rz(-2.5682243) q[0];
sx q[0];
rz(-2.0290802) q[0];
rz(1.6587229) q[1];
sx q[1];
rz(-0.590938) q[1];
sx q[1];
rz(0.15731752) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293662) q[0];
sx q[0];
rz(-0.88878814) q[0];
sx q[0];
rz(-1.3444882) q[0];
x q[1];
rz(2.7685086) q[2];
sx q[2];
rz(-1.6519128) q[2];
sx q[2];
rz(-0.099992601) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29780218) q[1];
sx q[1];
rz(-2.4239967) q[1];
sx q[1];
rz(-1.9466496) q[1];
x q[2];
rz(-0.3213082) q[3];
sx q[3];
rz(-2.9071593) q[3];
sx q[3];
rz(0.050416273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4979672) q[2];
sx q[2];
rz(-2.6622541) q[2];
sx q[2];
rz(0.41110006) q[2];
rz(2.5655365) q[3];
sx q[3];
rz(-2.1274302) q[3];
sx q[3];
rz(0.99803734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75329798) q[0];
sx q[0];
rz(-2.1385215) q[0];
sx q[0];
rz(0.74111795) q[0];
rz(1.6916212) q[1];
sx q[1];
rz(-1.3131817) q[1];
sx q[1];
rz(-0.57201874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93942552) q[0];
sx q[0];
rz(-1.8644445) q[0];
sx q[0];
rz(-2.9333326) q[0];
rz(-pi) q[1];
rz(2.4837844) q[2];
sx q[2];
rz(-0.33051046) q[2];
sx q[2];
rz(-0.1275488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3254889) q[1];
sx q[1];
rz(-1.1469335) q[1];
sx q[1];
rz(-0.44386379) q[1];
x q[2];
rz(2.1628863) q[3];
sx q[3];
rz(-2.7436069) q[3];
sx q[3];
rz(-0.0090816895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5570306) q[2];
sx q[2];
rz(-1.9494373) q[2];
sx q[2];
rz(-0.7977879) q[2];
rz(-1.8320463) q[3];
sx q[3];
rz(-1.2582425) q[3];
sx q[3];
rz(-1.7358739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13084594) q[0];
sx q[0];
rz(-0.89711419) q[0];
sx q[0];
rz(0.36706269) q[0];
rz(1.2182073) q[1];
sx q[1];
rz(-1.7299078) q[1];
sx q[1];
rz(0.22148111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9629134) q[0];
sx q[0];
rz(-1.7811462) q[0];
sx q[0];
rz(-1.5014929) q[0];
rz(-2.7985953) q[2];
sx q[2];
rz(-0.38831899) q[2];
sx q[2];
rz(1.266154) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74879157) q[1];
sx q[1];
rz(-1.1587156) q[1];
sx q[1];
rz(-1.2296299) q[1];
rz(-0.4226004) q[3];
sx q[3];
rz(-2.5739193) q[3];
sx q[3];
rz(-1.6740396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9941142) q[2];
sx q[2];
rz(-1.9693815) q[2];
sx q[2];
rz(1.8032903) q[2];
rz(-2.639468) q[3];
sx q[3];
rz(-0.10052557) q[3];
sx q[3];
rz(0.24371915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62622708) q[0];
sx q[0];
rz(-0.60972649) q[0];
sx q[0];
rz(1.6060265) q[0];
rz(-3.0961127) q[1];
sx q[1];
rz(-1.7605503) q[1];
sx q[1];
rz(-0.46612003) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461994) q[0];
sx q[0];
rz(-1.3211622) q[0];
sx q[0];
rz(-0.14139463) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.999249) q[2];
sx q[2];
rz(-0.68926297) q[2];
sx q[2];
rz(1.9826629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.032253232) q[1];
sx q[1];
rz(-1.6620599) q[1];
sx q[1];
rz(-1.5534459) q[1];
x q[2];
rz(-2.4908972) q[3];
sx q[3];
rz(-1.078372) q[3];
sx q[3];
rz(-0.7975815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2443627) q[2];
sx q[2];
rz(-1.4655317) q[2];
sx q[2];
rz(-2.3587295) q[2];
rz(-3.1252981) q[3];
sx q[3];
rz(-1.8330845) q[3];
sx q[3];
rz(2.079594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8472854) q[0];
sx q[0];
rz(-0.48536244) q[0];
sx q[0];
rz(2.7622727) q[0];
rz(-0.30917057) q[1];
sx q[1];
rz(-1.199017) q[1];
sx q[1];
rz(2.9177623) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081793) q[0];
sx q[0];
rz(-0.88841618) q[0];
sx q[0];
rz(0.40485415) q[0];
rz(-pi) q[1];
rz(1.3794704) q[2];
sx q[2];
rz(-1.680964) q[2];
sx q[2];
rz(-0.7502816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8221082) q[1];
sx q[1];
rz(-0.15076877) q[1];
sx q[1];
rz(2.3586078) q[1];
rz(-pi) q[2];
rz(0.70586127) q[3];
sx q[3];
rz(-2.0076723) q[3];
sx q[3];
rz(-1.0918416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0763187) q[2];
sx q[2];
rz(-1.167401) q[2];
sx q[2];
rz(0.25263146) q[2];
rz(-2.1624055) q[3];
sx q[3];
rz(-0.88117176) q[3];
sx q[3];
rz(2.7948936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.613649) q[0];
sx q[0];
rz(-2.376463) q[0];
sx q[0];
rz(1.2921523) q[0];
rz(-2.6076803) q[1];
sx q[1];
rz(-1.4521867) q[1];
sx q[1];
rz(-2.7814878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8468854) q[0];
sx q[0];
rz(-0.36723235) q[0];
sx q[0];
rz(-1.7647554) q[0];
rz(-1.2525642) q[2];
sx q[2];
rz(-1.1063853) q[2];
sx q[2];
rz(-2.0385252) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3715263) q[1];
sx q[1];
rz(-1.0213823) q[1];
sx q[1];
rz(-1.2868164) q[1];
rz(0.033161226) q[3];
sx q[3];
rz(-1.4680742) q[3];
sx q[3];
rz(-0.20513137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1503633) q[2];
sx q[2];
rz(-1.5183307) q[2];
sx q[2];
rz(-2.4156127) q[2];
rz(-0.43068543) q[3];
sx q[3];
rz(-2.3613598) q[3];
sx q[3];
rz(-0.45346692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3579269) q[0];
sx q[0];
rz(-1.4877321) q[0];
sx q[0];
rz(0.7861535) q[0];
rz(2.9642726) q[1];
sx q[1];
rz(-1.0847849) q[1];
sx q[1];
rz(-2.1254983) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02064146) q[0];
sx q[0];
rz(-0.30123152) q[0];
sx q[0];
rz(1.2094638) q[0];
x q[1];
rz(-2.5943726) q[2];
sx q[2];
rz(-0.24697082) q[2];
sx q[2];
rz(2.9647567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0018277) q[1];
sx q[1];
rz(-1.5878126) q[1];
sx q[1];
rz(-2.4775392) q[1];
x q[2];
rz(2.7402596) q[3];
sx q[3];
rz(-0.79425967) q[3];
sx q[3];
rz(-3.1242736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5440386) q[2];
sx q[2];
rz(-2.1860213) q[2];
sx q[2];
rz(-2.3717144) q[2];
rz(2.9652413) q[3];
sx q[3];
rz(-1.6917112) q[3];
sx q[3];
rz(2.8281853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65799323) q[0];
sx q[0];
rz(-3.0577116) q[0];
sx q[0];
rz(-1.0461079) q[0];
rz(1.5484035) q[1];
sx q[1];
rz(-1.4175339) q[1];
sx q[1];
rz(-1.9200602) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51325996) q[0];
sx q[0];
rz(-1.2157389) q[0];
sx q[0];
rz(-3.0299788) q[0];
rz(-pi) q[1];
rz(2.0030973) q[2];
sx q[2];
rz(-2.7932248) q[2];
sx q[2];
rz(1.5569725) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56961774) q[1];
sx q[1];
rz(-1.7844266) q[1];
sx q[1];
rz(-0.20628449) q[1];
x q[2];
rz(1.5658978) q[3];
sx q[3];
rz(-1.5272445) q[3];
sx q[3];
rz(2.5033308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.017642411) q[2];
sx q[2];
rz(-1.7240179) q[2];
sx q[2];
rz(1.1851912) q[2];
rz(2.8017398) q[3];
sx q[3];
rz(-2.161882) q[3];
sx q[3];
rz(3.0237696) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7793943) q[0];
sx q[0];
rz(-1.8147991) q[0];
sx q[0];
rz(-0.69778824) q[0];
rz(2.4367874) q[1];
sx q[1];
rz(-2.0185399) q[1];
sx q[1];
rz(-2.8885081) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4721203) q[0];
sx q[0];
rz(-1.4968431) q[0];
sx q[0];
rz(1.3373313) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0519876) q[2];
sx q[2];
rz(-0.98271433) q[2];
sx q[2];
rz(1.1747509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.733787) q[1];
sx q[1];
rz(-1.1152667) q[1];
sx q[1];
rz(-0.84750073) q[1];
rz(-pi) q[2];
rz(-1.3138103) q[3];
sx q[3];
rz(-0.8567613) q[3];
sx q[3];
rz(-1.0257126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7207328) q[2];
sx q[2];
rz(-1.3497817) q[2];
sx q[2];
rz(2.4534658) q[2];
rz(2.1963035) q[3];
sx q[3];
rz(-1.1752081) q[3];
sx q[3];
rz(-2.2139886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5505493) q[0];
sx q[0];
rz(-1.3591546) q[0];
sx q[0];
rz(-1.2426283) q[0];
rz(0.91724829) q[1];
sx q[1];
rz(-1.4731673) q[1];
sx q[1];
rz(0.57631667) q[1];
rz(-2.5653432) q[2];
sx q[2];
rz(-1.8787619) q[2];
sx q[2];
rz(0.93870434) q[2];
rz(0.096002738) q[3];
sx q[3];
rz(-1.9954372) q[3];
sx q[3];
rz(0.33653997) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
