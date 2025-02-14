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
rz(-2.0353844) q[0];
sx q[0];
rz(-0.68618965) q[0];
sx q[0];
rz(1.9880779) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(0.44841132) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8657239) q[0];
sx q[0];
rz(-0.65058904) q[0];
sx q[0];
rz(-3.0537729) q[0];
rz(-pi) q[1];
x q[1];
rz(2.885778) q[2];
sx q[2];
rz(-1.668715) q[2];
sx q[2];
rz(0.73543392) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4223215) q[1];
sx q[1];
rz(-1.0485253) q[1];
sx q[1];
rz(0.85776718) q[1];
rz(-0.087440447) q[3];
sx q[3];
rz(-2.6387847) q[3];
sx q[3];
rz(1.3399762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-2.5173729) q[2];
rz(0.37912399) q[3];
sx q[3];
rz(-0.99025327) q[3];
sx q[3];
rz(-0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666331) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(1.0354743) q[0];
rz(1.8677208) q[1];
sx q[1];
rz(-2.404411) q[1];
sx q[1];
rz(2.6587291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.567822) q[0];
sx q[0];
rz(-1.5250686) q[0];
sx q[0];
rz(3.0450495) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4522885) q[2];
sx q[2];
rz(-1.7359455) q[2];
sx q[2];
rz(1.120795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.711593) q[1];
sx q[1];
rz(-1.5801799) q[1];
sx q[1];
rz(1.1389772) q[1];
rz(-1.489352) q[3];
sx q[3];
rz(-1.2508498) q[3];
sx q[3];
rz(2.4575352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.024293385) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(1.8710322) q[2];
rz(-2.471586) q[3];
sx q[3];
rz(-2.5195401) q[3];
sx q[3];
rz(-2.2445934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4033177) q[0];
sx q[0];
rz(-2.6955695) q[0];
sx q[0];
rz(2.4025412) q[0];
rz(-0.96145472) q[1];
sx q[1];
rz(-0.32197222) q[1];
sx q[1];
rz(2.3846073) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3630545) q[0];
sx q[0];
rz(-0.86150384) q[0];
sx q[0];
rz(2.7947133) q[0];
x q[1];
rz(1.1951094) q[2];
sx q[2];
rz(-0.91478148) q[2];
sx q[2];
rz(3.0041848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.095159141) q[1];
sx q[1];
rz(-0.94645247) q[1];
sx q[1];
rz(1.4051564) q[1];
rz(-3.0717172) q[3];
sx q[3];
rz(-1.2638076) q[3];
sx q[3];
rz(-1.2941293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5028533) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(-0.70510954) q[2];
rz(2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(2.7307935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.131677) q[0];
sx q[0];
rz(-0.83808815) q[0];
sx q[0];
rz(-1.9158844) q[0];
rz(1.907584) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(-0.066224901) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28540885) q[0];
sx q[0];
rz(-1.3254306) q[0];
sx q[0];
rz(-0.092459926) q[0];
x q[1];
rz(-2.8429864) q[2];
sx q[2];
rz(-2.3176607) q[2];
sx q[2];
rz(-0.9048942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5410616) q[1];
sx q[1];
rz(-2.5503073) q[1];
sx q[1];
rz(-0.97452428) q[1];
rz(2.8230571) q[3];
sx q[3];
rz(-0.38357601) q[3];
sx q[3];
rz(0.84795241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0951198) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(2.7009916) q[2];
rz(1.0062086) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.7202268) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(-1.5059858) q[0];
rz(-1.5274564) q[1];
sx q[1];
rz(-1.9215877) q[1];
sx q[1];
rz(-1.45586) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9855749) q[0];
sx q[0];
rz(-0.7472207) q[0];
sx q[0];
rz(1.78563) q[0];
rz(-2.9142889) q[2];
sx q[2];
rz(-1.4571726) q[2];
sx q[2];
rz(2.794968) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.33259049) q[1];
sx q[1];
rz(-2.2223516) q[1];
sx q[1];
rz(-1.5452191) q[1];
rz(-pi) q[2];
rz(2.9221228) q[3];
sx q[3];
rz(-2.5220036) q[3];
sx q[3];
rz(2.0989024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8684034) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(-0.91147649) q[2];
rz(-0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8580496) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-0.13033303) q[0];
rz(2.6538972) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(-0.74388751) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.795563) q[0];
sx q[0];
rz(-1.6206121) q[0];
sx q[0];
rz(1.4883785) q[0];
rz(-pi) q[1];
rz(1.7623474) q[2];
sx q[2];
rz(-0.63427502) q[2];
sx q[2];
rz(-2.3630362) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6106657) q[1];
sx q[1];
rz(-2.1082889) q[1];
sx q[1];
rz(-2.8100138) q[1];
rz(-pi) q[2];
rz(1.3906995) q[3];
sx q[3];
rz(-2.0493747) q[3];
sx q[3];
rz(0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1193715) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(-0.96088299) q[2];
rz(2.7458701) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457086) q[0];
sx q[0];
rz(-1.1154024) q[0];
sx q[0];
rz(-2.0945666) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(2.7488757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1395124) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(2.9154791) q[0];
rz(-pi) q[1];
rz(-0.057633295) q[2];
sx q[2];
rz(-0.95256348) q[2];
sx q[2];
rz(-2.1098441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9153292) q[1];
sx q[1];
rz(-1.8460423) q[1];
sx q[1];
rz(-1.54285) q[1];
rz(-pi) q[2];
rz(-2.4684577) q[3];
sx q[3];
rz(-0.41134787) q[3];
sx q[3];
rz(-0.92837161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.087223209) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(-1.8325904) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(2.8712809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.8488309) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(-2.8344179) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-0.38506404) q[1];
sx q[1];
rz(-0.90075341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9884315) q[0];
sx q[0];
rz(-1.5539196) q[0];
sx q[0];
rz(2.8794226) q[0];
rz(-pi) q[1];
rz(-2.2239752) q[2];
sx q[2];
rz(-1.2003044) q[2];
sx q[2];
rz(-2.8965184) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.72462725) q[1];
sx q[1];
rz(-2.6748383) q[1];
sx q[1];
rz(2.3326567) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6314478) q[3];
sx q[3];
rz(-1.6854146) q[3];
sx q[3];
rz(1.7607911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8073392) q[2];
sx q[2];
rz(-3.0901577) q[2];
sx q[2];
rz(-0.61484289) q[2];
rz(-2.0452512) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(-2.443327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480302) q[0];
sx q[0];
rz(-0.63069558) q[0];
sx q[0];
rz(-0.50931859) q[0];
rz(-2.9604984) q[1];
sx q[1];
rz(-1.7362005) q[1];
sx q[1];
rz(2.1720355) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80705331) q[0];
sx q[0];
rz(-1.7248123) q[0];
sx q[0];
rz(-2.3898983) q[0];
rz(-pi) q[1];
rz(-1.4793899) q[2];
sx q[2];
rz(-2.6326523) q[2];
sx q[2];
rz(1.7569052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0694094) q[1];
sx q[1];
rz(-1.2320215) q[1];
sx q[1];
rz(1.6852001) q[1];
rz(-pi) q[2];
rz(1.182231) q[3];
sx q[3];
rz(-0.7668743) q[3];
sx q[3];
rz(0.23393133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8934882) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-2.060037) q[2];
rz(-0.23165101) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(2.933568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.76081) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(-0.64714062) q[0];
rz(0.45516792) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(2.3796577) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3477701) q[0];
sx q[0];
rz(-0.88632727) q[0];
sx q[0];
rz(1.1948299) q[0];
rz(-0.64024957) q[2];
sx q[2];
rz(-0.92776042) q[2];
sx q[2];
rz(1.4852448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7358491) q[1];
sx q[1];
rz(-0.65913504) q[1];
sx q[1];
rz(0.05593158) q[1];
x q[2];
rz(-1.5555218) q[3];
sx q[3];
rz(-0.46746436) q[3];
sx q[3];
rz(-1.338442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.053293856) q[2];
sx q[2];
rz(-0.21685728) q[2];
sx q[2];
rz(2.2376412) q[2];
rz(1.6186742) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(-1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.224613) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(-0.12679535) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(2.0218718) q[2];
sx q[2];
rz(-1.3430165) q[2];
sx q[2];
rz(-2.7271885) q[2];
rz(0.42999646) q[3];
sx q[3];
rz(-1.5681793) q[3];
sx q[3];
rz(1.5344686) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
