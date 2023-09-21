OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.598782) q[0];
sx q[0];
rz(-3.071975) q[0];
sx q[0];
rz(1.1845864) q[0];
rz(-0.047702567) q[2];
sx q[2];
rz(-2.45445) q[2];
sx q[2];
rz(1.7575761) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0114853) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(-2.6891522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59395091) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(-2.5719197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(-0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37061781) q[0];
sx q[0];
rz(-0.3807225) q[0];
sx q[0];
rz(-2.5230663) q[0];
rz(-pi) q[1];
rz(-1.2334521) q[2];
sx q[2];
rz(-1.0609027) q[2];
sx q[2];
rz(-2.9166729) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6795923) q[1];
sx q[1];
rz(-0.781171) q[1];
sx q[1];
rz(0.52266927) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46918418) q[3];
sx q[3];
rz(-2.5622534) q[3];
sx q[3];
rz(-2.7964696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-0.16361374) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(-1.7680426) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4931902) q[0];
sx q[0];
rz(-1.5957386) q[0];
sx q[0];
rz(-0.52669749) q[0];
rz(0.85767304) q[2];
sx q[2];
rz(-1.3366941) q[2];
sx q[2];
rz(2.8785994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.219017) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(0.69505691) q[1];
rz(-pi) q[2];
rz(0.80188607) q[3];
sx q[3];
rz(-2.826414) q[3];
sx q[3];
rz(-1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0107161) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(0.61292928) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-0.83909488) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(0.61027169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82204098) q[0];
sx q[0];
rz(-3.0317231) q[0];
sx q[0];
rz(1.6739474) q[0];
rz(-pi) q[1];
rz(-2.7856366) q[2];
sx q[2];
rz(-1.2095371) q[2];
sx q[2];
rz(-1.4796096) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0415045) q[1];
sx q[1];
rz(-1.8309635) q[1];
sx q[1];
rz(0.86160223) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9917594) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.5031682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(-2.9042517) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(1.3132494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816313) q[0];
sx q[0];
rz(-2.2035723) q[0];
sx q[0];
rz(-2.7347793) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1336742) q[2];
sx q[2];
rz(-0.067069947) q[2];
sx q[2];
rz(-1.1941393) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4876731) q[1];
sx q[1];
rz(-0.54033414) q[1];
sx q[1];
rz(0.78507702) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89574121) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(-0.18273396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076684549) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(-3.0867807) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(0.48318133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85855243) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(1.900308) q[0];
rz(1.5051571) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(1.0077196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.253787) q[1];
sx q[1];
rz(-1.0711545) q[1];
sx q[1];
rz(-2.3274724) q[1];
rz(-pi) q[2];
rz(1.9032848) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(0.40963848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(0.25209299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.945767) q[0];
sx q[0];
rz(-0.44647631) q[0];
sx q[0];
rz(2.8004942) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78105314) q[2];
sx q[2];
rz(-2.1456246) q[2];
sx q[2];
rz(-0.99036723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.568087) q[1];
sx q[1];
rz(-1.9316257) q[1];
sx q[1];
rz(-0.87330694) q[1];
x q[2];
rz(1.1105359) q[3];
sx q[3];
rz(-1.9122951) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.38348848) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7082108) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.4138387) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.3716912) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28384128) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(1.7218504) q[0];
rz(-pi) q[1];
rz(0.76099446) q[2];
sx q[2];
rz(-1.4174263) q[2];
sx q[2];
rz(-0.93190565) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.614537) q[1];
sx q[1];
rz(-1.0053047) q[1];
sx q[1];
rz(1.2177699) q[1];
x q[2];
rz(-1.4189818) q[3];
sx q[3];
rz(-1.2137128) q[3];
sx q[3];
rz(-0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82683212) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(-0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9234377) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.6802457) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(-0.27059069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9070248) q[0];
sx q[0];
rz(-1.3493291) q[0];
sx q[0];
rz(-1.7937167) q[0];
rz(1.2241237) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(-2.5172362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2968263) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(-1.5674577) q[1];
rz(-1.7646592) q[3];
sx q[3];
rz(-2.6541775) q[3];
sx q[3];
rz(-2.2262115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(-1.489893) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437317) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(2.399209) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091632) q[0];
sx q[0];
rz(-2.2262648) q[0];
sx q[0];
rz(2.3497033) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(-1.0647578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91126373) q[1];
sx q[1];
rz(-1.2250369) q[1];
sx q[1];
rz(-1.5589578) q[1];
rz(-pi) q[2];
rz(2.9418482) q[3];
sx q[3];
rz(-0.14530694) q[3];
sx q[3];
rz(2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8515274) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(0.32158357) q[2];
sx q[2];
rz(-2.2939773) q[2];
sx q[2];
rz(1.1800223) q[2];
rz(0.7704173) q[3];
sx q[3];
rz(-0.19376783) q[3];
sx q[3];
rz(2.7289058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];