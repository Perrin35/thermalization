OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(-2.8470319) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.784214) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(-1.635301) q[0];
x q[1];
rz(-0.047702567) q[2];
sx q[2];
rz(-2.45445) q[2];
sx q[2];
rz(1.7575761) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4946343) q[1];
sx q[1];
rz(-2.0152433) q[1];
sx q[1];
rz(1.7723945) q[1];
rz(-pi) q[2];
rz(-1.7105605) q[3];
sx q[3];
rz(-2.1602109) q[3];
sx q[3];
rz(-2.0624269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(8*pi/11) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(0.37047085) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1169491) q[0];
sx q[0];
rz(-1.2632217) q[0];
sx q[0];
rz(-1.3427539) q[0];
rz(-0.53441647) q[2];
sx q[2];
rz(-0.60306163) q[2];
sx q[2];
rz(2.7433928) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49711455) q[1];
sx q[1];
rz(-1.2116355) q[1];
sx q[1];
rz(-2.4317846) q[1];
rz(-pi) q[2];
rz(2.6133735) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(-1.6268829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.52477437) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(-1.3735501) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4931902) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(-0.52669749) q[0];
rz(0.85767304) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(-2.8785994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92257567) q[1];
sx q[1];
rz(-0.95239988) q[1];
sx q[1];
rz(-0.69505691) q[1];
rz(1.340629) q[3];
sx q[3];
rz(-1.3535415) q[3];
sx q[3];
rz(-2.1118856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(-2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(-1.5377195) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(0.61027169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.423324) q[0];
sx q[0];
rz(-1.4615131) q[0];
sx q[0];
rz(-3.1302343) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9539815) q[2];
sx q[2];
rz(-1.9028579) q[2];
sx q[2];
rz(0.039475723) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8205386) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(1.1827521) q[1];
rz(-0.9917594) q[3];
sx q[3];
rz(-1.990253) q[3];
sx q[3];
rz(1.5031682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(1.3132494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6812692) q[0];
sx q[0];
rz(-1.2460421) q[0];
sx q[0];
rz(-0.89694174) q[0];
rz(-3.1057642) q[2];
sx q[2];
rz(-1.6275068) q[2];
sx q[2];
rz(0.63024516) q[2];
rz(-pi) q[3];
x q[3];
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
rz(-0.58532183) q[3];
sx q[3];
rz(-0.96671852) q[3];
sx q[3];
rz(-0.6795336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(-2.4364566) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.054811906) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(2.6584113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9478215) q[0];
sx q[0];
rz(-1.3894488) q[0];
sx q[0];
rz(-3.0788455) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1081084) q[2];
sx q[2];
rz(-1.5414642) q[2];
sx q[2];
rz(0.5043475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15258458) q[1];
sx q[1];
rz(-2.2629988) q[1];
sx q[1];
rz(-2.2425376) q[1];
rz(-0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(1.9432817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9873535) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(-0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(-2.5700991) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(2.8894997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1958256) q[0];
sx q[0];
rz(-0.44647631) q[0];
sx q[0];
rz(-2.8004942) q[0];
x q[1];
rz(0.7438296) q[2];
sx q[2];
rz(-0.93223909) q[2];
sx q[2];
rz(-0.078787412) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5735057) q[1];
sx q[1];
rz(-1.9316257) q[1];
sx q[1];
rz(-0.87330694) q[1];
rz(-2.7639286) q[3];
sx q[3];
rz(-2.0026243) q[3];
sx q[3];
rz(0.8197195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(0.90448109) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(-0.65902567) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(1.4138387) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.3716912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1641952) q[0];
sx q[0];
rz(-1.482555) q[0];
sx q[0];
rz(0.95055687) q[0];
rz(-pi) q[1];
rz(2.3805982) q[2];
sx q[2];
rz(-1.7241663) q[2];
sx q[2];
rz(-0.93190565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9904069) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(2.5469261) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7226108) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(0.23432215) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-2.871002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9070248) q[0];
sx q[0];
rz(-1.7922635) q[0];
sx q[0];
rz(1.7937167) q[0];
rz(-pi) q[1];
rz(-0.19240304) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(-0.88142384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8671659) q[1];
sx q[1];
rz(-1.5741036) q[1];
sx q[1];
rz(3.0043688) q[1];
rz(1.7646592) q[3];
sx q[3];
rz(-2.6541775) q[3];
sx q[3];
rz(-0.91538115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(-1.6516997) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(-0.94296304) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(2.399209) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6615804) q[0];
sx q[0];
rz(-0.98012692) q[0];
sx q[0];
rz(-2.3175879) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43138357) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(-2.0768349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2303289) q[1];
sx q[1];
rz(-1.9165557) q[1];
sx q[1];
rz(-1.5589578) q[1];
rz(-pi) q[2];
rz(2.9418482) q[3];
sx q[3];
rz(-0.14530694) q[3];
sx q[3];
rz(-0.47606459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(2.6081086) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(0.82143299) q[2];
sx q[2];
rz(-1.8100304) q[2];
sx q[2];
rz(2.9678154) q[2];
rz(-0.13989862) q[3];
sx q[3];
rz(-1.7053053) q[3];
sx q[3];
rz(-2.7444621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];