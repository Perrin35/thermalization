OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30259351) q[0];
sx q[0];
rz(-1.6380881) q[0];
sx q[0];
rz(-1.6291314) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2201266) q[2];
sx q[2];
rz(-1.7979243) q[2];
sx q[2];
rz(-0.35958689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.623917) q[1];
sx q[1];
rz(-0.78849925) q[1];
sx q[1];
rz(2.4967525) q[1];
x q[2];
rz(-0.37791032) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(-2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.4651441) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11008308) q[0];
sx q[0];
rz(-0.84658909) q[0];
sx q[0];
rz(-1.5674885) q[0];
rz(-pi) q[1];
rz(-1.0826153) q[2];
sx q[2];
rz(-1.3435875) q[2];
sx q[2];
rz(0.28085923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84491731) q[1];
sx q[1];
rz(-1.0122609) q[1];
sx q[1];
rz(1.0122453) q[1];
x q[2];
rz(-1.2259237) q[3];
sx q[3];
rz(-2.2861835) q[3];
sx q[3];
rz(-1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-2.6254568) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(-0.80054545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0101937) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(1.2244768) q[0];
rz(-pi) q[1];
rz(-2.4750701) q[2];
sx q[2];
rz(-2.1608673) q[2];
sx q[2];
rz(-0.18596622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2785981) q[1];
sx q[1];
rz(-1.1578214) q[1];
sx q[1];
rz(-1.5652656) q[1];
x q[2];
rz(1.9242026) q[3];
sx q[3];
rz(-2.4089775) q[3];
sx q[3];
rz(-0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.7819972) q[2];
rz(2.1740186) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(1.4250071) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-2.676679) q[0];
rz(2.7930296) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0264498) q[0];
sx q[0];
rz(-2.0949674) q[0];
sx q[0];
rz(-1.8379704) q[0];
x q[1];
rz(1.9374574) q[2];
sx q[2];
rz(-1.3385834) q[2];
sx q[2];
rz(2.1249078) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76772056) q[1];
sx q[1];
rz(-1.4704736) q[1];
sx q[1];
rz(-2.9994681) q[1];
rz(-pi) q[2];
rz(-0.61120175) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-0.51182169) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27914771) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(-0.43235835) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3571346) q[0];
sx q[0];
rz(-1.1103837) q[0];
sx q[0];
rz(2.5255192) q[0];
rz(-1.8129187) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(-2.9850933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2213858) q[1];
sx q[1];
rz(-1.1077987) q[1];
sx q[1];
rz(1.7358001) q[1];
rz(3.0675689) q[3];
sx q[3];
rz(-0.97462666) q[3];
sx q[3];
rz(-2.1285469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.3283407) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(2.7691675) q[0];
rz(1.8107481) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(0.16528027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344791) q[0];
sx q[0];
rz(-1.665859) q[0];
sx q[0];
rz(1.6001742) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.084947305) q[2];
sx q[2];
rz(-0.68945486) q[2];
sx q[2];
rz(-1.8817608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2416934) q[1];
sx q[1];
rz(-0.69679931) q[1];
sx q[1];
rz(-0.65710575) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9061756) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(-0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2281987) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.6983263) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(0.73648891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6155646) q[0];
sx q[0];
rz(-0.051858735) q[0];
sx q[0];
rz(1.9066914) q[0];
x q[1];
rz(2.8924106) q[2];
sx q[2];
rz(-1.2417925) q[2];
sx q[2];
rz(0.56275425) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(2.3751395) q[1];
rz(-pi) q[2];
rz(1.1398846) q[3];
sx q[3];
rz(-1.1093372) q[3];
sx q[3];
rz(2.0737089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0050469) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(2.2198026) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1944273) q[0];
sx q[0];
rz(-0.88473407) q[0];
sx q[0];
rz(1.3961193) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8554101) q[2];
sx q[2];
rz(-0.85368644) q[2];
sx q[2];
rz(-1.4746427) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39107716) q[1];
sx q[1];
rz(-2.8120496) q[1];
sx q[1];
rz(-2.4604172) q[1];
rz(-pi) q[2];
rz(0.86254085) q[3];
sx q[3];
rz(-1.7077912) q[3];
sx q[3];
rz(-0.18602895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.1716589) q[2];
rz(1.3575859) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.6660447) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(2.4618861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8863227) q[0];
sx q[0];
rz(-0.5605883) q[0];
sx q[0];
rz(-0.38829304) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7858511) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(0.099345318) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45021536) q[1];
sx q[1];
rz(-0.9141578) q[1];
sx q[1];
rz(2.0037829) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1867832) q[3];
sx q[3];
rz(-2.1534352) q[3];
sx q[3];
rz(-0.59347502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(-2.4662468) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-0.85062406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.1766599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51417527) q[0];
sx q[0];
rz(-1.5222933) q[0];
sx q[0];
rz(1.4072627) q[0];
rz(-pi) q[1];
rz(0.87877019) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(-2.9825485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51917045) q[1];
sx q[1];
rz(-1.6455489) q[1];
sx q[1];
rz(1.012411) q[1];
rz(-1.6586187) q[3];
sx q[3];
rz(-1.0613958) q[3];
sx q[3];
rz(-2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-2.0142377) q[2];
rz(-0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.42416278) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-2.7453616) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(1.2926119) q[2];
sx q[2];
rz(-2.2879911) q[2];
sx q[2];
rz(0.0030980274) q[2];
rz(2.7908294) q[3];
sx q[3];
rz(-1.5784932) q[3];
sx q[3];
rz(-2.2898883) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];