OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(2.4584682) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0177512) q[0];
sx q[0];
rz(-0.089028247) q[0];
sx q[0];
rz(-0.71319367) q[0];
rz(-pi) q[1];
x q[1];
rz(1.921466) q[2];
sx q[2];
rz(-1.3436683) q[2];
sx q[2];
rz(0.35958689) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4343623) q[1];
sx q[1];
rz(-2.011236) q[1];
sx q[1];
rz(0.67727725) q[1];
rz(-0.37791032) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(0.4370673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(2.0430298) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.4651441) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4629048) q[0];
sx q[0];
rz(-1.5683187) q[0];
sx q[0];
rz(0.72420995) q[0];
x q[1];
rz(2.0589774) q[2];
sx q[2];
rz(-1.3435875) q[2];
sx q[2];
rz(-2.8607334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0959024) q[1];
sx q[1];
rz(-1.1046788) q[1];
sx q[1];
rz(2.5065266) q[1];
rz(-0.74554262) q[3];
sx q[3];
rz(-1.8288444) q[3];
sx q[3];
rz(-2.6170956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(3.0318276) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.4573147) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(2.6254568) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-0.80054545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24414177) q[0];
sx q[0];
rz(-1.7139009) q[0];
sx q[0];
rz(-1.98154) q[0];
x q[1];
rz(-0.86513743) q[2];
sx q[2];
rz(-2.1102998) q[2];
sx q[2];
rz(1.3441966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4360106) q[1];
sx q[1];
rz(-1.5657305) q[1];
sx q[1];
rz(-0.41298053) q[1];
rz(-pi) q[2];
rz(-1.21739) q[3];
sx q[3];
rz(-0.73261515) q[3];
sx q[3];
rz(0.64980799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67733726) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99252218) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(-0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(2.0565313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5497919) q[0];
sx q[0];
rz(-1.8013957) q[0];
sx q[0];
rz(2.6016298) q[0];
x q[1];
rz(-0.98767878) q[2];
sx q[2];
rz(-0.43118011) q[2];
sx q[2];
rz(1.093986) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19238732) q[1];
sx q[1];
rz(-0.17377033) q[1];
sx q[1];
rz(-2.5237571) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5303909) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(-1.9577648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(-0.26369035) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27914771) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(-2.7092343) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-0.98168215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9156993) q[0];
sx q[0];
rz(-2.3908273) q[0];
sx q[0];
rz(2.4322926) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.328674) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(2.9850933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.86451929) q[1];
sx q[1];
rz(-0.48950567) q[1];
sx q[1];
rz(-0.31788748) q[1];
rz(0.074023789) q[3];
sx q[3];
rz(-0.97462666) q[3];
sx q[3];
rz(2.1285469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-2.6521818) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(-1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(0.16528027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0342456) q[0];
sx q[0];
rz(-0.099485569) q[0];
sx q[0];
rz(-0.29883595) q[0];
rz(2.4539102) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(0.24535594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11152553) q[1];
sx q[1];
rz(-1.0377874) q[1];
sx q[1];
rz(2.0433321) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0323456) q[3];
sx q[3];
rz(-0.49406067) q[3];
sx q[3];
rz(0.1915313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.4432663) q[2];
rz(0.36744395) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(0.2120367) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-1.0914047) q[0];
sx q[0];
rz(-2.7923287) q[0];
rz(-2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-0.73648891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9518785) q[0];
sx q[0];
rz(-1.5218381) q[0];
sx q[0];
rz(-3.1244856) q[0];
rz(-pi) q[1];
rz(1.9095124) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(0.92600694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-2.3751395) q[1];
x q[2];
rz(0.50076671) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-0.70481833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625967) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(0.92179006) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9471654) q[0];
sx q[0];
rz(-2.2568586) q[0];
sx q[0];
rz(-1.7454733) q[0];
x q[1];
rz(1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(-1.2459754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31729749) q[1];
sx q[1];
rz(-1.3166787) q[1];
sx q[1];
rz(-1.3586678) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7796302) q[3];
sx q[3];
rz(-0.71912557) q[3];
sx q[3];
rz(-1.5428839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2395997) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.6931036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(0.67970651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.7734818) q[0];
sx q[0];
rz(0.52635877) q[0];
rz(-pi) q[1];
rz(-1.7858511) q[2];
sx q[2];
rz(-2.4348223) q[2];
sx q[2];
rz(3.0422473) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3956086) q[1];
sx q[1];
rz(-1.2320227) q[1];
sx q[1];
rz(-2.4376274) q[1];
rz(-2.4217442) q[3];
sx q[3];
rz(-2.3206629) q[3];
sx q[3];
rz(1.6380701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-2.2311907) q[2];
rz(-2.4662468) q[3];
sx q[3];
rz(-2.2048435) q[3];
sx q[3];
rz(-2.2909686) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(0.7146548) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.1766599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0929716) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(-3.0924348) q[0];
rz(-2.2628225) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(0.15904418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1365876) q[1];
sx q[1];
rz(-1.0141546) q[1];
sx q[1];
rz(0.088076061) q[1];
x q[2];
rz(1.482974) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(-1.0305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(1.127355) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-2.2704411) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(0.30504967) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(-1.5626004) q[3];
sx q[3];
rz(-1.9215487) q[3];
sx q[3];
rz(2.4196845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
