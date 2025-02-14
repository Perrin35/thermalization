OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(0.77603618) q[0];
rz(0.70932055) q[1];
sx q[1];
rz(4.6588916) q[1];
sx q[1];
rz(8.8257596) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8255492) q[0];
sx q[0];
rz(-0.49015309) q[0];
sx q[0];
rz(1.3087628) q[0];
x q[1];
rz(-2.2602735) q[2];
sx q[2];
rz(-1.7909087) q[2];
sx q[2];
rz(0.51267363) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6769961) q[1];
sx q[1];
rz(-1.076816) q[1];
sx q[1];
rz(1.9736273) q[1];
x q[2];
rz(2.3612411) q[3];
sx q[3];
rz(-1.2557185) q[3];
sx q[3];
rz(-1.9596069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0029995) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(0.19621672) q[2];
rz(-2.0972706) q[3];
sx q[3];
rz(-0.82572562) q[3];
sx q[3];
rz(2.830982) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0394548) q[0];
sx q[0];
rz(-0.65334833) q[0];
sx q[0];
rz(0.85773221) q[0];
rz(0.54025447) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(-0.78261715) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40561134) q[0];
sx q[0];
rz(-2.1850039) q[0];
sx q[0];
rz(-3.0354397) q[0];
x q[1];
rz(1.247632) q[2];
sx q[2];
rz(-2.4944381) q[2];
sx q[2];
rz(-2.0314856) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23797135) q[1];
sx q[1];
rz(-0.97717932) q[1];
sx q[1];
rz(1.3376544) q[1];
rz(-pi) q[2];
rz(0.55763839) q[3];
sx q[3];
rz(-2.1774315) q[3];
sx q[3];
rz(2.5733657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1648078) q[2];
sx q[2];
rz(-2.3471577) q[2];
sx q[2];
rz(0.59107333) q[2];
rz(2.1137386) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22981055) q[0];
sx q[0];
rz(-2.6214143) q[0];
sx q[0];
rz(1.0562563) q[0];
rz(-0.99110574) q[1];
sx q[1];
rz(-1.7678363) q[1];
sx q[1];
rz(-0.80449218) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4933518) q[0];
sx q[0];
rz(-1.3051148) q[0];
sx q[0];
rz(-0.65346395) q[0];
rz(-pi) q[1];
rz(-0.095011906) q[2];
sx q[2];
rz(-2.1974051) q[2];
sx q[2];
rz(-1.0209521) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3884405) q[1];
sx q[1];
rz(-2.8049251) q[1];
sx q[1];
rz(2.3364288) q[1];
rz(-0.25193416) q[3];
sx q[3];
rz(-2.478699) q[3];
sx q[3];
rz(0.54352647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.450401) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(-2.5737393) q[2];
rz(-1.9495226) q[3];
sx q[3];
rz(-2.6046533) q[3];
sx q[3];
rz(-2.9862064) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464226) q[0];
sx q[0];
rz(-1.3986873) q[0];
sx q[0];
rz(-1.1871185) q[0];
rz(-2.8861956) q[1];
sx q[1];
rz(-1.7385769) q[1];
sx q[1];
rz(-2.752221) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5449937) q[0];
sx q[0];
rz(-1.7905718) q[0];
sx q[0];
rz(1.4334428) q[0];
rz(-1.537048) q[2];
sx q[2];
rz(-2.4141888) q[2];
sx q[2];
rz(-2.6836065) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2823209) q[1];
sx q[1];
rz(-2.2208557) q[1];
sx q[1];
rz(2.9031624) q[1];
rz(-pi) q[2];
rz(-1.3741174) q[3];
sx q[3];
rz(-1.4405319) q[3];
sx q[3];
rz(-1.2127339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7369467) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(-1.7741989) q[2];
rz(0.5395475) q[3];
sx q[3];
rz(-1.5522141) q[3];
sx q[3];
rz(1.8168943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0668199) q[0];
sx q[0];
rz(-0.19817752) q[0];
sx q[0];
rz(-0.5603801) q[0];
rz(-0.21891521) q[1];
sx q[1];
rz(-2.0603265) q[1];
sx q[1];
rz(-2.970649) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6534916) q[0];
sx q[0];
rz(-0.91918901) q[0];
sx q[0];
rz(-0.96820684) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4291228) q[2];
sx q[2];
rz(-1.7729086) q[2];
sx q[2];
rz(-2.0703482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39300181) q[1];
sx q[1];
rz(-2.1868863) q[1];
sx q[1];
rz(-2.2115117) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6690977) q[3];
sx q[3];
rz(-2.519033) q[3];
sx q[3];
rz(-3.0873201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35974744) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(-2.183059) q[2];
rz(-2.9790699) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(-1.5490279) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.573134) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(-0.74813133) q[0];
rz(1.1185147) q[1];
sx q[1];
rz(-0.41901127) q[1];
sx q[1];
rz(-2.8245139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1823163) q[0];
sx q[0];
rz(-1.3249363) q[0];
sx q[0];
rz(1.5005174) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9002764) q[2];
sx q[2];
rz(-2.9394627) q[2];
sx q[2];
rz(-1.8692819) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78817716) q[1];
sx q[1];
rz(-1.6420351) q[1];
sx q[1];
rz(2.337238) q[1];
rz(-3.1012332) q[3];
sx q[3];
rz(-2.3642614) q[3];
sx q[3];
rz(-0.07829994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3743484) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(-1.413215) q[2];
rz(3.0610541) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(-2.9197689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0291075) q[0];
sx q[0];
rz(-1.8928098) q[0];
sx q[0];
rz(1.2741733) q[0];
rz(1.796465) q[1];
sx q[1];
rz(-2.3990264) q[1];
sx q[1];
rz(1.8003731) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91914908) q[0];
sx q[0];
rz(-2.0608927) q[0];
sx q[0];
rz(-0.63700139) q[0];
rz(2.6599081) q[2];
sx q[2];
rz(-1.8284214) q[2];
sx q[2];
rz(-0.017901808) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0643031) q[1];
sx q[1];
rz(-1.3796564) q[1];
sx q[1];
rz(0.15644507) q[1];
rz(-pi) q[2];
rz(-3.0242434) q[3];
sx q[3];
rz(-0.92399358) q[3];
sx q[3];
rz(2.5600159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.93366569) q[2];
sx q[2];
rz(-1.1873446) q[2];
sx q[2];
rz(0.46736091) q[2];
rz(-0.76006877) q[3];
sx q[3];
rz(-1.6642539) q[3];
sx q[3];
rz(-1.8925586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6811328) q[0];
sx q[0];
rz(-2.12119) q[0];
sx q[0];
rz(1.4469752) q[0];
rz(2.815822) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(1.6206585) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50624079) q[0];
sx q[0];
rz(-1.304848) q[0];
sx q[0];
rz(-2.6145934) q[0];
rz(-pi) q[1];
x q[1];
rz(0.09163945) q[2];
sx q[2];
rz(-1.5711725) q[2];
sx q[2];
rz(2.3918652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.28315) q[1];
sx q[1];
rz(-1.2179188) q[1];
sx q[1];
rz(-2.9550885) q[1];
x q[2];
rz(-2.9321947) q[3];
sx q[3];
rz(-0.87762524) q[3];
sx q[3];
rz(1.5537167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96523607) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(2.526324) q[2];
rz(1.2728914) q[3];
sx q[3];
rz(-2.8764909) q[3];
sx q[3];
rz(-2.7191775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0203005) q[0];
sx q[0];
rz(-1.8676119) q[0];
sx q[0];
rz(2.2913388) q[0];
rz(2.0203159) q[1];
sx q[1];
rz(-1.774615) q[1];
sx q[1];
rz(-0.77176315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5284755) q[0];
sx q[0];
rz(-0.69418797) q[0];
sx q[0];
rz(-1.608169) q[0];
rz(-pi) q[1];
rz(-2.775101) q[2];
sx q[2];
rz(-2.9482609) q[2];
sx q[2];
rz(1.1761348) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38852019) q[1];
sx q[1];
rz(-2.8917312) q[1];
sx q[1];
rz(1.7953963) q[1];
x q[2];
rz(1.2466627) q[3];
sx q[3];
rz(-1.3132902) q[3];
sx q[3];
rz(-0.83453629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0261592) q[2];
sx q[2];
rz(-1.5966281) q[2];
sx q[2];
rz(-2.3010632) q[2];
rz(0.080951512) q[3];
sx q[3];
rz(-0.5778802) q[3];
sx q[3];
rz(2.8988083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42846546) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(-1.1567098) q[0];
rz(-2.1139862) q[1];
sx q[1];
rz(-1.7637858) q[1];
sx q[1];
rz(0.81046945) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5549703) q[0];
sx q[0];
rz(-1.3548618) q[0];
sx q[0];
rz(-0.76540375) q[0];
rz(0.54525156) q[2];
sx q[2];
rz(-1.8601189) q[2];
sx q[2];
rz(-2.3985753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6582074) q[1];
sx q[1];
rz(-1.864363) q[1];
sx q[1];
rz(-2.844643) q[1];
rz(-pi) q[2];
rz(-0.23328554) q[3];
sx q[3];
rz(-2.5676227) q[3];
sx q[3];
rz(-0.026653224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.345574) q[2];
sx q[2];
rz(-1.8859665) q[2];
sx q[2];
rz(1.5239117) q[2];
rz(-1.1003305) q[3];
sx q[3];
rz(-1.1805781) q[3];
sx q[3];
rz(-2.9617917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1179467) q[0];
sx q[0];
rz(-1.4143586) q[0];
sx q[0];
rz(2.216862) q[0];
rz(-3.0991411) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(-1.4412075) q[2];
sx q[2];
rz(-2.7322506) q[2];
sx q[2];
rz(2.6553497) q[2];
rz(0.073033606) q[3];
sx q[3];
rz(-2.6831153) q[3];
sx q[3];
rz(-1.4014184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
