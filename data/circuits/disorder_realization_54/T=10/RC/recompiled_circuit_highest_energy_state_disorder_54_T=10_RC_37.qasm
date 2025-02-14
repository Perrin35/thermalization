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
rz(-0.71159166) q[0];
sx q[0];
rz(1.1049668) q[0];
sx q[0];
rz(8.2850716) q[0];
rz(-0.79144129) q[1];
sx q[1];
rz(-1.0410407) q[1];
sx q[1];
rz(2.0120373) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4576042) q[0];
sx q[0];
rz(-0.86375551) q[0];
sx q[0];
rz(-0.23914214) q[0];
rz(-0.66313498) q[2];
sx q[2];
rz(-1.5929993) q[2];
sx q[2];
rz(0.93142366) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12295402) q[1];
sx q[1];
rz(-1.8644511) q[1];
sx q[1];
rz(-0.83033009) q[1];
x q[2];
rz(-3.0409052) q[3];
sx q[3];
rz(-2.8071981) q[3];
sx q[3];
rz(-1.3833801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1019885) q[2];
sx q[2];
rz(-0.81059376) q[2];
sx q[2];
rz(-0.90768901) q[2];
rz(-3.0229819) q[3];
sx q[3];
rz(-1.5396996) q[3];
sx q[3];
rz(-2.8906726) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1949961) q[0];
sx q[0];
rz(-2.4793766) q[0];
sx q[0];
rz(0.10435852) q[0];
rz(-1.1907499) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(0.45713919) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32923082) q[0];
sx q[0];
rz(-0.19596772) q[0];
sx q[0];
rz(1.7295444) q[0];
rz(-0.89294101) q[2];
sx q[2];
rz(-1.2051851) q[2];
sx q[2];
rz(-1.9302238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7434843) q[1];
sx q[1];
rz(-2.8116075) q[1];
sx q[1];
rz(-1.8673926) q[1];
x q[2];
rz(-0.13808226) q[3];
sx q[3];
rz(-0.90681091) q[3];
sx q[3];
rz(2.9241204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5856058) q[2];
sx q[2];
rz(-1.9827236) q[2];
sx q[2];
rz(-2.9846094) q[2];
rz(2.6202776) q[3];
sx q[3];
rz(-0.83955228) q[3];
sx q[3];
rz(-0.014001525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(-0.34573063) q[0];
rz(1.455447) q[1];
sx q[1];
rz(-1.2724178) q[1];
sx q[1];
rz(1.5706496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60590505) q[0];
sx q[0];
rz(-0.55258026) q[0];
sx q[0];
rz(-1.4009757) q[0];
rz(-pi) q[1];
rz(-1.5758031) q[2];
sx q[2];
rz(-0.62267674) q[2];
sx q[2];
rz(-1.1345991) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45528835) q[1];
sx q[1];
rz(-1.521811) q[1];
sx q[1];
rz(-0.21360417) q[1];
rz(0.9216347) q[3];
sx q[3];
rz(-1.7329669) q[3];
sx q[3];
rz(-2.0258486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7741144) q[2];
sx q[2];
rz(-0.228129) q[2];
sx q[2];
rz(2.1841614) q[2];
rz(-1.1227603) q[3];
sx q[3];
rz(-1.9072396) q[3];
sx q[3];
rz(1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6289309) q[0];
sx q[0];
rz(-1.4224195) q[0];
sx q[0];
rz(1.5976394) q[0];
rz(1.3900025) q[1];
sx q[1];
rz(-1.3163687) q[1];
sx q[1];
rz(1.66473) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2643971) q[0];
sx q[0];
rz(-1.0478643) q[0];
sx q[0];
rz(1.8967129) q[0];
rz(1.3324758) q[2];
sx q[2];
rz(-1.4891948) q[2];
sx q[2];
rz(-1.9072744) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.99914306) q[1];
sx q[1];
rz(-1.8987185) q[1];
sx q[1];
rz(-2.1254005) q[1];
rz(-pi) q[2];
rz(-0.32962004) q[3];
sx q[3];
rz(-1.8330036) q[3];
sx q[3];
rz(-0.58673687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59820286) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(-2.1447935) q[2];
rz(-1.8761084) q[3];
sx q[3];
rz(-0.71804738) q[3];
sx q[3];
rz(-2.8640981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0811512) q[0];
sx q[0];
rz(-1.569898) q[0];
sx q[0];
rz(2.0695709) q[0];
rz(-2.187166) q[1];
sx q[1];
rz(-1.1773033) q[1];
sx q[1];
rz(1.4929474) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7537751) q[0];
sx q[0];
rz(-2.3507581) q[0];
sx q[0];
rz(-2.8641939) q[0];
rz(-pi) q[1];
rz(-0.49333879) q[2];
sx q[2];
rz(-2.3559921) q[2];
sx q[2];
rz(1.0935022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49661885) q[1];
sx q[1];
rz(-2.7826834) q[1];
sx q[1];
rz(-0.41618698) q[1];
x q[2];
rz(2.8209854) q[3];
sx q[3];
rz(-1.4757475) q[3];
sx q[3];
rz(-2.7966201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2299049) q[2];
sx q[2];
rz(-2.7547084) q[2];
sx q[2];
rz(-1.1502728) q[2];
rz(-0.041042717) q[3];
sx q[3];
rz(-1.5951472) q[3];
sx q[3];
rz(-2.7400147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7683485) q[0];
sx q[0];
rz(-1.1078438) q[0];
sx q[0];
rz(1.9819697) q[0];
rz(-1.6805964) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.9083091) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1430968) q[0];
sx q[0];
rz(-0.84718207) q[0];
sx q[0];
rz(1.1529403) q[0];
x q[1];
rz(-0.057282863) q[2];
sx q[2];
rz(-1.5174688) q[2];
sx q[2];
rz(1.48207) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.199904) q[1];
sx q[1];
rz(-0.91477981) q[1];
sx q[1];
rz(2.0787048) q[1];
rz(-2.6150675) q[3];
sx q[3];
rz(-1.3428146) q[3];
sx q[3];
rz(2.990852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.636574) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(1.150307) q[2];
rz(-0.73927528) q[3];
sx q[3];
rz(-1.247568) q[3];
sx q[3];
rz(2.6967743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4742541) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(2.6928103) q[0];
rz(1.5397286) q[1];
sx q[1];
rz(-1.1069143) q[1];
sx q[1];
rz(1.9805699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5509508) q[0];
sx q[0];
rz(-1.30473) q[0];
sx q[0];
rz(-2.7177285) q[0];
x q[1];
rz(0.49893219) q[2];
sx q[2];
rz(-1.1321486) q[2];
sx q[2];
rz(2.2427151) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55884493) q[1];
sx q[1];
rz(-0.44932355) q[1];
sx q[1];
rz(-0.72968633) q[1];
rz(-1.1584985) q[3];
sx q[3];
rz(-2.5479043) q[3];
sx q[3];
rz(0.75324654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(-1.1807582) q[2];
rz(-2.0598038) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(-2.7150174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54843724) q[0];
sx q[0];
rz(-1.8226382) q[0];
sx q[0];
rz(-0.68761188) q[0];
rz(1.1697191) q[1];
sx q[1];
rz(-0.97427383) q[1];
sx q[1];
rz(-2.7489472) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8075645) q[0];
sx q[0];
rz(-2.0453718) q[0];
sx q[0];
rz(2.9827098) q[0];
rz(-0.31198172) q[2];
sx q[2];
rz(-1.7034966) q[2];
sx q[2];
rz(3.0614982) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1117797) q[1];
sx q[1];
rz(-1.4832442) q[1];
sx q[1];
rz(-2.665641) q[1];
rz(-pi) q[2];
rz(-2.9526934) q[3];
sx q[3];
rz(-1.4903063) q[3];
sx q[3];
rz(2.1565425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9795867) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(-1.6458192) q[2];
rz(1.6188072) q[3];
sx q[3];
rz(-0.77097547) q[3];
sx q[3];
rz(-2.5405367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9252121) q[0];
sx q[0];
rz(-0.38327152) q[0];
sx q[0];
rz(1.8076757) q[0];
rz(-1.1434309) q[1];
sx q[1];
rz(-0.83088487) q[1];
sx q[1];
rz(2.3480031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0045354768) q[0];
sx q[0];
rz(-1.6414149) q[0];
sx q[0];
rz(0.15839346) q[0];
x q[1];
rz(-1.6657532) q[2];
sx q[2];
rz(-2.5296202) q[2];
sx q[2];
rz(2.8245935) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5836583) q[1];
sx q[1];
rz(-1.0685295) q[1];
sx q[1];
rz(-1.3468379) q[1];
rz(-0.10754866) q[3];
sx q[3];
rz(-1.1226153) q[3];
sx q[3];
rz(-1.5654237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6952343) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(1.0515155) q[2];
rz(2.4255883) q[3];
sx q[3];
rz(-2.1456238) q[3];
sx q[3];
rz(2.9314465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42651549) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(-2.9606384) q[0];
rz(1.1894233) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(-0.98035556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5073779) q[0];
sx q[0];
rz(-1.6537154) q[0];
sx q[0];
rz(-2.6392127) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6858079) q[2];
sx q[2];
rz(-1.1768627) q[2];
sx q[2];
rz(-1.822871) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7804523) q[1];
sx q[1];
rz(-1.7516416) q[1];
sx q[1];
rz(-2.3922582) q[1];
x q[2];
rz(1.3602363) q[3];
sx q[3];
rz(-1.1835775) q[3];
sx q[3];
rz(0.59084597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1159726) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(1.8170961) q[2];
rz(-2.4888511) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(-3.1033707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7272335) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(0.7242135) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(0.090567055) q[2];
sx q[2];
rz(-1.9048077) q[2];
sx q[2];
rz(2.122369) q[2];
rz(1.4606838) q[3];
sx q[3];
rz(-2.3398945) q[3];
sx q[3];
rz(-1.0490976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
