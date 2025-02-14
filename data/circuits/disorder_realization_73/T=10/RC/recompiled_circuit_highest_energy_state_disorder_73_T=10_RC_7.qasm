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
rz(2.4740648) q[0];
sx q[0];
rz(-0.83905554) q[0];
sx q[0];
rz(0.43247217) q[0];
rz(0.67983812) q[1];
sx q[1];
rz(-1.4601269) q[1];
sx q[1];
rz(2.3753994) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6941152) q[0];
sx q[0];
rz(-1.2060529) q[0];
sx q[0];
rz(-3.0459131) q[0];
rz(-1.7518429) q[2];
sx q[2];
rz(-1.6012601) q[2];
sx q[2];
rz(1.9150645) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39734617) q[1];
sx q[1];
rz(-1.0434217) q[1];
sx q[1];
rz(2.9832193) q[1];
rz(-pi) q[2];
rz(-3.0187383) q[3];
sx q[3];
rz(-1.4335504) q[3];
sx q[3];
rz(-0.48396971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5863775) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(-0.66184723) q[2];
rz(-1.2343179) q[3];
sx q[3];
rz(-1.581278) q[3];
sx q[3];
rz(-0.1525277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9041651) q[0];
sx q[0];
rz(-2.4894042) q[0];
sx q[0];
rz(2.7320614) q[0];
rz(2.883059) q[1];
sx q[1];
rz(-1.9347609) q[1];
sx q[1];
rz(2.5915367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6790228) q[0];
sx q[0];
rz(-1.4840811) q[0];
sx q[0];
rz(1.0957415) q[0];
rz(-1.0070877) q[2];
sx q[2];
rz(-1.1620143) q[2];
sx q[2];
rz(-0.43540149) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.066795863) q[1];
sx q[1];
rz(-0.48802146) q[1];
sx q[1];
rz(2.0171432) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82050555) q[3];
sx q[3];
rz(-1.3119999) q[3];
sx q[3];
rz(0.25322278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6253769) q[2];
sx q[2];
rz(-1.7552525) q[2];
sx q[2];
rz(-0.59164444) q[2];
rz(2.6243465) q[3];
sx q[3];
rz(-1.1040684) q[3];
sx q[3];
rz(-2.5669317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2194694) q[0];
sx q[0];
rz(-1.8485494) q[0];
sx q[0];
rz(-1.0009694) q[0];
rz(-1.7132022) q[1];
sx q[1];
rz(-0.75894558) q[1];
sx q[1];
rz(3.0975814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3428033) q[0];
sx q[0];
rz(-0.90315032) q[0];
sx q[0];
rz(-2.0602047) q[0];
x q[1];
rz(0.050878164) q[2];
sx q[2];
rz(-0.83453611) q[2];
sx q[2];
rz(0.60896192) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2340369) q[1];
sx q[1];
rz(-2.2436972) q[1];
sx q[1];
rz(2.2220045) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1009665) q[3];
sx q[3];
rz(-1.3990132) q[3];
sx q[3];
rz(-2.9295827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22237805) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(0.10981336) q[2];
rz(2.5114656) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(2.8540247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6269864) q[0];
sx q[0];
rz(-2.0872748) q[0];
sx q[0];
rz(2.0475673) q[0];
rz(-0.47359723) q[1];
sx q[1];
rz(-2.4351599) q[1];
sx q[1];
rz(-0.72023645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0648616) q[0];
sx q[0];
rz(-2.0429049) q[0];
sx q[0];
rz(-1.9143299) q[0];
rz(-pi) q[1];
rz(-1.6178151) q[2];
sx q[2];
rz(-1.8787307) q[2];
sx q[2];
rz(1.5297254) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1853237) q[1];
sx q[1];
rz(-2.2898623) q[1];
sx q[1];
rz(1.1683502) q[1];
rz(-pi) q[2];
rz(-1.9071753) q[3];
sx q[3];
rz(-2.3569538) q[3];
sx q[3];
rz(1.4799081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56095162) q[2];
sx q[2];
rz(-1.6683942) q[2];
sx q[2];
rz(-0.15216039) q[2];
rz(-1.8469319) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(-1.7744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691567) q[0];
sx q[0];
rz(-0.56273383) q[0];
sx q[0];
rz(-0.38134545) q[0];
rz(2.5533679) q[1];
sx q[1];
rz(-1.4071848) q[1];
sx q[1];
rz(3.0773967) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682024) q[0];
sx q[0];
rz(-1.9406835) q[0];
sx q[0];
rz(-1.3805519) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1663127) q[2];
sx q[2];
rz(-2.1475386) q[2];
sx q[2];
rz(-2.8167644) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.209268) q[1];
sx q[1];
rz(-1.2863683) q[1];
sx q[1];
rz(0.80764215) q[1];
x q[2];
rz(-1.3649356) q[3];
sx q[3];
rz(-1.4511709) q[3];
sx q[3];
rz(-2.4158166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.75307816) q[2];
sx q[2];
rz(-1.3103176) q[2];
sx q[2];
rz(-2.348796) q[2];
rz(0.13475284) q[3];
sx q[3];
rz(-0.56551352) q[3];
sx q[3];
rz(-1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(0.95293522) q[0];
rz(1.8903525) q[1];
sx q[1];
rz(-1.4616936) q[1];
sx q[1];
rz(-0.23426637) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.983646) q[0];
sx q[0];
rz(-1.0704213) q[0];
sx q[0];
rz(1.7364794) q[0];
rz(2.6855747) q[2];
sx q[2];
rz(-1.2022744) q[2];
sx q[2];
rz(0.8816661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3418048) q[1];
sx q[1];
rz(-1.4569286) q[1];
sx q[1];
rz(2.1569685) q[1];
x q[2];
rz(1.4552507) q[3];
sx q[3];
rz(-2.503814) q[3];
sx q[3];
rz(-1.6751643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17144063) q[2];
sx q[2];
rz(-0.60005391) q[2];
sx q[2];
rz(-0.94685143) q[2];
rz(1.6998477) q[3];
sx q[3];
rz(-1.3274679) q[3];
sx q[3];
rz(3.0202151) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22987188) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(0.3748931) q[0];
rz(-2.3986469) q[1];
sx q[1];
rz(-0.82287794) q[1];
sx q[1];
rz(-0.65623647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0699871) q[0];
sx q[0];
rz(-0.90815571) q[0];
sx q[0];
rz(1.9020984) q[0];
rz(-2.1383275) q[2];
sx q[2];
rz(-2.1250688) q[2];
sx q[2];
rz(0.90855937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9083169) q[1];
sx q[1];
rz(-1.9413949) q[1];
sx q[1];
rz(-2.6909687) q[1];
rz(-pi) q[2];
rz(2.8493622) q[3];
sx q[3];
rz(-2.5295895) q[3];
sx q[3];
rz(-2.0296869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52834213) q[2];
sx q[2];
rz(-2.2684147) q[2];
sx q[2];
rz(-0.64995107) q[2];
rz(-2.2925099) q[3];
sx q[3];
rz(-0.63747469) q[3];
sx q[3];
rz(1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7717302) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(0.318203) q[0];
rz(-2.8089306) q[1];
sx q[1];
rz(-2.5883364) q[1];
sx q[1];
rz(-2.7645848) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72412598) q[0];
sx q[0];
rz(-1.2340349) q[0];
sx q[0];
rz(0.048857852) q[0];
x q[1];
rz(1.6416574) q[2];
sx q[2];
rz(-0.46816816) q[2];
sx q[2];
rz(-0.62781321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.01025) q[1];
sx q[1];
rz(-2.2325062) q[1];
sx q[1];
rz(2.0493839) q[1];
rz(-pi) q[2];
rz(2.0539862) q[3];
sx q[3];
rz(-1.6738664) q[3];
sx q[3];
rz(-0.3576935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.13410021) q[2];
sx q[2];
rz(-1.8635805) q[2];
sx q[2];
rz(2.210614) q[2];
rz(1.5382986) q[3];
sx q[3];
rz(-1.337991) q[3];
sx q[3];
rz(1.3163756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304921) q[0];
sx q[0];
rz(-2.7340524) q[0];
sx q[0];
rz(-2.7000309) q[0];
rz(2.8764797) q[1];
sx q[1];
rz(-1.6997063) q[1];
sx q[1];
rz(1.5107059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2870474) q[0];
sx q[0];
rz(-1.1397188) q[0];
sx q[0];
rz(-1.8148941) q[0];
x q[1];
rz(-2.8042107) q[2];
sx q[2];
rz(-1.0058912) q[2];
sx q[2];
rz(-0.039410465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35263044) q[1];
sx q[1];
rz(-1.4711385) q[1];
sx q[1];
rz(-1.8046735) q[1];
rz(-pi) q[2];
rz(-2.2208414) q[3];
sx q[3];
rz(-1.8691321) q[3];
sx q[3];
rz(0.32978299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6049217) q[2];
sx q[2];
rz(-1.2769545) q[2];
sx q[2];
rz(-0.46860179) q[2];
rz(2.129715) q[3];
sx q[3];
rz(-1.7185271) q[3];
sx q[3];
rz(-1.0812149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45379797) q[0];
sx q[0];
rz(-0.84306222) q[0];
sx q[0];
rz(-1.4367562) q[0];
rz(0.94802481) q[1];
sx q[1];
rz(-1.3399622) q[1];
sx q[1];
rz(0.2737793) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9655436) q[0];
sx q[0];
rz(-1.4469441) q[0];
sx q[0];
rz(-1.6424422) q[0];
x q[1];
rz(-0.003652775) q[2];
sx q[2];
rz(-2.4953702) q[2];
sx q[2];
rz(-2.5346859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.35008365) q[1];
sx q[1];
rz(-0.71301642) q[1];
sx q[1];
rz(-0.75292315) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0300566) q[3];
sx q[3];
rz(-0.74990898) q[3];
sx q[3];
rz(0.53242079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56129366) q[2];
sx q[2];
rz(-1.0729921) q[2];
sx q[2];
rz(2.6936074) q[2];
rz(0.63730803) q[3];
sx q[3];
rz(-2.1154805) q[3];
sx q[3];
rz(0.89504939) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7359903) q[0];
sx q[0];
rz(-1.525282) q[0];
sx q[0];
rz(-1.4646144) q[0];
rz(-0.51392705) q[1];
sx q[1];
rz(-2.3873139) q[1];
sx q[1];
rz(0.25968459) q[1];
rz(2.6988637) q[2];
sx q[2];
rz(-1.6479658) q[2];
sx q[2];
rz(-2.992291) q[2];
rz(-3.1410029) q[3];
sx q[3];
rz(-1.1564217) q[3];
sx q[3];
rz(0.52903928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
