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
rz(1.1178782) q[0];
sx q[0];
rz(4.2269308) q[0];
sx q[0];
rz(9.7836054) q[0];
rz(-1.8817512) q[1];
sx q[1];
rz(-2.0460195) q[1];
sx q[1];
rz(-2.1020558) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7494237) q[0];
sx q[0];
rz(-0.75558973) q[0];
sx q[0];
rz(-2.6120793) q[0];
rz(-pi) q[1];
rz(-1.5284003) q[2];
sx q[2];
rz(-1.543352) q[2];
sx q[2];
rz(1.3224441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1778587) q[1];
sx q[1];
rz(-1.8699282) q[1];
sx q[1];
rz(-1.9480716) q[1];
rz(-pi) q[2];
rz(-1.667439) q[3];
sx q[3];
rz(-0.90296871) q[3];
sx q[3];
rz(-0.63233318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.871668) q[2];
sx q[2];
rz(-0.85447997) q[2];
sx q[2];
rz(-2.975115) q[2];
rz(3.0050333) q[3];
sx q[3];
rz(-2.2763493) q[3];
sx q[3];
rz(-1.6816007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.0110375) q[0];
sx q[0];
rz(-2.8189973) q[0];
sx q[0];
rz(0.36001298) q[0];
rz(0.75045466) q[1];
sx q[1];
rz(-2.2513335) q[1];
sx q[1];
rz(-0.043936122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.982807) q[0];
sx q[0];
rz(-1.7097397) q[0];
sx q[0];
rz(-1.7278633) q[0];
rz(-1.5055799) q[2];
sx q[2];
rz(-2.4963317) q[2];
sx q[2];
rz(1.1087904) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3219731) q[1];
sx q[1];
rz(-0.83998806) q[1];
sx q[1];
rz(-1.2801484) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93505164) q[3];
sx q[3];
rz(-2.8160444) q[3];
sx q[3];
rz(2.9071484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7875849) q[2];
sx q[2];
rz(-0.44469357) q[2];
sx q[2];
rz(1.0178052) q[2];
rz(2.9338037) q[3];
sx q[3];
rz(-2.5354247) q[3];
sx q[3];
rz(-0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(1.4491693) q[0];
sx q[0];
rz(-2.848931) q[0];
sx q[0];
rz(-1.0798651) q[0];
rz(2.1366468) q[1];
sx q[1];
rz(-2.4506073) q[1];
sx q[1];
rz(-1.2452004) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024187758) q[0];
sx q[0];
rz(-1.4862446) q[0];
sx q[0];
rz(-0.3106433) q[0];
rz(-0.4246906) q[2];
sx q[2];
rz(-1.7199368) q[2];
sx q[2];
rz(-1.5737267) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.019466467) q[1];
sx q[1];
rz(-1.2858741) q[1];
sx q[1];
rz(2.1287005) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1916189) q[3];
sx q[3];
rz(-1.432918) q[3];
sx q[3];
rz(2.8153489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9006742) q[2];
sx q[2];
rz(-1.2718879) q[2];
sx q[2];
rz(-2.7980878) q[2];
rz(3.1268696) q[3];
sx q[3];
rz(-0.69535178) q[3];
sx q[3];
rz(-1.2408313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13046509) q[0];
sx q[0];
rz(-1.0437597) q[0];
sx q[0];
rz(-0.68487942) q[0];
rz(0.20826805) q[1];
sx q[1];
rz(-1.8332278) q[1];
sx q[1];
rz(-2.3484255) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.255729) q[0];
sx q[0];
rz(-1.0287713) q[0];
sx q[0];
rz(-3.0266808) q[0];
rz(2.9689413) q[2];
sx q[2];
rz(-1.358629) q[2];
sx q[2];
rz(-1.6546951) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85911688) q[1];
sx q[1];
rz(-2.1819677) q[1];
sx q[1];
rz(3.1219367) q[1];
rz(-pi) q[2];
rz(-1.2007522) q[3];
sx q[3];
rz(-1.368282) q[3];
sx q[3];
rz(-2.5807192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76306474) q[2];
sx q[2];
rz(-1.6918162) q[2];
sx q[2];
rz(-2.6587963) q[2];
rz(1.5040846) q[3];
sx q[3];
rz(-2.5129694) q[3];
sx q[3];
rz(-2.8211236) q[3];
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
rz(-pi/2) q[0];
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
rz(0.046431635) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(-3.0401373) q[0];
rz(0.21508148) q[1];
sx q[1];
rz(-1.2098034) q[1];
sx q[1];
rz(1.2717517) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1908784) q[0];
sx q[0];
rz(-2.029065) q[0];
sx q[0];
rz(0.31744297) q[0];
rz(-pi) q[1];
rz(0.59753363) q[2];
sx q[2];
rz(-1.1206822) q[2];
sx q[2];
rz(-2.2240153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7600619) q[1];
sx q[1];
rz(-2.7193128) q[1];
sx q[1];
rz(-1.1997919) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9157294) q[3];
sx q[3];
rz(-1.3702379) q[3];
sx q[3];
rz(-0.95074703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7868906) q[2];
sx q[2];
rz(-0.44822794) q[2];
sx q[2];
rz(-0.1498214) q[2];
rz(-0.24770728) q[3];
sx q[3];
rz(-2.7713573) q[3];
sx q[3];
rz(0.80055922) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.596452) q[0];
sx q[0];
rz(-2.5189724) q[0];
sx q[0];
rz(-1.2860292) q[0];
rz(-0.81895685) q[1];
sx q[1];
rz(-0.30636925) q[1];
sx q[1];
rz(-0.18316306) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64252999) q[0];
sx q[0];
rz(-1.0181435) q[0];
sx q[0];
rz(-1.8091317) q[0];
rz(-pi) q[1];
rz(-1.6346857) q[2];
sx q[2];
rz(-1.3861292) q[2];
sx q[2];
rz(-1.9353362) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0350041) q[1];
sx q[1];
rz(-2.2640806) q[1];
sx q[1];
rz(0.50030746) q[1];
x q[2];
rz(2.5991393) q[3];
sx q[3];
rz(-0.56672719) q[3];
sx q[3];
rz(-2.7572214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5554819) q[2];
sx q[2];
rz(-1.4376983) q[2];
sx q[2];
rz(2.6044031) q[2];
rz(-1.4026583) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(-2.5130443) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1527767) q[0];
sx q[0];
rz(-2.7429136) q[0];
sx q[0];
rz(-0.11845778) q[0];
rz(-1.1064103) q[1];
sx q[1];
rz(-1.9170008) q[1];
sx q[1];
rz(-2.4647663) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9119916) q[0];
sx q[0];
rz(-1.6458479) q[0];
sx q[0];
rz(-1.5154535) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65420349) q[2];
sx q[2];
rz(-1.6779644) q[2];
sx q[2];
rz(0.079710641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21444328) q[1];
sx q[1];
rz(-2.8942558) q[1];
sx q[1];
rz(2.3995968) q[1];
rz(1.5854225) q[3];
sx q[3];
rz(-1.5717603) q[3];
sx q[3];
rz(0.87982501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38830882) q[2];
sx q[2];
rz(-0.94677418) q[2];
sx q[2];
rz(-0.050966144) q[2];
rz(-2.986159) q[3];
sx q[3];
rz(-1.4916462) q[3];
sx q[3];
rz(1.0203993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8841356) q[0];
sx q[0];
rz(-1.1549042) q[0];
sx q[0];
rz(-2.1575523) q[0];
rz(2.6718196) q[1];
sx q[1];
rz(-1.4642986) q[1];
sx q[1];
rz(0.22107548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8741989) q[0];
sx q[0];
rz(-1.4765863) q[0];
sx q[0];
rz(-1.5896958) q[0];
x q[1];
rz(-1.1936155) q[2];
sx q[2];
rz(-2.1862891) q[2];
sx q[2];
rz(0.35843682) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4692336) q[1];
sx q[1];
rz(-1.2738528) q[1];
sx q[1];
rz(0.23632144) q[1];
rz(2.2443767) q[3];
sx q[3];
rz(-2.4468337) q[3];
sx q[3];
rz(-1.0396797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2880963) q[2];
sx q[2];
rz(-0.6820389) q[2];
sx q[2];
rz(-1.3300396) q[2];
rz(2.5174777) q[3];
sx q[3];
rz(-2.1835486) q[3];
sx q[3];
rz(-2.7522855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5617111) q[0];
sx q[0];
rz(-1.4917253) q[0];
sx q[0];
rz(-2.0910182) q[0];
rz(-2.1613278) q[1];
sx q[1];
rz(-1.0824208) q[1];
sx q[1];
rz(-1.6260653) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1001871) q[0];
sx q[0];
rz(-0.55313191) q[0];
sx q[0];
rz(-2.83908) q[0];
rz(0.43584906) q[2];
sx q[2];
rz(-2.6245891) q[2];
sx q[2];
rz(2.2877501) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5055464) q[1];
sx q[1];
rz(-2.6719173) q[1];
sx q[1];
rz(-2.4081025) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.100176) q[3];
sx q[3];
rz(-1.017316) q[3];
sx q[3];
rz(-3.0867004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3773697) q[2];
sx q[2];
rz(-1.2219967) q[2];
sx q[2];
rz(-0.73966217) q[2];
rz(2.9861279) q[3];
sx q[3];
rz(-2.7630617) q[3];
sx q[3];
rz(0.19708656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4231606) q[0];
sx q[0];
rz(-2.6078556) q[0];
sx q[0];
rz(1.1370132) q[0];
rz(-2.5088189) q[1];
sx q[1];
rz(-0.67291617) q[1];
sx q[1];
rz(-0.64579642) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.150313) q[0];
sx q[0];
rz(-1.8984573) q[0];
sx q[0];
rz(1.3929266) q[0];
rz(0.77984758) q[2];
sx q[2];
rz(-1.7371685) q[2];
sx q[2];
rz(0.8188687) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.551522) q[1];
sx q[1];
rz(-2.3209718) q[1];
sx q[1];
rz(3.0469826) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0539551) q[3];
sx q[3];
rz(-2.5541383) q[3];
sx q[3];
rz(-1.7980391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48468963) q[2];
sx q[2];
rz(-2.778896) q[2];
sx q[2];
rz(-0.038385782) q[2];
rz(3.0302327) q[3];
sx q[3];
rz(-0.42698082) q[3];
sx q[3];
rz(-0.67489433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69042324) q[0];
sx q[0];
rz(-1.7043865) q[0];
sx q[0];
rz(-1.7519328) q[0];
rz(1.2870862) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(2.3310326) q[2];
sx q[2];
rz(-2.3702224) q[2];
sx q[2];
rz(0.67077256) q[2];
rz(1.7467563) q[3];
sx q[3];
rz(-0.8192438) q[3];
sx q[3];
rz(1.4024709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
