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
rz(2.3025371) q[0];
sx q[0];
rz(8.9923058) q[0];
rz(-2.4617545) q[1];
sx q[1];
rz(-1.6814657) q[1];
sx q[1];
rz(0.7661933) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1846518) q[0];
sx q[0];
rz(-0.37654018) q[0];
sx q[0];
rz(1.3256289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030969663) q[2];
sx q[2];
rz(-1.7517579) q[2];
sx q[2];
rz(-0.33869264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39734617) q[1];
sx q[1];
rz(-2.098171) q[1];
sx q[1];
rz(2.9832193) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0187383) q[3];
sx q[3];
rz(-1.7080423) q[3];
sx q[3];
rz(-2.6576229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55521518) q[2];
sx q[2];
rz(-2.5368498) q[2];
sx q[2];
rz(-0.66184723) q[2];
rz(1.2343179) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(2.9890649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2374275) q[0];
sx q[0];
rz(-0.65218848) q[0];
sx q[0];
rz(0.40953127) q[0];
rz(-0.25853363) q[1];
sx q[1];
rz(-1.9347609) q[1];
sx q[1];
rz(-0.55005598) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15274238) q[0];
sx q[0];
rz(-2.0439195) q[0];
sx q[0];
rz(-3.0441441) q[0];
rz(-pi) q[1];
rz(-2.2520355) q[2];
sx q[2];
rz(-2.4585138) q[2];
sx q[2];
rz(1.6967333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9039144) q[1];
sx q[1];
rz(-1.774607) q[1];
sx q[1];
rz(1.1242178) q[1];
x q[2];
rz(2.7943479) q[3];
sx q[3];
rz(-0.85111516) q[3];
sx q[3];
rz(1.589845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5162158) q[2];
sx q[2];
rz(-1.7552525) q[2];
sx q[2];
rz(-0.59164444) q[2];
rz(-0.51724616) q[3];
sx q[3];
rz(-2.0375242) q[3];
sx q[3];
rz(-0.57466093) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2194694) q[0];
sx q[0];
rz(-1.8485494) q[0];
sx q[0];
rz(1.0009694) q[0];
rz(-1.4283904) q[1];
sx q[1];
rz(-0.75894558) q[1];
sx q[1];
rz(-3.0975814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.509217) q[0];
sx q[0];
rz(-2.3365743) q[0];
sx q[0];
rz(2.6039327) q[0];
rz(2.3077009) q[2];
sx q[2];
rz(-1.6084889) q[2];
sx q[2];
rz(2.1455763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3487425) q[1];
sx q[1];
rz(-0.89952911) q[1];
sx q[1];
rz(0.6502186) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3408362) q[3];
sx q[3];
rz(-2.965117) q[3];
sx q[3];
rz(-3.120144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9192146) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(-0.10981336) q[2];
rz(0.63012704) q[3];
sx q[3];
rz(-2.5870242) q[3];
sx q[3];
rz(-0.28756791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.0940254) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(2.4213562) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41021529) q[0];
sx q[0];
rz(-0.57613724) q[0];
sx q[0];
rz(2.558484) q[0];
rz(-pi) q[1];
rz(-0.30825395) q[2];
sx q[2];
rz(-1.6156019) q[2];
sx q[2];
rz(-3.0862612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1853237) q[1];
sx q[1];
rz(-2.2898623) q[1];
sx q[1];
rz(-1.1683502) q[1];
rz(-0.31836001) q[3];
sx q[3];
rz(-2.3009319) q[3];
sx q[3];
rz(2.1206253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56095162) q[2];
sx q[2];
rz(-1.6683942) q[2];
sx q[2];
rz(-0.15216039) q[2];
rz(-1.2946607) q[3];
sx q[3];
rz(-1.4363689) q[3];
sx q[3];
rz(1.367182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072435943) q[0];
sx q[0];
rz(-2.5788588) q[0];
sx q[0];
rz(0.38134545) q[0];
rz(2.5533679) q[1];
sx q[1];
rz(-1.7344079) q[1];
sx q[1];
rz(0.064195976) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31606693) q[0];
sx q[0];
rz(-2.7276606) q[0];
sx q[0];
rz(-0.45376208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5975294) q[2];
sx q[2];
rz(-0.69101101) q[2];
sx q[2];
rz(2.8008583) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.209268) q[1];
sx q[1];
rz(-1.2863683) q[1];
sx q[1];
rz(-0.80764215) q[1];
rz(-pi) q[2];
rz(-2.1023687) q[3];
sx q[3];
rz(-2.9039249) q[3];
sx q[3];
rz(-0.32581926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75307816) q[2];
sx q[2];
rz(-1.831275) q[2];
sx q[2];
rz(-2.348796) q[2];
rz(-3.0068398) q[3];
sx q[3];
rz(-0.56551352) q[3];
sx q[3];
rz(1.6511065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7029) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(0.95293522) q[0];
rz(-1.8903525) q[1];
sx q[1];
rz(-1.6798991) q[1];
sx q[1];
rz(2.9073263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.983646) q[0];
sx q[0];
rz(-1.0704213) q[0];
sx q[0];
rz(1.7364794) q[0];
rz(-2.6855747) q[2];
sx q[2];
rz(-1.2022744) q[2];
sx q[2];
rz(-0.8816661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.059522337) q[1];
sx q[1];
rz(-0.59585427) q[1];
sx q[1];
rz(1.3669307) q[1];
x q[2];
rz(-1.686342) q[3];
sx q[3];
rz(-2.503814) q[3];
sx q[3];
rz(-1.6751643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17144063) q[2];
sx q[2];
rz(-0.60005391) q[2];
sx q[2];
rz(0.94685143) q[2];
rz(-1.4417449) q[3];
sx q[3];
rz(-1.3274679) q[3];
sx q[3];
rz(3.0202151) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9117208) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(0.3748931) q[0];
rz(-0.74294576) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(-0.65623647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0716055) q[0];
sx q[0];
rz(-2.2334369) q[0];
sx q[0];
rz(1.2394942) q[0];
rz(-pi) q[1];
rz(2.1383275) q[2];
sx q[2];
rz(-1.0165239) q[2];
sx q[2];
rz(-2.2330333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69493964) q[1];
sx q[1];
rz(-2.5663554) q[1];
sx q[1];
rz(0.72845643) q[1];
x q[2];
rz(2.5498059) q[3];
sx q[3];
rz(-1.4045219) q[3];
sx q[3];
rz(-2.4412654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52834213) q[2];
sx q[2];
rz(-0.87317792) q[2];
sx q[2];
rz(0.64995107) q[2];
rz(-0.84908271) q[3];
sx q[3];
rz(-2.504118) q[3];
sx q[3];
rz(1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.36986247) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(2.8233897) q[0];
rz(-0.33266208) q[1];
sx q[1];
rz(-2.5883364) q[1];
sx q[1];
rz(2.7645848) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72412598) q[0];
sx q[0];
rz(-1.9075577) q[0];
sx q[0];
rz(-0.048857852) q[0];
rz(1.4999352) q[2];
sx q[2];
rz(-2.6734245) q[2];
sx q[2];
rz(2.5137794) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.711447) q[1];
sx q[1];
rz(-0.79497577) q[1];
sx q[1];
rz(-2.6076016) q[1];
rz(-pi) q[2];
rz(1.0876065) q[3];
sx q[3];
rz(-1.4677262) q[3];
sx q[3];
rz(-0.3576935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0074924) q[2];
sx q[2];
rz(-1.2780122) q[2];
sx q[2];
rz(-2.210614) q[2];
rz(-1.603294) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61304921) q[0];
sx q[0];
rz(-2.7340524) q[0];
sx q[0];
rz(2.7000309) q[0];
rz(0.26511296) q[1];
sx q[1];
rz(-1.6997063) q[1];
sx q[1];
rz(1.6308867) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.961542) q[0];
sx q[0];
rz(-1.7921711) q[0];
sx q[0];
rz(2.6989797) q[0];
rz(-pi) q[1];
rz(-2.8042107) q[2];
sx q[2];
rz(-2.1357015) q[2];
sx q[2];
rz(0.039410465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35263044) q[1];
sx q[1];
rz(-1.6704541) q[1];
sx q[1];
rz(1.3369192) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7729553) q[3];
sx q[3];
rz(-0.95392272) q[3];
sx q[3];
rz(-1.6807238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.6049217) q[2];
sx q[2];
rz(-1.8646381) q[2];
sx q[2];
rz(-0.46860179) q[2];
rz(-2.129715) q[3];
sx q[3];
rz(-1.7185271) q[3];
sx q[3];
rz(-2.0603777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6877947) q[0];
sx q[0];
rz(-0.84306222) q[0];
sx q[0];
rz(-1.7048365) q[0];
rz(0.94802481) q[1];
sx q[1];
rz(-1.8016305) q[1];
sx q[1];
rz(2.8678133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38588148) q[0];
sx q[0];
rz(-1.4997002) q[0];
sx q[0];
rz(0.12416755) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5735515) q[2];
sx q[2];
rz(-0.92457891) q[2];
sx q[2];
rz(-0.60233145) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54160454) q[1];
sx q[1];
rz(-2.0683824) q[1];
sx q[1];
rz(1.0367839) q[1];
rz(2.2446452) q[3];
sx q[3];
rz(-1.2123143) q[3];
sx q[3];
rz(1.6892387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56129366) q[2];
sx q[2];
rz(-2.0686006) q[2];
sx q[2];
rz(0.44798526) q[2];
rz(0.63730803) q[3];
sx q[3];
rz(-2.1154805) q[3];
sx q[3];
rz(0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4056024) q[0];
sx q[0];
rz(-1.525282) q[0];
sx q[0];
rz(-1.4646144) q[0];
rz(2.6276656) q[1];
sx q[1];
rz(-2.3873139) q[1];
sx q[1];
rz(0.25968459) q[1];
rz(0.17856715) q[2];
sx q[2];
rz(-0.44896508) q[2];
sx q[2];
rz(1.8812897) q[2];
rz(-1.5694554) q[3];
sx q[3];
rz(-2.7272177) q[3];
sx q[3];
rz(0.5275744) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
