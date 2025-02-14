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
rz(-2.7091205) q[0];
rz(-2.4617545) q[1];
sx q[1];
rz(-1.6814657) q[1];
sx q[1];
rz(0.7661933) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1575398) q[0];
sx q[0];
rz(-1.6601642) q[0];
sx q[0];
rz(-1.9370701) q[0];
rz(1.7384479) q[2];
sx q[2];
rz(-2.9580287) q[2];
sx q[2];
rz(-0.50915424) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.39734617) q[1];
sx q[1];
rz(-1.0434217) q[1];
sx q[1];
rz(-0.1583734) q[1];
rz(1.7090714) q[3];
sx q[3];
rz(-1.449103) q[3];
sx q[3];
rz(-1.0699348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5863775) q[2];
sx q[2];
rz(-2.5368498) q[2];
sx q[2];
rz(-2.4797454) q[2];
rz(1.2343179) q[3];
sx q[3];
rz(-1.581278) q[3];
sx q[3];
rz(0.1525277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9041651) q[0];
sx q[0];
rz(-2.4894042) q[0];
sx q[0];
rz(-2.7320614) q[0];
rz(2.883059) q[1];
sx q[1];
rz(-1.2068318) q[1];
sx q[1];
rz(-2.5915367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15274238) q[0];
sx q[0];
rz(-1.0976731) q[0];
sx q[0];
rz(0.097448601) q[0];
rz(-2.2520355) q[2];
sx q[2];
rz(-2.4585138) q[2];
sx q[2];
rz(-1.4448593) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42973516) q[1];
sx q[1];
rz(-2.0074872) q[1];
sx q[1];
rz(-0.22526422) q[1];
rz(2.7943479) q[3];
sx q[3];
rz(-0.85111516) q[3];
sx q[3];
rz(-1.5517476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6253769) q[2];
sx q[2];
rz(-1.3863401) q[2];
sx q[2];
rz(-0.59164444) q[2];
rz(2.6243465) q[3];
sx q[3];
rz(-1.1040684) q[3];
sx q[3];
rz(0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2194694) q[0];
sx q[0];
rz(-1.8485494) q[0];
sx q[0];
rz(-2.1406232) q[0];
rz(1.7132022) q[1];
sx q[1];
rz(-0.75894558) q[1];
sx q[1];
rz(-3.0975814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45346584) q[0];
sx q[0];
rz(-1.1926873) q[0];
sx q[0];
rz(-2.4124958) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6268544) q[2];
sx q[2];
rz(-0.73768697) q[2];
sx q[2];
rz(0.53327582) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7928501) q[1];
sx q[1];
rz(-0.89952911) q[1];
sx q[1];
rz(2.491374) q[1];
rz(-pi) q[2];
rz(1.3408362) q[3];
sx q[3];
rz(-0.17647568) q[3];
sx q[3];
rz(-0.0214487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9192146) q[2];
sx q[2];
rz(-1.5488449) q[2];
sx q[2];
rz(0.10981336) q[2];
rz(-0.63012704) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(2.8540247) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.0940254) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(2.4213562) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7313774) q[0];
sx q[0];
rz(-2.5654554) q[0];
sx q[0];
rz(0.58310862) q[0];
x q[1];
rz(0.30825395) q[2];
sx q[2];
rz(-1.6156019) q[2];
sx q[2];
rz(-0.055331478) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1853237) q[1];
sx q[1];
rz(-0.8517304) q[1];
sx q[1];
rz(1.1683502) q[1];
x q[2];
rz(-0.31836001) q[3];
sx q[3];
rz(-2.3009319) q[3];
sx q[3];
rz(2.1206253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.580641) q[2];
sx q[2];
rz(-1.6683942) q[2];
sx q[2];
rz(2.9894323) q[2];
rz(1.8469319) q[3];
sx q[3];
rz(-1.4363689) q[3];
sx q[3];
rz(1.367182) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072435943) q[0];
sx q[0];
rz(-0.56273383) q[0];
sx q[0];
rz(2.7602472) q[0];
rz(-2.5533679) q[1];
sx q[1];
rz(-1.7344079) q[1];
sx q[1];
rz(3.0773967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682024) q[0];
sx q[0];
rz(-1.2009092) q[0];
sx q[0];
rz(-1.7610407) q[0];
rz(0.54406329) q[2];
sx q[2];
rz(-2.4505816) q[2];
sx q[2];
rz(-0.34073433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9323247) q[1];
sx q[1];
rz(-1.8552244) q[1];
sx q[1];
rz(0.80764215) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.776657) q[3];
sx q[3];
rz(-1.6904217) q[3];
sx q[3];
rz(-2.4158166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3885145) q[2];
sx q[2];
rz(-1.831275) q[2];
sx q[2];
rz(-0.7927967) q[2];
rz(0.13475284) q[3];
sx q[3];
rz(-2.5760791) q[3];
sx q[3];
rz(1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43869269) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(2.1886574) q[0];
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
rz(0.4933245) q[0];
sx q[0];
rz(-2.6167194) q[0];
sx q[0];
rz(2.8486445) q[0];
x q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3418048) q[1];
sx q[1];
rz(-1.4569286) q[1];
sx q[1];
rz(-0.98462414) q[1];
x q[2];
rz(-1.4552507) q[3];
sx q[3];
rz(-0.63777861) q[3];
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
rz(-1.8141247) q[3];
sx q[3];
rz(0.12137752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.22987188) q[0];
sx q[0];
rz(-1.0026824) q[0];
sx q[0];
rz(-2.7666996) q[0];
rz(-0.74294576) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(2.4853562) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5617758) q[0];
sx q[0];
rz(-0.72951376) q[0];
sx q[0];
rz(-0.39493409) q[0];
x q[1];
rz(-2.5083582) q[2];
sx q[2];
rz(-1.0960964) q[2];
sx q[2];
rz(-2.1556319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.446653) q[1];
sx q[1];
rz(-0.57523721) q[1];
sx q[1];
rz(2.4131362) q[1];
x q[2];
rz(0.2922304) q[3];
sx q[3];
rz(-0.61200313) q[3];
sx q[3];
rz(-2.0296869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52834213) q[2];
sx q[2];
rz(-2.2684147) q[2];
sx q[2];
rz(-0.64995107) q[2];
rz(-2.2925099) q[3];
sx q[3];
rz(-2.504118) q[3];
sx q[3];
rz(-1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36986247) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(-2.8233897) q[0];
rz(-0.33266208) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(0.37700787) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57721389) q[0];
sx q[0];
rz(-2.8014392) q[0];
sx q[0];
rz(-1.4321912) q[0];
rz(-pi) q[1];
rz(-3.105806) q[2];
sx q[2];
rz(-2.0376959) q[2];
sx q[2];
rz(2.4344079) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7480638) q[1];
sx q[1];
rz(-1.1989582) q[1];
sx q[1];
rz(2.4213874) q[1];
x q[2];
rz(3.0253103) q[3];
sx q[3];
rz(-1.0903889) q[3];
sx q[3];
rz(-1.2670328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13410021) q[2];
sx q[2];
rz(-1.2780122) q[2];
sx q[2];
rz(-2.210614) q[2];
rz(1.603294) q[3];
sx q[3];
rz(-1.8036017) q[3];
sx q[3];
rz(-1.8252171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304921) q[0];
sx q[0];
rz(-0.40754023) q[0];
sx q[0];
rz(2.7000309) q[0];
rz(0.26511296) q[1];
sx q[1];
rz(-1.4418863) q[1];
sx q[1];
rz(1.5107059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8245659) q[0];
sx q[0];
rz(-0.49158934) q[0];
sx q[0];
rz(2.6577709) q[0];
rz(1.0894906) q[2];
sx q[2];
rz(-2.4931452) q[2];
sx q[2];
rz(0.61948739) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5278219) q[1];
sx q[1];
rz(-0.25386679) q[1];
sx q[1];
rz(1.1634747) q[1];
x q[2];
rz(2.7729553) q[3];
sx q[3];
rz(-0.95392272) q[3];
sx q[3];
rz(-1.4608689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.536671) q[2];
sx q[2];
rz(-1.2769545) q[2];
sx q[2];
rz(2.6729909) q[2];
rz(-1.0118777) q[3];
sx q[3];
rz(-1.4230655) q[3];
sx q[3];
rz(-2.0603777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6877947) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(-1.4367562) q[0];
rz(2.1935678) q[1];
sx q[1];
rz(-1.3399622) q[1];
sx q[1];
rz(-0.2737793) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.439246) q[0];
sx q[0];
rz(-2.9986023) q[0];
sx q[0];
rz(0.52185328) q[0];
rz(2.4953734) q[2];
sx q[2];
rz(-1.5729959) q[2];
sx q[2];
rz(-0.96680582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5999881) q[1];
sx q[1];
rz(-1.0732103) q[1];
sx q[1];
rz(-1.0367839) q[1];
rz(-1.0300566) q[3];
sx q[3];
rz(-2.3916837) q[3];
sx q[3];
rz(-0.53242079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.580299) q[2];
sx q[2];
rz(-2.0686006) q[2];
sx q[2];
rz(0.44798526) q[2];
rz(-0.63730803) q[3];
sx q[3];
rz(-1.0261122) q[3];
sx q[3];
rz(-2.2465433) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4056024) q[0];
sx q[0];
rz(-1.525282) q[0];
sx q[0];
rz(-1.4646144) q[0];
rz(-2.6276656) q[1];
sx q[1];
rz(-0.75427873) q[1];
sx q[1];
rz(-2.8819081) q[1];
rz(0.44272896) q[2];
sx q[2];
rz(-1.4936269) q[2];
sx q[2];
rz(0.14930162) q[2];
rz(-1.5721372) q[3];
sx q[3];
rz(-0.41437498) q[3];
sx q[3];
rz(-2.6140183) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
