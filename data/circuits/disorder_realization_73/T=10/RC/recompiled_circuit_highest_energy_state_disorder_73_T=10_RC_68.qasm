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
rz(0.67983812) q[1];
sx q[1];
rz(-1.4601269) q[1];
sx q[1];
rz(2.3753994) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9840529) q[0];
sx q[0];
rz(-1.6601642) q[0];
sx q[0];
rz(-1.9370701) q[0];
x q[1];
rz(-1.7384479) q[2];
sx q[2];
rz(-0.18356398) q[2];
sx q[2];
rz(-0.50915424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39734617) q[1];
sx q[1];
rz(-1.0434217) q[1];
sx q[1];
rz(-0.1583734) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12285439) q[3];
sx q[3];
rz(-1.7080423) q[3];
sx q[3];
rz(-0.48396971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5863775) q[2];
sx q[2];
rz(-2.5368498) q[2];
sx q[2];
rz(-2.4797454) q[2];
rz(1.9072748) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(0.1525277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.9041651) q[0];
sx q[0];
rz(-0.65218848) q[0];
sx q[0];
rz(-0.40953127) q[0];
rz(-2.883059) q[1];
sx q[1];
rz(-1.2068318) q[1];
sx q[1];
rz(2.5915367) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0830005) q[0];
sx q[0];
rz(-2.6592836) q[0];
sx q[0];
rz(1.7586209) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0070877) q[2];
sx q[2];
rz(-1.1620143) q[2];
sx q[2];
rz(-2.7061912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7118575) q[1];
sx q[1];
rz(-1.1341055) q[1];
sx q[1];
rz(-0.22526422) q[1];
rz(-0.82050555) q[3];
sx q[3];
rz(-1.8295928) q[3];
sx q[3];
rz(-0.25322278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.5162158) q[2];
sx q[2];
rz(-1.3863401) q[2];
sx q[2];
rz(2.5499482) q[2];
rz(-2.6243465) q[3];
sx q[3];
rz(-1.1040684) q[3];
sx q[3];
rz(2.5669317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2194694) q[0];
sx q[0];
rz(-1.2930433) q[0];
sx q[0];
rz(-2.1406232) q[0];
rz(1.7132022) q[1];
sx q[1];
rz(-2.3826471) q[1];
sx q[1];
rz(3.0975814) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.509217) q[0];
sx q[0];
rz(-0.80501834) q[0];
sx q[0];
rz(-0.53765996) q[0];
rz(-0.050878164) q[2];
sx q[2];
rz(-0.83453611) q[2];
sx q[2];
rz(2.5326307) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3487425) q[1];
sx q[1];
rz(-0.89952911) q[1];
sx q[1];
rz(-2.491374) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.040626133) q[3];
sx q[3];
rz(-1.3990132) q[3];
sx q[3];
rz(-0.21200997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9192146) q[2];
sx q[2];
rz(-1.5488449) q[2];
sx q[2];
rz(-0.10981336) q[2];
rz(0.63012704) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(0.28756791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6269864) q[0];
sx q[0];
rz(-1.0543178) q[0];
sx q[0];
rz(-2.0475673) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(-0.72023645) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.076731) q[0];
sx q[0];
rz(-2.0429049) q[0];
sx q[0];
rz(1.9143299) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14671756) q[2];
sx q[2];
rz(-0.31139049) q[2];
sx q[2];
rz(-1.7658833) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.95626895) q[1];
sx q[1];
rz(-0.8517304) q[1];
sx q[1];
rz(1.1683502) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81497808) q[3];
sx q[3];
rz(-1.8061823) q[3];
sx q[3];
rz(0.3334563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.580641) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(-2.9894323) q[2];
rz(-1.2946607) q[3];
sx q[3];
rz(-1.4363689) q[3];
sx q[3];
rz(-1.7744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072435943) q[0];
sx q[0];
rz(-0.56273383) q[0];
sx q[0];
rz(-2.7602472) q[0];
rz(-2.5533679) q[1];
sx q[1];
rz(-1.7344079) q[1];
sx q[1];
rz(-0.064195976) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4669111) q[0];
sx q[0];
rz(-1.7480326) q[0];
sx q[0];
rz(-2.7655274) q[0];
rz(2.5975294) q[2];
sx q[2];
rz(-0.69101101) q[2];
sx q[2];
rz(-0.34073433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.076157454) q[1];
sx q[1];
rz(-2.3373051) q[1];
sx q[1];
rz(1.1706513) q[1];
x q[2];
rz(-1.776657) q[3];
sx q[3];
rz(-1.4511709) q[3];
sx q[3];
rz(-0.72577602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.75307816) q[2];
sx q[2];
rz(-1.831275) q[2];
sx q[2];
rz(0.7927967) q[2];
rz(3.0068398) q[3];
sx q[3];
rz(-0.56551352) q[3];
sx q[3];
rz(-1.6511065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43869269) q[0];
sx q[0];
rz(-1.5543992) q[0];
sx q[0];
rz(0.95293522) q[0];
rz(-1.8903525) q[1];
sx q[1];
rz(-1.4616936) q[1];
sx q[1];
rz(0.23426637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.983646) q[0];
sx q[0];
rz(-1.0704213) q[0];
sx q[0];
rz(1.4051132) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45601797) q[2];
sx q[2];
rz(-1.9393183) q[2];
sx q[2];
rz(-2.2599266) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7997879) q[1];
sx q[1];
rz(-1.4569286) q[1];
sx q[1];
rz(2.1569685) q[1];
rz(-pi) q[2];
rz(2.2053777) q[3];
sx q[3];
rz(-1.6394947) q[3];
sx q[3];
rz(0.19734624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.970152) q[2];
sx q[2];
rz(-2.5415387) q[2];
sx q[2];
rz(0.94685143) q[2];
rz(-1.6998477) q[3];
sx q[3];
rz(-1.8141247) q[3];
sx q[3];
rz(-0.12137752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22987188) q[0];
sx q[0];
rz(-1.0026824) q[0];
sx q[0];
rz(-2.7666996) q[0];
rz(-2.3986469) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(0.65623647) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0716055) q[0];
sx q[0];
rz(-0.90815571) q[0];
sx q[0];
rz(-1.9020984) q[0];
rz(-0.71509174) q[2];
sx q[2];
rz(-0.77132853) q[2];
sx q[2];
rz(0.027931438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6306047) q[1];
sx q[1];
rz(-1.1527677) q[1];
sx q[1];
rz(-1.9782778) q[1];
rz(-pi) q[2];
rz(2.8493622) q[3];
sx q[3];
rz(-2.5295895) q[3];
sx q[3];
rz(-2.0296869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52834213) q[2];
sx q[2];
rz(-2.2684147) q[2];
sx q[2];
rz(2.4916416) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36986247) q[0];
sx q[0];
rz(-1.5700392) q[0];
sx q[0];
rz(-0.318203) q[0];
rz(-2.8089306) q[1];
sx q[1];
rz(-2.5883364) q[1];
sx q[1];
rz(0.37700787) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72412598) q[0];
sx q[0];
rz(-1.2340349) q[0];
sx q[0];
rz(-0.048857852) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.105806) q[2];
sx q[2];
rz(-2.0376959) q[2];
sx q[2];
rz(-0.70718471) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3935289) q[1];
sx q[1];
rz(-1.1989582) q[1];
sx q[1];
rz(-2.4213874) q[1];
x q[2];
rz(-1.351736) q[3];
sx q[3];
rz(-2.6483834) q[3];
sx q[3];
rz(1.0194609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0074924) q[2];
sx q[2];
rz(-1.8635805) q[2];
sx q[2];
rz(0.93097869) q[2];
rz(1.5382986) q[3];
sx q[3];
rz(-1.8036017) q[3];
sx q[3];
rz(-1.3163756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5285434) q[0];
sx q[0];
rz(-0.40754023) q[0];
sx q[0];
rz(0.44156179) q[0];
rz(0.26511296) q[1];
sx q[1];
rz(-1.4418863) q[1];
sx q[1];
rz(-1.6308867) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8245659) q[0];
sx q[0];
rz(-2.6500033) q[0];
sx q[0];
rz(-0.48382171) q[0];
rz(2.8042107) q[2];
sx q[2];
rz(-2.1357015) q[2];
sx q[2];
rz(3.1021822) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6137707) q[1];
sx q[1];
rz(-2.8877259) q[1];
sx q[1];
rz(-1.978118) q[1];
rz(0.36863736) q[3];
sx q[3];
rz(-0.95392272) q[3];
sx q[3];
rz(-1.6807238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.6049217) q[2];
sx q[2];
rz(-1.2769545) q[2];
sx q[2];
rz(-2.6729909) q[2];
rz(-2.129715) q[3];
sx q[3];
rz(-1.4230655) q[3];
sx q[3];
rz(-1.0812149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.45379797) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(1.7048365) q[0];
rz(-2.1935678) q[1];
sx q[1];
rz(-1.3399622) q[1];
sx q[1];
rz(0.2737793) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38588148) q[0];
sx q[0];
rz(-1.4997002) q[0];
sx q[0];
rz(-3.0174251) q[0];
rz(1.5680412) q[2];
sx q[2];
rz(-0.92457891) q[2];
sx q[2];
rz(0.60233145) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54160454) q[1];
sx q[1];
rz(-1.0732103) q[1];
sx q[1];
rz(-2.1048087) q[1];
rz(-pi) q[2];
rz(1.0300566) q[3];
sx q[3];
rz(-2.3916837) q[3];
sx q[3];
rz(0.53242079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56129366) q[2];
sx q[2];
rz(-1.0729921) q[2];
sx q[2];
rz(-2.6936074) q[2];
rz(-2.5042846) q[3];
sx q[3];
rz(-1.0261122) q[3];
sx q[3];
rz(-0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4056024) q[0];
sx q[0];
rz(-1.6163106) q[0];
sx q[0];
rz(1.6769782) q[0];
rz(-0.51392705) q[1];
sx q[1];
rz(-2.3873139) q[1];
sx q[1];
rz(0.25968459) q[1];
rz(1.6561618) q[2];
sx q[2];
rz(-2.0121147) q[2];
sx q[2];
rz(1.683563) q[2];
rz(1.985171) q[3];
sx q[3];
rz(-1.5713362) q[3];
sx q[3];
rz(2.0995981) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
