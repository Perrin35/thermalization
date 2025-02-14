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
rz(-2.3753994) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1846518) q[0];
sx q[0];
rz(-0.37654018) q[0];
sx q[0];
rz(1.3256289) q[0];
x q[1];
rz(1.4031448) q[2];
sx q[2];
rz(-2.9580287) q[2];
sx q[2];
rz(0.50915424) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0483472) q[1];
sx q[1];
rz(-1.4340869) q[1];
sx q[1];
rz(-1.0379278) q[1];
rz(0.84505524) q[3];
sx q[3];
rz(-2.9576506) q[3];
sx q[3];
rz(1.2181653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55521518) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(0.66184723) q[2];
rz(-1.9072748) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(-0.1525277) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-1.2068318) q[1];
sx q[1];
rz(-2.5915367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0830005) q[0];
sx q[0];
rz(-0.48230902) q[0];
sx q[0];
rz(-1.7586209) q[0];
x q[1];
rz(-2.134505) q[2];
sx q[2];
rz(-1.1620143) q[2];
sx q[2];
rz(0.43540149) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0747968) q[1];
sx q[1];
rz(-2.6535712) q[1];
sx q[1];
rz(2.0171432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7943479) q[3];
sx q[3];
rz(-0.85111516) q[3];
sx q[3];
rz(-1.5517476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.5162158) q[2];
sx q[2];
rz(-1.7552525) q[2];
sx q[2];
rz(-2.5499482) q[2];
rz(-0.51724616) q[3];
sx q[3];
rz(-2.0375242) q[3];
sx q[3];
rz(-0.57466093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.9221233) q[0];
sx q[0];
rz(-1.8485494) q[0];
sx q[0];
rz(-1.0009694) q[0];
rz(1.7132022) q[1];
sx q[1];
rz(-2.3826471) q[1];
sx q[1];
rz(-0.044011291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.509217) q[0];
sx q[0];
rz(-0.80501834) q[0];
sx q[0];
rz(-2.6039327) q[0];
rz(-pi) q[1];
rz(2.3077009) q[2];
sx q[2];
rz(-1.6084889) q[2];
sx q[2];
rz(2.1455763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3487425) q[1];
sx q[1];
rz(-0.89952911) q[1];
sx q[1];
rz(-0.6502186) q[1];
x q[2];
rz(-0.040626133) q[3];
sx q[3];
rz(-1.3990132) q[3];
sx q[3];
rz(-0.21200997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22237805) q[2];
sx q[2];
rz(-1.5488449) q[2];
sx q[2];
rz(-3.0317793) q[2];
rz(0.63012704) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(0.28756791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5146062) q[0];
sx q[0];
rz(-1.0543178) q[0];
sx q[0];
rz(-1.0940254) q[0];
rz(0.47359723) q[1];
sx q[1];
rz(-0.70643276) q[1];
sx q[1];
rz(2.4213562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.076731) q[0];
sx q[0];
rz(-1.0986878) q[0];
sx q[0];
rz(1.9143299) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6178151) q[2];
sx q[2];
rz(-1.2628619) q[2];
sx q[2];
rz(-1.5297254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2536972) q[1];
sx q[1];
rz(-1.869939) q[1];
sx q[1];
rz(0.76051416) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9071753) q[3];
sx q[3];
rz(-0.78463882) q[3];
sx q[3];
rz(1.6616846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.56095162) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(-0.15216039) q[2];
rz(1.8469319) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(-1.367182) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072435943) q[0];
sx q[0];
rz(-2.5788588) q[0];
sx q[0];
rz(-2.7602472) q[0];
rz(0.58822477) q[1];
sx q[1];
rz(-1.7344079) q[1];
sx q[1];
rz(-0.064195976) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9682024) q[0];
sx q[0];
rz(-1.9406835) q[0];
sx q[0];
rz(1.3805519) q[0];
rz(2.5257686) q[2];
sx q[2];
rz(-1.2346137) q[2];
sx q[2];
rz(-1.4752963) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5178495) q[1];
sx q[1];
rz(-0.84539834) q[1];
sx q[1];
rz(-0.38442594) q[1];
rz(-pi) q[2];
rz(3.0194123) q[3];
sx q[3];
rz(-1.3664277) q[3];
sx q[3];
rz(2.271657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75307816) q[2];
sx q[2];
rz(-1.3103176) q[2];
sx q[2];
rz(2.348796) q[2];
rz(-3.0068398) q[3];
sx q[3];
rz(-0.56551352) q[3];
sx q[3];
rz(-1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-2.7029) q[0];
sx q[0];
rz(-1.5543992) q[0];
sx q[0];
rz(0.95293522) q[0];
rz(1.2512402) q[1];
sx q[1];
rz(-1.6798991) q[1];
sx q[1];
rz(-0.23426637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3327985) q[0];
sx q[0];
rz(-1.7160128) q[0];
sx q[0];
rz(0.50620075) q[0];
rz(1.1646005) q[2];
sx q[2];
rz(-1.147454) q[2];
sx q[2];
rz(-2.6273531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3418048) q[1];
sx q[1];
rz(-1.684664) q[1];
sx q[1];
rz(-0.98462414) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2053777) q[3];
sx q[3];
rz(-1.6394947) q[3];
sx q[3];
rz(2.9442464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.17144063) q[2];
sx q[2];
rz(-2.5415387) q[2];
sx q[2];
rz(-0.94685143) q[2];
rz(1.6998477) q[3];
sx q[3];
rz(-1.3274679) q[3];
sx q[3];
rz(-0.12137752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117208) q[0];
sx q[0];
rz(-2.1389102) q[0];
sx q[0];
rz(2.7666996) q[0];
rz(0.74294576) q[1];
sx q[1];
rz(-2.3187147) q[1];
sx q[1];
rz(0.65623647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8493238) q[0];
sx q[0];
rz(-1.3114623) q[0];
sx q[0];
rz(-0.68993129) q[0];
rz(-2.4265009) q[2];
sx q[2];
rz(-2.3702641) q[2];
sx q[2];
rz(0.027931438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.446653) q[1];
sx q[1];
rz(-2.5663554) q[1];
sx q[1];
rz(-0.72845643) q[1];
rz(-2.5498059) q[3];
sx q[3];
rz(-1.7370708) q[3];
sx q[3];
rz(0.70032728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52834213) q[2];
sx q[2];
rz(-0.87317792) q[2];
sx q[2];
rz(2.4916416) q[2];
rz(2.2925099) q[3];
sx q[3];
rz(-2.504118) q[3];
sx q[3];
rz(1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36986247) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(-0.318203) q[0];
rz(-0.33266208) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(0.37700787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5643788) q[0];
sx q[0];
rz(-2.8014392) q[0];
sx q[0];
rz(-1.7094014) q[0];
rz(-pi) q[1];
rz(-1.1036393) q[2];
sx q[2];
rz(-1.6027513) q[2];
sx q[2];
rz(-0.87972537) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3935289) q[1];
sx q[1];
rz(-1.1989582) q[1];
sx q[1];
rz(0.7202052) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11628233) q[3];
sx q[3];
rz(-2.0512037) q[3];
sx q[3];
rz(-1.2670328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0074924) q[2];
sx q[2];
rz(-1.8635805) q[2];
sx q[2];
rz(-2.210614) q[2];
rz(-1.5382986) q[3];
sx q[3];
rz(-1.337991) q[3];
sx q[3];
rz(-1.3163756) q[3];
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
rz(-0.61304921) q[0];
sx q[0];
rz(-2.7340524) q[0];
sx q[0];
rz(2.7000309) q[0];
rz(2.8764797) q[1];
sx q[1];
rz(-1.6997063) q[1];
sx q[1];
rz(1.5107059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8545452) q[0];
sx q[0];
rz(-2.0018739) q[0];
sx q[0];
rz(-1.8148941) q[0];
x q[1];
rz(-1.0894906) q[2];
sx q[2];
rz(-2.4931452) q[2];
sx q[2];
rz(-0.61948739) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1944674) q[1];
sx q[1];
rz(-1.3381011) q[1];
sx q[1];
rz(0.10242771) q[1];
rz(-0.36863736) q[3];
sx q[3];
rz(-0.95392272) q[3];
sx q[3];
rz(-1.4608689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.536671) q[2];
sx q[2];
rz(-1.8646381) q[2];
sx q[2];
rz(-2.6729909) q[2];
rz(2.129715) q[3];
sx q[3];
rz(-1.4230655) q[3];
sx q[3];
rz(-2.0603777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6877947) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(-1.7048365) q[0];
rz(0.94802481) q[1];
sx q[1];
rz(-1.8016305) q[1];
sx q[1];
rz(2.8678133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.439246) q[0];
sx q[0];
rz(-0.14299031) q[0];
sx q[0];
rz(0.52185328) q[0];
x q[1];
rz(-1.5735515) q[2];
sx q[2];
rz(-0.92457891) q[2];
sx q[2];
rz(0.60233145) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3042766) q[1];
sx q[1];
rz(-1.1070861) q[1];
sx q[1];
rz(-2.5786693) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6945029) q[3];
sx q[3];
rz(-0.94674094) q[3];
sx q[3];
rz(0.15472356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.580299) q[2];
sx q[2];
rz(-2.0686006) q[2];
sx q[2];
rz(-0.44798526) q[2];
rz(-0.63730803) q[3];
sx q[3];
rz(-2.1154805) q[3];
sx q[3];
rz(-0.89504939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-1.6561618) q[2];
sx q[2];
rz(-1.1294779) q[2];
sx q[2];
rz(-1.4580296) q[2];
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
