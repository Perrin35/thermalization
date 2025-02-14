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
rz(-0.66752783) q[0];
sx q[0];
rz(-2.3025371) q[0];
sx q[0];
rz(2.7091205) q[0];
rz(-2.4617545) q[1];
sx q[1];
rz(-1.6814657) q[1];
sx q[1];
rz(0.7661933) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9569409) q[0];
sx q[0];
rz(-0.37654018) q[0];
sx q[0];
rz(-1.8159637) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030969663) q[2];
sx q[2];
rz(-1.3898347) q[2];
sx q[2];
rz(0.33869264) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0932454) q[1];
sx q[1];
rz(-1.4340869) q[1];
sx q[1];
rz(-1.0379278) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2965374) q[3];
sx q[3];
rz(-2.9576506) q[3];
sx q[3];
rz(-1.9234274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5863775) q[2];
sx q[2];
rz(-0.60474288) q[2];
sx q[2];
rz(2.4797454) q[2];
rz(1.2343179) q[3];
sx q[3];
rz(-1.5603147) q[3];
sx q[3];
rz(-0.1525277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9041651) q[0];
sx q[0];
rz(-0.65218848) q[0];
sx q[0];
rz(2.7320614) q[0];
rz(-0.25853363) q[1];
sx q[1];
rz(-1.9347609) q[1];
sx q[1];
rz(-0.55005598) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9888503) q[0];
sx q[0];
rz(-2.0439195) q[0];
sx q[0];
rz(3.0441441) q[0];
rz(0.47357719) q[2];
sx q[2];
rz(-2.0832295) q[2];
sx q[2];
rz(-2.2523675) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.066795863) q[1];
sx q[1];
rz(-0.48802146) q[1];
sx q[1];
rz(-1.1244494) q[1];
x q[2];
rz(2.7943479) q[3];
sx q[3];
rz(-0.85111516) q[3];
sx q[3];
rz(1.589845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.5162158) q[2];
sx q[2];
rz(-1.3863401) q[2];
sx q[2];
rz(-2.5499482) q[2];
rz(2.6243465) q[3];
sx q[3];
rz(-2.0375242) q[3];
sx q[3];
rz(2.5669317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2194694) q[0];
sx q[0];
rz(-1.2930433) q[0];
sx q[0];
rz(-1.0009694) q[0];
rz(-1.7132022) q[1];
sx q[1];
rz(-2.3826471) q[1];
sx q[1];
rz(0.044011291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79878935) q[0];
sx q[0];
rz(-2.2384423) q[0];
sx q[0];
rz(-2.0602047) q[0];
x q[1];
rz(3.0907145) q[2];
sx q[2];
rz(-0.83453611) q[2];
sx q[2];
rz(-0.60896192) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2340369) q[1];
sx q[1];
rz(-2.2436972) q[1];
sx q[1];
rz(-0.91958811) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3988741) q[3];
sx q[3];
rz(-1.5307685) q[3];
sx q[3];
rz(-1.775858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9192146) q[2];
sx q[2];
rz(-1.5927477) q[2];
sx q[2];
rz(-3.0317793) q[2];
rz(0.63012704) q[3];
sx q[3];
rz(-0.55456847) q[3];
sx q[3];
rz(-2.8540247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(2.0475673) q[0];
rz(-2.6679954) q[1];
sx q[1];
rz(-2.4351599) q[1];
sx q[1];
rz(-2.4213562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4862625) q[0];
sx q[0];
rz(-1.8754601) q[0];
sx q[0];
rz(2.644674) q[0];
rz(-2.9948751) q[2];
sx q[2];
rz(-2.8302022) q[2];
sx q[2];
rz(-1.7658833) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88789543) q[1];
sx q[1];
rz(-1.869939) q[1];
sx q[1];
rz(2.3810785) q[1];
rz(-pi) q[2];
rz(-0.81497808) q[3];
sx q[3];
rz(-1.3354104) q[3];
sx q[3];
rz(2.8081364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.580641) q[2];
sx q[2];
rz(-1.4731984) q[2];
sx q[2];
rz(-2.9894323) q[2];
rz(-1.8469319) q[3];
sx q[3];
rz(-1.7052238) q[3];
sx q[3];
rz(-1.7744106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-3.0691567) q[0];
sx q[0];
rz(-0.56273383) q[0];
sx q[0];
rz(2.7602472) q[0];
rz(2.5533679) q[1];
sx q[1];
rz(-1.4071848) q[1];
sx q[1];
rz(3.0773967) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9682024) q[0];
sx q[0];
rz(-1.2009092) q[0];
sx q[0];
rz(-1.7610407) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5257686) q[2];
sx q[2];
rz(-1.9069789) q[2];
sx q[2];
rz(1.4752963) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.076157454) q[1];
sx q[1];
rz(-2.3373051) q[1];
sx q[1];
rz(-1.9709414) q[1];
x q[2];
rz(3.0194123) q[3];
sx q[3];
rz(-1.3664277) q[3];
sx q[3];
rz(2.271657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3885145) q[2];
sx q[2];
rz(-1.831275) q[2];
sx q[2];
rz(0.7927967) q[2];
rz(-3.0068398) q[3];
sx q[3];
rz(-0.56551352) q[3];
sx q[3];
rz(-1.4904862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7029) q[0];
sx q[0];
rz(-1.5871935) q[0];
sx q[0];
rz(-0.95293522) q[0];
rz(-1.8903525) q[1];
sx q[1];
rz(-1.4616936) q[1];
sx q[1];
rz(0.23426637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8087942) q[0];
sx q[0];
rz(-1.4255798) q[0];
sx q[0];
rz(-0.50620075) q[0];
x q[1];
rz(-0.71990196) q[2];
sx q[2];
rz(-0.57804075) q[2];
sx q[2];
rz(1.8190011) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0820703) q[1];
sx q[1];
rz(-0.59585427) q[1];
sx q[1];
rz(-1.3669307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.686342) q[3];
sx q[3];
rz(-0.63777861) q[3];
sx q[3];
rz(1.4664283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17144063) q[2];
sx q[2];
rz(-0.60005391) q[2];
sx q[2];
rz(2.1947412) q[2];
rz(-1.4417449) q[3];
sx q[3];
rz(-1.8141247) q[3];
sx q[3];
rz(-3.0202151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22987188) q[0];
sx q[0];
rz(-1.0026824) q[0];
sx q[0];
rz(-0.3748931) q[0];
rz(-2.3986469) q[1];
sx q[1];
rz(-0.82287794) q[1];
sx q[1];
rz(2.4853562) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29226881) q[0];
sx q[0];
rz(-1.3114623) q[0];
sx q[0];
rz(0.68993129) q[0];
rz(-pi) q[1];
rz(-2.4265009) q[2];
sx q[2];
rz(-2.3702641) q[2];
sx q[2];
rz(-3.1136612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23327574) q[1];
sx q[1];
rz(-1.2001977) q[1];
sx q[1];
rz(2.6909687) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8493622) q[3];
sx q[3];
rz(-2.5295895) q[3];
sx q[3];
rz(-1.1119058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6132505) q[2];
sx q[2];
rz(-0.87317792) q[2];
sx q[2];
rz(-2.4916416) q[2];
rz(2.2925099) q[3];
sx q[3];
rz(-2.504118) q[3];
sx q[3];
rz(1.953663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717302) q[0];
sx q[0];
rz(-1.5715535) q[0];
sx q[0];
rz(0.318203) q[0];
rz(-2.8089306) q[1];
sx q[1];
rz(-0.55325621) q[1];
sx q[1];
rz(-0.37700787) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4174667) q[0];
sx q[0];
rz(-1.9075577) q[0];
sx q[0];
rz(-3.0927348) q[0];
rz(0.035786619) q[2];
sx q[2];
rz(-1.1038968) q[2];
sx q[2];
rz(-2.4344079) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3935289) q[1];
sx q[1];
rz(-1.1989582) q[1];
sx q[1];
rz(-0.7202052) q[1];
x q[2];
rz(1.351736) q[3];
sx q[3];
rz(-2.6483834) q[3];
sx q[3];
rz(-1.0194609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13410021) q[2];
sx q[2];
rz(-1.2780122) q[2];
sx q[2];
rz(-0.93097869) q[2];
rz(1.603294) q[3];
sx q[3];
rz(-1.337991) q[3];
sx q[3];
rz(-1.3163756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61304921) q[0];
sx q[0];
rz(-2.7340524) q[0];
sx q[0];
rz(0.44156179) q[0];
rz(-0.26511296) q[1];
sx q[1];
rz(-1.4418863) q[1];
sx q[1];
rz(1.6308867) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2870474) q[0];
sx q[0];
rz(-1.1397188) q[0];
sx q[0];
rz(-1.3266985) q[0];
x q[1];
rz(-0.33738193) q[2];
sx q[2];
rz(-2.1357015) q[2];
sx q[2];
rz(-0.039410465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5278219) q[1];
sx q[1];
rz(-0.25386679) q[1];
sx q[1];
rz(-1.1634747) q[1];
rz(-2.7729553) q[3];
sx q[3];
rz(-2.1876699) q[3];
sx q[3];
rz(-1.4608689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.6049217) q[2];
sx q[2];
rz(-1.2769545) q[2];
sx q[2];
rz(2.6729909) q[2];
rz(2.129715) q[3];
sx q[3];
rz(-1.7185271) q[3];
sx q[3];
rz(-1.0812149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.45379797) q[0];
sx q[0];
rz(-2.2985304) q[0];
sx q[0];
rz(-1.7048365) q[0];
rz(-2.1935678) q[1];
sx q[1];
rz(-1.8016305) q[1];
sx q[1];
rz(-0.2737793) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1760491) q[0];
sx q[0];
rz(-1.4469441) q[0];
sx q[0];
rz(-1.6424422) q[0];
rz(-pi) q[1];
x q[1];
rz(0.003652775) q[2];
sx q[2];
rz(-2.4953702) q[2];
sx q[2];
rz(2.5346859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.791509) q[1];
sx q[1];
rz(-2.4285762) q[1];
sx q[1];
rz(-2.3886695) q[1];
rz(-pi) q[2];
rz(-0.44708973) q[3];
sx q[3];
rz(-0.94674094) q[3];
sx q[3];
rz(2.9868691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.56129366) q[2];
sx q[2];
rz(-2.0686006) q[2];
sx q[2];
rz(0.44798526) q[2];
rz(0.63730803) q[3];
sx q[3];
rz(-1.0261122) q[3];
sx q[3];
rz(2.2465433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7359903) q[0];
sx q[0];
rz(-1.6163106) q[0];
sx q[0];
rz(1.6769782) q[0];
rz(2.6276656) q[1];
sx q[1];
rz(-2.3873139) q[1];
sx q[1];
rz(0.25968459) q[1];
rz(-1.6561618) q[2];
sx q[2];
rz(-1.1294779) q[2];
sx q[2];
rz(-1.4580296) q[2];
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
