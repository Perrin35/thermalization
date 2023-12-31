OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(5.8435506) q[0];
sx q[0];
rz(6.2018659) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(0.78279701) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047526) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(0.063637861) q[0];
rz(-pi) q[1];
rz(1.9185669) q[2];
sx q[2];
rz(-0.43453056) q[2];
sx q[2];
rz(0.44473106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0040972) q[1];
sx q[1];
rz(-1.4427408) q[1];
sx q[1];
rz(-2.0069684) q[1];
rz(-pi) q[2];
rz(0.26773914) q[3];
sx q[3];
rz(-0.32666884) q[3];
sx q[3];
rz(0.79145811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.2540934) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-0.23981747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46562425) q[0];
sx q[0];
rz(-1.8166607) q[0];
sx q[0];
rz(1.7322391) q[0];
rz(2.9302836) q[2];
sx q[2];
rz(-1.276187) q[2];
sx q[2];
rz(-0.78795563) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51057306) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(1.6454979) q[1];
x q[2];
rz(-0.014157045) q[3];
sx q[3];
rz(-1.5870968) q[3];
sx q[3];
rz(2.1099159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(-1.3298539) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.8168861) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-2.4734316) q[0];
rz(1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0720172) q[0];
sx q[0];
rz(-1.6117275) q[0];
sx q[0];
rz(1.8191562) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5553603) q[2];
sx q[2];
rz(-1.4484324) q[2];
sx q[2];
rz(-2.8107779) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.26820688) q[1];
sx q[1];
rz(-1.8516487) q[1];
sx q[1];
rz(-2.1594949) q[1];
rz(-1.7245674) q[3];
sx q[3];
rz(-2.4745998) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9540017) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-1.1266174) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.6548086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0006205) q[0];
sx q[0];
rz(-1.9019433) q[0];
sx q[0];
rz(1.8676057) q[0];
x q[1];
rz(-1.2232259) q[2];
sx q[2];
rz(-0.98993694) q[2];
sx q[2];
rz(2.6207404) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1176016) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(-2.8531122) q[1];
rz(-1.5781457) q[3];
sx q[3];
rz(-2.3971933) q[3];
sx q[3];
rz(3.0878382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(-1.5530855) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-0.016074093) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8668629) q[0];
sx q[0];
rz(-0.99039536) q[0];
sx q[0];
rz(-2.4972563) q[0];
rz(-pi) q[1];
rz(-1.1400914) q[2];
sx q[2];
rz(-1.5632196) q[2];
sx q[2];
rz(-1.8097144) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4453455) q[1];
sx q[1];
rz(-2.1228585) q[1];
sx q[1];
rz(-0.50260431) q[1];
rz(-pi) q[2];
rz(0.011868422) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(3.0167937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1092704) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(-2.0444929) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9428403) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-2.2348256) q[0];
rz(-0.81470195) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46366102) q[0];
sx q[0];
rz(-1.7860798) q[0];
sx q[0];
rz(-2.2023375) q[0];
x q[1];
rz(1.1282863) q[2];
sx q[2];
rz(-1.3820717) q[2];
sx q[2];
rz(0.21274266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9310301) q[1];
sx q[1];
rz(-1.5157053) q[1];
sx q[1];
rz(-2.1767666) q[1];
x q[2];
rz(1.6812427) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(-2.9810867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99888745) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(-0.21329221) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9796824) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(2.5612223) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(-2.4408128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.472815) q[0];
sx q[0];
rz(-2.6758709) q[0];
sx q[0];
rz(1.3186243) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3406616) q[2];
sx q[2];
rz(-0.33764687) q[2];
sx q[2];
rz(1.9549191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3157981) q[1];
sx q[1];
rz(-2.1343532) q[1];
sx q[1];
rz(-3.0347996) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3015162) q[3];
sx q[3];
rz(-1.5080161) q[3];
sx q[3];
rz(1.9907794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(3.11943) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-0.40063342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.78684029) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(2.2498806) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69670024) q[0];
sx q[0];
rz(-0.80388821) q[0];
sx q[0];
rz(-2.9219887) q[0];
x q[1];
rz(-2.0217998) q[2];
sx q[2];
rz(-2.1037256) q[2];
sx q[2];
rz(-1.9110796) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6139639) q[1];
sx q[1];
rz(-2.2288537) q[1];
sx q[1];
rz(-2.0975454) q[1];
rz(0.91120054) q[3];
sx q[3];
rz(-1.3174787) q[3];
sx q[3];
rz(0.60318702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6531758) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(-0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.7766215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46699539) q[0];
sx q[0];
rz(-2.31384) q[0];
sx q[0];
rz(1.5145495) q[0];
rz(-pi) q[1];
x q[1];
rz(0.083655595) q[2];
sx q[2];
rz(-1.9803279) q[2];
sx q[2];
rz(1.4510029) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2862128) q[1];
sx q[1];
rz(-0.78821048) q[1];
sx q[1];
rz(0.49459313) q[1];
rz(-pi) q[2];
rz(-2.9421259) q[3];
sx q[3];
rz(-1.8623127) q[3];
sx q[3];
rz(-2.9150972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(-0.03406295) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.05474) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.7419086) q[1];
sx q[1];
rz(-2.1733984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23046872) q[0];
sx q[0];
rz(-1.9096806) q[0];
sx q[0];
rz(-1.1662657) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41843157) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(0.07428169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26607516) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(2.7060899) q[1];
rz(-pi) q[2];
rz(-1.2486542) q[3];
sx q[3];
rz(-1.5710186) q[3];
sx q[3];
rz(-0.9957046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(3.1402804) q[2];
rz(-2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447727) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(0.36623476) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(1.8832072) q[2];
sx q[2];
rz(-0.84597833) q[2];
sx q[2];
rz(-0.99920263) q[2];
rz(-1.0127388) q[3];
sx q[3];
rz(-1.2663208) q[3];
sx q[3];
rz(-2.3940621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
