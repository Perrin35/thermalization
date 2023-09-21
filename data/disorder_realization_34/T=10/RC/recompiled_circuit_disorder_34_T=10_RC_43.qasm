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
rz(-0.62277943) q[0];
sx q[0];
rz(-2.8840027) q[0];
sx q[0];
rz(-1.326484) q[0];
rz(-1.2230258) q[2];
sx q[2];
rz(-2.7070621) q[2];
sx q[2];
rz(-0.44473106) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4408735) q[1];
sx q[1];
rz(-0.45342017) q[1];
sx q[1];
rz(-1.2749626) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4814032) q[3];
sx q[3];
rz(-1.2561744) q[3];
sx q[3];
rz(0.50953007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.89447442) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2540934) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(-0.37503606) q[1];
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
rz(-1.4093536) q[0];
x q[1];
rz(-1.8717143) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(-2.2965455) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6310196) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(1.4960947) q[1];
rz(-1.5544942) q[3];
sx q[3];
rz(-1.5849515) q[3];
sx q[3];
rz(0.53888884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5043162) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(-1.8117388) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(0.66816107) q[0];
rz(-1.4913303) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(-1.0659165) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33879963) q[0];
sx q[0];
rz(-2.8899513) q[0];
sx q[0];
rz(-1.405707) q[0];
rz(2.9124444) q[2];
sx q[2];
rz(-2.5742968) q[2];
sx q[2];
rz(-1.0457525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.26820688) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(-2.1594949) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90956456) q[3];
sx q[3];
rz(-1.4759016) q[3];
sx q[3];
rz(0.50300099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(-2.7815946) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(0.70670635) q[0];
rz(-1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.6548086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0006205) q[0];
sx q[0];
rz(-1.9019433) q[0];
sx q[0];
rz(1.2739869) q[0];
rz(-pi) q[1];
rz(-0.47866486) q[2];
sx q[2];
rz(-2.475111) q[2];
sx q[2];
rz(-2.037231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1292124) q[1];
sx q[1];
rz(-2.817569) q[1];
sx q[1];
rz(0.48735168) q[1];
rz(-pi) q[2];
rz(-2.3151822) q[3];
sx q[3];
rz(-1.5658169) q[3];
sx q[3];
rz(-1.5116364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-2.13307) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(3.0116459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0356045) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.460357) q[0];
rz(-1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(3.1255186) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92656266) q[0];
sx q[0];
rz(-0.83850551) q[0];
sx q[0];
rz(-0.82920427) q[0];
x q[1];
rz(-3.1332544) q[2];
sx q[2];
rz(-2.001488) q[2];
sx q[2];
rz(2.8991933) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2542418) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(2.2345047) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36230476) q[3];
sx q[3];
rz(-1.575003) q[3];
sx q[3];
rz(1.4570953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1092704) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(1.2247359) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8794494) q[0];
sx q[0];
rz(-0.95603846) q[0];
sx q[0];
rz(2.8770147) q[0];
rz(-pi) q[1];
rz(0.20829006) q[2];
sx q[2];
rz(-1.1366833) q[2];
sx q[2];
rz(1.8722033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21056255) q[1];
sx q[1];
rz(-1.5157053) q[1];
sx q[1];
rz(2.1767666) q[1];
rz(-pi) q[2];
x q[2];
rz(1.46035) q[3];
sx q[3];
rz(-1.9697646) q[3];
sx q[3];
rz(2.9810867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(-0.93969807) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(0.70077983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0133007) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(2.0237472) q[0];
rz(-pi) q[1];
rz(-1.6875661) q[2];
sx q[2];
rz(-1.2532557) q[2];
sx q[2];
rz(-1.5460528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.339387) q[1];
sx q[1];
rz(-1.6610258) q[1];
sx q[1];
rz(-1.0046563) q[1];
x q[2];
rz(0.084214597) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(0.47615151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3547524) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(1.7250852) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(2.2498806) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4211593) q[0];
sx q[0];
rz(-1.4132858) q[0];
sx q[0];
rz(-0.79173761) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6363836) q[2];
sx q[2];
rz(-0.6837662) q[2];
sx q[2];
rz(-2.6725339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2878694) q[1];
sx q[1];
rz(-2.3239377) q[1];
sx q[1];
rz(-2.5649648) q[1];
rz(2.8250152) q[3];
sx q[3];
rz(-0.93571767) q[3];
sx q[3];
rz(-2.3659335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(0.78424224) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(0.23323664) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.3649712) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7576335) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(3.0804759) q[0];
rz(-pi) q[1];
rz(-3.0579371) q[2];
sx q[2];
rz(-1.9803279) q[2];
sx q[2];
rz(1.4510029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.939146) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(1.1254315) q[1];
rz(-pi) q[2];
rz(-1.2737208) q[3];
sx q[3];
rz(-1.7617412) q[3];
sx q[3];
rz(1.2862658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0329131) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(-1.9101248) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(-1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(0.96819425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68010846) q[0];
sx q[0];
rz(-0.52163863) q[0];
sx q[0];
rz(-2.3011544) q[0];
rz(0.88629006) q[2];
sx q[2];
rz(-0.61294014) q[2];
sx q[2];
rz(-0.71000368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4149975) q[1];
sx q[1];
rz(-1.8907049) q[1];
sx q[1];
rz(-1.725561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8929385) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-2.145888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(-3.1402804) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(-1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.39682) q[0];
sx q[0];
rz(-1.1068494) q[0];
sx q[0];
rz(-1.1882991) q[0];
rz(0.36623476) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-1.8832072) q[2];
sx q[2];
rz(-2.2956143) q[2];
sx q[2];
rz(2.14239) q[2];
rz(-2.1288539) q[3];
sx q[3];
rz(-1.8752718) q[3];
sx q[3];
rz(0.74753052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
