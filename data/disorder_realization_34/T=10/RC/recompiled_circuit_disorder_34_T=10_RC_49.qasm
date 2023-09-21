OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570219) q[0];
sx q[0];
rz(-1.5091358) q[0];
sx q[0];
rz(-1.8210568) q[0];
x q[1];
rz(1.159367) q[2];
sx q[2];
rz(-1.7147658) q[2];
sx q[2];
rz(1.6978482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.767747) q[1];
sx q[1];
rz(-2.0031553) q[1];
sx q[1];
rz(-0.14111622) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8738535) q[3];
sx q[3];
rz(-2.8149238) q[3];
sx q[3];
rz(-0.79145811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2471182) q[2];
sx q[2];
rz(-2.1366182) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.2540934) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(-2.7665566) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(2.9017752) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12407423) q[0];
sx q[0];
rz(-0.29323175) q[0];
sx q[0];
rz(-2.571884) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8717143) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(2.2965455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61566831) q[1];
sx q[1];
rz(-0.79155542) q[1];
sx q[1];
rz(3.0676003) q[1];
x q[2];
rz(-1.5544942) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-0.53888884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6372765) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(-1.3298539) q[2];
rz(-1.3416393) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2770237) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(-1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(-2.0756762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0695755) q[0];
sx q[0];
rz(-1.6117275) q[0];
sx q[0];
rz(-1.3224365) q[0];
x q[1];
rz(-1.714528) q[2];
sx q[2];
rz(-1.020069) q[2];
sx q[2];
rz(-1.3155754) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2324893) q[1];
sx q[1];
rz(-2.4965582) q[1];
sx q[1];
rz(2.049936) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4170253) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.2325341) q[2];
sx q[2];
rz(3.1090453) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(-0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(-0.70670635) q[0];
rz(1.9354405) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.6548086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.810881) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(2.7964554) q[0];
rz(0.47866486) q[2];
sx q[2];
rz(-2.475111) q[2];
sx q[2];
rz(2.037231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50258892) q[1];
sx q[1];
rz(-1.8559615) q[1];
sx q[1];
rz(1.7267978) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1348226) q[3];
sx q[3];
rz(-2.3151708) q[3];
sx q[3];
rz(-0.063746728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(2.13307) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0356045) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(-1.6812356) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(-0.016074093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094729312) q[0];
sx q[0];
rz(-2.097058) q[0];
sx q[0];
rz(-2.2577283) q[0];
rz(-pi) q[1];
rz(1.1400914) q[2];
sx q[2];
rz(-1.5783731) q[2];
sx q[2];
rz(-1.8097144) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88735089) q[1];
sx q[1];
rz(-2.4130531) q[1];
sx q[1];
rz(-0.90708797) q[1];
rz(-3.1297242) q[3];
sx q[3];
rz(-0.36232811) q[3];
sx q[3];
rz(0.12479898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1092704) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(2.5197022) q[2];
rz(-1.0970998) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(-0.81470195) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46366102) q[0];
sx q[0];
rz(-1.7860798) q[0];
sx q[0];
rz(2.2023375) q[0];
rz(-pi) q[1];
rz(-1.1512418) q[2];
sx q[2];
rz(-2.662979) q[2];
sx q[2];
rz(1.4065557) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.280937) q[1];
sx q[1];
rz(-2.5334364) q[1];
sx q[1];
rz(-1.6673191) q[1];
x q[2];
rz(1.46035) q[3];
sx q[3];
rz(-1.1718281) q[3];
sx q[3];
rz(0.16050592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1427052) q[2];
sx q[2];
rz(-0.83054101) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(2.9283004) q[3];
sx q[3];
rz(-2.8010938) q[3];
sx q[3];
rz(-1.3903769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.6886728) q[1];
sx q[1];
rz(2.4408128) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66877767) q[0];
sx q[0];
rz(-0.46572177) q[0];
sx q[0];
rz(-1.8229683) q[0];
rz(-pi) q[1];
rz(2.8220196) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(-0.061352913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.62774819) q[1];
sx q[1];
rz(-0.57250896) q[1];
sx q[1];
rz(1.4036914) q[1];
rz(0.084214597) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(-2.6654411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(3.11943) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3547524) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(-1.4165075) q[0];
rz(-1.7658866) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1336846) q[0];
sx q[0];
rz(-2.3500751) q[0];
sx q[0];
rz(-1.3484811) q[0];
x q[1];
rz(-2.5052091) q[2];
sx q[2];
rz(-0.6837662) q[2];
sx q[2];
rz(-0.46905876) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8430431) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(-0.72960735) q[1];
rz(0.91120054) q[3];
sx q[3];
rz(-1.3174787) q[3];
sx q[3];
rz(0.60318702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(-0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4760251) q[0];
sx q[0];
rz(-1.6044171) q[0];
sx q[0];
rz(0.72203565) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
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
rz(-0.82775263) q[0];
sx q[0];
rz(1.6270431) q[0];
rz(3.0579371) q[2];
sx q[2];
rz(-1.1612648) q[2];
sx q[2];
rz(1.4510029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.939146) q[1];
sx q[1];
rz(-2.2447963) q[1];
sx q[1];
rz(-2.0161611) q[1];
rz(-pi) q[2];
rz(-0.19946675) q[3];
sx q[3];
rz(-1.27928) q[3];
sx q[3];
rz(0.22649543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1086796) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(-1.9101248) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4817081) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(-0.36614059) q[0];
rz(-pi) q[1];
rz(-0.88629006) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(2.431589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8755175) q[1];
sx q[1];
rz(-2.7873758) q[1];
sx q[1];
rz(-0.43550272) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8929385) q[3];
sx q[3];
rz(-1.5710186) q[3];
sx q[3];
rz(-0.9957046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(1.1408268) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.39682) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(-0.74950779) q[2];
sx q[2];
rz(-1.3386249) q[2];
sx q[2];
rz(-2.7809536) q[2];
rz(1.0127388) q[3];
sx q[3];
rz(-1.8752718) q[3];
sx q[3];
rz(0.74753052) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
