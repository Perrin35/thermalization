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
rz(-0.44218818) q[0];
sx q[0];
rz(-1.1531885) q[0];
sx q[0];
rz(-0.12864223) q[0];
rz(1.6755942) q[1];
sx q[1];
rz(-2.0425551) q[1];
sx q[1];
rz(-2.0597982) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2483467) q[0];
sx q[0];
rz(-1.5061646) q[0];
sx q[0];
rz(1.6016225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8063534) q[2];
sx q[2];
rz(-1.2150107) q[2];
sx q[2];
rz(0.29976732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1638559) q[1];
sx q[1];
rz(-0.1739991) q[1];
sx q[1];
rz(-1.15088) q[1];
rz(-1.8954738) q[3];
sx q[3];
rz(-0.28936181) q[3];
sx q[3];
rz(1.0288256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.094540207) q[2];
sx q[2];
rz(-0.69539842) q[2];
sx q[2];
rz(-2.409234) q[2];
rz(-3.0426466) q[3];
sx q[3];
rz(-2.1061335) q[3];
sx q[3];
rz(-1.8100479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.052213) q[0];
sx q[0];
rz(-0.78240028) q[0];
sx q[0];
rz(0.041444929) q[0];
rz(2.4502358) q[1];
sx q[1];
rz(-1.8212916) q[1];
sx q[1];
rz(2.2533805) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0109581) q[0];
sx q[0];
rz(-0.97773363) q[0];
sx q[0];
rz(-0.61390604) q[0];
x q[1];
rz(2.5609547) q[2];
sx q[2];
rz(-2.3742445) q[2];
sx q[2];
rz(-1.5754171) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6710818) q[1];
sx q[1];
rz(-0.24674079) q[1];
sx q[1];
rz(-2.9637439) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.290894) q[3];
sx q[3];
rz(-2.0655819) q[3];
sx q[3];
rz(-1.6542357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1516252) q[2];
sx q[2];
rz(-1.7513559) q[2];
sx q[2];
rz(-1.3803231) q[2];
rz(3.0917998) q[3];
sx q[3];
rz(-2.4536665) q[3];
sx q[3];
rz(-0.5602347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3675073) q[0];
sx q[0];
rz(-2.9893576) q[0];
sx q[0];
rz(1.2138858) q[0];
rz(-1.8269352) q[1];
sx q[1];
rz(-1.0379125) q[1];
sx q[1];
rz(1.5710057) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0308098) q[0];
sx q[0];
rz(-2.6834496) q[0];
sx q[0];
rz(-0.58266298) q[0];
rz(0.41061398) q[2];
sx q[2];
rz(-1.4272235) q[2];
sx q[2];
rz(-2.5559354) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6398869) q[1];
sx q[1];
rz(-0.16777953) q[1];
sx q[1];
rz(1.3659992) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5475471) q[3];
sx q[3];
rz(-1.0540773) q[3];
sx q[3];
rz(1.7917716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35459241) q[2];
sx q[2];
rz(-0.74942333) q[2];
sx q[2];
rz(2.836239) q[2];
rz(-2.2954156) q[3];
sx q[3];
rz(-1.927522) q[3];
sx q[3];
rz(2.2727216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884035) q[0];
sx q[0];
rz(-1.0193634) q[0];
sx q[0];
rz(0.97002059) q[0];
rz(2.8963529) q[1];
sx q[1];
rz(-2.0742564) q[1];
sx q[1];
rz(1.6397569) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0787112) q[0];
sx q[0];
rz(-0.20096603) q[0];
sx q[0];
rz(-2.4522792) q[0];
rz(-pi) q[1];
rz(1.7443329) q[2];
sx q[2];
rz(-1.3385941) q[2];
sx q[2];
rz(-2.6628138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0027155) q[1];
sx q[1];
rz(-1.381523) q[1];
sx q[1];
rz(1.7258964) q[1];
rz(-pi) q[2];
rz(2.7795622) q[3];
sx q[3];
rz(-0.76043441) q[3];
sx q[3];
rz(2.8575051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60845) q[2];
sx q[2];
rz(-1.5773982) q[2];
sx q[2];
rz(-1.2987785) q[2];
rz(-0.53457824) q[3];
sx q[3];
rz(-0.60408533) q[3];
sx q[3];
rz(0.66876137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0534441) q[0];
sx q[0];
rz(-1.7867333) q[0];
sx q[0];
rz(0.97766367) q[0];
rz(2.4560302) q[1];
sx q[1];
rz(-2.3473163) q[1];
sx q[1];
rz(2.175144) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9622877) q[0];
sx q[0];
rz(-2.0317269) q[0];
sx q[0];
rz(-2.3910752) q[0];
rz(-0.62018779) q[2];
sx q[2];
rz(-2.475707) q[2];
sx q[2];
rz(2.263139) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9271665) q[1];
sx q[1];
rz(-0.71950235) q[1];
sx q[1];
rz(-1.7822766) q[1];
rz(2.590476) q[3];
sx q[3];
rz(-1.8752021) q[3];
sx q[3];
rz(1.9929043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55296772) q[2];
sx q[2];
rz(-1.3543465) q[2];
sx q[2];
rz(0.28263131) q[2];
rz(-2.1770554) q[3];
sx q[3];
rz(-1.5610361) q[3];
sx q[3];
rz(-0.36981043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033217) q[0];
sx q[0];
rz(-2.7300457) q[0];
sx q[0];
rz(2.0630398) q[0];
rz(2.318553) q[1];
sx q[1];
rz(-1.1230725) q[1];
sx q[1];
rz(-2.3413234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1861411) q[0];
sx q[0];
rz(-1.5793043) q[0];
sx q[0];
rz(1.5415989) q[0];
rz(-0.79485017) q[2];
sx q[2];
rz(-2.7140116) q[2];
sx q[2];
rz(2.7045369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9046968) q[1];
sx q[1];
rz(-2.176732) q[1];
sx q[1];
rz(0.69463457) q[1];
x q[2];
rz(0.54366248) q[3];
sx q[3];
rz(-2.398536) q[3];
sx q[3];
rz(-2.0280968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0018953) q[2];
sx q[2];
rz(-1.8997833) q[2];
sx q[2];
rz(-0.099420698) q[2];
rz(2.5771778) q[3];
sx q[3];
rz(-1.5921009) q[3];
sx q[3];
rz(2.2216589) q[3];
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
rz(2.1426706) q[0];
sx q[0];
rz(-2.8746334) q[0];
sx q[0];
rz(3.1166742) q[0];
rz(2.0492367) q[1];
sx q[1];
rz(-1.9905636) q[1];
sx q[1];
rz(-0.54164642) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9077907) q[0];
sx q[0];
rz(-0.72392091) q[0];
sx q[0];
rz(-1.0095327) q[0];
rz(0.32633467) q[2];
sx q[2];
rz(-1.0251364) q[2];
sx q[2];
rz(1.5991925) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6691963) q[1];
sx q[1];
rz(-1.7150567) q[1];
sx q[1];
rz(-1.902182) q[1];
x q[2];
rz(0.17766134) q[3];
sx q[3];
rz(-0.97148669) q[3];
sx q[3];
rz(0.29479542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3148552) q[2];
sx q[2];
rz(-1.9172226) q[2];
sx q[2];
rz(-1.1401736) q[2];
rz(1.4402116) q[3];
sx q[3];
rz(-0.39339104) q[3];
sx q[3];
rz(-2.9668729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.2575271) q[0];
sx q[0];
rz(-1.9767569) q[0];
sx q[0];
rz(2.2648532) q[0];
rz(0.17403099) q[1];
sx q[1];
rz(-2.048025) q[1];
sx q[1];
rz(1.3678975) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68326658) q[0];
sx q[0];
rz(-1.84671) q[0];
sx q[0];
rz(-2.1349195) q[0];
x q[1];
rz(1.2064916) q[2];
sx q[2];
rz(-1.3576531) q[2];
sx q[2];
rz(-0.079696004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0292082) q[1];
sx q[1];
rz(-1.9416326) q[1];
sx q[1];
rz(0.98354323) q[1];
rz(-pi) q[2];
rz(1.01497) q[3];
sx q[3];
rz(-2.9983302) q[3];
sx q[3];
rz(-1.48326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6690663) q[2];
sx q[2];
rz(-1.8146699) q[2];
sx q[2];
rz(-0.18829045) q[2];
rz(-1.3663728) q[3];
sx q[3];
rz(-1.7732737) q[3];
sx q[3];
rz(-1.202047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.61757225) q[0];
sx q[0];
rz(-3.0078648) q[0];
sx q[0];
rz(0.66620859) q[0];
rz(-2.0062402) q[1];
sx q[1];
rz(-2.1609781) q[1];
sx q[1];
rz(-1.1599783) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85765094) q[0];
sx q[0];
rz(-2.2277624) q[0];
sx q[0];
rz(2.3131671) q[0];
rz(0.23187821) q[2];
sx q[2];
rz(-2.3740951) q[2];
sx q[2];
rz(2.9434443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.682232) q[1];
sx q[1];
rz(-2.5256093) q[1];
sx q[1];
rz(1.5840696) q[1];
x q[2];
rz(-2.1362392) q[3];
sx q[3];
rz(-1.3015038) q[3];
sx q[3];
rz(-2.5776742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74405115) q[2];
sx q[2];
rz(-2.5752189) q[2];
sx q[2];
rz(-2.9446824) q[2];
rz(0.87013733) q[3];
sx q[3];
rz(-1.0715485) q[3];
sx q[3];
rz(-1.2675233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63721913) q[0];
sx q[0];
rz(-1.0727896) q[0];
sx q[0];
rz(-3.1043501) q[0];
rz(1.0875018) q[1];
sx q[1];
rz(-1.0481513) q[1];
sx q[1];
rz(0.16669272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.113362) q[0];
sx q[0];
rz(-2.2489002) q[0];
sx q[0];
rz(-1.9166758) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6564381) q[2];
sx q[2];
rz(-1.8388766) q[2];
sx q[2];
rz(1.175665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5350668) q[1];
sx q[1];
rz(-1.6589612) q[1];
sx q[1];
rz(-2.8463581) q[1];
x q[2];
rz(-2.9965472) q[3];
sx q[3];
rz(-1.6393472) q[3];
sx q[3];
rz(1.9653494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4166261) q[2];
sx q[2];
rz(-0.61369696) q[2];
sx q[2];
rz(1.7566768) q[2];
rz(0.40063217) q[3];
sx q[3];
rz(-1.2847565) q[3];
sx q[3];
rz(-1.0111151) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2597802) q[0];
sx q[0];
rz(-1.1008982) q[0];
sx q[0];
rz(-0.63006054) q[0];
rz(-0.13314816) q[1];
sx q[1];
rz(-1.1590191) q[1];
sx q[1];
rz(-1.0551183) q[1];
rz(-1.3721977) q[2];
sx q[2];
rz(-0.19459859) q[2];
sx q[2];
rz(2.9698402) q[2];
rz(-1.6937485) q[3];
sx q[3];
rz(-0.80514859) q[3];
sx q[3];
rz(-2.2283016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
