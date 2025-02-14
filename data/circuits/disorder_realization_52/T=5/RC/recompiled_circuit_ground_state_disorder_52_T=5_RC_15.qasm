OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0527394) q[0];
sx q[0];
rz(-2.7272447) q[0];
sx q[0];
rz(-2.156884) q[0];
rz(-0.86496487) q[1];
sx q[1];
rz(-0.84671658) q[1];
sx q[1];
rz(-2.4196978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74367172) q[0];
sx q[0];
rz(-1.1211485) q[0];
sx q[0];
rz(-2.0364773) q[0];
rz(-pi) q[1];
rz(1.9198011) q[2];
sx q[2];
rz(-0.86599444) q[2];
sx q[2];
rz(0.26217957) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2014247) q[1];
sx q[1];
rz(-0.076751953) q[1];
sx q[1];
rz(1.1809262) q[1];
rz(-pi) q[2];
rz(-0.15318449) q[3];
sx q[3];
rz(-1.1878052) q[3];
sx q[3];
rz(-2.9944978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7559173) q[2];
sx q[2];
rz(-2.6188681) q[2];
sx q[2];
rz(-0.44504607) q[2];
rz(-1.2612777) q[3];
sx q[3];
rz(-1.6877561) q[3];
sx q[3];
rz(-1.5104347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181013) q[0];
sx q[0];
rz(-2.0913251) q[0];
sx q[0];
rz(2.4871248) q[0];
rz(2.0596313) q[1];
sx q[1];
rz(-1.315821) q[1];
sx q[1];
rz(1.4644324) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.678955) q[0];
sx q[0];
rz(-1.5206095) q[0];
sx q[0];
rz(2.0431014) q[0];
rz(-pi) q[1];
rz(-2.1747443) q[2];
sx q[2];
rz(-1.6308846) q[2];
sx q[2];
rz(-2.0397358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6738715) q[1];
sx q[1];
rz(-0.74127156) q[1];
sx q[1];
rz(-0.060972496) q[1];
rz(2.5188279) q[3];
sx q[3];
rz(-0.55985427) q[3];
sx q[3];
rz(0.44677681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4945041) q[2];
sx q[2];
rz(-3.0454128) q[2];
sx q[2];
rz(0.75867009) q[2];
rz(-1.8132973) q[3];
sx q[3];
rz(-2.0601065) q[3];
sx q[3];
rz(-1.7648511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3674304) q[0];
sx q[0];
rz(-1.4752911) q[0];
sx q[0];
rz(0.037516315) q[0];
rz(1.0388177) q[1];
sx q[1];
rz(-2.807834) q[1];
sx q[1];
rz(2.2822101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276703) q[0];
sx q[0];
rz(-1.3827208) q[0];
sx q[0];
rz(-2.5187224) q[0];
rz(-3.0065303) q[2];
sx q[2];
rz(-2.3762868) q[2];
sx q[2];
rz(1.0926334) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1990314) q[1];
sx q[1];
rz(-1.9936526) q[1];
sx q[1];
rz(-2.7884952) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5435197) q[3];
sx q[3];
rz(-1.9864286) q[3];
sx q[3];
rz(1.9126836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77217707) q[2];
sx q[2];
rz(-2.206216) q[2];
sx q[2];
rz(-2.3773362) q[2];
rz(-2.1956826) q[3];
sx q[3];
rz(-1.9352244) q[3];
sx q[3];
rz(-0.8980155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1979444) q[0];
sx q[0];
rz(-2.2479489) q[0];
sx q[0];
rz(1.8290895) q[0];
rz(-0.30458826) q[1];
sx q[1];
rz(-2.1423788) q[1];
sx q[1];
rz(2.8056858) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48913025) q[0];
sx q[0];
rz(-1.5565598) q[0];
sx q[0];
rz(1.4984309) q[0];
rz(0.44563771) q[2];
sx q[2];
rz(-1.012371) q[2];
sx q[2];
rz(1.4134917) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2780032) q[1];
sx q[1];
rz(-1.9801323) q[1];
sx q[1];
rz(2.5990942) q[1];
x q[2];
rz(-2.920067) q[3];
sx q[3];
rz(-1.5153839) q[3];
sx q[3];
rz(1.6694091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37561068) q[2];
sx q[2];
rz(-2.7178662) q[2];
sx q[2];
rz(1.83164) q[2];
rz(1.8464108) q[3];
sx q[3];
rz(-2.3056307) q[3];
sx q[3];
rz(2.0415993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.497371) q[0];
sx q[0];
rz(-2.8692684) q[0];
sx q[0];
rz(2.3261133) q[0];
rz(0.71715912) q[1];
sx q[1];
rz(-1.6970789) q[1];
sx q[1];
rz(-1.1517634) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718147) q[0];
sx q[0];
rz(-1.0564532) q[0];
sx q[0];
rz(-1.3199947) q[0];
rz(-pi) q[1];
rz(-2.8377811) q[2];
sx q[2];
rz(-1.3388435) q[2];
sx q[2];
rz(1.0396921) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5260701) q[1];
sx q[1];
rz(-2.263219) q[1];
sx q[1];
rz(-2.7500765) q[1];
rz(-pi) q[2];
rz(0.11687704) q[3];
sx q[3];
rz(-0.3038376) q[3];
sx q[3];
rz(-0.60910329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.069245) q[2];
sx q[2];
rz(-2.3657511) q[2];
sx q[2];
rz(-1.5838712) q[2];
rz(-0.97258687) q[3];
sx q[3];
rz(-1.2310622) q[3];
sx q[3];
rz(1.7044273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4438181) q[0];
sx q[0];
rz(-2.8688718) q[0];
sx q[0];
rz(2.4716603) q[0];
rz(2.2386235) q[1];
sx q[1];
rz(-2.4882856) q[1];
sx q[1];
rz(3.0526551) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020821559) q[0];
sx q[0];
rz(-2.1820314) q[0];
sx q[0];
rz(-2.8957248) q[0];
rz(-2.495359) q[2];
sx q[2];
rz(-1.3696732) q[2];
sx q[2];
rz(1.0758019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6429421) q[1];
sx q[1];
rz(-1.6410259) q[1];
sx q[1];
rz(1.8204017) q[1];
x q[2];
rz(0.1286147) q[3];
sx q[3];
rz(-2.0484784) q[3];
sx q[3];
rz(-1.1123808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4778121) q[2];
sx q[2];
rz(-2.0974443) q[2];
sx q[2];
rz(0.021154724) q[2];
rz(0.92709213) q[3];
sx q[3];
rz(-1.6095716) q[3];
sx q[3];
rz(-2.482614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-3.0878736) q[0];
sx q[0];
rz(-1.7807732) q[0];
sx q[0];
rz(1.1370283) q[0];
rz(-2.673705) q[1];
sx q[1];
rz(-1.7158022) q[1];
sx q[1];
rz(-0.55380026) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0040941) q[0];
sx q[0];
rz(-1.8589673) q[0];
sx q[0];
rz(2.1779446) q[0];
rz(-pi) q[1];
rz(-2.2442859) q[2];
sx q[2];
rz(-2.6061432) q[2];
sx q[2];
rz(-0.4430807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97275954) q[1];
sx q[1];
rz(-1.0565896) q[1];
sx q[1];
rz(2.7943816) q[1];
x q[2];
rz(-1.6662237) q[3];
sx q[3];
rz(-0.57668873) q[3];
sx q[3];
rz(2.5043023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0566473) q[2];
sx q[2];
rz(-2.2979996) q[2];
sx q[2];
rz(-2.4455369) q[2];
rz(2.7737235) q[3];
sx q[3];
rz(-1.1080247) q[3];
sx q[3];
rz(-2.8554816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2680161) q[0];
sx q[0];
rz(-1.3883256) q[0];
sx q[0];
rz(2.3263113) q[0];
rz(2.6878327) q[1];
sx q[1];
rz(-2.2397857) q[1];
sx q[1];
rz(0.59660965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94039791) q[0];
sx q[0];
rz(-1.0642645) q[0];
sx q[0];
rz(-0.80123637) q[0];
rz(-pi) q[1];
rz(-1.6758133) q[2];
sx q[2];
rz(-0.74286443) q[2];
sx q[2];
rz(0.32203963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4547465) q[1];
sx q[1];
rz(-1.0300714) q[1];
sx q[1];
rz(-0.35130934) q[1];
rz(-2.7030029) q[3];
sx q[3];
rz(-2.1206754) q[3];
sx q[3];
rz(-0.97766961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10670184) q[2];
sx q[2];
rz(-1.4406349) q[2];
sx q[2];
rz(1.2197257) q[2];
rz(1.8969511) q[3];
sx q[3];
rz(-1.5361667) q[3];
sx q[3];
rz(2.9419148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91584665) q[0];
sx q[0];
rz(-0.93996489) q[0];
sx q[0];
rz(1.655727) q[0];
rz(2.0303717) q[1];
sx q[1];
rz(-1.8467555) q[1];
sx q[1];
rz(1.750754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6538453) q[0];
sx q[0];
rz(-1.5731205) q[0];
sx q[0];
rz(2.8043158) q[0];
rz(1.6860854) q[2];
sx q[2];
rz(-1.2284245) q[2];
sx q[2];
rz(1.1631249) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2272033) q[1];
sx q[1];
rz(-1.7787284) q[1];
sx q[1];
rz(-0.035391594) q[1];
x q[2];
rz(2.7928915) q[3];
sx q[3];
rz(-1.9174272) q[3];
sx q[3];
rz(1.5894401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2755462) q[2];
sx q[2];
rz(-1.6547497) q[2];
sx q[2];
rz(-0.1869959) q[2];
rz(-0.20728076) q[3];
sx q[3];
rz(-0.63396251) q[3];
sx q[3];
rz(1.0435587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2279376) q[0];
sx q[0];
rz(-0.74988237) q[0];
sx q[0];
rz(0.0023181152) q[0];
rz(-1.2222611) q[1];
sx q[1];
rz(-1.5536676) q[1];
sx q[1];
rz(-1.4171756) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.173439) q[0];
sx q[0];
rz(-1.7002956) q[0];
sx q[0];
rz(-2.9347675) q[0];
rz(-1.8882636) q[2];
sx q[2];
rz(-1.0523044) q[2];
sx q[2];
rz(2.7310892) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0884888) q[1];
sx q[1];
rz(-1.7758867) q[1];
sx q[1];
rz(0.049909485) q[1];
rz(-pi) q[2];
rz(2.1571092) q[3];
sx q[3];
rz(-1.7147736) q[3];
sx q[3];
rz(1.7881623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8861387) q[2];
sx q[2];
rz(-1.6243434) q[2];
sx q[2];
rz(-0.21190602) q[2];
rz(-2.741277) q[3];
sx q[3];
rz(-0.90462697) q[3];
sx q[3];
rz(-1.8109842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0935852) q[0];
sx q[0];
rz(-1.891991) q[0];
sx q[0];
rz(1.4461507) q[0];
rz(-2.6493337) q[1];
sx q[1];
rz(-1.6620363) q[1];
sx q[1];
rz(-1.5580039) q[1];
rz(1.3561495) q[2];
sx q[2];
rz(-1.3673906) q[2];
sx q[2];
rz(-2.2833952) q[2];
rz(2.0782804) q[3];
sx q[3];
rz(-1.1125269) q[3];
sx q[3];
rz(-0.0081882523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
