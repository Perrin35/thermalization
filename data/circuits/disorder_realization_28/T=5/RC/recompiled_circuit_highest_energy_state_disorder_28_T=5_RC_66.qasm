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
rz(-2.2313843) q[0];
sx q[0];
rz(-2.3910523) q[0];
sx q[0];
rz(2.5810177) q[0];
rz(-1.203546) q[1];
sx q[1];
rz(-0.37625852) q[1];
sx q[1];
rz(-2.9484152) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5420509) q[0];
sx q[0];
rz(-2.1452869) q[0];
sx q[0];
rz(2.1035139) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4273663) q[2];
sx q[2];
rz(-0.71038112) q[2];
sx q[2];
rz(-1.8153035) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3456125) q[1];
sx q[1];
rz(-0.66540816) q[1];
sx q[1];
rz(-1.3434975) q[1];
rz(-pi) q[2];
rz(1.3985996) q[3];
sx q[3];
rz(-1.2463257) q[3];
sx q[3];
rz(-0.29968047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2304077) q[2];
sx q[2];
rz(-0.62778968) q[2];
sx q[2];
rz(-0.5578624) q[2];
rz(-1.9134391) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(3.0465928) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8519583) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(-1.5767545) q[0];
rz(1.1467038) q[1];
sx q[1];
rz(-1.4253989) q[1];
sx q[1];
rz(-2.1099405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1250336) q[0];
sx q[0];
rz(-1.8047338) q[0];
sx q[0];
rz(2.1159322) q[0];
rz(-pi) q[1];
rz(-1.8081153) q[2];
sx q[2];
rz(-1.9734335) q[2];
sx q[2];
rz(1.7466842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.116392) q[1];
sx q[1];
rz(-0.78577367) q[1];
sx q[1];
rz(-1.0791732) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7835763) q[3];
sx q[3];
rz(-2.6081134) q[3];
sx q[3];
rz(1.250035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8138294) q[2];
sx q[2];
rz(-0.29080614) q[2];
sx q[2];
rz(-1.7247058) q[2];
rz(-0.82301569) q[3];
sx q[3];
rz(-1.6133285) q[3];
sx q[3];
rz(2.4970162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7773892) q[0];
sx q[0];
rz(-1.1193898) q[0];
sx q[0];
rz(0.35225824) q[0];
rz(0.38905713) q[1];
sx q[1];
rz(-1.16951) q[1];
sx q[1];
rz(-2.9152117) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1539696) q[0];
sx q[0];
rz(-0.26096086) q[0];
sx q[0];
rz(1.8960192) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1445072) q[2];
sx q[2];
rz(-1.2752394) q[2];
sx q[2];
rz(-1.3546582) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4714053) q[1];
sx q[1];
rz(-0.67882631) q[1];
sx q[1];
rz(0.77969867) q[1];
rz(-pi) q[2];
x q[2];
rz(2.705615) q[3];
sx q[3];
rz(-1.5020026) q[3];
sx q[3];
rz(-1.1747557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9809197) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(0.36160198) q[2];
rz(2.8412039) q[3];
sx q[3];
rz(-1.3114248) q[3];
sx q[3];
rz(0.96295199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2196197) q[0];
sx q[0];
rz(-1.0924871) q[0];
sx q[0];
rz(-3.0203982) q[0];
rz(1.1990064) q[1];
sx q[1];
rz(-2.0620748) q[1];
sx q[1];
rz(-1.8669063) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5124487) q[0];
sx q[0];
rz(-1.7928018) q[0];
sx q[0];
rz(-1.6187526) q[0];
rz(-pi) q[1];
rz(0.2054487) q[2];
sx q[2];
rz(-1.5866071) q[2];
sx q[2];
rz(0.39218536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9642942) q[1];
sx q[1];
rz(-1.7907447) q[1];
sx q[1];
rz(0.58348685) q[1];
rz(1.2043456) q[3];
sx q[3];
rz(-0.38903686) q[3];
sx q[3];
rz(1.4938482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.092209665) q[2];
sx q[2];
rz(-1.4468687) q[2];
sx q[2];
rz(-1.3775187) q[2];
rz(-2.0130017) q[3];
sx q[3];
rz(-2.1994574) q[3];
sx q[3];
rz(0.82730627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89151299) q[0];
sx q[0];
rz(-2.2703607) q[0];
sx q[0];
rz(-2.0569892) q[0];
rz(0.67101038) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(-0.074140851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0947572) q[0];
sx q[0];
rz(-2.1291385) q[0];
sx q[0];
rz(-1.895491) q[0];
x q[1];
rz(0.83168516) q[2];
sx q[2];
rz(-2.7234969) q[2];
sx q[2];
rz(1.8036133) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74846327) q[1];
sx q[1];
rz(-0.65915758) q[1];
sx q[1];
rz(-2.0214969) q[1];
rz(2.0228593) q[3];
sx q[3];
rz(-0.30395711) q[3];
sx q[3];
rz(3.0641276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.025617754) q[2];
sx q[2];
rz(-2.3927549) q[2];
sx q[2];
rz(-0.96308127) q[2];
rz(-1.3747619) q[3];
sx q[3];
rz(-2.0120967) q[3];
sx q[3];
rz(-2.6095384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(1.9328203) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(1.8836841) q[0];
rz(-2.3014297) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(-2.449583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0449311) q[0];
sx q[0];
rz(-2.7122444) q[0];
sx q[0];
rz(-2.2675687) q[0];
rz(0.18628405) q[2];
sx q[2];
rz(-0.77464235) q[2];
sx q[2];
rz(1.8587405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.55992014) q[1];
sx q[1];
rz(-2.2894303) q[1];
sx q[1];
rz(2.2581359) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63661544) q[3];
sx q[3];
rz(-2.0354384) q[3];
sx q[3];
rz(2.1615504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.96364) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(0.5160416) q[2];
rz(1.3456723) q[3];
sx q[3];
rz(-2.7005152) q[3];
sx q[3];
rz(2.1399982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26821414) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(-1.0299261) q[0];
rz(-0.65451199) q[1];
sx q[1];
rz(-1.3348568) q[1];
sx q[1];
rz(-0.17394224) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2692341) q[0];
sx q[0];
rz(-0.96251153) q[0];
sx q[0];
rz(0.66156265) q[0];
x q[1];
rz(2.6268466) q[2];
sx q[2];
rz(-0.87640773) q[2];
sx q[2];
rz(1.0389164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0566114) q[1];
sx q[1];
rz(-2.1826996) q[1];
sx q[1];
rz(-2.8065744) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44486041) q[3];
sx q[3];
rz(-0.56021094) q[3];
sx q[3];
rz(-1.0154533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8339463) q[2];
sx q[2];
rz(-2.6070194) q[2];
sx q[2];
rz(1.8024811) q[2];
rz(-2.3329959) q[3];
sx q[3];
rz(-1.1491038) q[3];
sx q[3];
rz(-1.2861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9209552) q[0];
sx q[0];
rz(-2.0957102) q[0];
sx q[0];
rz(1.9761696) q[0];
rz(-0.19365817) q[1];
sx q[1];
rz(-2.3083189) q[1];
sx q[1];
rz(-1.1509034) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347558) q[0];
sx q[0];
rz(-0.92972219) q[0];
sx q[0];
rz(1.2527554) q[0];
rz(-pi) q[1];
rz(-0.80457689) q[2];
sx q[2];
rz(-1.6940306) q[2];
sx q[2];
rz(-0.096913902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3079466) q[1];
sx q[1];
rz(-2.4599791) q[1];
sx q[1];
rz(0.8591842) q[1];
rz(-pi) q[2];
rz(0.70827534) q[3];
sx q[3];
rz(-2.357956) q[3];
sx q[3];
rz(-0.9140789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9264441) q[2];
sx q[2];
rz(-2.1200659) q[2];
sx q[2];
rz(1.273217) q[2];
rz(-2.6269954) q[3];
sx q[3];
rz(-2.8223346) q[3];
sx q[3];
rz(-2.4014421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1505245) q[0];
sx q[0];
rz(-0.06572289) q[0];
sx q[0];
rz(-2.7142628) q[0];
rz(1.7006251) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(2.2429121) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60481614) q[0];
sx q[0];
rz(-1.1425848) q[0];
sx q[0];
rz(0.34128071) q[0];
x q[1];
rz(2.1601474) q[2];
sx q[2];
rz(-1.9716096) q[2];
sx q[2];
rz(-0.0055731853) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8898892) q[1];
sx q[1];
rz(-2.8475175) q[1];
sx q[1];
rz(-2.2132323) q[1];
rz(-pi) q[2];
rz(-0.69435223) q[3];
sx q[3];
rz(-2.4789841) q[3];
sx q[3];
rz(2.724805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.906189) q[2];
sx q[2];
rz(-2.0960505) q[2];
sx q[2];
rz(-1.7380627) q[2];
rz(-0.67052001) q[3];
sx q[3];
rz(-1.5961921) q[3];
sx q[3];
rz(1.8286573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32187605) q[0];
sx q[0];
rz(-2.1724367) q[0];
sx q[0];
rz(-0.72625351) q[0];
rz(0.67604524) q[1];
sx q[1];
rz(-1.36422) q[1];
sx q[1];
rz(0.77686754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5345316) q[0];
sx q[0];
rz(-2.2103469) q[0];
sx q[0];
rz(0.056900815) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9954722) q[2];
sx q[2];
rz(-1.5263128) q[2];
sx q[2];
rz(-1.5715949) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8269413) q[1];
sx q[1];
rz(-2.6562641) q[1];
sx q[1];
rz(-0.86599533) q[1];
x q[2];
rz(-0.59263521) q[3];
sx q[3];
rz(-1.0909972) q[3];
sx q[3];
rz(2.0247839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1399416) q[2];
sx q[2];
rz(-2.9652014) q[2];
sx q[2];
rz(1.4886935) q[2];
rz(2.0148924) q[3];
sx q[3];
rz(-1.1688346) q[3];
sx q[3];
rz(0.53762976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5269932) q[0];
sx q[0];
rz(-0.64170964) q[0];
sx q[0];
rz(1.5668305) q[0];
rz(-2.3552786) q[1];
sx q[1];
rz(-2.96824) q[1];
sx q[1];
rz(2.7460964) q[1];
rz(2.4509571) q[2];
sx q[2];
rz(-0.80123676) q[2];
sx q[2];
rz(-1.5581808) q[2];
rz(0.24675225) q[3];
sx q[3];
rz(-0.79083957) q[3];
sx q[3];
rz(-1.4141465) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
