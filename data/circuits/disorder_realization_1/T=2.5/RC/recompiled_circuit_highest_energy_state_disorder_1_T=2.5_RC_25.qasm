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
rz(-2.3020307) q[0];
sx q[0];
rz(-0.65930128) q[0];
sx q[0];
rz(-1.0681485) q[0];
rz(0.90574342) q[1];
sx q[1];
rz(-1.4248166) q[1];
sx q[1];
rz(-0.6583156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6092583) q[0];
sx q[0];
rz(-0.89920767) q[0];
sx q[0];
rz(0.39686578) q[0];
rz(-pi) q[1];
rz(0.80950244) q[2];
sx q[2];
rz(-0.90478071) q[2];
sx q[2];
rz(1.8272546) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5707455) q[1];
sx q[1];
rz(-0.98624228) q[1];
sx q[1];
rz(-0.61573124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0327037) q[3];
sx q[3];
rz(-1.7600087) q[3];
sx q[3];
rz(3.12836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1540404) q[2];
sx q[2];
rz(-1.1215569) q[2];
sx q[2];
rz(1.7169607) q[2];
rz(-0.79036653) q[3];
sx q[3];
rz(-2.3873886) q[3];
sx q[3];
rz(-2.8823901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454999) q[0];
sx q[0];
rz(-0.51829618) q[0];
sx q[0];
rz(-1.8081283) q[0];
rz(1.2893527) q[1];
sx q[1];
rz(-0.63356304) q[1];
sx q[1];
rz(3.0234171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2545276) q[0];
sx q[0];
rz(-0.62572563) q[0];
sx q[0];
rz(1.5805946) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9508544) q[2];
sx q[2];
rz(-0.59761492) q[2];
sx q[2];
rz(0.80035287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66023405) q[1];
sx q[1];
rz(-0.80029857) q[1];
sx q[1];
rz(-2.1651398) q[1];
x q[2];
rz(0.66690655) q[3];
sx q[3];
rz(-0.42508891) q[3];
sx q[3];
rz(2.5626884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3893343) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(-2.3040859) q[2];
rz(0.62344712) q[3];
sx q[3];
rz(-2.0313171) q[3];
sx q[3];
rz(-0.71201223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74152827) q[0];
sx q[0];
rz(-0.15151227) q[0];
sx q[0];
rz(-0.22756273) q[0];
rz(2.172982) q[1];
sx q[1];
rz(-2.249735) q[1];
sx q[1];
rz(1.4898941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.299376) q[0];
sx q[0];
rz(-0.56952667) q[0];
sx q[0];
rz(2.1511908) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.030185791) q[2];
sx q[2];
rz(-1.7329114) q[2];
sx q[2];
rz(-2.8679304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19514485) q[1];
sx q[1];
rz(-2.01457) q[1];
sx q[1];
rz(1.5002314) q[1];
rz(-pi) q[2];
rz(-2.3132626) q[3];
sx q[3];
rz(-1.7975494) q[3];
sx q[3];
rz(1.5080719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0554793) q[2];
sx q[2];
rz(-1.5986634) q[2];
sx q[2];
rz(-1.0184658) q[2];
rz(2.5633519) q[3];
sx q[3];
rz(-0.48931229) q[3];
sx q[3];
rz(-1.1967292) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2932435) q[0];
sx q[0];
rz(-1.542955) q[0];
sx q[0];
rz(1.4620713) q[0];
rz(2.2734185) q[1];
sx q[1];
rz(-1.236092) q[1];
sx q[1];
rz(2.6624534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86159443) q[0];
sx q[0];
rz(-0.71528212) q[0];
sx q[0];
rz(-2.7378553) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62792471) q[2];
sx q[2];
rz(-2.6158545) q[2];
sx q[2];
rz(3.1319654) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16461611) q[1];
sx q[1];
rz(-2.6880126) q[1];
sx q[1];
rz(0.11670392) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3684351) q[3];
sx q[3];
rz(-0.36423238) q[3];
sx q[3];
rz(2.3868167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9118871) q[2];
sx q[2];
rz(-0.90596002) q[2];
sx q[2];
rz(0.68944302) q[2];
rz(-2.7202969) q[3];
sx q[3];
rz(-0.7553941) q[3];
sx q[3];
rz(-1.9740483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01920779) q[0];
sx q[0];
rz(-2.8715219) q[0];
sx q[0];
rz(-2.5558359) q[0];
rz(2.6981804) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(-2.401039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6950615) q[0];
sx q[0];
rz(-0.90264812) q[0];
sx q[0];
rz(-0.88094385) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9890063) q[2];
sx q[2];
rz(-0.64179342) q[2];
sx q[2];
rz(-0.60371232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49879229) q[1];
sx q[1];
rz(-1.8329002) q[1];
sx q[1];
rz(-2.8828055) q[1];
rz(-pi) q[2];
rz(-1.3199453) q[3];
sx q[3];
rz(-1.2870868) q[3];
sx q[3];
rz(-1.2267867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.022078557) q[2];
sx q[2];
rz(-1.7064648) q[2];
sx q[2];
rz(-2.4100927) q[2];
rz(0.59705192) q[3];
sx q[3];
rz(-0.67295939) q[3];
sx q[3];
rz(1.2520242) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5090094) q[0];
sx q[0];
rz(-2.3409797) q[0];
sx q[0];
rz(1.5923694) q[0];
rz(1.0410615) q[1];
sx q[1];
rz(-2.5544303) q[1];
sx q[1];
rz(2.7802173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.310285) q[0];
sx q[0];
rz(-0.84670369) q[0];
sx q[0];
rz(-1.4115439) q[0];
rz(-pi) q[1];
rz(-1.560426) q[2];
sx q[2];
rz(-0.25744312) q[2];
sx q[2];
rz(-2.4598222) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2126218) q[1];
sx q[1];
rz(-2.5849301) q[1];
sx q[1];
rz(-2.5257439) q[1];
rz(-pi) q[2];
rz(-2.1704223) q[3];
sx q[3];
rz(-1.4145734) q[3];
sx q[3];
rz(1.0167817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8063712) q[2];
sx q[2];
rz(-2.1101895) q[2];
sx q[2];
rz(1.6356989) q[2];
rz(-0.38883543) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(0.65765643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1049221) q[0];
sx q[0];
rz(-2.4668283) q[0];
sx q[0];
rz(2.2210333) q[0];
rz(-0.29832828) q[1];
sx q[1];
rz(-1.0519271) q[1];
sx q[1];
rz(-1.9783609) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6434403) q[0];
sx q[0];
rz(-1.5609682) q[0];
sx q[0];
rz(0.32443829) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1569491) q[2];
sx q[2];
rz(-2.317642) q[2];
sx q[2];
rz(-2.3904843) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6400078) q[1];
sx q[1];
rz(-0.54646508) q[1];
sx q[1];
rz(-3.1390921) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13917346) q[3];
sx q[3];
rz(-1.8669898) q[3];
sx q[3];
rz(0.20352669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7982911) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(3.0291271) q[2];
rz(-2.8602142) q[3];
sx q[3];
rz(-2.3795542) q[3];
sx q[3];
rz(-1.5828097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.52044582) q[0];
sx q[0];
rz(-2.1805094) q[0];
sx q[0];
rz(3.0302826) q[0];
rz(-2.3217907) q[1];
sx q[1];
rz(-2.3078121) q[1];
sx q[1];
rz(-0.05096635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726377) q[0];
sx q[0];
rz(-1.2737455) q[0];
sx q[0];
rz(0.1764651) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62306632) q[2];
sx q[2];
rz(-1.754481) q[2];
sx q[2];
rz(2.3206835) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6283577) q[1];
sx q[1];
rz(-1.7019148) q[1];
sx q[1];
rz(0.031946957) q[1];
rz(-pi) q[2];
rz(0.9370798) q[3];
sx q[3];
rz(-1.9670319) q[3];
sx q[3];
rz(-1.2806438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0593947) q[2];
sx q[2];
rz(-2.0323362) q[2];
sx q[2];
rz(1.8471897) q[2];
rz(-2.1829677) q[3];
sx q[3];
rz(-2.7599823) q[3];
sx q[3];
rz(0.81058782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0628919) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(-0.94619757) q[0];
rz(1.9323438) q[1];
sx q[1];
rz(-0.47018662) q[1];
sx q[1];
rz(-1.5399923) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14808753) q[0];
sx q[0];
rz(-0.74310857) q[0];
sx q[0];
rz(-2.00032) q[0];
x q[1];
rz(-2.7230324) q[2];
sx q[2];
rz(-0.53024923) q[2];
sx q[2];
rz(0.53539941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9609449) q[1];
sx q[1];
rz(-2.1259318) q[1];
sx q[1];
rz(-1.8349232) q[1];
rz(-pi) q[2];
rz(2.8084635) q[3];
sx q[3];
rz(-2.5605943) q[3];
sx q[3];
rz(2.6521386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20751791) q[2];
sx q[2];
rz(-1.3413651) q[2];
sx q[2];
rz(1.6259469) q[2];
rz(2.3385284) q[3];
sx q[3];
rz(-2.8290437) q[3];
sx q[3];
rz(1.9934742) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.76252) q[0];
sx q[0];
rz(-1.2305434) q[0];
sx q[0];
rz(1.8050964) q[0];
rz(1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(-1.9660827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208082) q[0];
sx q[0];
rz(-2.3649923) q[0];
sx q[0];
rz(-0.34790502) q[0];
rz(2.3606775) q[2];
sx q[2];
rz(-1.7405207) q[2];
sx q[2];
rz(-3.0220424) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.43688135) q[1];
sx q[1];
rz(-0.14848868) q[1];
sx q[1];
rz(0.36098934) q[1];
rz(-pi) q[2];
rz(-0.43607278) q[3];
sx q[3];
rz(-0.64358866) q[3];
sx q[3];
rz(1.3042579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0846348) q[2];
sx q[2];
rz(-1.4738169) q[2];
sx q[2];
rz(-1.3062306) q[2];
rz(0.42825395) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(0.38124198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.320095) q[0];
sx q[0];
rz(-1.7392673) q[0];
sx q[0];
rz(2.0489954) q[0];
rz(-0.26662695) q[1];
sx q[1];
rz(-1.0933924) q[1];
sx q[1];
rz(-0.46270121) q[1];
rz(2.5746018) q[2];
sx q[2];
rz(-1.2556354) q[2];
sx q[2];
rz(-3.1047594) q[2];
rz(1.7206031) q[3];
sx q[3];
rz(-1.6394184) q[3];
sx q[3];
rz(1.7240789) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
