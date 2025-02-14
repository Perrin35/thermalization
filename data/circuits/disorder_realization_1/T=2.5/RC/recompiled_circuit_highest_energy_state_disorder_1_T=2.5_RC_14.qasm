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
rz(0.83956194) q[0];
sx q[0];
rz(-2.4822914) q[0];
sx q[0];
rz(1.0681485) q[0];
rz(0.90574342) q[1];
sx q[1];
rz(-1.4248166) q[1];
sx q[1];
rz(-0.6583156) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8480331) q[0];
sx q[0];
rz(-1.8782037) q[0];
sx q[0];
rz(0.859476) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4211123) q[2];
sx q[2];
rz(-2.176365) q[2];
sx q[2];
rz(-2.8228113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66246819) q[1];
sx q[1];
rz(-0.82193804) q[1];
sx q[1];
rz(2.2884018) q[1];
rz(-pi) q[2];
rz(-1.9283965) q[3];
sx q[3];
rz(-2.5743067) q[3];
sx q[3];
rz(1.889313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1540404) q[2];
sx q[2];
rz(-1.1215569) q[2];
sx q[2];
rz(1.424632) q[2];
rz(0.79036653) q[3];
sx q[3];
rz(-2.3873886) q[3];
sx q[3];
rz(-0.2592026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454999) q[0];
sx q[0];
rz(-2.6232965) q[0];
sx q[0];
rz(-1.3334644) q[0];
rz(1.2893527) q[1];
sx q[1];
rz(-2.5080296) q[1];
sx q[1];
rz(0.11817558) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3242107) q[0];
sx q[0];
rz(-1.576535) q[0];
sx q[0];
rz(-2.1964992) q[0];
x q[1];
rz(1.0071272) q[2];
sx q[2];
rz(-1.360513) q[2];
sx q[2];
rz(-2.0521833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0314583) q[1];
sx q[1];
rz(-2.2074568) q[1];
sx q[1];
rz(-0.52326556) q[1];
rz(-pi) q[2];
rz(-2.7998522) q[3];
sx q[3];
rz(-1.3128508) q[3];
sx q[3];
rz(0.3697357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3893343) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(0.83750677) q[2];
rz(0.62344712) q[3];
sx q[3];
rz(-2.0313171) q[3];
sx q[3];
rz(-0.71201223) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74152827) q[0];
sx q[0];
rz(-2.9900804) q[0];
sx q[0];
rz(2.9140299) q[0];
rz(0.96861068) q[1];
sx q[1];
rz(-0.89185762) q[1];
sx q[1];
rz(-1.6516986) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1806396) q[0];
sx q[0];
rz(-1.102987) q[0];
sx q[0];
rz(-0.33766467) q[0];
x q[1];
rz(1.4086087) q[2];
sx q[2];
rz(-1.5410064) q[2];
sx q[2];
rz(-1.3020077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3453133) q[1];
sx q[1];
rz(-1.6345164) q[1];
sx q[1];
rz(0.44474067) q[1];
x q[2];
rz(-0.82833008) q[3];
sx q[3];
rz(-1.3440432) q[3];
sx q[3];
rz(-1.6335207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0554793) q[2];
sx q[2];
rz(-1.5429292) q[2];
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
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8483491) q[0];
sx q[0];
rz(-1.5986377) q[0];
sx q[0];
rz(-1.4620713) q[0];
rz(-2.2734185) q[1];
sx q[1];
rz(-1.9055007) q[1];
sx q[1];
rz(2.6624534) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3972319) q[0];
sx q[0];
rz(-1.3102089) q[0];
sx q[0];
rz(-0.67407383) q[0];
x q[1];
rz(-0.43897668) q[2];
sx q[2];
rz(-1.2715281) q[2];
sx q[2];
rz(2.1411095) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.034923542) q[1];
sx q[1];
rz(-2.0210631) q[1];
sx q[1];
rz(1.627497) q[1];
rz(-pi) q[2];
rz(-0.26637443) q[3];
sx q[3];
rz(-1.3193635) q[3];
sx q[3];
rz(3.0648674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9118871) q[2];
sx q[2];
rz(-0.90596002) q[2];
sx q[2];
rz(2.4521496) q[2];
rz(2.7202969) q[3];
sx q[3];
rz(-0.7553941) q[3];
sx q[3];
rz(-1.1675444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1223849) q[0];
sx q[0];
rz(-2.8715219) q[0];
sx q[0];
rz(0.58575678) q[0];
rz(-0.44341227) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(0.74055368) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5447841) q[0];
sx q[0];
rz(-2.0939026) q[0];
sx q[0];
rz(0.79687065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1699723) q[2];
sx q[2];
rz(-1.3252153) q[2];
sx q[2];
rz(-1.309091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1380746) q[1];
sx q[1];
rz(-1.3210396) q[1];
sx q[1];
rz(-1.8415007) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70526975) q[3];
sx q[3];
rz(-2.765145) q[3];
sx q[3];
rz(-1.173526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.022078557) q[2];
sx q[2];
rz(-1.7064648) q[2];
sx q[2];
rz(-0.73149991) q[2];
rz(-0.59705192) q[3];
sx q[3];
rz(-0.67295939) q[3];
sx q[3];
rz(1.8895684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6325833) q[0];
sx q[0];
rz(-0.80061299) q[0];
sx q[0];
rz(-1.5923694) q[0];
rz(-2.1005311) q[1];
sx q[1];
rz(-2.5544303) q[1];
sx q[1];
rz(2.7802173) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7750814) q[0];
sx q[0];
rz(-1.4517227) q[0];
sx q[0];
rz(0.7304169) q[0];
rz(3.1388624) q[2];
sx q[2];
rz(-1.3133674) q[2];
sx q[2];
rz(0.67104679) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2126218) q[1];
sx q[1];
rz(-0.55666258) q[1];
sx q[1];
rz(2.5257439) q[1];
rz(-pi) q[2];
rz(-0.97117037) q[3];
sx q[3];
rz(-1.7270192) q[3];
sx q[3];
rz(1.0167817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.3352215) q[2];
sx q[2];
rz(-2.1101895) q[2];
sx q[2];
rz(-1.5058937) q[2];
rz(-0.38883543) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(0.65765643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1049221) q[0];
sx q[0];
rz(-0.67476434) q[0];
sx q[0];
rz(2.2210333) q[0];
rz(2.8432644) q[1];
sx q[1];
rz(-2.0896656) q[1];
sx q[1];
rz(-1.1632318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0434303) q[0];
sx q[0];
rz(-0.32458186) q[0];
sx q[0];
rz(3.1107706) q[0];
x q[1];
rz(-0.40980659) q[2];
sx q[2];
rz(-0.834045) q[2];
sx q[2];
rz(-2.9644239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6400078) q[1];
sx q[1];
rz(-0.54646508) q[1];
sx q[1];
rz(3.1390921) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2718795) q[3];
sx q[3];
rz(-1.7038725) q[3];
sx q[3];
rz(1.7334593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3433015) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(3.0291271) q[2];
rz(0.28137842) q[3];
sx q[3];
rz(-0.76203841) q[3];
sx q[3];
rz(-1.5587829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52044582) q[0];
sx q[0];
rz(-0.96108323) q[0];
sx q[0];
rz(-3.0302826) q[0];
rz(-2.3217907) q[1];
sx q[1];
rz(-2.3078121) q[1];
sx q[1];
rz(-0.05096635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5539885) q[0];
sx q[0];
rz(-1.4021356) q[0];
sx q[0];
rz(1.8722562) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30822513) q[2];
sx q[2];
rz(-0.64611084) q[2];
sx q[2];
rz(-2.6406312) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.061739616) q[1];
sx q[1];
rz(-1.5391237) q[1];
sx q[1];
rz(1.701981) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2045129) q[3];
sx q[3];
rz(-1.1745608) q[3];
sx q[3];
rz(-1.8609488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0593947) q[2];
sx q[2];
rz(-2.0323362) q[2];
sx q[2];
rz(-1.294403) q[2];
rz(2.1829677) q[3];
sx q[3];
rz(-2.7599823) q[3];
sx q[3];
rz(2.3310048) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078700773) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(2.1953951) q[0];
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
rz(1.748016) q[0];
sx q[0];
rz(-1.8564176) q[0];
sx q[0];
rz(-0.87487166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6498383) q[2];
sx q[2];
rz(-1.3637614) q[2];
sx q[2];
rz(-1.4018261) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.18064776) q[1];
sx q[1];
rz(-1.0156609) q[1];
sx q[1];
rz(1.8349232) q[1];
x q[2];
rz(-0.55535613) q[3];
sx q[3];
rz(-1.7512519) q[3];
sx q[3];
rz(0.79977913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9340747) q[2];
sx q[2];
rz(-1.8002276) q[2];
sx q[2];
rz(1.6259469) q[2];
rz(2.3385284) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(-1.9934742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.76252) q[0];
sx q[0];
rz(-1.9110492) q[0];
sx q[0];
rz(1.3364963) q[0];
rz(-1.9876009) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(-1.9660827) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82078445) q[0];
sx q[0];
rz(-0.77660038) q[0];
sx q[0];
rz(2.7936876) q[0];
x q[1];
rz(2.9027862) q[2];
sx q[2];
rz(-2.3462834) q[2];
sx q[2];
rz(1.8591566) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4912632) q[1];
sx q[1];
rz(-1.5185188) q[1];
sx q[1];
rz(3.0025474) q[1];
rz(-pi) q[2];
rz(2.7055199) q[3];
sx q[3];
rz(-0.64358866) q[3];
sx q[3];
rz(-1.8373348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0846348) q[2];
sx q[2];
rz(-1.4738169) q[2];
sx q[2];
rz(1.8353621) q[2];
rz(-2.7133387) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(0.38124198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.320095) q[0];
sx q[0];
rz(-1.7392673) q[0];
sx q[0];
rz(2.0489954) q[0];
rz(0.26662695) q[1];
sx q[1];
rz(-2.0482002) q[1];
sx q[1];
rz(2.6788914) q[1];
rz(0.54556432) q[2];
sx q[2];
rz(-2.5013899) q[2];
sx q[2];
rz(-1.987006) q[2];
rz(3.0721957) q[3];
sx q[3];
rz(-1.7202478) q[3];
sx q[3];
rz(-2.9779609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
