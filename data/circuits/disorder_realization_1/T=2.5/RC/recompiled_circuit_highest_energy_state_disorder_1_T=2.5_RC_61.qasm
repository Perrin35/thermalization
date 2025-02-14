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
rz(3.8008939) q[0];
sx q[0];
rz(10.492926) q[0];
rz(0.90574342) q[1];
sx q[1];
rz(-1.4248166) q[1];
sx q[1];
rz(2.4832771) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6092583) q[0];
sx q[0];
rz(-0.89920767) q[0];
sx q[0];
rz(-0.39686578) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82633965) q[2];
sx q[2];
rz(-0.99747083) q[2];
sx q[2];
rz(-2.353015) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5708471) q[1];
sx q[1];
rz(-0.98624228) q[1];
sx q[1];
rz(0.61573124) q[1];
rz(-pi) q[2];
rz(1.0327037) q[3];
sx q[3];
rz(-1.7600087) q[3];
sx q[3];
rz(0.013232649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1540404) q[2];
sx q[2];
rz(-1.1215569) q[2];
sx q[2];
rz(-1.424632) q[2];
rz(-2.3512261) q[3];
sx q[3];
rz(-2.3873886) q[3];
sx q[3];
rz(2.8823901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454999) q[0];
sx q[0];
rz(-0.51829618) q[0];
sx q[0];
rz(-1.8081283) q[0];
rz(-1.2893527) q[1];
sx q[1];
rz(-0.63356304) q[1];
sx q[1];
rz(-3.0234171) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2545276) q[0];
sx q[0];
rz(-2.515867) q[0];
sx q[0];
rz(-1.560998) q[0];
rz(-pi) q[1];
rz(1.1907383) q[2];
sx q[2];
rz(-2.5439777) q[2];
sx q[2];
rz(-2.3412398) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0314583) q[1];
sx q[1];
rz(-2.2074568) q[1];
sx q[1];
rz(-0.52326556) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66690655) q[3];
sx q[3];
rz(-2.7165037) q[3];
sx q[3];
rz(-2.5626884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3893343) q[2];
sx q[2];
rz(-2.3576333) q[2];
sx q[2];
rz(-2.3040859) q[2];
rz(2.5181455) q[3];
sx q[3];
rz(-2.0313171) q[3];
sx q[3];
rz(0.71201223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74152827) q[0];
sx q[0];
rz(-0.15151227) q[0];
sx q[0];
rz(-2.9140299) q[0];
rz(-2.172982) q[1];
sx q[1];
rz(-2.249735) q[1];
sx q[1];
rz(1.6516986) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422166) q[0];
sx q[0];
rz(-2.572066) q[0];
sx q[0];
rz(0.99040188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4086087) q[2];
sx q[2];
rz(-1.5410064) q[2];
sx q[2];
rz(-1.839585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7962793) q[1];
sx q[1];
rz(-1.6345164) q[1];
sx q[1];
rz(-2.696852) q[1];
rz(0.82833008) q[3];
sx q[3];
rz(-1.7975494) q[3];
sx q[3];
rz(1.5080719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0861133) q[2];
sx q[2];
rz(-1.5429292) q[2];
sx q[2];
rz(-1.0184658) q[2];
rz(2.5633519) q[3];
sx q[3];
rz(-2.6522804) q[3];
sx q[3];
rz(1.1967292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8483491) q[0];
sx q[0];
rz(-1.5986377) q[0];
sx q[0];
rz(1.4620713) q[0];
rz(0.86817414) q[1];
sx q[1];
rz(-1.9055007) q[1];
sx q[1];
rz(-0.47913924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2799982) q[0];
sx q[0];
rz(-0.71528212) q[0];
sx q[0];
rz(0.4037373) q[0];
rz(-2.702616) q[2];
sx q[2];
rz(-1.8700646) q[2];
sx q[2];
rz(2.1411095) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1066691) q[1];
sx q[1];
rz(-2.0210631) q[1];
sx q[1];
rz(1.5140956) q[1];
rz(-pi) q[2];
rz(1.3105743) q[3];
sx q[3];
rz(-1.3129915) q[3];
sx q[3];
rz(1.5797404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9118871) q[2];
sx q[2];
rz(-0.90596002) q[2];
sx q[2];
rz(0.68944302) q[2];
rz(0.42129579) q[3];
sx q[3];
rz(-0.7553941) q[3];
sx q[3];
rz(-1.9740483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1223849) q[0];
sx q[0];
rz(-0.27007073) q[0];
sx q[0];
rz(-0.58575678) q[0];
rz(-2.6981804) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(-0.74055368) q[1];
rz(-pi) q[2];
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
rz(-1.9890063) q[2];
sx q[2];
rz(-0.64179342) q[2];
sx q[2];
rz(-0.60371232) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6428004) q[1];
sx q[1];
rz(-1.3086924) q[1];
sx q[1];
rz(0.2587872) q[1];
x q[2];
rz(1.8216474) q[3];
sx q[3];
rz(-1.2870868) q[3];
sx q[3];
rz(1.914806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.022078557) q[2];
sx q[2];
rz(-1.4351279) q[2];
sx q[2];
rz(-2.4100927) q[2];
rz(-2.5445407) q[3];
sx q[3];
rz(-0.67295939) q[3];
sx q[3];
rz(1.2520242) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5090094) q[0];
sx q[0];
rz(-0.80061299) q[0];
sx q[0];
rz(-1.5923694) q[0];
rz(-1.0410615) q[1];
sx q[1];
rz(-0.58716232) q[1];
sx q[1];
rz(2.7802173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7750814) q[0];
sx q[0];
rz(-1.4517227) q[0];
sx q[0];
rz(0.7304169) q[0];
rz(-pi) q[1];
rz(-3.1388624) q[2];
sx q[2];
rz(-1.3133674) q[2];
sx q[2];
rz(-0.67104679) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6238072) q[1];
sx q[1];
rz(-2.0167162) q[1];
sx q[1];
rz(-1.2257027) q[1];
x q[2];
rz(-1.2986211) q[3];
sx q[3];
rz(-2.5243773) q[3];
sx q[3];
rz(-0.77780741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8063712) q[2];
sx q[2];
rz(-1.0314032) q[2];
sx q[2];
rz(1.5058937) q[2];
rz(-0.38883543) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(-2.4839362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(0.036670551) q[0];
sx q[0];
rz(-2.4668283) q[0];
sx q[0];
rz(-0.92055935) q[0];
rz(2.8432644) q[1];
sx q[1];
rz(-2.0896656) q[1];
sx q[1];
rz(1.9783609) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0759493) q[0];
sx q[0];
rz(-1.8952184) q[0];
sx q[0];
rz(-1.5811654) q[0];
rz(-2.3506868) q[2];
sx q[2];
rz(-1.8703572) q[2];
sx q[2];
rz(-1.1096481) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.498658) q[1];
sx q[1];
rz(-1.0243332) q[1];
sx q[1];
rz(-1.5723173) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9974532) q[3];
sx q[3];
rz(-2.8152044) q[3];
sx q[3];
rz(0.24392621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7982911) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(0.11246559) q[2];
rz(2.8602142) q[3];
sx q[3];
rz(-2.3795542) q[3];
sx q[3];
rz(1.5828097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.52044582) q[0];
sx q[0];
rz(-0.96108323) q[0];
sx q[0];
rz(-3.0302826) q[0];
rz(2.3217907) q[1];
sx q[1];
rz(-2.3078121) q[1];
sx q[1];
rz(0.05096635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5539885) q[0];
sx q[0];
rz(-1.7394571) q[0];
sx q[0];
rz(1.8722562) q[0];
rz(-0.62306632) q[2];
sx q[2];
rz(-1.754481) q[2];
sx q[2];
rz(-2.3206835) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8680893) q[1];
sx q[1];
rz(-3.00666) q[1];
sx q[1];
rz(-1.3331624) q[1];
rz(-0.47886465) q[3];
sx q[3];
rz(-0.99289809) q[3];
sx q[3];
rz(-3.1277872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0593947) q[2];
sx q[2];
rz(-1.1092564) q[2];
sx q[2];
rz(1.8471897) q[2];
rz(0.95862499) q[3];
sx q[3];
rz(-2.7599823) q[3];
sx q[3];
rz(0.81058782) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628919) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(0.94619757) q[0];
rz(1.9323438) q[1];
sx q[1];
rz(-2.671406) q[1];
sx q[1];
rz(1.5399923) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.748016) q[0];
sx q[0];
rz(-1.8564176) q[0];
sx q[0];
rz(2.266721) q[0];
rz(-pi) q[1];
rz(-1.8047137) q[2];
sx q[2];
rz(-1.0904462) q[2];
sx q[2];
rz(0.059305517) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8481503) q[1];
sx q[1];
rz(-0.60876011) q[1];
sx q[1];
rz(-0.3984299) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3592992) q[3];
sx q[3];
rz(-2.1161079) q[3];
sx q[3];
rz(-2.2596667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20751791) q[2];
sx q[2];
rz(-1.8002276) q[2];
sx q[2];
rz(1.5156457) q[2];
rz(0.80306426) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(1.9934742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37907264) q[0];
sx q[0];
rz(-1.2305434) q[0];
sx q[0];
rz(1.3364963) q[0];
rz(1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(1.1755099) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7911691) q[0];
sx q[0];
rz(-2.2901112) q[0];
sx q[0];
rz(-1.8940303) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3606775) q[2];
sx q[2];
rz(-1.7405207) q[2];
sx q[2];
rz(-0.11955027) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43688135) q[1];
sx q[1];
rz(-0.14848868) q[1];
sx q[1];
rz(-2.7806033) q[1];
rz(-pi) q[2];
rz(-1.8776348) q[3];
sx q[3];
rz(-2.1458907) q[3];
sx q[3];
rz(-1.3098615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0846348) q[2];
sx q[2];
rz(-1.6677758) q[2];
sx q[2];
rz(1.8353621) q[2];
rz(-2.7133387) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(-2.7603507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8214977) q[0];
sx q[0];
rz(-1.7392673) q[0];
sx q[0];
rz(2.0489954) q[0];
rz(-2.8749657) q[1];
sx q[1];
rz(-2.0482002) q[1];
sx q[1];
rz(2.6788914) q[1];
rz(0.54556432) q[2];
sx q[2];
rz(-2.5013899) q[2];
sx q[2];
rz(-1.987006) q[2];
rz(1.1392351) q[3];
sx q[3];
rz(-0.16466864) q[3];
sx q[3];
rz(2.8684657) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
