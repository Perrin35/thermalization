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
rz(2.0734442) q[0];
rz(0.90574342) q[1];
sx q[1];
rz(-1.4248166) q[1];
sx q[1];
rz(2.4832771) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93958116) q[0];
sx q[0];
rz(-2.3774703) q[0];
sx q[0];
rz(2.0234213) q[0];
rz(-pi) q[1];
rz(0.7204804) q[2];
sx q[2];
rz(-2.176365) q[2];
sx q[2];
rz(0.31878135) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.66246819) q[1];
sx q[1];
rz(-2.3196546) q[1];
sx q[1];
rz(2.2884018) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-0.98755223) q[2];
sx q[2];
rz(-1.1215569) q[2];
sx q[2];
rz(1.424632) q[2];
rz(2.3512261) q[3];
sx q[3];
rz(-0.75420403) q[3];
sx q[3];
rz(-0.2592026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49609274) q[0];
sx q[0];
rz(-2.6232965) q[0];
sx q[0];
rz(1.8081283) q[0];
rz(1.2893527) q[1];
sx q[1];
rz(-2.5080296) q[1];
sx q[1];
rz(-3.0234171) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2545276) q[0];
sx q[0];
rz(-0.62572563) q[0];
sx q[0];
rz(1.5805946) q[0];
rz(-1.1907383) q[2];
sx q[2];
rz(-2.5439777) q[2];
sx q[2];
rz(2.3412398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4813586) q[1];
sx q[1];
rz(-2.3412941) q[1];
sx q[1];
rz(-2.1651398) q[1];
x q[2];
rz(-0.34174049) q[3];
sx q[3];
rz(-1.8287418) q[3];
sx q[3];
rz(-2.771857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3893343) q[2];
sx q[2];
rz(-2.3576333) q[2];
sx q[2];
rz(0.83750677) q[2];
rz(-2.5181455) q[3];
sx q[3];
rz(-1.1102755) q[3];
sx q[3];
rz(0.71201223) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74152827) q[0];
sx q[0];
rz(-0.15151227) q[0];
sx q[0];
rz(2.9140299) q[0];
rz(2.172982) q[1];
sx q[1];
rz(-2.249735) q[1];
sx q[1];
rz(1.4898941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8422166) q[0];
sx q[0];
rz(-0.56952667) q[0];
sx q[0];
rz(2.1511908) q[0];
x q[1];
rz(0.030185791) q[2];
sx q[2];
rz(-1.7329114) q[2];
sx q[2];
rz(-0.27366226) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7832853) q[1];
sx q[1];
rz(-0.44898012) q[1];
sx q[1];
rz(-2.994356) q[1];
rz(-1.8996542) q[3];
sx q[3];
rz(-0.76995459) q[3];
sx q[3];
rz(-0.17754031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0861133) q[2];
sx q[2];
rz(-1.5986634) q[2];
sx q[2];
rz(-1.0184658) q[2];
rz(-2.5633519) q[3];
sx q[3];
rz(-2.6522804) q[3];
sx q[3];
rz(-1.1967292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2932435) q[0];
sx q[0];
rz(-1.5986377) q[0];
sx q[0];
rz(-1.6795213) q[0];
rz(2.2734185) q[1];
sx q[1];
rz(-1.236092) q[1];
sx q[1];
rz(2.6624534) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7443607) q[0];
sx q[0];
rz(-1.3102089) q[0];
sx q[0];
rz(0.67407383) q[0];
x q[1];
rz(-2.5136679) q[2];
sx q[2];
rz(-2.6158545) q[2];
sx q[2];
rz(0.0096273011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1066691) q[1];
sx q[1];
rz(-2.0210631) q[1];
sx q[1];
rz(-1.5140956) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3684351) q[3];
sx q[3];
rz(-2.7773603) q[3];
sx q[3];
rz(-2.3868167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9118871) q[2];
sx q[2];
rz(-2.2356326) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01920779) q[0];
sx q[0];
rz(-2.8715219) q[0];
sx q[0];
rz(-0.58575678) q[0];
rz(-2.6981804) q[1];
sx q[1];
rz(-1.5250008) q[1];
sx q[1];
rz(-2.401039) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6218005) q[0];
sx q[0];
rz(-0.92042887) q[0];
sx q[0];
rz(-2.4629794) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1699723) q[2];
sx q[2];
rz(-1.8163774) q[2];
sx q[2];
rz(-1.309091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8466928) q[1];
sx q[1];
rz(-2.7753662) q[1];
sx q[1];
rz(-0.80898714) q[1];
x q[2];
rz(-1.8216474) q[3];
sx q[3];
rz(-1.2870868) q[3];
sx q[3];
rz(1.2267867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1195141) q[2];
sx q[2];
rz(-1.7064648) q[2];
sx q[2];
rz(2.4100927) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5090094) q[0];
sx q[0];
rz(-0.80061299) q[0];
sx q[0];
rz(1.5923694) q[0];
rz(-1.0410615) q[1];
sx q[1];
rz(-0.58716232) q[1];
sx q[1];
rz(-0.36137533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3665112) q[0];
sx q[0];
rz(-1.68987) q[0];
sx q[0];
rz(0.7304169) q[0];
rz(1.560426) q[2];
sx q[2];
rz(-2.8841495) q[2];
sx q[2];
rz(0.68177047) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92897086) q[1];
sx q[1];
rz(-0.55666258) q[1];
sx q[1];
rz(-0.61584872) q[1];
rz(1.8429716) q[3];
sx q[3];
rz(-2.5243773) q[3];
sx q[3];
rz(2.3637852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8063712) q[2];
sx q[2];
rz(-1.0314032) q[2];
sx q[2];
rz(-1.6356989) q[2];
rz(-0.38883543) q[3];
sx q[3];
rz(-2.5917228) q[3];
sx q[3];
rz(2.4839362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036670551) q[0];
sx q[0];
rz(-0.67476434) q[0];
sx q[0];
rz(-0.92055935) q[0];
rz(-2.8432644) q[1];
sx q[1];
rz(-1.0519271) q[1];
sx q[1];
rz(-1.1632318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6434403) q[0];
sx q[0];
rz(-1.5806245) q[0];
sx q[0];
rz(0.32443829) q[0];
rz(-pi) q[1];
rz(2.3506868) q[2];
sx q[2];
rz(-1.8703572) q[2];
sx q[2];
rz(-2.0319446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5015848) q[1];
sx q[1];
rz(-0.54646508) q[1];
sx q[1];
rz(-0.0025005977) q[1];
rz(-pi) q[2];
rz(-1.1441394) q[3];
sx q[3];
rz(-0.32638829) q[3];
sx q[3];
rz(-2.8976664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7982911) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(0.11246559) q[2];
rz(2.8602142) q[3];
sx q[3];
rz(-2.3795542) q[3];
sx q[3];
rz(-1.5587829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.52044582) q[0];
sx q[0];
rz(-0.96108323) q[0];
sx q[0];
rz(0.1113101) q[0];
rz(0.81980199) q[1];
sx q[1];
rz(-0.83378053) q[1];
sx q[1];
rz(0.05096635) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5876041) q[0];
sx q[0];
rz(-1.7394571) q[0];
sx q[0];
rz(1.8722562) q[0];
rz(-pi) q[1];
rz(-1.3459008) q[2];
sx q[2];
rz(-2.1818265) q[2];
sx q[2];
rz(2.261206) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.513235) q[1];
sx q[1];
rz(-1.4396778) q[1];
sx q[1];
rz(0.031946957) q[1];
rz(-pi) q[2];
rz(2.662728) q[3];
sx q[3];
rz(-2.1486946) q[3];
sx q[3];
rz(-0.013805429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.082197949) q[2];
sx q[2];
rz(-1.1092564) q[2];
sx q[2];
rz(1.294403) q[2];
rz(2.1829677) q[3];
sx q[3];
rz(-0.38161033) q[3];
sx q[3];
rz(-2.3310048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078700773) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(-0.94619757) q[0];
rz(-1.2092489) q[1];
sx q[1];
rz(-0.47018662) q[1];
sx q[1];
rz(1.6016003) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7332184) q[0];
sx q[0];
rz(-0.90825496) q[0];
sx q[0];
rz(-2.7761534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8047137) q[2];
sx q[2];
rz(-2.0511464) q[2];
sx q[2];
rz(-0.059305517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6098574) q[1];
sx q[1];
rz(-1.3470728) q[1];
sx q[1];
rz(-2.5705264) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3592992) q[3];
sx q[3];
rz(-1.0254847) q[3];
sx q[3];
rz(-0.88192597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20751791) q[2];
sx q[2];
rz(-1.8002276) q[2];
sx q[2];
rz(-1.6259469) q[2];
rz(-0.80306426) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(-1.9934742) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37907264) q[0];
sx q[0];
rz(-1.9110492) q[0];
sx q[0];
rz(-1.8050964) q[0];
rz(-1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(-1.1755099) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0031457) q[0];
sx q[0];
rz(-1.3295191) q[0];
sx q[0];
rz(-2.3958415) q[0];
rz(2.3606775) q[2];
sx q[2];
rz(-1.7405207) q[2];
sx q[2];
rz(0.11955027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4912632) q[1];
sx q[1];
rz(-1.6230738) q[1];
sx q[1];
rz(0.13904528) q[1];
rz(-0.59713364) q[3];
sx q[3];
rz(-1.8270511) q[3];
sx q[3];
rz(-0.090274752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.056957873) q[2];
sx q[2];
rz(-1.4738169) q[2];
sx q[2];
rz(-1.8353621) q[2];
rz(-0.42825395) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(-0.38124198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.9396176) q[2];
sx q[2];
rz(-2.106728) q[2];
sx q[2];
rz(1.8024924) q[2];
rz(-2.0023575) q[3];
sx q[3];
rz(-0.16466864) q[3];
sx q[3];
rz(2.8684657) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
