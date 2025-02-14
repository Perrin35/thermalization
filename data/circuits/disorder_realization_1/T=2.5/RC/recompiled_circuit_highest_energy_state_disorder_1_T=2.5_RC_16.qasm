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
rz(-0.6583156) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6092583) q[0];
sx q[0];
rz(-2.242385) q[0];
sx q[0];
rz(2.7447269) q[0];
rz(-pi) q[1];
rz(-2.3320902) q[2];
sx q[2];
rz(-0.90478071) q[2];
sx q[2];
rz(1.8272546) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66246819) q[1];
sx q[1];
rz(-2.3196546) q[1];
sx q[1];
rz(0.85319086) q[1];
rz(-pi) q[2];
rz(-2.9221651) q[3];
sx q[3];
rz(-1.0433222) q[3];
sx q[3];
rz(-1.4722401) q[3];
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
rz(2.3512261) q[3];
sx q[3];
rz(-0.75420403) q[3];
sx q[3];
rz(2.8823901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49609274) q[0];
sx q[0];
rz(-0.51829618) q[0];
sx q[0];
rz(-1.3334644) q[0];
rz(1.2893527) q[1];
sx q[1];
rz(-2.5080296) q[1];
sx q[1];
rz(0.11817558) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8991535) q[0];
sx q[0];
rz(-0.94510539) q[0];
sx q[0];
rz(3.1345128) q[0];
rz(2.8942621) q[2];
sx q[2];
rz(-1.0209903) q[2];
sx q[2];
rz(-2.791385) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3505248) q[1];
sx q[1];
rz(-1.1573029) q[1];
sx q[1];
rz(0.86422635) q[1];
rz(2.7998522) q[3];
sx q[3];
rz(-1.3128508) q[3];
sx q[3];
rz(-0.3697357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3893343) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(-2.3040859) q[2];
rz(-0.62344712) q[3];
sx q[3];
rz(-2.0313171) q[3];
sx q[3];
rz(-2.4295804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.299376) q[0];
sx q[0];
rz(-0.56952667) q[0];
sx q[0];
rz(2.1511908) q[0];
x q[1];
rz(-3.1114069) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9464478) q[1];
sx q[1];
rz(-1.1270226) q[1];
sx q[1];
rz(1.5002314) q[1];
rz(-pi) q[2];
rz(-1.8996542) q[3];
sx q[3];
rz(-2.3716381) q[3];
sx q[3];
rz(0.17754031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0554793) q[2];
sx q[2];
rz(-1.5986634) q[2];
sx q[2];
rz(1.0184658) q[2];
rz(-0.57824072) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2932435) q[0];
sx q[0];
rz(-1.542955) q[0];
sx q[0];
rz(-1.6795213) q[0];
rz(2.2734185) q[1];
sx q[1];
rz(-1.9055007) q[1];
sx q[1];
rz(0.47913924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3972319) q[0];
sx q[0];
rz(-1.8313837) q[0];
sx q[0];
rz(0.67407383) q[0];
rz(-pi) q[1];
rz(0.43897668) q[2];
sx q[2];
rz(-1.8700646) q[2];
sx q[2];
rz(-1.0004832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16461611) q[1];
sx q[1];
rz(-0.45358001) q[1];
sx q[1];
rz(-3.0248887) q[1];
rz(2.3684351) q[3];
sx q[3];
rz(-0.36423238) q[3];
sx q[3];
rz(0.75477597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.22970557) q[2];
sx q[2];
rz(-0.90596002) q[2];
sx q[2];
rz(-0.68944302) q[2];
rz(2.7202969) q[3];
sx q[3];
rz(-2.3861986) q[3];
sx q[3];
rz(1.1675444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1223849) q[0];
sx q[0];
rz(-2.8715219) q[0];
sx q[0];
rz(-0.58575678) q[0];
rz(-0.44341227) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(-2.401039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51979216) q[0];
sx q[0];
rz(-2.2211638) q[0];
sx q[0];
rz(2.4629794) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1525864) q[2];
sx q[2];
rz(-0.64179342) q[2];
sx q[2];
rz(-0.60371232) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2948998) q[1];
sx q[1];
rz(-0.36622643) q[1];
sx q[1];
rz(-0.80898714) q[1];
rz(1.8216474) q[3];
sx q[3];
rz(-1.8545058) q[3];
sx q[3];
rz(1.2267867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1195141) q[2];
sx q[2];
rz(-1.4351279) q[2];
sx q[2];
rz(-2.4100927) q[2];
rz(2.5445407) q[3];
sx q[3];
rz(-2.4686333) q[3];
sx q[3];
rz(1.2520242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6325833) q[0];
sx q[0];
rz(-2.3409797) q[0];
sx q[0];
rz(-1.5492232) q[0];
rz(2.1005311) q[1];
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
rz(0.072424732) q[0];
sx q[0];
rz(-0.7382881) q[0];
sx q[0];
rz(2.9641572) q[0];
x q[1];
rz(-3.1388624) q[2];
sx q[2];
rz(-1.8282253) q[2];
sx q[2];
rz(-2.4705459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10080316) q[1];
sx q[1];
rz(-1.2606422) q[1];
sx q[1];
rz(0.47001229) q[1];
rz(-pi) q[2];
rz(-0.18852461) q[3];
sx q[3];
rz(-0.97947165) q[3];
sx q[3];
rz(-2.6935379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8063712) q[2];
sx q[2];
rz(-1.0314032) q[2];
sx q[2];
rz(-1.5058937) q[2];
rz(0.38883543) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(2.4839362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.036670551) q[0];
sx q[0];
rz(-0.67476434) q[0];
sx q[0];
rz(-0.92055935) q[0];
rz(0.29832828) q[1];
sx q[1];
rz(-2.0896656) q[1];
sx q[1];
rz(-1.9783609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0434303) q[0];
sx q[0];
rz(-0.32458186) q[0];
sx q[0];
rz(-3.1107706) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79090581) q[2];
sx q[2];
rz(-1.8703572) q[2];
sx q[2];
rz(-1.1096481) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.498658) q[1];
sx q[1];
rz(-2.1172595) q[1];
sx q[1];
rz(-1.5692753) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(1.3433015) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(-3.0291271) q[2];
rz(-0.28137842) q[3];
sx q[3];
rz(-2.3795542) q[3];
sx q[3];
rz(1.5828097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6211468) q[0];
sx q[0];
rz(-0.96108323) q[0];
sx q[0];
rz(3.0302826) q[0];
rz(0.81980199) q[1];
sx q[1];
rz(-0.83378053) q[1];
sx q[1];
rz(0.05096635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5539885) q[0];
sx q[0];
rz(-1.7394571) q[0];
sx q[0];
rz(-1.8722562) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3459008) q[2];
sx q[2];
rz(-2.1818265) q[2];
sx q[2];
rz(-2.261206) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8680893) q[1];
sx q[1];
rz(-0.13493262) q[1];
sx q[1];
rz(-1.8084303) q[1];
rz(2.662728) q[3];
sx q[3];
rz(-2.1486946) q[3];
sx q[3];
rz(-0.013805429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0593947) q[2];
sx q[2];
rz(-1.1092564) q[2];
sx q[2];
rz(1.8471897) q[2];
rz(2.1829677) q[3];
sx q[3];
rz(-0.38161033) q[3];
sx q[3];
rz(-2.3310048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0628919) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(-2.1953951) q[0];
rz(1.9323438) q[1];
sx q[1];
rz(-0.47018662) q[1];
sx q[1];
rz(1.6016003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9935051) q[0];
sx q[0];
rz(-0.74310857) q[0];
sx q[0];
rz(2.00032) q[0];
rz(-pi) q[1];
rz(-0.49175434) q[2];
sx q[2];
rz(-1.3637614) q[2];
sx q[2];
rz(-1.7397665) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-2.8084635) q[3];
sx q[3];
rz(-0.58099834) q[3];
sx q[3];
rz(2.6521386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20751791) q[2];
sx q[2];
rz(-1.3413651) q[2];
sx q[2];
rz(1.6259469) q[2];
rz(-0.80306426) q[3];
sx q[3];
rz(-2.8290437) q[3];
sx q[3];
rz(-1.1481185) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37907264) q[0];
sx q[0];
rz(-1.9110492) q[0];
sx q[0];
rz(-1.8050964) q[0];
rz(1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(1.1755099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0031457) q[0];
sx q[0];
rz(-1.8120736) q[0];
sx q[0];
rz(0.74575119) q[0];
rz(-pi) q[1];
rz(1.8075502) q[2];
sx q[2];
rz(-2.3375698) q[2];
sx q[2];
rz(-1.6171111) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0693722) q[1];
sx q[1];
rz(-1.4319422) q[1];
sx q[1];
rz(-1.6235823) q[1];
rz(-pi) q[2];
rz(-2.544459) q[3];
sx q[3];
rz(-1.8270511) q[3];
sx q[3];
rz(0.090274752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0846348) q[2];
sx q[2];
rz(-1.6677758) q[2];
sx q[2];
rz(-1.8353621) q[2];
rz(-2.7133387) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(-2.7603507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-1.320095) q[0];
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
rz(0.069396917) q[3];
sx q[3];
rz(-1.4213448) q[3];
sx q[3];
rz(0.16363174) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
