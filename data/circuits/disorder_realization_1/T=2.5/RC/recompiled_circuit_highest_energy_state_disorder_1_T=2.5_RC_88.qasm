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
rz(-2.2358492) q[1];
sx q[1];
rz(-1.716776) q[1];
sx q[1];
rz(-2.4832771) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8480331) q[0];
sx q[0];
rz(-1.2633889) q[0];
sx q[0];
rz(-2.2821167) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7204804) q[2];
sx q[2];
rz(-2.176365) q[2];
sx q[2];
rz(-0.31878135) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4791245) q[1];
sx q[1];
rz(-0.82193804) q[1];
sx q[1];
rz(2.2884018) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0327037) q[3];
sx q[3];
rz(-1.3815839) q[3];
sx q[3];
rz(-3.12836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1540404) q[2];
sx q[2];
rz(-2.0200358) q[2];
sx q[2];
rz(1.7169607) q[2];
rz(0.79036653) q[3];
sx q[3];
rz(-0.75420403) q[3];
sx q[3];
rz(0.2592026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(1.8081283) q[0];
rz(1.85224) q[1];
sx q[1];
rz(-2.5080296) q[1];
sx q[1];
rz(-0.11817558) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2424392) q[0];
sx q[0];
rz(-0.94510539) q[0];
sx q[0];
rz(-3.1345128) q[0];
rz(-pi) q[1];
rz(-1.9508544) q[2];
sx q[2];
rz(-0.59761492) q[2];
sx q[2];
rz(2.3412398) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11013438) q[1];
sx q[1];
rz(-2.2074568) q[1];
sx q[1];
rz(2.6183271) q[1];
rz(2.7998522) q[3];
sx q[3];
rz(-1.8287418) q[3];
sx q[3];
rz(0.3697357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3893343) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(-2.3040859) q[2];
rz(0.62344712) q[3];
sx q[3];
rz(-1.1102755) q[3];
sx q[3];
rz(0.71201223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74152827) q[0];
sx q[0];
rz(-2.9900804) q[0];
sx q[0];
rz(-0.22756273) q[0];
rz(0.96861068) q[1];
sx q[1];
rz(-2.249735) q[1];
sx q[1];
rz(-1.4898941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9609531) q[0];
sx q[0];
rz(-2.0386057) q[0];
sx q[0];
rz(2.803928) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7329839) q[2];
sx q[2];
rz(-1.6005862) q[2];
sx q[2];
rz(-1.839585) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35830733) q[1];
sx q[1];
rz(-0.44898012) q[1];
sx q[1];
rz(2.994356) q[1];
rz(2.3132626) q[3];
sx q[3];
rz(-1.7975494) q[3];
sx q[3];
rz(-1.5080719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0554793) q[2];
sx q[2];
rz(-1.5986634) q[2];
sx q[2];
rz(-1.0184658) q[2];
rz(0.57824072) q[3];
sx q[3];
rz(-2.6522804) q[3];
sx q[3];
rz(1.9448634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2932435) q[0];
sx q[0];
rz(-1.542955) q[0];
sx q[0];
rz(-1.6795213) q[0];
rz(0.86817414) q[1];
sx q[1];
rz(-1.236092) q[1];
sx q[1];
rz(0.47913924) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7650267) q[0];
sx q[0];
rz(-0.92341316) q[0];
sx q[0];
rz(1.8996973) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2422945) q[2];
sx q[2];
rz(-1.9890033) q[2];
sx q[2];
rz(2.433726) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5111749) q[1];
sx q[1];
rz(-1.6218405) q[1];
sx q[1];
rz(-0.45089727) q[1];
rz(1.3105743) q[3];
sx q[3];
rz(-1.3129915) q[3];
sx q[3];
rz(1.5797404) q[3];
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
rz(-0.68944302) q[2];
rz(0.42129579) q[3];
sx q[3];
rz(-2.3861986) q[3];
sx q[3];
rz(1.9740483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01920779) q[0];
sx q[0];
rz(-2.8715219) q[0];
sx q[0];
rz(0.58575678) q[0];
rz(2.6981804) q[1];
sx q[1];
rz(-1.6165918) q[1];
sx q[1];
rz(-2.401039) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4465312) q[0];
sx q[0];
rz(-0.90264812) q[0];
sx q[0];
rz(-0.88094385) q[0];
rz(-pi) q[1];
rz(1.1525864) q[2];
sx q[2];
rz(-2.4997992) q[2];
sx q[2];
rz(-2.5378803) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49879229) q[1];
sx q[1];
rz(-1.3086924) q[1];
sx q[1];
rz(2.8828055) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4363229) q[3];
sx q[3];
rz(-2.765145) q[3];
sx q[3];
rz(-1.173526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1195141) q[2];
sx q[2];
rz(-1.7064648) q[2];
sx q[2];
rz(-0.73149991) q[2];
rz(-0.59705192) q[3];
sx q[3];
rz(-2.4686333) q[3];
sx q[3];
rz(1.2520242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(2.1005311) q[1];
sx q[1];
rz(-0.58716232) q[1];
sx q[1];
rz(2.7802173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8313077) q[0];
sx q[0];
rz(-0.84670369) q[0];
sx q[0];
rz(-1.4115439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.560426) q[2];
sx q[2];
rz(-2.8841495) q[2];
sx q[2];
rz(-0.68177047) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0407895) q[1];
sx q[1];
rz(-1.2606422) q[1];
sx q[1];
rz(2.6715804) q[1];
rz(-pi) q[2];
rz(0.18852461) q[3];
sx q[3];
rz(-0.97947165) q[3];
sx q[3];
rz(-0.44805474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8063712) q[2];
sx q[2];
rz(-2.1101895) q[2];
sx q[2];
rz(1.5058937) q[2];
rz(-0.38883543) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(0.65765643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036670551) q[0];
sx q[0];
rz(-0.67476434) q[0];
sx q[0];
rz(0.92055935) q[0];
rz(2.8432644) q[1];
sx q[1];
rz(-2.0896656) q[1];
sx q[1];
rz(1.9783609) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0981623) q[0];
sx q[0];
rz(-0.32458186) q[0];
sx q[0];
rz(-3.1107706) q[0];
rz(-0.79090581) q[2];
sx q[2];
rz(-1.2712354) q[2];
sx q[2];
rz(-1.1096481) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0702447) q[1];
sx q[1];
rz(-1.5694968) q[1];
sx q[1];
rz(-2.595129) q[1];
x q[2];
rz(1.8697132) q[3];
sx q[3];
rz(-1.7038725) q[3];
sx q[3];
rz(-1.4081334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3433015) q[2];
sx q[2];
rz(-1.9140665) q[2];
sx q[2];
rz(3.0291271) q[2];
rz(-2.8602142) q[3];
sx q[3];
rz(-0.76203841) q[3];
sx q[3];
rz(1.5828097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
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
rz(-3.0906263) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726377) q[0];
sx q[0];
rz(-1.8678471) q[0];
sx q[0];
rz(-0.1764651) q[0];
rz(2.5185263) q[2];
sx q[2];
rz(-1.754481) q[2];
sx q[2];
rz(-2.3206835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6283577) q[1];
sx q[1];
rz(-1.7019148) q[1];
sx q[1];
rz(3.1096457) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2045129) q[3];
sx q[3];
rz(-1.9670319) q[3];
sx q[3];
rz(-1.8609488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0593947) q[2];
sx q[2];
rz(-1.1092564) q[2];
sx q[2];
rz(-1.294403) q[2];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078700773) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(-0.94619757) q[0];
rz(-1.9323438) q[1];
sx q[1];
rz(-2.671406) q[1];
sx q[1];
rz(1.6016003) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14808753) q[0];
sx q[0];
rz(-0.74310857) q[0];
sx q[0];
rz(-2.00032) q[0];
x q[1];
rz(-1.3368789) q[2];
sx q[2];
rz(-2.0511464) q[2];
sx q[2];
rz(-3.0822871) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29344236) q[1];
sx q[1];
rz(-2.5328325) q[1];
sx q[1];
rz(-2.7431627) q[1];
rz(0.55535613) q[3];
sx q[3];
rz(-1.7512519) q[3];
sx q[3];
rz(2.3418135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20751791) q[2];
sx q[2];
rz(-1.8002276) q[2];
sx q[2];
rz(-1.5156457) q[2];
rz(-2.3385284) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(-1.1481185) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37907264) q[0];
sx q[0];
rz(-1.2305434) q[0];
sx q[0];
rz(1.8050964) q[0];
rz(1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(-1.9660827) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208082) q[0];
sx q[0];
rz(-2.3649923) q[0];
sx q[0];
rz(-0.34790502) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3340424) q[2];
sx q[2];
rz(-0.8040229) q[2];
sx q[2];
rz(1.5244816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4912632) q[1];
sx q[1];
rz(-1.5185188) q[1];
sx q[1];
rz(3.0025474) q[1];
rz(-pi) q[2];
rz(2.7055199) q[3];
sx q[3];
rz(-2.498004) q[3];
sx q[3];
rz(-1.3042579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.056957873) q[2];
sx q[2];
rz(-1.4738169) q[2];
sx q[2];
rz(1.3062306) q[2];
rz(2.7133387) q[3];
sx q[3];
rz(-1.8249325) q[3];
sx q[3];
rz(-0.38124198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8214977) q[0];
sx q[0];
rz(-1.7392673) q[0];
sx q[0];
rz(2.0489954) q[0];
rz(-0.26662695) q[1];
sx q[1];
rz(-1.0933924) q[1];
sx q[1];
rz(-0.46270121) q[1];
rz(0.54556432) q[2];
sx q[2];
rz(-2.5013899) q[2];
sx q[2];
rz(-1.987006) q[2];
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
