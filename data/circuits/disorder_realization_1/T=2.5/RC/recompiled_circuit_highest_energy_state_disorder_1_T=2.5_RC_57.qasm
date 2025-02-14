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
rz(-0.6583156) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29355958) q[0];
sx q[0];
rz(-1.8782037) q[0];
sx q[0];
rz(0.859476) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7204804) q[2];
sx q[2];
rz(-0.96522766) q[2];
sx q[2];
rz(-0.31878135) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5708471) q[1];
sx q[1];
rz(-2.1553504) q[1];
sx q[1];
rz(-0.61573124) q[1];
rz(0.2194276) q[3];
sx q[3];
rz(-1.0433222) q[3];
sx q[3];
rz(-1.4722401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98755223) q[2];
sx q[2];
rz(-1.1215569) q[2];
sx q[2];
rz(-1.7169607) q[2];
rz(-0.79036653) q[3];
sx q[3];
rz(-0.75420403) q[3];
sx q[3];
rz(2.8823901) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454999) q[0];
sx q[0];
rz(-0.51829618) q[0];
sx q[0];
rz(-1.3334644) q[0];
rz(-1.2893527) q[1];
sx q[1];
rz(-0.63356304) q[1];
sx q[1];
rz(0.11817558) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2424392) q[0];
sx q[0];
rz(-2.1964873) q[0];
sx q[0];
rz(0.0070798955) q[0];
rz(-1.0071272) q[2];
sx q[2];
rz(-1.360513) q[2];
sx q[2];
rz(2.0521833) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11013438) q[1];
sx q[1];
rz(-0.93413583) q[1];
sx q[1];
rz(-2.6183271) q[1];
rz(-pi) q[2];
rz(-0.66690655) q[3];
sx q[3];
rz(-2.7165037) q[3];
sx q[3];
rz(-0.57890427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3893343) q[2];
sx q[2];
rz(-0.78395939) q[2];
sx q[2];
rz(-0.83750677) q[2];
rz(0.62344712) q[3];
sx q[3];
rz(-1.1102755) q[3];
sx q[3];
rz(0.71201223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74152827) q[0];
sx q[0];
rz(-0.15151227) q[0];
sx q[0];
rz(-0.22756273) q[0];
rz(-0.96861068) q[1];
sx q[1];
rz(-2.249735) q[1];
sx q[1];
rz(-1.6516986) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.299376) q[0];
sx q[0];
rz(-0.56952667) q[0];
sx q[0];
rz(-0.99040188) q[0];
rz(-pi) q[1];
rz(1.7532811) q[2];
sx q[2];
rz(-2.9767155) q[2];
sx q[2];
rz(0.08872513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.35830733) q[1];
sx q[1];
rz(-0.44898012) q[1];
sx q[1];
rz(0.14723666) q[1];
rz(-pi) q[2];
rz(-1.2419385) q[3];
sx q[3];
rz(-2.3716381) q[3];
sx q[3];
rz(2.9640523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0861133) q[2];
sx q[2];
rz(-1.5429292) q[2];
sx q[2];
rz(-1.0184658) q[2];
rz(0.57824072) q[3];
sx q[3];
rz(-0.48931229) q[3];
sx q[3];
rz(-1.9448634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
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
rz(0.86817414) q[1];
sx q[1];
rz(-1.9055007) q[1];
sx q[1];
rz(2.6624534) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3972319) q[0];
sx q[0];
rz(-1.3102089) q[0];
sx q[0];
rz(0.67407383) q[0];
rz(-pi) q[1];
rz(-2.702616) q[2];
sx q[2];
rz(-1.8700646) q[2];
sx q[2];
rz(-1.0004832) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9769765) q[1];
sx q[1];
rz(-0.45358001) q[1];
sx q[1];
rz(3.0248887) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3684351) q[3];
sx q[3];
rz(-0.36423238) q[3];
sx q[3];
rz(0.75477597) q[3];
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
rz(-0.68944302) q[2];
rz(-0.42129579) q[3];
sx q[3];
rz(-2.3861986) q[3];
sx q[3];
rz(1.1675444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01920779) q[0];
sx q[0];
rz(-2.8715219) q[0];
sx q[0];
rz(-0.58575678) q[0];
rz(0.44341227) q[1];
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
rz(-2.6218005) q[0];
sx q[0];
rz(-2.2211638) q[0];
sx q[0];
rz(0.67861329) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9716203) q[2];
sx q[2];
rz(-1.8163774) q[2];
sx q[2];
rz(1.8325017) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2948998) q[1];
sx q[1];
rz(-0.36622643) q[1];
sx q[1];
rz(-0.80898714) q[1];
x q[2];
rz(-2.8492227) q[3];
sx q[3];
rz(-1.3301759) q[3];
sx q[3];
rz(-0.27240348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1195141) q[2];
sx q[2];
rz(-1.7064648) q[2];
sx q[2];
rz(-0.73149991) q[2];
rz(-0.59705192) q[3];
sx q[3];
rz(-2.4686333) q[3];
sx q[3];
rz(-1.8895684) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6325833) q[0];
sx q[0];
rz(-2.3409797) q[0];
sx q[0];
rz(1.5492232) q[0];
rz(2.1005311) q[1];
sx q[1];
rz(-2.5544303) q[1];
sx q[1];
rz(-2.7802173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072424732) q[0];
sx q[0];
rz(-2.4033045) q[0];
sx q[0];
rz(-2.9641572) q[0];
x q[1];
rz(-1.5811666) q[2];
sx q[2];
rz(-0.25744312) q[2];
sx q[2];
rz(2.4598222) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10080316) q[1];
sx q[1];
rz(-1.2606422) q[1];
sx q[1];
rz(-2.6715804) q[1];
x q[2];
rz(0.18852461) q[3];
sx q[3];
rz(-2.162121) q[3];
sx q[3];
rz(0.44805474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3352215) q[2];
sx q[2];
rz(-2.1101895) q[2];
sx q[2];
rz(-1.6356989) q[2];
rz(0.38883543) q[3];
sx q[3];
rz(-0.54986984) q[3];
sx q[3];
rz(2.4839362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036670551) q[0];
sx q[0];
rz(-0.67476434) q[0];
sx q[0];
rz(0.92055935) q[0];
rz(-0.29832828) q[1];
sx q[1];
rz(-2.0896656) q[1];
sx q[1];
rz(-1.1632318) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0656433) q[0];
sx q[0];
rz(-1.8952184) q[0];
sx q[0];
rz(1.5811654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7317861) q[2];
sx q[2];
rz(-0.834045) q[2];
sx q[2];
rz(-0.17716874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6400078) q[1];
sx q[1];
rz(-0.54646508) q[1];
sx q[1];
rz(0.0025005977) q[1];
rz(-pi) q[2];
rz(1.8697132) q[3];
sx q[3];
rz(-1.7038725) q[3];
sx q[3];
rz(1.7334593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3433015) q[2];
sx q[2];
rz(-1.2275262) q[2];
sx q[2];
rz(-3.0291271) q[2];
rz(-2.8602142) q[3];
sx q[3];
rz(-0.76203841) q[3];
sx q[3];
rz(-1.5587829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726377) q[0];
sx q[0];
rz(-1.8678471) q[0];
sx q[0];
rz(-0.1764651) q[0];
rz(1.3459008) q[2];
sx q[2];
rz(-0.95976613) q[2];
sx q[2];
rz(2.261206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2735034) q[1];
sx q[1];
rz(-3.00666) q[1];
sx q[1];
rz(-1.3331624) q[1];
rz(2.662728) q[3];
sx q[3];
rz(-0.99289809) q[3];
sx q[3];
rz(-3.1277872) q[3];
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
rz(-2.1829677) q[3];
sx q[3];
rz(-0.38161033) q[3];
sx q[3];
rz(2.3310048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078700773) q[0];
sx q[0];
rz(-1.8475516) q[0];
sx q[0];
rz(0.94619757) q[0];
rz(1.2092489) q[1];
sx q[1];
rz(-2.671406) q[1];
sx q[1];
rz(-1.5399923) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3935767) q[0];
sx q[0];
rz(-1.8564176) q[0];
sx q[0];
rz(-2.266721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41856022) q[2];
sx q[2];
rz(-2.6113434) q[2];
sx q[2];
rz(-2.6061932) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6098574) q[1];
sx q[1];
rz(-1.3470728) q[1];
sx q[1];
rz(-0.57106627) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8084635) q[3];
sx q[3];
rz(-2.5605943) q[3];
sx q[3];
rz(-0.48945409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9340747) q[2];
sx q[2];
rz(-1.8002276) q[2];
sx q[2];
rz(-1.5156457) q[2];
rz(-0.80306426) q[3];
sx q[3];
rz(-0.31254891) q[3];
sx q[3];
rz(1.1481185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37907264) q[0];
sx q[0];
rz(-1.9110492) q[0];
sx q[0];
rz(1.3364963) q[0];
rz(-1.1539917) q[1];
sx q[1];
rz(-2.3497252) q[1];
sx q[1];
rz(-1.1755099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35042351) q[0];
sx q[0];
rz(-0.8514815) q[0];
sx q[0];
rz(1.2475623) q[0];
x q[1];
rz(-1.8075502) q[2];
sx q[2];
rz(-2.3375698) q[2];
sx q[2];
rz(-1.5244816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6503295) q[1];
sx q[1];
rz(-1.6230738) q[1];
sx q[1];
rz(-3.0025474) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8776348) q[3];
sx q[3];
rz(-2.1458907) q[3];
sx q[3];
rz(1.8317312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.056957873) q[2];
sx q[2];
rz(-1.6677758) q[2];
sx q[2];
rz(1.3062306) q[2];
rz(-2.7133387) q[3];
sx q[3];
rz(-1.3166602) q[3];
sx q[3];
rz(2.7603507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(0.56699085) q[2];
sx q[2];
rz(-1.8859572) q[2];
sx q[2];
rz(0.036833298) q[2];
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
