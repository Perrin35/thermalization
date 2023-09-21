OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
rz(-0.31495467) q[1];
sx q[1];
rz(-2.0575674) q[1];
sx q[1];
rz(-1.6853583) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0460912) q[0];
sx q[0];
rz(-1.4399733) q[0];
sx q[0];
rz(0.019898947) q[0];
rz(0.37402447) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(-2.6467269) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.734182) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(0.68024866) q[1];
rz(-0.16430328) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(-2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2125856) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(0.47098413) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1176227) q[0];
sx q[0];
rz(-1.6256486) q[0];
sx q[0];
rz(-1.0883925) q[0];
x q[1];
rz(-1.2330301) q[2];
sx q[2];
rz(-0.62440364) q[2];
sx q[2];
rz(1.0242467) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4091332) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(-0.061846102) q[1];
x q[2];
rz(3.0456411) q[3];
sx q[3];
rz(-0.6978242) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.058078893) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(0.2557959) q[2];
rz(1.6563709) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(2.9799057) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3210179) q[0];
sx q[0];
rz(-1.6512172) q[0];
sx q[0];
rz(-3.1119425) q[0];
x q[1];
rz(0.36053948) q[2];
sx q[2];
rz(-1.5224824) q[2];
sx q[2];
rz(-0.53256065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1767133) q[1];
sx q[1];
rz(-2.5039154) q[1];
sx q[1];
rz(1.8646851) q[1];
x q[2];
rz(-0.93033354) q[3];
sx q[3];
rz(-0.6593245) q[3];
sx q[3];
rz(1.7742771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(0.20047323) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4865206) q[0];
sx q[0];
rz(-2.1694896) q[0];
sx q[0];
rz(0.18874164) q[0];
x q[1];
rz(-0.38509102) q[2];
sx q[2];
rz(-1.4355806) q[2];
sx q[2];
rz(1.1831634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7728459) q[1];
sx q[1];
rz(-2.6987942) q[1];
sx q[1];
rz(0.55261353) q[1];
rz(-pi) q[2];
rz(3.0365385) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(0.74329174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071881) q[0];
sx q[0];
rz(-0.92538639) q[0];
sx q[0];
rz(2.8651644) q[0];
x q[1];
rz(-2.6402316) q[2];
sx q[2];
rz(-1.6622346) q[2];
sx q[2];
rz(-2.3102592) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6965461) q[1];
sx q[1];
rz(-1.7940709) q[1];
sx q[1];
rz(-2.9162507) q[1];
rz(-pi) q[2];
rz(0.27637847) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(2.7793022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-1.0553029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366093) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(1.2654632) q[0];
rz(2.6199117) q[2];
sx q[2];
rz(-0.98419596) q[2];
sx q[2];
rz(-2.4751543) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1736974) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(2.150279) q[1];
x q[2];
rz(0.25552337) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(-1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-2.5678182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.505578) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(-1.3789603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4622719) q[2];
sx q[2];
rz(-0.2013686) q[2];
sx q[2];
rz(0.078171922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3473914) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(2.3458523) q[1];
x q[2];
rz(-1.0049099) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(-2.3435081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-0.77073628) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(-1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.27281) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(-1.8458813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0062795) q[0];
sx q[0];
rz(-1.2426002) q[0];
sx q[0];
rz(-0.4853863) q[0];
rz(1.9240575) q[2];
sx q[2];
rz(-1.8262987) q[2];
sx q[2];
rz(3.0614292) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2669275) q[1];
sx q[1];
rz(-1.7942567) q[1];
sx q[1];
rz(2.2063971) q[1];
rz(3.1136884) q[3];
sx q[3];
rz(-1.1301665) q[3];
sx q[3];
rz(-0.84907109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(-0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-2.303404) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121376) q[0];
sx q[0];
rz(-1.9337618) q[0];
sx q[0];
rz(-1.8118993) q[0];
rz(-pi) q[1];
rz(1.6423177) q[2];
sx q[2];
rz(-0.75746775) q[2];
sx q[2];
rz(-2.1069991) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6088241) q[1];
sx q[1];
rz(-2.3012487) q[1];
sx q[1];
rz(0.26809147) q[1];
x q[2];
rz(-2.8832199) q[3];
sx q[3];
rz(-1.6515454) q[3];
sx q[3];
rz(-1.0563869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(1.2667123) q[2];
rz(2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(-0.094749711) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-0.231803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7451413) q[0];
sx q[0];
rz(-1.6661577) q[0];
sx q[0];
rz(-2.0800637) q[0];
rz(-pi) q[1];
rz(-2.899029) q[2];
sx q[2];
rz(-1.3360268) q[2];
sx q[2];
rz(0.64005062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.787549) q[1];
sx q[1];
rz(-1.429306) q[1];
sx q[1];
rz(-0.267412) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2139123) q[3];
sx q[3];
rz(-2.247346) q[3];
sx q[3];
rz(0.032534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(0.96881962) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(-1.8112524) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
rz(1.590329) q[3];
sx q[3];
rz(-1.3828779) q[3];
sx q[3];
rz(-1.8120017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];