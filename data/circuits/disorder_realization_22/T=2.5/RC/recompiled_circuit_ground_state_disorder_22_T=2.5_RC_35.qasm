OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9353256) q[0];
sx q[0];
rz(4.7630035) q[0];
sx q[0];
rz(8.0061316) q[0];
rz(6.3267646) q[1];
sx q[1];
rz(5.0001231) q[1];
sx q[1];
rz(12.725886) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.537764) q[0];
sx q[0];
rz(-1.7366752) q[0];
sx q[0];
rz(-0.028847522) q[0];
x q[1];
rz(2.5297727) q[2];
sx q[2];
rz(-1.7828336) q[2];
sx q[2];
rz(3.0481047) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2251896) q[1];
sx q[1];
rz(-1.708751) q[1];
sx q[1];
rz(-1.0911343) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1662935) q[3];
sx q[3];
rz(-2.3784436) q[3];
sx q[3];
rz(2.7119702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.688545) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(0.090252074) q[2];
rz(3.0669751) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(0.41372764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8892882) q[0];
sx q[0];
rz(-1.1587208) q[0];
sx q[0];
rz(1.4537551) q[0];
rz(-2.3989035) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(0.76505605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.830421) q[0];
sx q[0];
rz(-0.97872916) q[0];
sx q[0];
rz(1.4903846) q[0];
rz(0.86642577) q[2];
sx q[2];
rz(-1.0110628) q[2];
sx q[2];
rz(2.100796) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2606352) q[1];
sx q[1];
rz(-2.2876014) q[1];
sx q[1];
rz(-2.6997801) q[1];
rz(-1.5734777) q[3];
sx q[3];
rz(-2.1187927) q[3];
sx q[3];
rz(0.23220095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9146156) q[2];
sx q[2];
rz(-2.1299629) q[2];
sx q[2];
rz(-1.6870618) q[2];
rz(-2.4948273) q[3];
sx q[3];
rz(-2.5819467) q[3];
sx q[3];
rz(-1.8618934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0525368) q[0];
sx q[0];
rz(-1.4707969) q[0];
sx q[0];
rz(-3.0964417) q[0];
rz(1.9106983) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(-2.0848354) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741183) q[0];
sx q[0];
rz(-1.2454659) q[0];
sx q[0];
rz(2.1528457) q[0];
rz(2.4963794) q[2];
sx q[2];
rz(-1.9786547) q[2];
sx q[2];
rz(2.4347664) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9581117) q[1];
sx q[1];
rz(-0.90881244) q[1];
sx q[1];
rz(-0.32420154) q[1];
rz(2.908024) q[3];
sx q[3];
rz(-2.2456777) q[3];
sx q[3];
rz(2.6949331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7354273) q[2];
sx q[2];
rz(-2.1801517) q[2];
sx q[2];
rz(0.65262922) q[2];
rz(-1.2930219) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(-0.098171083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80766135) q[0];
sx q[0];
rz(-2.6787391) q[0];
sx q[0];
rz(1.7474784) q[0];
rz(1.8380503) q[1];
sx q[1];
rz(-1.1640254) q[1];
sx q[1];
rz(2.0533144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1980285) q[0];
sx q[0];
rz(-1.8450415) q[0];
sx q[0];
rz(0.79601566) q[0];
x q[1];
rz(0.59594229) q[2];
sx q[2];
rz(-2.1998027) q[2];
sx q[2];
rz(-0.12550917) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40068059) q[1];
sx q[1];
rz(-2.0517382) q[1];
sx q[1];
rz(0.6948284) q[1];
rz(-0.49934629) q[3];
sx q[3];
rz(-1.5245228) q[3];
sx q[3];
rz(2.3114227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3950222) q[2];
sx q[2];
rz(-1.3777379) q[2];
sx q[2];
rz(2.6611888) q[2];
rz(0.26614842) q[3];
sx q[3];
rz(-1.1689309) q[3];
sx q[3];
rz(-0.84421414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10923037) q[0];
sx q[0];
rz(-2.5205595) q[0];
sx q[0];
rz(1.9367628) q[0];
rz(3.0989528) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(0.94243324) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401817) q[0];
sx q[0];
rz(-1.4352006) q[0];
sx q[0];
rz(-1.1419161) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93330002) q[2];
sx q[2];
rz(-1.0663053) q[2];
sx q[2];
rz(1.2802893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.98676315) q[1];
sx q[1];
rz(-1.0752605) q[1];
sx q[1];
rz(-2.7096728) q[1];
rz(-2.5025401) q[3];
sx q[3];
rz(-1.6585232) q[3];
sx q[3];
rz(1.4301585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.11613906) q[2];
sx q[2];
rz(-2.4567273) q[2];
sx q[2];
rz(1.2241414) q[2];
rz(0.72743607) q[3];
sx q[3];
rz(-1.4801315) q[3];
sx q[3];
rz(-0.050189655) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3392451) q[0];
sx q[0];
rz(-0.50528637) q[0];
sx q[0];
rz(-0.25830609) q[0];
rz(0.91336617) q[1];
sx q[1];
rz(-0.86562997) q[1];
sx q[1];
rz(0.9679274) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0570404) q[0];
sx q[0];
rz(-1.0121317) q[0];
sx q[0];
rz(-2.0303805) q[0];
x q[1];
rz(2.2931678) q[2];
sx q[2];
rz(-1.8471187) q[2];
sx q[2];
rz(1.487731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.63686759) q[1];
sx q[1];
rz(-0.47885413) q[1];
sx q[1];
rz(-0.07696159) q[1];
x q[2];
rz(2.4960732) q[3];
sx q[3];
rz(-1.1214167) q[3];
sx q[3];
rz(-1.9674894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4919058) q[2];
sx q[2];
rz(-0.68444362) q[2];
sx q[2];
rz(0.41927949) q[2];
rz(-0.78420091) q[3];
sx q[3];
rz(-1.0201642) q[3];
sx q[3];
rz(-1.7694582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.075031) q[0];
sx q[0];
rz(-3.0176268) q[0];
sx q[0];
rz(0.01874622) q[0];
rz(0.092451267) q[1];
sx q[1];
rz(-1.8094614) q[1];
sx q[1];
rz(0.094955347) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2128747) q[0];
sx q[0];
rz(-1.2480019) q[0];
sx q[0];
rz(2.5514437) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58231852) q[2];
sx q[2];
rz(-2.4597801) q[2];
sx q[2];
rz(-0.63722975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.590789) q[1];
sx q[1];
rz(-1.992464) q[1];
sx q[1];
rz(1.8487537) q[1];
rz(1.9390887) q[3];
sx q[3];
rz(-2.7164258) q[3];
sx q[3];
rz(1.7452546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9065325) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(2.2825784) q[2];
rz(-2.5281995) q[3];
sx q[3];
rz(-0.20039138) q[3];
sx q[3];
rz(1.332351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5253946) q[0];
sx q[0];
rz(-1.8680251) q[0];
sx q[0];
rz(-1.4816351) q[0];
rz(2.0741277) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(-1.8225089) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049371) q[0];
sx q[0];
rz(-0.67064136) q[0];
sx q[0];
rz(-1.4192307) q[0];
rz(-pi) q[1];
x q[1];
rz(1.943687) q[2];
sx q[2];
rz(-2.0569062) q[2];
sx q[2];
rz(-1.1463469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47506491) q[1];
sx q[1];
rz(-1.804978) q[1];
sx q[1];
rz(0.98824309) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6982001) q[3];
sx q[3];
rz(-1.9249328) q[3];
sx q[3];
rz(2.4417801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2655502) q[2];
sx q[2];
rz(-0.37028131) q[2];
sx q[2];
rz(3.1004356) q[2];
rz(1.1228849) q[3];
sx q[3];
rz(-1.4833781) q[3];
sx q[3];
rz(2.6179204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48162833) q[0];
sx q[0];
rz(-0.97844231) q[0];
sx q[0];
rz(1.6690669) q[0];
rz(2.4166079) q[1];
sx q[1];
rz(-1.1712733) q[1];
sx q[1];
rz(2.5670126) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32550463) q[0];
sx q[0];
rz(-0.34127745) q[0];
sx q[0];
rz(0.14320548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99733229) q[2];
sx q[2];
rz(-2.6509299) q[2];
sx q[2];
rz(0.034324797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1657527) q[1];
sx q[1];
rz(-0.71893349) q[1];
sx q[1];
rz(-0.045054212) q[1];
rz(1.3909886) q[3];
sx q[3];
rz(-2.4459908) q[3];
sx q[3];
rz(2.4814388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4943115) q[2];
sx q[2];
rz(-2.8936671) q[2];
sx q[2];
rz(-0.9606804) q[2];
rz(-2.658355) q[3];
sx q[3];
rz(-1.5471231) q[3];
sx q[3];
rz(-1.5794992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305337) q[0];
sx q[0];
rz(-2.647825) q[0];
sx q[0];
rz(0.44179398) q[0];
rz(2.4590625) q[1];
sx q[1];
rz(-1.4492757) q[1];
sx q[1];
rz(2.0880879) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614612) q[0];
sx q[0];
rz(-2.2765358) q[0];
sx q[0];
rz(-2.8012026) q[0];
x q[1];
rz(1.7948661) q[2];
sx q[2];
rz(-2.7599979) q[2];
sx q[2];
rz(3.0004108) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9113831) q[1];
sx q[1];
rz(-2.1911096) q[1];
sx q[1];
rz(2.7432454) q[1];
x q[2];
rz(-2.7476596) q[3];
sx q[3];
rz(-1.6389168) q[3];
sx q[3];
rz(1.3501271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0423476) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(-2.919) q[2];
rz(-1.7706722) q[3];
sx q[3];
rz(-1.418908) q[3];
sx q[3];
rz(-2.5194949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.901578) q[0];
sx q[0];
rz(-2.1299025) q[0];
sx q[0];
rz(1.0255751) q[0];
rz(-1.1657794) q[1];
sx q[1];
rz(-2.5318601) q[1];
sx q[1];
rz(-0.21808521) q[1];
rz(0.72543421) q[2];
sx q[2];
rz(-0.67580961) q[2];
sx q[2];
rz(0.11563822) q[2];
rz(-1.5693657) q[3];
sx q[3];
rz(-0.32162255) q[3];
sx q[3];
rz(-2.376775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
