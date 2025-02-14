OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.206267) q[0];
sx q[0];
rz(-1.6214108) q[0];
sx q[0];
rz(-1.7229463) q[0];
rz(-3.0980134) q[1];
sx q[1];
rz(-1.8585304) q[1];
sx q[1];
rz(-0.15951523) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60382868) q[0];
sx q[0];
rz(-1.4049174) q[0];
sx q[0];
rz(-3.1127451) q[0];
rz(-pi) q[1];
rz(2.5297727) q[2];
sx q[2];
rz(-1.7828336) q[2];
sx q[2];
rz(-0.093487948) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.91640305) q[1];
sx q[1];
rz(-1.708751) q[1];
sx q[1];
rz(2.0504584) q[1];
x q[2];
rz(2.2920424) q[3];
sx q[3];
rz(-1.8462984) q[3];
sx q[3];
rz(-1.7003789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45304766) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(0.090252074) q[2];
rz(-0.074617535) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(-2.727865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8892882) q[0];
sx q[0];
rz(-1.1587208) q[0];
sx q[0];
rz(-1.4537551) q[0];
rz(-0.74268913) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(-0.76505605) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1677645) q[0];
sx q[0];
rz(-2.544738) q[0];
sx q[0];
rz(3.0227227) q[0];
x q[1];
rz(2.3397331) q[2];
sx q[2];
rz(-0.86881402) q[2];
sx q[2];
rz(-3.1131256) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7530427) q[1];
sx q[1];
rz(-1.8990128) q[1];
sx q[1];
rz(-0.80372827) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0043930769) q[3];
sx q[3];
rz(-2.5935904) q[3];
sx q[3];
rz(0.22705423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9146156) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(1.6870618) q[2];
rz(0.64676532) q[3];
sx q[3];
rz(-0.55964595) q[3];
sx q[3];
rz(1.8618934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890559) q[0];
sx q[0];
rz(-1.6707957) q[0];
sx q[0];
rz(-3.0964417) q[0];
rz(-1.9106983) q[1];
sx q[1];
rz(-1.3094333) q[1];
sx q[1];
rz(-2.0848354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089356747) q[0];
sx q[0];
rz(-1.0228906) q[0];
sx q[0];
rz(2.7578081) q[0];
rz(-2.0665523) q[2];
sx q[2];
rz(-0.98598749) q[2];
sx q[2];
rz(-2.5677447) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31683982) q[1];
sx q[1];
rz(-0.72622967) q[1];
sx q[1];
rz(-1.1827501) q[1];
rz(-pi) q[2];
rz(2.908024) q[3];
sx q[3];
rz(-0.89591494) q[3];
sx q[3];
rz(-2.6949331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4061654) q[2];
sx q[2];
rz(-0.96144095) q[2];
sx q[2];
rz(0.65262922) q[2];
rz(1.2930219) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(-3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3339313) q[0];
sx q[0];
rz(-2.6787391) q[0];
sx q[0];
rz(-1.3941143) q[0];
rz(-1.8380503) q[1];
sx q[1];
rz(-1.9775672) q[1];
sx q[1];
rz(-1.0882783) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1980285) q[0];
sx q[0];
rz(-1.2965512) q[0];
sx q[0];
rz(2.345577) q[0];
x q[1];
rz(2.5456504) q[2];
sx q[2];
rz(-0.94178994) q[2];
sx q[2];
rz(3.0160835) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6034607) q[1];
sx q[1];
rz(-2.1744121) q[1];
sx q[1];
rz(0.97410283) q[1];
x q[2];
rz(2.6422464) q[3];
sx q[3];
rz(-1.6170698) q[3];
sx q[3];
rz(0.83016992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3950222) q[2];
sx q[2];
rz(-1.7638548) q[2];
sx q[2];
rz(-2.6611888) q[2];
rz(-0.26614842) q[3];
sx q[3];
rz(-1.9726617) q[3];
sx q[3];
rz(2.2973785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10923037) q[0];
sx q[0];
rz(-2.5205595) q[0];
sx q[0];
rz(1.2048298) q[0];
rz(3.0989528) q[1];
sx q[1];
rz(-1.7667021) q[1];
sx q[1];
rz(2.1991594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.248837) q[0];
sx q[0];
rz(-1.1461094) q[0];
sx q[0];
rz(-0.14890944) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3187014) q[2];
sx q[2];
rz(-2.3510755) q[2];
sx q[2];
rz(-2.2728777) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.79980274) q[1];
sx q[1];
rz(-1.9479349) q[1];
sx q[1];
rz(2.1076635) q[1];
rz(-pi) q[2];
rz(-0.63905255) q[3];
sx q[3];
rz(-1.6585232) q[3];
sx q[3];
rz(-1.4301585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11613906) q[2];
sx q[2];
rz(-0.68486539) q[2];
sx q[2];
rz(-1.2241414) q[2];
rz(2.4141566) q[3];
sx q[3];
rz(-1.4801315) q[3];
sx q[3];
rz(-3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8023476) q[0];
sx q[0];
rz(-0.50528637) q[0];
sx q[0];
rz(-0.25830609) q[0];
rz(0.91336617) q[1];
sx q[1];
rz(-0.86562997) q[1];
sx q[1];
rz(-2.1736653) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3334078) q[0];
sx q[0];
rz(-0.70751941) q[0];
sx q[0];
rz(0.61715166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2931678) q[2];
sx q[2];
rz(-1.8471187) q[2];
sx q[2];
rz(1.6538617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63686759) q[1];
sx q[1];
rz(-2.6627385) q[1];
sx q[1];
rz(3.0646311) q[1];
rz(0.67575309) q[3];
sx q[3];
rz(-2.3738513) q[3];
sx q[3];
rz(-0.91986419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4919058) q[2];
sx q[2];
rz(-0.68444362) q[2];
sx q[2];
rz(-2.7223132) q[2];
rz(2.3573917) q[3];
sx q[3];
rz(-2.1214285) q[3];
sx q[3];
rz(-1.3721344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066561617) q[0];
sx q[0];
rz(-3.0176268) q[0];
sx q[0];
rz(-3.1228464) q[0];
rz(0.092451267) q[1];
sx q[1];
rz(-1.8094614) q[1];
sx q[1];
rz(-3.0466373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92871794) q[0];
sx q[0];
rz(-1.2480019) q[0];
sx q[0];
rz(-0.59014894) q[0];
rz(2.5592741) q[2];
sx q[2];
rz(-0.68181255) q[2];
sx q[2];
rz(-2.5043629) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5508037) q[1];
sx q[1];
rz(-1.992464) q[1];
sx q[1];
rz(-1.292839) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9800013) q[3];
sx q[3];
rz(-1.1757903) q[3];
sx q[3];
rz(2.145951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9065325) q[2];
sx q[2];
rz(-1.9573213) q[2];
sx q[2];
rz(0.85901421) q[2];
rz(-2.5281995) q[3];
sx q[3];
rz(-2.9412013) q[3];
sx q[3];
rz(-1.332351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61619806) q[0];
sx q[0];
rz(-1.2735676) q[0];
sx q[0];
rz(1.4816351) q[0];
rz(2.0741277) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(-1.8225089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5123915) q[0];
sx q[0];
rz(-0.90921697) q[0];
sx q[0];
rz(-3.0223836) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60361344) q[2];
sx q[2];
rz(-0.60342646) q[2];
sx q[2];
rz(-0.44925424) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7072152) q[1];
sx q[1];
rz(-0.62276749) q[1];
sx q[1];
rz(1.9799443) q[1];
rz(2.784801) q[3];
sx q[3];
rz(-1.4513375) q[3];
sx q[3];
rz(-0.91537634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2655502) q[2];
sx q[2];
rz(-2.7713113) q[2];
sx q[2];
rz(-3.1004356) q[2];
rz(2.0187078) q[3];
sx q[3];
rz(-1.4833781) q[3];
sx q[3];
rz(0.52367228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48162833) q[0];
sx q[0];
rz(-0.97844231) q[0];
sx q[0];
rz(-1.4725257) q[0];
rz(-2.4166079) q[1];
sx q[1];
rz(-1.1712733) q[1];
sx q[1];
rz(-2.5670126) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6642484) q[0];
sx q[0];
rz(-1.2331523) q[0];
sx q[0];
rz(1.6214423) q[0];
rz(-1.1489611) q[2];
sx q[2];
rz(-1.8293194) q[2];
sx q[2];
rz(2.0542893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2255937) q[1];
sx q[1];
rz(-2.288842) q[1];
sx q[1];
rz(1.6101933) q[1];
rz(-pi) q[2];
rz(-2.2584142) q[3];
sx q[3];
rz(-1.6856582) q[3];
sx q[3];
rz(-2.3695947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64728111) q[2];
sx q[2];
rz(-0.24792555) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.1305337) q[0];
sx q[0];
rz(-0.49376765) q[0];
sx q[0];
rz(-2.6997987) q[0];
rz(-0.68253016) q[1];
sx q[1];
rz(-1.4492757) q[1];
sx q[1];
rz(2.0880879) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2164648) q[0];
sx q[0];
rz(-1.3138714) q[0];
sx q[0];
rz(-0.83576699) q[0];
rz(-pi) q[1];
rz(1.7948661) q[2];
sx q[2];
rz(-2.7599979) q[2];
sx q[2];
rz(-0.14118186) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5610841) q[1];
sx q[1];
rz(-1.891916) q[1];
sx q[1];
rz(-2.2301484) q[1];
x q[2];
rz(2.7476596) q[3];
sx q[3];
rz(-1.6389168) q[3];
sx q[3];
rz(-1.3501271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0992451) q[2];
sx q[2];
rz(-1.8727563) q[2];
sx q[2];
rz(-0.22259268) q[2];
rz(1.7706722) q[3];
sx q[3];
rz(-1.418908) q[3];
sx q[3];
rz(2.5194949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2400146) q[0];
sx q[0];
rz(-1.0116901) q[0];
sx q[0];
rz(-2.1160175) q[0];
rz(1.9758132) q[1];
sx q[1];
rz(-2.5318601) q[1];
sx q[1];
rz(-0.21808521) q[1];
rz(2.6012602) q[2];
sx q[2];
rz(-1.9987543) q[2];
sx q[2];
rz(1.0812154) q[2];
rz(-1.5722269) q[3];
sx q[3];
rz(-2.8199701) q[3];
sx q[3];
rz(0.7648177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
