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
rz(0.043579276) q[1];
sx q[1];
rz(-1.2830623) q[1];
sx q[1];
rz(-2.9820774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.537764) q[0];
sx q[0];
rz(-1.7366752) q[0];
sx q[0];
rz(3.1127451) q[0];
rz(-pi) q[1];
rz(1.8279499) q[2];
sx q[2];
rz(-2.1669877) q[2];
sx q[2];
rz(1.330708) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.415787) q[1];
sx q[1];
rz(-1.09607) q[1];
sx q[1];
rz(-0.15523796) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2920424) q[3];
sx q[3];
rz(-1.2952943) q[3];
sx q[3];
rz(1.7003789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45304766) q[2];
sx q[2];
rz(-1.9828372) q[2];
sx q[2];
rz(3.0513406) q[2];
rz(0.074617535) q[3];
sx q[3];
rz(-1.4679694) q[3];
sx q[3];
rz(-0.41372764) q[3];
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
sx q[0];
rz(-pi) q[0];
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
rz(-2.3989035) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(0.76505605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30456736) q[0];
sx q[0];
rz(-1.6374987) q[0];
sx q[0];
rz(0.59356687) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68816784) q[2];
sx q[2];
rz(-2.151647) q[2];
sx q[2];
rz(0.95385393) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6366406) q[1];
sx q[1];
rz(-2.3206298) q[1];
sx q[1];
rz(2.0269462) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5935947) q[3];
sx q[3];
rz(-1.573085) q[3];
sx q[3];
rz(1.3399923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9146156) q[2];
sx q[2];
rz(-1.0116297) q[2];
sx q[2];
rz(1.6870618) q[2];
rz(2.4948273) q[3];
sx q[3];
rz(-0.55964595) q[3];
sx q[3];
rz(1.2796992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890559) q[0];
sx q[0];
rz(-1.4707969) q[0];
sx q[0];
rz(-0.045150969) q[0];
rz(-1.2308944) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(-2.0848354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74879828) q[0];
sx q[0];
rz(-2.4840889) q[0];
sx q[0];
rz(2.1211521) q[0];
rz(-0.62303876) q[2];
sx q[2];
rz(-2.3941561) q[2];
sx q[2];
rz(-1.7926271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18348098) q[1];
sx q[1];
rz(-2.2327802) q[1];
sx q[1];
rz(0.32420154) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8523434) q[3];
sx q[3];
rz(-2.4334638) q[3];
sx q[3];
rz(2.331108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7354273) q[2];
sx q[2];
rz(-2.1801517) q[2];
sx q[2];
rz(-0.65262922) q[2];
rz(-1.2930219) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80766135) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(1.3941143) q[0];
rz(-1.8380503) q[1];
sx q[1];
rz(-1.1640254) q[1];
sx q[1];
rz(-2.0533144) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9435642) q[0];
sx q[0];
rz(-1.8450415) q[0];
sx q[0];
rz(-2.345577) q[0];
rz(-pi) q[1];
rz(2.5456504) q[2];
sx q[2];
rz(-2.1998027) q[2];
sx q[2];
rz(0.12550917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5381319) q[1];
sx q[1];
rz(-2.1744121) q[1];
sx q[1];
rz(2.1674898) q[1];
rz(0.096403555) q[3];
sx q[3];
rz(-0.50130522) q[3];
sx q[3];
rz(-0.65602428) q[3];
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
rz(2.8754442) q[3];
sx q[3];
rz(-1.9726617) q[3];
sx q[3];
rz(2.2973785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(3.0323623) q[0];
sx q[0];
rz(-0.6210331) q[0];
sx q[0];
rz(-1.9367628) q[0];
rz(-3.0989528) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(2.1991594) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5980127) q[0];
sx q[0];
rz(-2.6930598) q[0];
sx q[0];
rz(1.2537812) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3187014) q[2];
sx q[2];
rz(-2.3510755) q[2];
sx q[2];
rz(-0.86871494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98676315) q[1];
sx q[1];
rz(-2.0663321) q[1];
sx q[1];
rz(-2.7096728) q[1];
x q[2];
rz(2.995184) q[3];
sx q[3];
rz(-0.64420986) q[3];
sx q[3];
rz(-0.25800426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0254536) q[2];
sx q[2];
rz(-0.68486539) q[2];
sx q[2];
rz(1.2241414) q[2];
rz(0.72743607) q[3];
sx q[3];
rz(-1.4801315) q[3];
sx q[3];
rz(3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8023476) q[0];
sx q[0];
rz(-2.6363063) q[0];
sx q[0];
rz(0.25830609) q[0];
rz(-2.2282265) q[1];
sx q[1];
rz(-2.2759627) q[1];
sx q[1];
rz(-0.9679274) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0845522) q[0];
sx q[0];
rz(-2.1294609) q[0];
sx q[0];
rz(-1.1112122) q[0];
rz(-pi) q[1];
rz(0.36138205) q[2];
sx q[2];
rz(-2.2602096) q[2];
sx q[2];
rz(2.9887082) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63686759) q[1];
sx q[1];
rz(-0.47885413) q[1];
sx q[1];
rz(0.07696159) q[1];
rz(-0.67575309) q[3];
sx q[3];
rz(-2.3738513) q[3];
sx q[3];
rz(-2.2217285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6496868) q[2];
sx q[2];
rz(-0.68444362) q[2];
sx q[2];
rz(2.7223132) q[2];
rz(-2.3573917) q[3];
sx q[3];
rz(-1.0201642) q[3];
sx q[3];
rz(-1.3721344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.075031) q[0];
sx q[0];
rz(-0.12396585) q[0];
sx q[0];
rz(-3.1228464) q[0];
rz(-3.0491414) q[1];
sx q[1];
rz(-1.3321313) q[1];
sx q[1];
rz(-0.094955347) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088647) q[0];
sx q[0];
rz(-1.0148541) q[0];
sx q[0];
rz(-1.9535319) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58231852) q[2];
sx q[2];
rz(-2.4597801) q[2];
sx q[2];
rz(-2.5043629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.590789) q[1];
sx q[1];
rz(-1.992464) q[1];
sx q[1];
rz(1.8487537) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9800013) q[3];
sx q[3];
rz(-1.9658024) q[3];
sx q[3];
rz(0.99564161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9065325) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(-2.2825784) q[2];
rz(-0.61339316) q[3];
sx q[3];
rz(-0.20039138) q[3];
sx q[3];
rz(-1.332351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61619806) q[0];
sx q[0];
rz(-1.8680251) q[0];
sx q[0];
rz(1.4816351) q[0];
rz(-2.0741277) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(-1.3190837) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5123915) q[0];
sx q[0];
rz(-0.90921697) q[0];
sx q[0];
rz(-0.11920905) q[0];
x q[1];
rz(-0.51609938) q[2];
sx q[2];
rz(-1.8987499) q[2];
sx q[2];
rz(-2.536377) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6665277) q[1];
sx q[1];
rz(-1.804978) q[1];
sx q[1];
rz(-0.98824309) q[1];
rz(-pi) q[2];
rz(-2.784801) q[3];
sx q[3];
rz(-1.4513375) q[3];
sx q[3];
rz(-2.2262163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2655502) q[2];
sx q[2];
rz(-0.37028131) q[2];
sx q[2];
rz(-3.1004356) q[2];
rz(-2.0187078) q[3];
sx q[3];
rz(-1.6582146) q[3];
sx q[3];
rz(-2.6179204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48162833) q[0];
sx q[0];
rz(-0.97844231) q[0];
sx q[0];
rz(-1.4725257) q[0];
rz(-0.72498471) q[1];
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
rz(-2.8003152) q[0];
sx q[0];
rz(-0.14320548) q[0];
rz(-pi) q[1];
rz(0.99733229) q[2];
sx q[2];
rz(-0.49066273) q[2];
sx q[2];
rz(3.1072679) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.37112999) q[1];
sx q[1];
rz(-1.5411301) q[1];
sx q[1];
rz(2.4231623) q[1];
x q[2];
rz(-1.7506041) q[3];
sx q[3];
rz(-2.4459908) q[3];
sx q[3];
rz(2.4814388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4943115) q[2];
sx q[2];
rz(-2.8936671) q[2];
sx q[2];
rz(-0.9606804) q[2];
rz(-2.658355) q[3];
sx q[3];
rz(-1.5471231) q[3];
sx q[3];
rz(1.5620935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305337) q[0];
sx q[0];
rz(-0.49376765) q[0];
sx q[0];
rz(-0.44179398) q[0];
rz(0.68253016) q[1];
sx q[1];
rz(-1.4492757) q[1];
sx q[1];
rz(1.0535047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58013142) q[0];
sx q[0];
rz(-0.86505689) q[0];
sx q[0];
rz(2.8012026) q[0];
rz(3.0526673) q[2];
sx q[2];
rz(-1.94238) q[2];
sx q[2];
rz(0.099612923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9113831) q[1];
sx q[1];
rz(-2.1911096) q[1];
sx q[1];
rz(-0.39834724) q[1];
x q[2];
rz(-2.9656762) q[3];
sx q[3];
rz(-2.7421167) q[3];
sx q[3];
rz(2.7586058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0992451) q[2];
sx q[2];
rz(-1.2688364) q[2];
sx q[2];
rz(-0.22259268) q[2];
rz(1.3709204) q[3];
sx q[3];
rz(-1.418908) q[3];
sx q[3];
rz(0.62209779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(0.54033242) q[2];
sx q[2];
rz(-1.1428383) q[2];
sx q[2];
rz(-2.0603772) q[2];
rz(-1.8924186) q[3];
sx q[3];
rz(-1.5712486) q[3];
sx q[3];
rz(-0.80733588) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
