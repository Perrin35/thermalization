OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3736149) q[0];
sx q[0];
rz(-0.64426214) q[0];
sx q[0];
rz(2.2777519) q[0];
rz(4.7148352) q[1];
sx q[1];
rz(5.0532053) q[1];
sx q[1];
rz(2.4911575) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7614131) q[0];
sx q[0];
rz(-0.54523796) q[0];
sx q[0];
rz(-2.7392967) q[0];
rz(-0.09133275) q[2];
sx q[2];
rz(-0.31158456) q[2];
sx q[2];
rz(-2.950719) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0155458) q[1];
sx q[1];
rz(-2.0728557) q[1];
sx q[1];
rz(2.253336) q[1];
rz(1.4635383) q[3];
sx q[3];
rz(-1.4441731) q[3];
sx q[3];
rz(-2.7789314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9170561) q[2];
sx q[2];
rz(-1.7660331) q[2];
sx q[2];
rz(-1.3517693) q[2];
rz(-0.81870643) q[3];
sx q[3];
rz(-1.6823781) q[3];
sx q[3];
rz(1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025892) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(-2.5658521) q[0];
rz(-0.88306824) q[1];
sx q[1];
rz(-1.6022316) q[1];
sx q[1];
rz(-1.3188804) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91399804) q[0];
sx q[0];
rz(-0.84216252) q[0];
sx q[0];
rz(-1.9084068) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1167745) q[2];
sx q[2];
rz(-0.78486004) q[2];
sx q[2];
rz(-1.0944686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0106338) q[1];
sx q[1];
rz(-1.9234227) q[1];
sx q[1];
rz(-2.3214798) q[1];
rz(-1.3926199) q[3];
sx q[3];
rz(-1.2548073) q[3];
sx q[3];
rz(-0.23366486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64112249) q[2];
sx q[2];
rz(-1.0176071) q[2];
sx q[2];
rz(-2.9024331) q[2];
rz(-0.33254361) q[3];
sx q[3];
rz(-1.4556689) q[3];
sx q[3];
rz(2.0502162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34514937) q[0];
sx q[0];
rz(-0.74757663) q[0];
sx q[0];
rz(-0.24459608) q[0];
rz(-2.7469475) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(1.2044725) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0255677) q[0];
sx q[0];
rz(-3.126156) q[0];
sx q[0];
rz(-0.77657338) q[0];
rz(-2.7543805) q[2];
sx q[2];
rz(-1.1128328) q[2];
sx q[2];
rz(2.2004623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2358346) q[1];
sx q[1];
rz(-1.8595962) q[1];
sx q[1];
rz(-2.1101084) q[1];
rz(1.0706903) q[3];
sx q[3];
rz(-1.2104038) q[3];
sx q[3];
rz(0.57756027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7977153) q[2];
sx q[2];
rz(-2.158973) q[2];
sx q[2];
rz(0.24125153) q[2];
rz(1.3681715) q[3];
sx q[3];
rz(-1.7521114) q[3];
sx q[3];
rz(0.97197896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57206804) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(-2.0401814) q[0];
rz(1.6481579) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(-1.000584) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7977622) q[0];
sx q[0];
rz(-1.5008108) q[0];
sx q[0];
rz(-1.153341) q[0];
rz(-2.1270833) q[2];
sx q[2];
rz(-1.5624701) q[2];
sx q[2];
rz(0.82288344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.31187181) q[1];
sx q[1];
rz(-1.5964566) q[1];
sx q[1];
rz(3.0405634) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8832753) q[3];
sx q[3];
rz(-0.44033209) q[3];
sx q[3];
rz(0.60763121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3558041) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(2.918952) q[2];
rz(-1.6629793) q[3];
sx q[3];
rz(-2.7361349) q[3];
sx q[3];
rz(1.2610029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3123689) q[0];
sx q[0];
rz(-1.1382505) q[0];
sx q[0];
rz(-2.9421222) q[0];
rz(1.624674) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(-2.6768501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366531) q[0];
sx q[0];
rz(-0.31855983) q[0];
sx q[0];
rz(-2.0037988) q[0];
rz(-0.22146341) q[2];
sx q[2];
rz(-1.5128947) q[2];
sx q[2];
rz(2.1277949) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6256959) q[1];
sx q[1];
rz(-0.40696884) q[1];
sx q[1];
rz(-0.76677983) q[1];
rz(-0.0031849234) q[3];
sx q[3];
rz(-2.4441458) q[3];
sx q[3];
rz(1.7980597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0988079) q[2];
sx q[2];
rz(-1.3036737) q[2];
sx q[2];
rz(0.25351563) q[2];
rz(-0.86067307) q[3];
sx q[3];
rz(-2.6200675) q[3];
sx q[3];
rz(1.1389987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(3.0694224) q[0];
sx q[0];
rz(-1.5266029) q[0];
sx q[0];
rz(1.3602863) q[0];
rz(2.5379429) q[1];
sx q[1];
rz(-2.344548) q[1];
sx q[1];
rz(2.6545677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384772) q[0];
sx q[0];
rz(-0.91775187) q[0];
sx q[0];
rz(1.6621291) q[0];
x q[1];
rz(-0.73871805) q[2];
sx q[2];
rz(-1.6411348) q[2];
sx q[2];
rz(2.3073151) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54283318) q[1];
sx q[1];
rz(-0.39879754) q[1];
sx q[1];
rz(2.0098575) q[1];
rz(1.9726874) q[3];
sx q[3];
rz(-1.3580267) q[3];
sx q[3];
rz(-0.63829845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6765678) q[2];
sx q[2];
rz(-2.1239231) q[2];
sx q[2];
rz(-0.49099311) q[2];
rz(2.6077014) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(-3.0965613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8462867) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(1.8701766) q[0];
rz(-2.2170587) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(-1.0251934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3640396) q[0];
sx q[0];
rz(-1.784783) q[0];
sx q[0];
rz(0.63978934) q[0];
x q[1];
rz(1.2751638) q[2];
sx q[2];
rz(-0.64648333) q[2];
sx q[2];
rz(-1.531284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89910903) q[1];
sx q[1];
rz(-0.61990737) q[1];
sx q[1];
rz(2.1371045) q[1];
rz(-pi) q[2];
rz(0.2853395) q[3];
sx q[3];
rz(-0.72721823) q[3];
sx q[3];
rz(0.83608999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2123432) q[2];
sx q[2];
rz(-2.847147) q[2];
sx q[2];
rz(2.219131) q[2];
rz(2.755002) q[3];
sx q[3];
rz(-1.2873193) q[3];
sx q[3];
rz(1.5132025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5492822) q[0];
sx q[0];
rz(-3.1070502) q[0];
sx q[0];
rz(2.4336245) q[0];
rz(0.51900807) q[1];
sx q[1];
rz(-2.612807) q[1];
sx q[1];
rz(2.4755898) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1493269) q[0];
sx q[0];
rz(-1.0368057) q[0];
sx q[0];
rz(-3.0342681) q[0];
rz(3.0254499) q[2];
sx q[2];
rz(-1.2525038) q[2];
sx q[2];
rz(-1.7081941) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.076326536) q[1];
sx q[1];
rz(-2.2398723) q[1];
sx q[1];
rz(-1.5764023) q[1];
rz(-pi) q[2];
rz(1.610393) q[3];
sx q[3];
rz(-2.5259113) q[3];
sx q[3];
rz(2.1510268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4697326) q[2];
sx q[2];
rz(-1.9817151) q[2];
sx q[2];
rz(-2.0812017) q[2];
rz(1.4276069) q[3];
sx q[3];
rz(-1.5771882) q[3];
sx q[3];
rz(0.76712999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0107182) q[0];
sx q[0];
rz(-0.81881443) q[0];
sx q[0];
rz(-0.70097104) q[0];
rz(2.2531807) q[1];
sx q[1];
rz(-2.7263434) q[1];
sx q[1];
rz(1.7078687) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2727391) q[0];
sx q[0];
rz(-2.0827385) q[0];
sx q[0];
rz(-1.4466132) q[0];
rz(-pi) q[1];
rz(-1.5436243) q[2];
sx q[2];
rz(-1.1650231) q[2];
sx q[2];
rz(2.9494065) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.10987416) q[1];
sx q[1];
rz(-1.5919135) q[1];
sx q[1];
rz(-1.5892522) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7319331) q[3];
sx q[3];
rz(-0.83935347) q[3];
sx q[3];
rz(1.5956772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1514757) q[2];
sx q[2];
rz(-1.2274123) q[2];
sx q[2];
rz(-1.6647476) q[2];
rz(-2.9448523) q[3];
sx q[3];
rz(-1.1367831) q[3];
sx q[3];
rz(-0.93973947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5340586) q[0];
sx q[0];
rz(-1.5628096) q[0];
sx q[0];
rz(1.1210572) q[0];
rz(2.6634482) q[1];
sx q[1];
rz(-0.68111626) q[1];
sx q[1];
rz(0.71472439) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5276565) q[0];
sx q[0];
rz(-0.52866565) q[0];
sx q[0];
rz(0.40724004) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98866971) q[2];
sx q[2];
rz(-0.63245648) q[2];
sx q[2];
rz(2.0159231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0352382) q[1];
sx q[1];
rz(-2.0788621) q[1];
sx q[1];
rz(-0.50915995) q[1];
x q[2];
rz(-2.5190936) q[3];
sx q[3];
rz(-0.9676515) q[3];
sx q[3];
rz(0.12605099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5152682) q[2];
sx q[2];
rz(-1.8147261) q[2];
sx q[2];
rz(-1.2088306) q[2];
rz(1.6995466) q[3];
sx q[3];
rz(-2.2550826) q[3];
sx q[3];
rz(-3.076021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47990738) q[0];
sx q[0];
rz(-1.8330782) q[0];
sx q[0];
rz(-0.080396419) q[0];
rz(-0.54795625) q[1];
sx q[1];
rz(-0.68956551) q[1];
sx q[1];
rz(-1.7908304) q[1];
rz(-1.1153658) q[2];
sx q[2];
rz(-2.6276988) q[2];
sx q[2];
rz(-2.217271) q[2];
rz(0.73765909) q[3];
sx q[3];
rz(-0.59288965) q[3];
sx q[3];
rz(2.2137143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
