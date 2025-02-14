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
rz(-0.86384073) q[0];
rz(4.7148352) q[1];
sx q[1];
rz(5.0532053) q[1];
sx q[1];
rz(2.4911575) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38017958) q[0];
sx q[0];
rz(-0.54523796) q[0];
sx q[0];
rz(0.40229599) q[0];
x q[1];
rz(-3.0502599) q[2];
sx q[2];
rz(-2.8300081) q[2];
sx q[2];
rz(-2.950719) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0155458) q[1];
sx q[1];
rz(-2.0728557) q[1];
sx q[1];
rz(0.8882567) q[1];
x q[2];
rz(-1.6780544) q[3];
sx q[3];
rz(-1.4441731) q[3];
sx q[3];
rz(-2.7789314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9170561) q[2];
sx q[2];
rz(-1.3755596) q[2];
sx q[2];
rz(1.7898233) q[2];
rz(-0.81870643) q[3];
sx q[3];
rz(-1.6823781) q[3];
sx q[3];
rz(1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025892) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(-0.57574058) q[0];
rz(-0.88306824) q[1];
sx q[1];
rz(-1.539361) q[1];
sx q[1];
rz(-1.8227122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7144189) q[0];
sx q[0];
rz(-1.321081) q[0];
sx q[0];
rz(0.7575794) q[0];
x q[1];
rz(-2.3568866) q[2];
sx q[2];
rz(-1.5883351) q[2];
sx q[2];
rz(0.49388805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3471862) q[1];
sx q[1];
rz(-0.81450317) q[1];
sx q[1];
rz(-2.0655355) q[1];
x q[2];
rz(2.8208597) q[3];
sx q[3];
rz(-1.7400636) q[3];
sx q[3];
rz(-1.2812268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64112249) q[2];
sx q[2];
rz(-2.1239855) q[2];
sx q[2];
rz(-0.23915954) q[2];
rz(-2.809049) q[3];
sx q[3];
rz(-1.4556689) q[3];
sx q[3];
rz(1.0913764) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964433) q[0];
sx q[0];
rz(-0.74757663) q[0];
sx q[0];
rz(2.8969966) q[0];
rz(0.39464513) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(-1.9371202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0255677) q[0];
sx q[0];
rz(-0.015436643) q[0];
sx q[0];
rz(0.77657338) q[0];
x q[1];
rz(-2.2245046) q[2];
sx q[2];
rz(-2.5508891) q[2];
sx q[2];
rz(-1.4554254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50379163) q[1];
sx q[1];
rz(-2.0855013) q[1];
sx q[1];
rz(0.33333244) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0709023) q[3];
sx q[3];
rz(-1.9311889) q[3];
sx q[3];
rz(2.5640324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7977153) q[2];
sx q[2];
rz(-0.98261967) q[2];
sx q[2];
rz(0.24125153) q[2];
rz(1.7734211) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(-2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57206804) q[0];
sx q[0];
rz(-1.2867186) q[0];
sx q[0];
rz(1.1014112) q[0];
rz(-1.6481579) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(-2.1410087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9456317) q[0];
sx q[0];
rz(-1.1544265) q[0];
sx q[0];
rz(-3.065057) q[0];
x q[1];
rz(-1.0145093) q[2];
sx q[2];
rz(-1.5791225) q[2];
sx q[2];
rz(0.82288344) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1305589) q[1];
sx q[1];
rz(-3.0373664) q[1];
sx q[1];
rz(2.8924045) q[1];
x q[2];
rz(1.9922929) q[3];
sx q[3];
rz(-1.4393844) q[3];
sx q[3];
rz(0.67883713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3558041) q[2];
sx q[2];
rz(-1.4996837) q[2];
sx q[2];
rz(2.918952) q[2];
rz(-1.4786134) q[3];
sx q[3];
rz(-2.7361349) q[3];
sx q[3];
rz(1.8805898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3123689) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(-2.9421222) q[0];
rz(-1.624674) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(-0.46474251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7620649) q[0];
sx q[0];
rz(-1.4389973) q[0];
sx q[0];
rz(-1.2799311) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9201292) q[2];
sx q[2];
rz(-1.5128947) q[2];
sx q[2];
rz(2.1277949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4721663) q[1];
sx q[1];
rz(-1.8490044) q[1];
sx q[1];
rz(0.30097715) q[1];
x q[2];
rz(0.69744436) q[3];
sx q[3];
rz(-1.5728419) q[3];
sx q[3];
rz(-0.22482219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.042784721) q[2];
sx q[2];
rz(-1.837919) q[2];
sx q[2];
rz(-0.25351563) q[2];
rz(0.86067307) q[3];
sx q[3];
rz(-2.6200675) q[3];
sx q[3];
rz(-1.1389987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.48702494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384772) q[0];
sx q[0];
rz(-0.91775187) q[0];
sx q[0];
rz(-1.4794635) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4757882) q[2];
sx q[2];
rz(-2.3072647) q[2];
sx q[2];
rz(-2.468994) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.071515171) q[1];
sx q[1];
rz(-1.9299475) q[1];
sx q[1];
rz(0.17724322) q[1];
rz(-pi) q[2];
rz(0.2305687) q[3];
sx q[3];
rz(-1.9631223) q[3];
sx q[3];
rz(-2.2986064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46502486) q[2];
sx q[2];
rz(-1.0176696) q[2];
sx q[2];
rz(0.49099311) q[2];
rz(2.6077014) q[3];
sx q[3];
rz(-2.5940354) q[3];
sx q[3];
rz(-0.045031358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.295306) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(-1.2714161) q[0];
rz(0.92453399) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(2.1163993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7775531) q[0];
sx q[0];
rz(-1.784783) q[0];
sx q[0];
rz(-2.5018033) q[0];
rz(-pi) q[1];
rz(1.2751638) q[2];
sx q[2];
rz(-0.64648333) q[2];
sx q[2];
rz(-1.531284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89910903) q[1];
sx q[1];
rz(-0.61990737) q[1];
sx q[1];
rz(1.0044881) q[1];
rz(-2.4348172) q[3];
sx q[3];
rz(-1.7590344) q[3];
sx q[3];
rz(-2.1911603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2123432) q[2];
sx q[2];
rz(-0.29444567) q[2];
sx q[2];
rz(-0.92246169) q[2];
rz(-2.755002) q[3];
sx q[3];
rz(-1.2873193) q[3];
sx q[3];
rz(1.6283901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59231049) q[0];
sx q[0];
rz(-0.03454241) q[0];
sx q[0];
rz(-0.70796815) q[0];
rz(0.51900807) q[1];
sx q[1];
rz(-2.612807) q[1];
sx q[1];
rz(2.4755898) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7749044) q[0];
sx q[0];
rz(-1.4784593) q[0];
sx q[0];
rz(-2.1073186) q[0];
rz(1.2504891) q[2];
sx q[2];
rz(-1.4605117) q[2];
sx q[2];
rz(-3.0406893) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0652661) q[1];
sx q[1];
rz(-2.2398723) q[1];
sx q[1];
rz(1.5764023) q[1];
rz(-pi) q[2];
x q[2];
rz(2.186108) q[3];
sx q[3];
rz(-1.5479327) q[3];
sx q[3];
rz(-0.61256204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4697326) q[2];
sx q[2];
rz(-1.1598776) q[2];
sx q[2];
rz(1.0603909) q[2];
rz(1.7139858) q[3];
sx q[3];
rz(-1.5644045) q[3];
sx q[3];
rz(0.76712999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1308744) q[0];
sx q[0];
rz(-2.3227782) q[0];
sx q[0];
rz(0.70097104) q[0];
rz(-2.2531807) q[1];
sx q[1];
rz(-2.7263434) q[1];
sx q[1];
rz(1.433724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2727391) q[0];
sx q[0];
rz(-2.0827385) q[0];
sx q[0];
rz(-1.4466132) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5979684) q[2];
sx q[2];
rz(-1.1650231) q[2];
sx q[2];
rz(0.19218615) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82821694) q[1];
sx q[1];
rz(-3.113548) q[1];
sx q[1];
rz(-0.7181479) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5340586) q[0];
sx q[0];
rz(-1.5787831) q[0];
sx q[0];
rz(1.1210572) q[0];
rz(2.6634482) q[1];
sx q[1];
rz(-0.68111626) q[1];
sx q[1];
rz(-2.4268683) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5276565) q[0];
sx q[0];
rz(-2.612927) q[0];
sx q[0];
rz(-0.40724004) q[0];
rz(-2.7585538) q[2];
sx q[2];
rz(-2.0872119) q[2];
sx q[2];
rz(0.4412152) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1063544) q[1];
sx q[1];
rz(-2.0788621) q[1];
sx q[1];
rz(2.6324327) q[1];
x q[2];
rz(-0.86758763) q[3];
sx q[3];
rz(-1.0699268) q[3];
sx q[3];
rz(-1.3102368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62632442) q[2];
sx q[2];
rz(-1.8147261) q[2];
sx q[2];
rz(-1.2088306) q[2];
rz(1.442046) q[3];
sx q[3];
rz(-0.88651005) q[3];
sx q[3];
rz(-3.076021) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6616853) q[0];
sx q[0];
rz(-1.3085145) q[0];
sx q[0];
rz(3.0611962) q[0];
rz(0.54795625) q[1];
sx q[1];
rz(-2.4520271) q[1];
sx q[1];
rz(1.3507623) q[1];
rz(-2.0399848) q[2];
sx q[2];
rz(-1.3528578) q[2];
sx q[2];
rz(2.8982671) q[2];
rz(-2.6790621) q[3];
sx q[3];
rz(-1.1855385) q[3];
sx q[3];
rz(-0.0029467646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
