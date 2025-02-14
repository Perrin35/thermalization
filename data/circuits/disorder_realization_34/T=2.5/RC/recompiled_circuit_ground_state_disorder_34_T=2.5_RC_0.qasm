OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.9095705) q[0];
sx q[0];
rz(2.4973305) q[0];
sx q[0];
rz(7.147026) q[0];
rz(-1.5683501) q[1];
sx q[1];
rz(-1.22998) q[1];
sx q[1];
rz(-0.65043515) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7614131) q[0];
sx q[0];
rz(-2.5963547) q[0];
sx q[0];
rz(-2.7392967) q[0];
rz(-pi) q[1];
rz(-0.09133275) q[2];
sx q[2];
rz(-0.31158456) q[2];
sx q[2];
rz(0.1908737) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0697223) q[1];
sx q[1];
rz(-0.98492565) q[1];
sx q[1];
rz(0.61572325) q[1];
rz(-3.0142455) q[3];
sx q[3];
rz(-1.4644002) q[3];
sx q[3];
rz(1.947054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2245366) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025892) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(2.5658521) q[0];
rz(2.2585244) q[1];
sx q[1];
rz(-1.539361) q[1];
sx q[1];
rz(-1.8227122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42717375) q[0];
sx q[0];
rz(-1.321081) q[0];
sx q[0];
rz(0.7575794) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3568866) q[2];
sx q[2];
rz(-1.5883351) q[2];
sx q[2];
rz(-0.49388805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0134351) q[1];
sx q[1];
rz(-0.87608209) q[1];
sx q[1];
rz(-0.46626587) q[1];
rz(-pi) q[2];
rz(-0.32073297) q[3];
sx q[3];
rz(-1.401529) q[3];
sx q[3];
rz(1.2812268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64112249) q[2];
sx q[2];
rz(-1.0176071) q[2];
sx q[2];
rz(2.9024331) q[2];
rz(2.809049) q[3];
sx q[3];
rz(-1.6859237) q[3];
sx q[3];
rz(-2.0502162) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7964433) q[0];
sx q[0];
rz(-2.394016) q[0];
sx q[0];
rz(0.24459608) q[0];
rz(2.7469475) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(1.9371202) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0255677) q[0];
sx q[0];
rz(-3.126156) q[0];
sx q[0];
rz(2.3650193) q[0];
rz(-pi) q[1];
rz(-0.38721217) q[2];
sx q[2];
rz(-1.1128328) q[2];
sx q[2];
rz(-2.2004623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.637801) q[1];
sx q[1];
rz(-2.0855013) q[1];
sx q[1];
rz(0.33333244) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0709023) q[3];
sx q[3];
rz(-1.9311889) q[3];
sx q[3];
rz(0.57756027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7977153) q[2];
sx q[2];
rz(-0.98261967) q[2];
sx q[2];
rz(-0.24125153) q[2];
rz(1.3681715) q[3];
sx q[3];
rz(-1.7521114) q[3];
sx q[3];
rz(0.97197896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5695246) q[0];
sx q[0];
rz(-1.854874) q[0];
sx q[0];
rz(-2.0401814) q[0];
rz(-1.4934348) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(-1.000584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9456317) q[0];
sx q[0];
rz(-1.9871662) q[0];
sx q[0];
rz(0.076535688) q[0];
rz(-pi) q[1];
rz(-1.0145093) q[2];
sx q[2];
rz(-1.5624701) q[2];
sx q[2];
rz(-0.82288344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1305589) q[1];
sx q[1];
rz(-3.0373664) q[1];
sx q[1];
rz(-2.8924045) q[1];
rz(-pi) q[2];
rz(1.9922929) q[3];
sx q[3];
rz(-1.7022082) q[3];
sx q[3];
rz(-0.67883713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3558041) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(-0.22264063) q[2];
rz(-1.6629793) q[3];
sx q[3];
rz(-2.7361349) q[3];
sx q[3];
rz(1.2610029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.8292238) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(0.19947048) q[0];
rz(-1.5169187) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(0.46474251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6049395) q[0];
sx q[0];
rz(-2.8230328) q[0];
sx q[0];
rz(-1.1377938) q[0];
x q[1];
rz(2.88358) q[2];
sx q[2];
rz(-2.9128053) q[2];
sx q[2];
rz(2.8361671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3252651) q[1];
sx q[1];
rz(-1.8598598) q[1];
sx q[1];
rz(1.2802034) q[1];
rz(-1.5681276) q[3];
sx q[3];
rz(-0.87335372) q[3];
sx q[3];
rz(1.3476882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.042784721) q[2];
sx q[2];
rz(-1.837919) q[2];
sx q[2];
rz(-0.25351563) q[2];
rz(0.86067307) q[3];
sx q[3];
rz(-0.52152514) q[3];
sx q[3];
rz(-2.0025939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072170243) q[0];
sx q[0];
rz(-1.6149898) q[0];
sx q[0];
rz(1.3602863) q[0];
rz(2.5379429) q[1];
sx q[1];
rz(-0.79704469) q[1];
sx q[1];
rz(0.48702494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2880858) q[0];
sx q[0];
rz(-2.4831193) q[0];
sx q[0];
rz(-0.11866026) q[0];
rz(-pi) q[1];
rz(-1.4757882) q[2];
sx q[2];
rz(-2.3072647) q[2];
sx q[2];
rz(-2.468994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54283318) q[1];
sx q[1];
rz(-2.7427951) q[1];
sx q[1];
rz(-2.0098575) q[1];
x q[2];
rz(-1.9726874) q[3];
sx q[3];
rz(-1.783566) q[3];
sx q[3];
rz(-0.63829845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46502486) q[2];
sx q[2];
rz(-2.1239231) q[2];
sx q[2];
rz(0.49099311) q[2];
rz(2.6077014) q[3];
sx q[3];
rz(-2.5940354) q[3];
sx q[3];
rz(3.0965613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.295306) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(1.2714161) q[0];
rz(0.92453399) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(2.1163993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071200095) q[0];
sx q[0];
rz(-2.4717583) q[0];
sx q[0];
rz(-0.34908648) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2164257) q[2];
sx q[2];
rz(-2.1849298) q[2];
sx q[2];
rz(-1.9748109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5793663) q[1];
sx q[1];
rz(-2.0831897) q[1];
sx q[1];
rz(2.775869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8562532) q[3];
sx q[3];
rz(-2.4143744) q[3];
sx q[3];
rz(-0.83608999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9292494) q[2];
sx q[2];
rz(-2.847147) q[2];
sx q[2];
rz(-2.219131) q[2];
rz(-2.755002) q[3];
sx q[3];
rz(-1.8542733) q[3];
sx q[3];
rz(-1.6283901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5492822) q[0];
sx q[0];
rz(-0.03454241) q[0];
sx q[0];
rz(0.70796815) q[0];
rz(0.51900807) q[1];
sx q[1];
rz(-2.612807) q[1];
sx q[1];
rz(-0.66600287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7836664) q[0];
sx q[0];
rz(-0.54364288) q[0];
sx q[0];
rz(-1.3915791) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2504891) q[2];
sx q[2];
rz(-1.6810809) q[2];
sx q[2];
rz(-0.10090339) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6506001) q[1];
sx q[1];
rz(-1.5751936) q[1];
sx q[1];
rz(-0.66908361) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95548463) q[3];
sx q[3];
rz(-1.59366) q[3];
sx q[3];
rz(2.5290306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6718601) q[2];
sx q[2];
rz(-1.1598776) q[2];
sx q[2];
rz(-1.0603909) q[2];
rz(1.7139858) q[3];
sx q[3];
rz(-1.5644045) q[3];
sx q[3];
rz(-2.3744627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1308744) q[0];
sx q[0];
rz(-0.81881443) q[0];
sx q[0];
rz(0.70097104) q[0];
rz(0.88841191) q[1];
sx q[1];
rz(-0.41524926) q[1];
sx q[1];
rz(-1.433724) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86885356) q[0];
sx q[0];
rz(-1.0588542) q[0];
sx q[0];
rz(1.6949795) q[0];
rz(1.5979684) q[2];
sx q[2];
rz(-1.1650231) q[2];
sx q[2];
rz(2.9494065) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4605324) q[1];
sx q[1];
rz(-1.5523445) q[1];
sx q[1];
rz(-0.021120763) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9307487) q[3];
sx q[3];
rz(-2.1574855) q[3];
sx q[3];
rz(0.66419475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1514757) q[2];
sx q[2];
rz(-1.2274123) q[2];
sx q[2];
rz(-1.6647476) q[2];
rz(0.19674033) q[3];
sx q[3];
rz(-2.0048095) q[3];
sx q[3];
rz(-2.2018532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
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
rz(-2.4604764) q[1];
sx q[1];
rz(-0.71472439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1506468) q[0];
sx q[0];
rz(-2.0523221) q[0];
sx q[0];
rz(1.343438) q[0];
rz(2.1529229) q[2];
sx q[2];
rz(-2.5091362) q[2];
sx q[2];
rz(-1.1256696) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27036941) q[1];
sx q[1];
rz(-2.0107234) q[1];
sx q[1];
rz(-1.0031071) q[1];
rz(0.62249903) q[3];
sx q[3];
rz(-0.9676515) q[3];
sx q[3];
rz(-3.0155417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.62632442) q[2];
sx q[2];
rz(-1.3268665) q[2];
sx q[2];
rz(-1.2088306) q[2];
rz(1.6995466) q[3];
sx q[3];
rz(-2.2550826) q[3];
sx q[3];
rz(0.065571688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6616853) q[0];
sx q[0];
rz(-1.8330782) q[0];
sx q[0];
rz(-0.080396419) q[0];
rz(2.5936364) q[1];
sx q[1];
rz(-0.68956551) q[1];
sx q[1];
rz(-1.7908304) q[1];
rz(-1.1153658) q[2];
sx q[2];
rz(-2.6276988) q[2];
sx q[2];
rz(-2.217271) q[2];
rz(-2.4039336) q[3];
sx q[3];
rz(-0.59288965) q[3];
sx q[3];
rz(2.2137143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
