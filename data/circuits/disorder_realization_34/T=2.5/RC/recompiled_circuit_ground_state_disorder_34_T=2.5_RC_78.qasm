OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.7679778) q[0];
sx q[0];
rz(-2.4973305) q[0];
sx q[0];
rz(0.86384073) q[0];
rz(-1.5683501) q[1];
sx q[1];
rz(-1.22998) q[1];
sx q[1];
rz(-0.65043515) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5395376) q[0];
sx q[0];
rz(-1.7752753) q[0];
sx q[0];
rz(-2.6325167) q[0];
rz(-pi) q[1];
x q[1];
rz(1.600163) q[2];
sx q[2];
rz(-1.8810388) q[2];
sx q[2];
rz(-0.28679906) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9792773) q[1];
sx q[1];
rz(-2.3188496) q[1];
sx q[1];
rz(2.2869799) q[1];
rz(-pi) q[2];
rz(-0.1273472) q[3];
sx q[3];
rz(-1.6771924) q[3];
sx q[3];
rz(-1.1945386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2245366) q[2];
sx q[2];
rz(-1.7660331) q[2];
sx q[2];
rz(-1.7898233) q[2];
rz(2.3228862) q[3];
sx q[3];
rz(-1.6823781) q[3];
sx q[3];
rz(1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8390035) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(-0.57574058) q[0];
rz(-2.2585244) q[1];
sx q[1];
rz(-1.539361) q[1];
sx q[1];
rz(-1.3188804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42717375) q[0];
sx q[0];
rz(-1.321081) q[0];
sx q[0];
rz(2.3840133) q[0];
rz(-pi) q[1];
rz(-1.5460125) q[2];
sx q[2];
rz(-0.78624386) q[2];
sx q[2];
rz(-1.0593965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0134351) q[1];
sx q[1];
rz(-0.87608209) q[1];
sx q[1];
rz(0.46626587) q[1];
x q[2];
rz(2.6448375) q[3];
sx q[3];
rz(-2.7803034) q[3];
sx q[3];
rz(0.75888097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64112249) q[2];
sx q[2];
rz(-1.0176071) q[2];
sx q[2];
rz(0.23915954) q[2];
rz(-2.809049) q[3];
sx q[3];
rz(-1.6859237) q[3];
sx q[3];
rz(2.0502162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34514937) q[0];
sx q[0];
rz(-0.74757663) q[0];
sx q[0];
rz(2.8969966) q[0];
rz(0.39464513) q[1];
sx q[1];
rz(-0.60401812) q[1];
sx q[1];
rz(-1.9371202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.116025) q[0];
sx q[0];
rz(-3.126156) q[0];
sx q[0];
rz(-0.77657338) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38721217) q[2];
sx q[2];
rz(-2.0287598) q[2];
sx q[2];
rz(2.2004623) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10919103) q[1];
sx q[1];
rz(-0.60495678) q[1];
sx q[1];
rz(-2.095286) q[1];
rz(-pi) q[2];
rz(-2.0709023) q[3];
sx q[3];
rz(-1.2104038) q[3];
sx q[3];
rz(0.57756027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.34387732) q[2];
sx q[2];
rz(-2.158973) q[2];
sx q[2];
rz(0.24125153) q[2];
rz(-1.7734211) q[3];
sx q[3];
rz(-1.7521114) q[3];
sx q[3];
rz(0.97197896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5695246) q[0];
sx q[0];
rz(-1.2867186) q[0];
sx q[0];
rz(-1.1014112) q[0];
rz(1.6481579) q[1];
sx q[1];
rz(-2.8547574) q[1];
sx q[1];
rz(-2.1410087) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1959609) q[0];
sx q[0];
rz(-1.9871662) q[0];
sx q[0];
rz(0.076535688) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5865636) q[2];
sx q[2];
rz(-2.5852499) q[2];
sx q[2];
rz(-0.76130262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0110338) q[1];
sx q[1];
rz(-3.0373664) q[1];
sx q[1];
rz(-0.24918814) q[1];
x q[2];
rz(1.2583174) q[3];
sx q[3];
rz(-0.44033209) q[3];
sx q[3];
rz(-2.5339614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7857886) q[2];
sx q[2];
rz(-1.641909) q[2];
sx q[2];
rz(2.918952) q[2];
rz(1.4786134) q[3];
sx q[3];
rz(-0.40545774) q[3];
sx q[3];
rz(1.8805898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8292238) q[0];
sx q[0];
rz(-1.1382505) q[0];
sx q[0];
rz(-2.9421222) q[0];
rz(-1.5169187) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(-2.6768501) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7620649) q[0];
sx q[0];
rz(-1.4389973) q[0];
sx q[0];
rz(1.2799311) q[0];
rz(-pi) q[1];
rz(0.25801261) q[2];
sx q[2];
rz(-2.9128053) q[2];
sx q[2];
rz(0.30542557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5158968) q[1];
sx q[1];
rz(-0.40696884) q[1];
sx q[1];
rz(2.3748128) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0031849234) q[3];
sx q[3];
rz(-2.4441458) q[3];
sx q[3];
rz(-1.3435329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0988079) q[2];
sx q[2];
rz(-1.3036737) q[2];
sx q[2];
rz(-0.25351563) q[2];
rz(-0.86067307) q[3];
sx q[3];
rz(-2.6200675) q[3];
sx q[3];
rz(-2.0025939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072170243) q[0];
sx q[0];
rz(-1.6149898) q[0];
sx q[0];
rz(-1.3602863) q[0];
rz(0.6036497) q[1];
sx q[1];
rz(-0.79704469) q[1];
sx q[1];
rz(-0.48702494) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2880858) q[0];
sx q[0];
rz(-0.65847337) q[0];
sx q[0];
rz(-3.0229324) q[0];
x q[1];
rz(-1.6658044) q[2];
sx q[2];
rz(-2.3072647) q[2];
sx q[2];
rz(2.468994) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54283318) q[1];
sx q[1];
rz(-2.7427951) q[1];
sx q[1];
rz(-1.1317352) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.911024) q[3];
sx q[3];
rz(-1.9631223) q[3];
sx q[3];
rz(0.84298625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6765678) q[2];
sx q[2];
rz(-1.0176696) q[2];
sx q[2];
rz(0.49099311) q[2];
rz(-0.53389126) q[3];
sx q[3];
rz(-0.54755727) q[3];
sx q[3];
rz(0.045031358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8462867) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(-1.8701766) q[0];
rz(-2.2170587) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(2.1163993) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7780925) q[0];
sx q[0];
rz(-0.94787593) q[0];
sx q[0];
rz(1.3062551) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94550262) q[2];
sx q[2];
rz(-1.3943814) q[2];
sx q[2];
rz(2.8636065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2424836) q[1];
sx q[1];
rz(-2.5216853) q[1];
sx q[1];
rz(2.1371045) q[1];
rz(-pi) q[2];
rz(1.8162433) q[3];
sx q[3];
rz(-2.2625828) q[3];
sx q[3];
rz(0.46190767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2123432) q[2];
sx q[2];
rz(-2.847147) q[2];
sx q[2];
rz(0.92246169) q[2];
rz(-0.38659066) q[3];
sx q[3];
rz(-1.2873193) q[3];
sx q[3];
rz(1.5132025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59231049) q[0];
sx q[0];
rz(-0.03454241) q[0];
sx q[0];
rz(2.4336245) q[0];
rz(-2.6225846) q[1];
sx q[1];
rz(-0.52878562) q[1];
sx q[1];
rz(0.66600287) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749044) q[0];
sx q[0];
rz(-1.4784593) q[0];
sx q[0];
rz(-1.0342741) q[0];
rz(-1.8911036) q[2];
sx q[2];
rz(-1.6810809) q[2];
sx q[2];
rz(-0.10090339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0652661) q[1];
sx q[1];
rz(-2.2398723) q[1];
sx q[1];
rz(-1.5651903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5311997) q[3];
sx q[3];
rz(-2.5259113) q[3];
sx q[3];
rz(0.99056584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4697326) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.0107182) q[0];
sx q[0];
rz(-0.81881443) q[0];
sx q[0];
rz(2.4406216) q[0];
rz(-2.2531807) q[1];
sx q[1];
rz(-0.41524926) q[1];
sx q[1];
rz(-1.433724) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1183567) q[0];
sx q[0];
rz(-0.5254841) q[0];
sx q[0];
rz(0.21696399) q[0];
x q[1];
rz(1.5979684) q[2];
sx q[2];
rz(-1.1650231) q[2];
sx q[2];
rz(2.9494065) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6810602) q[1];
sx q[1];
rz(-1.5523445) q[1];
sx q[1];
rz(0.021120763) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4096595) q[3];
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
rz(1.4768451) q[2];
rz(0.19674033) q[3];
sx q[3];
rz(-2.0048095) q[3];
sx q[3];
rz(-2.2018532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5340586) q[0];
sx q[0];
rz(-1.5787831) q[0];
sx q[0];
rz(2.0205355) q[0];
rz(-2.6634482) q[1];
sx q[1];
rz(-0.68111626) q[1];
sx q[1];
rz(-0.71472439) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31340718) q[0];
sx q[0];
rz(-1.3696693) q[0];
sx q[0];
rz(0.49230663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7585538) q[2];
sx q[2];
rz(-1.0543807) q[2];
sx q[2];
rz(0.4412152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8712232) q[1];
sx q[1];
rz(-2.0107234) q[1];
sx q[1];
rz(2.1384856) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86831324) q[3];
sx q[3];
rz(-2.3037617) q[3];
sx q[3];
rz(0.7759076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5152682) q[2];
sx q[2];
rz(-1.8147261) q[2];
sx q[2];
rz(1.932762) q[2];
rz(-1.6995466) q[3];
sx q[3];
rz(-2.2550826) q[3];
sx q[3];
rz(3.076021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.47990738) q[0];
sx q[0];
rz(-1.3085145) q[0];
sx q[0];
rz(3.0611962) q[0];
rz(-0.54795625) q[1];
sx q[1];
rz(-0.68956551) q[1];
sx q[1];
rz(-1.7908304) q[1];
rz(-2.0399848) q[2];
sx q[2];
rz(-1.3528578) q[2];
sx q[2];
rz(2.8982671) q[2];
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
