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
rz(1.5732425) q[1];
sx q[1];
rz(-1.9116126) q[1];
sx q[1];
rz(0.65043515) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5395376) q[0];
sx q[0];
rz(-1.3663174) q[0];
sx q[0];
rz(0.50907593) q[0];
rz(-pi) q[1];
rz(-1.600163) q[2];
sx q[2];
rz(-1.8810388) q[2];
sx q[2];
rz(0.28679906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9792773) q[1];
sx q[1];
rz(-0.82274306) q[1];
sx q[1];
rz(2.2869799) q[1];
rz(-0.69919805) q[3];
sx q[3];
rz(-0.16575925) q[3];
sx q[3];
rz(2.0729271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9170561) q[2];
sx q[2];
rz(-1.3755596) q[2];
sx q[2];
rz(1.7898233) q[2];
rz(-2.3228862) q[3];
sx q[3];
rz(-1.6823781) q[3];
sx q[3];
rz(-1.5514577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3025892) q[0];
sx q[0];
rz(-2.4636457) q[0];
sx q[0];
rz(2.5658521) q[0];
rz(-2.2585244) q[1];
sx q[1];
rz(-1.6022316) q[1];
sx q[1];
rz(-1.8227122) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2275946) q[0];
sx q[0];
rz(-0.84216252) q[0];
sx q[0];
rz(1.2331859) q[0];
rz(-2.3568866) q[2];
sx q[2];
rz(-1.5883351) q[2];
sx q[2];
rz(0.49388805) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1309588) q[1];
sx q[1];
rz(-1.2181699) q[1];
sx q[1];
rz(2.3214798) q[1];
rz(-pi) q[2];
rz(-1.3926199) q[3];
sx q[3];
rz(-1.8867853) q[3];
sx q[3];
rz(-2.9079278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64112249) q[2];
sx q[2];
rz(-2.1239855) q[2];
sx q[2];
rz(0.23915954) q[2];
rz(0.33254361) q[3];
sx q[3];
rz(-1.4556689) q[3];
sx q[3];
rz(1.0913764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34514937) q[0];
sx q[0];
rz(-0.74757663) q[0];
sx q[0];
rz(-2.8969966) q[0];
rz(2.7469475) q[1];
sx q[1];
rz(-2.5375745) q[1];
sx q[1];
rz(-1.9371202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24893471) q[0];
sx q[0];
rz(-1.5818074) q[0];
sx q[0];
rz(1.5599773) q[0];
rz(2.7543805) q[2];
sx q[2];
rz(-1.1128328) q[2];
sx q[2];
rz(0.94113038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50379163) q[1];
sx q[1];
rz(-1.0560913) q[1];
sx q[1];
rz(-2.8082602) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0709023) q[3];
sx q[3];
rz(-1.9311889) q[3];
sx q[3];
rz(-2.5640324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7977153) q[2];
sx q[2];
rz(-0.98261967) q[2];
sx q[2];
rz(2.9003411) q[2];
rz(1.3681715) q[3];
sx q[3];
rz(-1.3894812) q[3];
sx q[3];
rz(2.1696137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57206804) q[0];
sx q[0];
rz(-1.2867186) q[0];
sx q[0];
rz(2.0401814) q[0];
rz(-1.6481579) q[1];
sx q[1];
rz(-0.28683528) q[1];
sx q[1];
rz(-2.1410087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1959609) q[0];
sx q[0];
rz(-1.1544265) q[0];
sx q[0];
rz(-0.076535688) q[0];
rz(-3.1317883) q[2];
sx q[2];
rz(-2.1270618) q[2];
sx q[2];
rz(-2.3988568) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0110338) q[1];
sx q[1];
rz(-0.10422626) q[1];
sx q[1];
rz(0.24918814) q[1];
x q[2];
rz(1.1492997) q[3];
sx q[3];
rz(-1.7022082) q[3];
sx q[3];
rz(0.67883713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7857886) q[2];
sx q[2];
rz(-1.641909) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8292238) q[0];
sx q[0];
rz(-2.0033422) q[0];
sx q[0];
rz(2.9421222) q[0];
rz(1.5169187) q[1];
sx q[1];
rz(-1.3360887) q[1];
sx q[1];
rz(2.6768501) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7620649) q[0];
sx q[0];
rz(-1.4389973) q[0];
sx q[0];
rz(1.2799311) q[0];
rz(0.22146341) q[2];
sx q[2];
rz(-1.628698) q[2];
sx q[2];
rz(2.1277949) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5158968) q[1];
sx q[1];
rz(-2.7346238) q[1];
sx q[1];
rz(-0.76677983) q[1];
x q[2];
rz(-3.1384077) q[3];
sx q[3];
rz(-0.69744686) q[3];
sx q[3];
rz(-1.3435329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0988079) q[2];
sx q[2];
rz(-1.837919) q[2];
sx q[2];
rz(-2.888077) q[2];
rz(-0.86067307) q[3];
sx q[3];
rz(-2.6200675) q[3];
sx q[3];
rz(-2.0025939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072170243) q[0];
sx q[0];
rz(-1.5266029) q[0];
sx q[0];
rz(1.3602863) q[0];
rz(0.6036497) q[1];
sx q[1];
rz(-0.79704469) q[1];
sx q[1];
rz(-0.48702494) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62327281) q[0];
sx q[0];
rz(-1.643299) q[0];
sx q[0];
rz(-2.4865315) q[0];
rz(-0.10425514) q[2];
sx q[2];
rz(-2.4001636) q[2];
sx q[2];
rz(-2.3280581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4364061) q[1];
sx q[1];
rz(-1.7366221) q[1];
sx q[1];
rz(-1.2064183) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0754082) q[3];
sx q[3];
rz(-0.45201354) q[3];
sx q[3];
rz(1.3936123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46502486) q[2];
sx q[2];
rz(-1.0176696) q[2];
sx q[2];
rz(2.6505995) q[2];
rz(2.6077014) q[3];
sx q[3];
rz(-2.5940354) q[3];
sx q[3];
rz(3.0965613) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8462867) q[0];
sx q[0];
rz(-2.5745109) q[0];
sx q[0];
rz(1.2714161) q[0];
rz(2.2170587) q[1];
sx q[1];
rz(-1.144616) q[1];
sx q[1];
rz(-2.1163993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36350016) q[0];
sx q[0];
rz(-2.1937167) q[0];
sx q[0];
rz(1.3062551) q[0];
rz(2.19609) q[2];
sx q[2];
rz(-1.7472113) q[2];
sx q[2];
rz(-2.8636065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5793663) q[1];
sx q[1];
rz(-1.0584029) q[1];
sx q[1];
rz(2.775869) q[1];
rz(-1.8162433) q[3];
sx q[3];
rz(-2.2625828) q[3];
sx q[3];
rz(-0.46190767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9292494) q[2];
sx q[2];
rz(-0.29444567) q[2];
sx q[2];
rz(2.219131) q[2];
rz(-2.755002) q[3];
sx q[3];
rz(-1.8542733) q[3];
sx q[3];
rz(-1.6283901) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(2.6225846) q[1];
sx q[1];
rz(-0.52878562) q[1];
sx q[1];
rz(-0.66600287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7836664) q[0];
sx q[0];
rz(-0.54364288) q[0];
sx q[0];
rz(1.3915791) q[0];
rz(-pi) q[1];
rz(-1.2326127) q[2];
sx q[2];
rz(-2.8034503) q[2];
sx q[2];
rz(-1.7902439) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0652661) q[1];
sx q[1];
rz(-0.90172036) q[1];
sx q[1];
rz(1.5651903) q[1];
x q[2];
rz(-1.5311997) q[3];
sx q[3];
rz(-2.5259113) q[3];
sx q[3];
rz(2.1510268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4697326) q[2];
sx q[2];
rz(-1.9817151) q[2];
sx q[2];
rz(1.0603909) q[2];
rz(1.4276069) q[3];
sx q[3];
rz(-1.5771882) q[3];
sx q[3];
rz(-2.3744627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(-0.41524926) q[1];
sx q[1];
rz(-1.433724) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64087039) q[0];
sx q[0];
rz(-1.4626012) q[0];
sx q[0];
rz(-0.51525028) q[0];
rz(0.063155576) q[2];
sx q[2];
rz(-0.40663162) q[2];
sx q[2];
rz(-2.8806608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4605324) q[1];
sx q[1];
rz(-1.5523445) q[1];
sx q[1];
rz(0.021120763) q[1];
rz(0.7319331) q[3];
sx q[3];
rz(-2.3022392) q[3];
sx q[3];
rz(-1.5956772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1514757) q[2];
sx q[2];
rz(-1.9141804) q[2];
sx q[2];
rz(1.4768451) q[2];
rz(-0.19674033) q[3];
sx q[3];
rz(-1.1367831) q[3];
sx q[3];
rz(0.93973947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5340586) q[0];
sx q[0];
rz(-1.5787831) q[0];
sx q[0];
rz(1.1210572) q[0];
rz(0.47814449) q[1];
sx q[1];
rz(-0.68111626) q[1];
sx q[1];
rz(-0.71472439) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31340718) q[0];
sx q[0];
rz(-1.7719233) q[0];
sx q[0];
rz(-2.649286) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38303886) q[2];
sx q[2];
rz(-1.0543807) q[2];
sx q[2];
rz(0.4412152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2523697) q[1];
sx q[1];
rz(-2.438527) q[1];
sx q[1];
rz(-0.85173793) q[1];
x q[2];
rz(-0.86831324) q[3];
sx q[3];
rz(-0.83783093) q[3];
sx q[3];
rz(0.7759076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.62632442) q[2];
sx q[2];
rz(-1.8147261) q[2];
sx q[2];
rz(1.932762) q[2];
rz(1.442046) q[3];
sx q[3];
rz(-0.88651005) q[3];
sx q[3];
rz(0.065571688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6616853) q[0];
sx q[0];
rz(-1.8330782) q[0];
sx q[0];
rz(-0.080396419) q[0];
rz(-0.54795625) q[1];
sx q[1];
rz(-0.68956551) q[1];
sx q[1];
rz(-1.7908304) q[1];
rz(2.0262268) q[2];
sx q[2];
rz(-2.6276988) q[2];
sx q[2];
rz(-2.217271) q[2];
rz(2.6790621) q[3];
sx q[3];
rz(-1.9560541) q[3];
sx q[3];
rz(3.1386459) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
