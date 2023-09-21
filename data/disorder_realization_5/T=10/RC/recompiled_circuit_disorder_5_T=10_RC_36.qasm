OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(-0.45698693) q[0];
rz(2.5198088) q[1];
sx q[1];
rz(3.8222651) q[1];
sx q[1];
rz(8.1487976) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9306643) q[0];
sx q[0];
rz(-1.7893447) q[0];
sx q[0];
rz(-1.4320089) q[0];
rz(-pi) q[1];
rz(-0.692042) q[2];
sx q[2];
rz(-3.0243052) q[2];
sx q[2];
rz(2.8534944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4133271) q[1];
sx q[1];
rz(-1.3682377) q[1];
sx q[1];
rz(-3.0800746) q[1];
x q[2];
rz(2.7492417) q[3];
sx q[3];
rz(-0.96071834) q[3];
sx q[3];
rz(-3.0685134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0322545) q[2];
sx q[2];
rz(-0.87301746) q[2];
sx q[2];
rz(0.16513744) q[2];
rz(1.5103546) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(-0.74222773) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6363268) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-3.0965565) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(0.80274686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6532324) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(-2.1064208) q[0];
rz(-pi) q[1];
rz(-1.0674092) q[2];
sx q[2];
rz(-1.479584) q[2];
sx q[2];
rz(0.78546055) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3853559) q[1];
sx q[1];
rz(-1.5192351) q[1];
sx q[1];
rz(-0.95183422) q[1];
rz(-pi) q[2];
rz(1.4071776) q[3];
sx q[3];
rz(-2.46582) q[3];
sx q[3];
rz(0.27392745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37614432) q[2];
sx q[2];
rz(-2.5120698) q[2];
sx q[2];
rz(1.7155898) q[2];
rz(-2.5358893) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2633318) q[0];
sx q[0];
rz(-0.86255951) q[0];
sx q[0];
rz(-1.0017298) q[0];
x q[1];
rz(0.3091829) q[2];
sx q[2];
rz(-2.2367466) q[2];
sx q[2];
rz(0.14112976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7741751) q[1];
sx q[1];
rz(-1.4111102) q[1];
sx q[1];
rz(-1.3819441) q[1];
x q[2];
rz(1.9698824) q[3];
sx q[3];
rz(-2.2200924) q[3];
sx q[3];
rz(2.2513575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4703935) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(-1.408067) q[2];
rz(-2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6614439) q[0];
sx q[0];
rz(-0.34879768) q[0];
sx q[0];
rz(-0.21425042) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-2.0206101) q[1];
sx q[1];
rz(-2.4750211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19569451) q[0];
sx q[0];
rz(-1.7576522) q[0];
sx q[0];
rz(2.9325783) q[0];
x q[1];
rz(0.33118403) q[2];
sx q[2];
rz(-2.3996155) q[2];
sx q[2];
rz(-2.3664898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9208593) q[1];
sx q[1];
rz(-1.024106) q[1];
sx q[1];
rz(2.2842555) q[1];
rz(-pi) q[2];
rz(1.7188844) q[3];
sx q[3];
rz(-1.423096) q[3];
sx q[3];
rz(-1.577391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(-1.9541698) q[3];
sx q[3];
rz(-1.238845) q[3];
sx q[3];
rz(-2.58113) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-2.3045325) q[0];
rz(0.95056668) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.1594835) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449074) q[0];
sx q[0];
rz(-0.7888182) q[0];
sx q[0];
rz(-2.6916654) q[0];
rz(2.1690364) q[2];
sx q[2];
rz(-1.6998569) q[2];
sx q[2];
rz(1.8051403) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9002478) q[1];
sx q[1];
rz(-0.7175788) q[1];
sx q[1];
rz(-0.89300968) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5043143) q[3];
sx q[3];
rz(-1.1253469) q[3];
sx q[3];
rz(0.74419903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(1.8699899) q[2];
rz(-1.050625) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(0.064855382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5746675) q[0];
sx q[0];
rz(-2.8820679) q[0];
sx q[0];
rz(-2.0264453) q[0];
rz(3.0976345) q[1];
sx q[1];
rz(-1.3479439) q[1];
sx q[1];
rz(0.10087068) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2642919) q[0];
sx q[0];
rz(-1.4192686) q[0];
sx q[0];
rz(-2.507188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2783373) q[2];
sx q[2];
rz(-2.0403701) q[2];
sx q[2];
rz(1.3172319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2821591) q[1];
sx q[1];
rz(-1.7267425) q[1];
sx q[1];
rz(-1.7755416) q[1];
rz(2.7363339) q[3];
sx q[3];
rz(-1.3558946) q[3];
sx q[3];
rz(0.98291558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4042525) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(0.92528382) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(1.293175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3826542) q[0];
sx q[0];
rz(-2.721334) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(1.5902279) q[1];
sx q[1];
rz(-2.1236877) q[1];
sx q[1];
rz(2.4538453) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15935414) q[0];
sx q[0];
rz(-1.237545) q[0];
sx q[0];
rz(-0.23157816) q[0];
rz(-pi) q[1];
rz(2.1805448) q[2];
sx q[2];
rz(-2.7597722) q[2];
sx q[2];
rz(-2.080999) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3878165) q[1];
sx q[1];
rz(-0.13145914) q[1];
sx q[1];
rz(-2.9629575) q[1];
x q[2];
rz(-1.5457821) q[3];
sx q[3];
rz(-2.5731312) q[3];
sx q[3];
rz(0.75077885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(2.2040099) q[2];
rz(-0.98085105) q[3];
sx q[3];
rz(-2.8068481) q[3];
sx q[3];
rz(-2.3576095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6136318) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(0.079285346) q[0];
rz(-2.0203363) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.942873) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51533651) q[0];
sx q[0];
rz(-1.6394098) q[0];
sx q[0];
rz(-0.4347883) q[0];
rz(-pi) q[1];
rz(1.7940117) q[2];
sx q[2];
rz(-2.7701021) q[2];
sx q[2];
rz(1.5944634) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.291154) q[1];
sx q[1];
rz(-1.6206695) q[1];
sx q[1];
rz(0.50317851) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6783887) q[3];
sx q[3];
rz(-1.8052881) q[3];
sx q[3];
rz(-2.6686058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2856059) q[2];
sx q[2];
rz(-0.14582835) q[2];
sx q[2];
rz(0.78906995) q[2];
rz(2.7834535) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(-1.8607128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(2.5532706) q[0];
rz(-0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(-0.9334329) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.371884) q[0];
sx q[0];
rz(-3.0196307) q[0];
sx q[0];
rz(0.50942771) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2869296) q[2];
sx q[2];
rz(-2.6551506) q[2];
sx q[2];
rz(1.9180627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0837299) q[1];
sx q[1];
rz(-1.5340367) q[1];
sx q[1];
rz(-2.4576393) q[1];
rz(-pi) q[2];
rz(-0.65255717) q[3];
sx q[3];
rz(-1.2823294) q[3];
sx q[3];
rz(2.1289338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-0.82010078) q[2];
sx q[2];
rz(2.771647) q[2];
rz(0.89921078) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(0.95329154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1870435) q[0];
sx q[0];
rz(-1.8148913) q[0];
sx q[0];
rz(0.6533587) q[0];
rz(-1.5006784) q[1];
sx q[1];
rz(-1.4705642) q[1];
sx q[1];
rz(0.18383372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18098022) q[0];
sx q[0];
rz(-1.4502589) q[0];
sx q[0];
rz(-1.6880077) q[0];
rz(-0.75283639) q[2];
sx q[2];
rz(-2.2194153) q[2];
sx q[2];
rz(-0.93968771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7067691) q[1];
sx q[1];
rz(-0.69908792) q[1];
sx q[1];
rz(-1.8301151) q[1];
rz(-0.023601836) q[3];
sx q[3];
rz(-0.58348237) q[3];
sx q[3];
rz(-2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9188149) q[2];
sx q[2];
rz(-2.1464525) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(-2.0683794) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.3674659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.68207537) q[0];
sx q[0];
rz(-2.2820602) q[0];
sx q[0];
rz(2.5862502) q[0];
rz(-0.36322414) q[1];
sx q[1];
rz(-1.8050615) q[1];
sx q[1];
rz(-0.25711679) q[1];
rz(1.9831052) q[2];
sx q[2];
rz(-1.9956797) q[2];
sx q[2];
rz(3.0394625) q[2];
rz(1.3115053) q[3];
sx q[3];
rz(-2.5355831) q[3];
sx q[3];
rz(-2.2021962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
