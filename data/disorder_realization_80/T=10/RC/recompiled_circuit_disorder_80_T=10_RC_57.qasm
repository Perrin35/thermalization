OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39188448) q[0];
sx q[0];
rz(-0.19667974) q[0];
sx q[0];
rz(1.1893907) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(-2.7639311) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084378622) q[0];
sx q[0];
rz(-0.14176653) q[0];
sx q[0];
rz(-2.111582) q[0];
rz(1.2871735) q[2];
sx q[2];
rz(-1.5429153) q[2];
sx q[2];
rz(1.6698128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9610112) q[1];
sx q[1];
rz(-2.2716224) q[1];
sx q[1];
rz(2.5518774) q[1];
rz(1.430106) q[3];
sx q[3];
rz(-1.1555539) q[3];
sx q[3];
rz(-1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1315786) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(-2.4689891) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(-1.8030362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23068962) q[0];
sx q[0];
rz(-2.2407273) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(-0.28796089) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32591336) q[0];
sx q[0];
rz(-2.8796112) q[0];
sx q[0];
rz(-2.322305) q[0];
rz(-pi) q[1];
rz(-2.1177883) q[2];
sx q[2];
rz(-2.4734801) q[2];
sx q[2];
rz(0.39272768) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8781232) q[1];
sx q[1];
rz(-0.97349226) q[1];
sx q[1];
rz(1.0695446) q[1];
rz(-pi) q[2];
rz(-0.1045707) q[3];
sx q[3];
rz(-0.35211709) q[3];
sx q[3];
rz(0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5471389) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(-1.6681558) q[2];
rz(-2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8163452) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(2.1898988) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0591473) q[0];
sx q[0];
rz(-1.7623616) q[0];
sx q[0];
rz(3.0491769) q[0];
rz(-0.2256514) q[2];
sx q[2];
rz(-1.0661085) q[2];
sx q[2];
rz(0.80863189) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3141331) q[1];
sx q[1];
rz(-2.3312807) q[1];
sx q[1];
rz(2.00287) q[1];
rz(2.2610407) q[3];
sx q[3];
rz(-0.75062245) q[3];
sx q[3];
rz(2.7504138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.13016985) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(0.34417957) q[2];
rz(-2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0043871) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.244506) q[0];
rz(3.0124774) q[1];
sx q[1];
rz(-1.7872417) q[1];
sx q[1];
rz(2.7688162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0817954) q[0];
sx q[0];
rz(-0.76515388) q[0];
sx q[0];
rz(-2.9761936) q[0];
x q[1];
rz(-1.9984841) q[2];
sx q[2];
rz(-2.3017985) q[2];
sx q[2];
rz(1.0243624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2671632) q[1];
sx q[1];
rz(-1.0871372) q[1];
sx q[1];
rz(2.8664385) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1552912) q[3];
sx q[3];
rz(-2.5875475) q[3];
sx q[3];
rz(0.98741764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.63561511) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(0.44057009) q[3];
sx q[3];
rz(-0.4959271) q[3];
sx q[3];
rz(2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.2878081) q[0];
rz(1.6632535) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(0.53422654) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6144845) q[0];
sx q[0];
rz(-1.8728349) q[0];
sx q[0];
rz(-1.8943647) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1557501) q[2];
sx q[2];
rz(-2.5702229) q[2];
sx q[2];
rz(2.5845598) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91671645) q[1];
sx q[1];
rz(-1.6460878) q[1];
sx q[1];
rz(-3.0776943) q[1];
x q[2];
rz(2.4228135) q[3];
sx q[3];
rz(-1.6703509) q[3];
sx q[3];
rz(2.0841141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(2.9980998) q[2];
rz(-1.7701373) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(0.21970704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7025529) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(-0.23705661) q[0];
rz(1.9006231) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(1.3751078) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9000589) q[0];
sx q[0];
rz(-2.7068479) q[0];
sx q[0];
rz(2.4777806) q[0];
rz(-pi) q[1];
rz(-1.6429971) q[2];
sx q[2];
rz(-2.318396) q[2];
sx q[2];
rz(-0.88897926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4543537) q[1];
sx q[1];
rz(-2.7308186) q[1];
sx q[1];
rz(0.06128581) q[1];
rz(-2.4414805) q[3];
sx q[3];
rz(-2.5317319) q[3];
sx q[3];
rz(-1.3175347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.53283006) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(-2.9928845) q[2];
rz(3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(2.9490411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(0.87316978) q[0];
rz(0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5859563) q[0];
sx q[0];
rz(-1.3480098) q[0];
sx q[0];
rz(-2.1369834) q[0];
rz(2.7883045) q[2];
sx q[2];
rz(-1.7751667) q[2];
sx q[2];
rz(-2.3452961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9675688) q[1];
sx q[1];
rz(-2.1597383) q[1];
sx q[1];
rz(-0.80935652) q[1];
x q[2];
rz(-2.300802) q[3];
sx q[3];
rz(-1.519324) q[3];
sx q[3];
rz(-1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7579047) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(-0.54523462) q[2];
rz(-2.7111354) q[3];
sx q[3];
rz(-1.1377708) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51171821) q[0];
sx q[0];
rz(-0.015462333) q[0];
sx q[0];
rz(-1.9301201) q[0];
rz(-0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57293939) q[0];
sx q[0];
rz(-1.4927673) q[0];
sx q[0];
rz(0.22139876) q[0];
rz(-pi) q[1];
rz(2.321645) q[2];
sx q[2];
rz(-2.2689399) q[2];
sx q[2];
rz(1.1993711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9981873) q[1];
sx q[1];
rz(-0.37969509) q[1];
sx q[1];
rz(-3.0642964) q[1];
x q[2];
rz(-2.4754727) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(0.70762779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(0.28042173) q[2];
rz(2.5366606) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1542926) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(0.61532414) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(0.5756793) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8262779) q[0];
sx q[0];
rz(-1.1018254) q[0];
sx q[0];
rz(0.6420352) q[0];
rz(-pi) q[1];
rz(-1.0506389) q[2];
sx q[2];
rz(-2.9165604) q[2];
sx q[2];
rz(0.88423836) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3153783) q[1];
sx q[1];
rz(-1.4634906) q[1];
sx q[1];
rz(-0.025749287) q[1];
x q[2];
rz(1.0300893) q[3];
sx q[3];
rz(-1.732015) q[3];
sx q[3];
rz(-2.2459523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82970396) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(3.1196307) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(-2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-2.3529265) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4983738) q[0];
sx q[0];
rz(-1.6341097) q[0];
sx q[0];
rz(-1.4988585) q[0];
rz(-pi) q[1];
x q[1];
rz(1*pi/15) q[2];
sx q[2];
rz(-2.3839715) q[2];
sx q[2];
rz(2.9486738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.22419588) q[1];
sx q[1];
rz(-2.2317413) q[1];
sx q[1];
rz(1.7978653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7006847) q[3];
sx q[3];
rz(-1.3135859) q[3];
sx q[3];
rz(2.5901026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6897631) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(-2.7815212) q[2];
rz(-1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2071028) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(0.33640484) q[2];
sx q[2];
rz(-1.9336666) q[2];
sx q[2];
rz(-1.578707) q[2];
rz(1.885407) q[3];
sx q[3];
rz(-2.0046069) q[3];
sx q[3];
rz(-0.14179695) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
