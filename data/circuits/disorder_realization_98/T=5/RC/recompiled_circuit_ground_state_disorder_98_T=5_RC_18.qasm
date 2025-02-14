OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7336361) q[0];
sx q[0];
rz(-0.1853369) q[0];
sx q[0];
rz(1.4987401) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(2.7979538) q[1];
sx q[1];
rz(9.8024896) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728017) q[0];
sx q[0];
rz(-2.8110162) q[0];
sx q[0];
rz(-2.3982993) q[0];
x q[1];
rz(-0.42141098) q[2];
sx q[2];
rz(-1.1666036) q[2];
sx q[2];
rz(-0.39502783) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6785665) q[1];
sx q[1];
rz(-2.4837079) q[1];
sx q[1];
rz(-0.055906217) q[1];
rz(1.6781647) q[3];
sx q[3];
rz(-2.8753548) q[3];
sx q[3];
rz(-3.0913946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7517884) q[2];
sx q[2];
rz(-1.3499539) q[2];
sx q[2];
rz(0.26503116) q[2];
rz(-0.42936471) q[3];
sx q[3];
rz(-1.8931484) q[3];
sx q[3];
rz(-3.140669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.65381831) q[0];
sx q[0];
rz(-0.14350292) q[0];
sx q[0];
rz(1.8408467) q[0];
rz(0.067151345) q[1];
sx q[1];
rz(-2.2323699) q[1];
sx q[1];
rz(2.2815509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9165827) q[0];
sx q[0];
rz(-2.2448178) q[0];
sx q[0];
rz(-2.3300578) q[0];
rz(0.70921398) q[2];
sx q[2];
rz(-1.0150036) q[2];
sx q[2];
rz(-1.9947987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7886161) q[1];
sx q[1];
rz(-1.1644272) q[1];
sx q[1];
rz(0.95384903) q[1];
rz(-pi) q[2];
rz(2.1116756) q[3];
sx q[3];
rz(-2.7949998) q[3];
sx q[3];
rz(-2.5427713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5106875) q[2];
sx q[2];
rz(-2.5248933) q[2];
sx q[2];
rz(2.5751233) q[2];
rz(-2.412292) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(-2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94839621) q[0];
sx q[0];
rz(-2.0553135) q[0];
sx q[0];
rz(-0.36619827) q[0];
rz(-1.5557479) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(0.050447024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18264601) q[0];
sx q[0];
rz(-3.052127) q[0];
sx q[0];
rz(1.0854118) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5620805) q[2];
sx q[2];
rz(-0.77634927) q[2];
sx q[2];
rz(0.001359847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.76169) q[1];
sx q[1];
rz(-0.27728785) q[1];
sx q[1];
rz(0.65245858) q[1];
x q[2];
rz(2.9962323) q[3];
sx q[3];
rz(-1.1499377) q[3];
sx q[3];
rz(0.03736729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(-1.3107497) q[2];
rz(-1.7835167) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.8686789) q[0];
sx q[0];
rz(-0.22628117) q[0];
sx q[0];
rz(-1.3847466) q[0];
rz(-1.2387431) q[1];
sx q[1];
rz(-1.2307931) q[1];
sx q[1];
rz(1.967427) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730302) q[0];
sx q[0];
rz(-1.1449887) q[0];
sx q[0];
rz(-2.9338475) q[0];
rz(-pi) q[1];
rz(-0.73100369) q[2];
sx q[2];
rz(-0.36492294) q[2];
sx q[2];
rz(2.475955) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2432855) q[1];
sx q[1];
rz(-1.2445868) q[1];
sx q[1];
rz(2.1065358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0405868) q[3];
sx q[3];
rz(-0.76602606) q[3];
sx q[3];
rz(1.9296196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7063286) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(0.038979385) q[2];
rz(-2.4987761) q[3];
sx q[3];
rz(-1.0182074) q[3];
sx q[3];
rz(-2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83127999) q[0];
sx q[0];
rz(-2.1273002) q[0];
sx q[0];
rz(-0.020462791) q[0];
rz(2.9488355) q[1];
sx q[1];
rz(-2.5242476) q[1];
sx q[1];
rz(-1.1654759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7503238) q[0];
sx q[0];
rz(-1.9708777) q[0];
sx q[0];
rz(-2.4720936) q[0];
rz(-pi) q[1];
rz(1.1212249) q[2];
sx q[2];
rz(-0.94749852) q[2];
sx q[2];
rz(2.4433608) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9679035) q[1];
sx q[1];
rz(-2.6375131) q[1];
sx q[1];
rz(-1.5289115) q[1];
rz(-pi) q[2];
rz(-3.0655954) q[3];
sx q[3];
rz(-1.381187) q[3];
sx q[3];
rz(-1.9698683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8274902) q[2];
sx q[2];
rz(-2.1496488) q[2];
sx q[2];
rz(2.6143383) q[2];
rz(-2.5984247) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(-2.7988722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602033) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(0.40670893) q[0];
rz(-1.1609062) q[1];
sx q[1];
rz(-0.91861594) q[1];
sx q[1];
rz(-1.0103753) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079491678) q[0];
sx q[0];
rz(-1.6915503) q[0];
sx q[0];
rz(-2.1202205) q[0];
rz(-pi) q[1];
rz(-0.40596227) q[2];
sx q[2];
rz(-0.34892198) q[2];
sx q[2];
rz(-2.5418848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1878309) q[1];
sx q[1];
rz(-1.7825025) q[1];
sx q[1];
rz(0.37734887) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80548894) q[3];
sx q[3];
rz(-1.7542766) q[3];
sx q[3];
rz(0.047540548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0668209) q[2];
sx q[2];
rz(-1.2819042) q[2];
sx q[2];
rz(-1.034896) q[2];
rz(-1.3567989) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9676301) q[0];
sx q[0];
rz(-1.6169463) q[0];
sx q[0];
rz(0.3279283) q[0];
rz(3.0294042) q[1];
sx q[1];
rz(-1.2022377) q[1];
sx q[1];
rz(-0.97253886) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5060726) q[0];
sx q[0];
rz(-1.0131665) q[0];
sx q[0];
rz(-0.41801674) q[0];
rz(-pi) q[1];
rz(2.0123291) q[2];
sx q[2];
rz(-1.8427094) q[2];
sx q[2];
rz(-0.11833469) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33148872) q[1];
sx q[1];
rz(-1.5520025) q[1];
sx q[1];
rz(-2.5522425) q[1];
rz(-2.5150033) q[3];
sx q[3];
rz(-1.8712776) q[3];
sx q[3];
rz(1.8412875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5358676) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(2.4884339) q[2];
rz(2.7086835) q[3];
sx q[3];
rz(-0.50986367) q[3];
sx q[3];
rz(0.21231095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9369649) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(2.9168108) q[0];
rz(2.3460491) q[1];
sx q[1];
rz(-1.8276151) q[1];
sx q[1];
rz(2.0733817) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.381425) q[0];
sx q[0];
rz(-1.0720709) q[0];
sx q[0];
rz(-2.4295761) q[0];
rz(-2.3722367) q[2];
sx q[2];
rz(-2.0343668) q[2];
sx q[2];
rz(-0.052841436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2676937) q[1];
sx q[1];
rz(-1.5963449) q[1];
sx q[1];
rz(2.4539095) q[1];
x q[2];
rz(1.8269038) q[3];
sx q[3];
rz(-2.912622) q[3];
sx q[3];
rz(-3.1310905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3855359) q[2];
sx q[2];
rz(-1.3827366) q[2];
sx q[2];
rz(-0.63560152) q[2];
rz(-0.33291891) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(2.6788768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8024837) q[0];
sx q[0];
rz(-3.1153296) q[0];
sx q[0];
rz(-3.0185757) q[0];
rz(-1.3975551) q[1];
sx q[1];
rz(-1.7536283) q[1];
sx q[1];
rz(-2.3618598) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1542669) q[0];
sx q[0];
rz(-0.92181081) q[0];
sx q[0];
rz(-0.34190468) q[0];
rz(-pi) q[1];
rz(0.68569195) q[2];
sx q[2];
rz(-1.6848411) q[2];
sx q[2];
rz(1.6673078) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9537264) q[1];
sx q[1];
rz(-2.7147016) q[1];
sx q[1];
rz(-1.0736476) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.443583) q[3];
sx q[3];
rz(-1.6774584) q[3];
sx q[3];
rz(-2.9677109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66703779) q[2];
sx q[2];
rz(-1.8755308) q[2];
sx q[2];
rz(-3.0511268) q[2];
rz(-2.7247834) q[3];
sx q[3];
rz(-0.79206812) q[3];
sx q[3];
rz(-2.1478103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6592634) q[0];
sx q[0];
rz(-1.8195131) q[0];
sx q[0];
rz(-3.0521159) q[0];
rz(-0.80884519) q[1];
sx q[1];
rz(-0.50061148) q[1];
sx q[1];
rz(-2.083875) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19729511) q[0];
sx q[0];
rz(-2.5677201) q[0];
sx q[0];
rz(-0.7446592) q[0];
x q[1];
rz(-1.5799149) q[2];
sx q[2];
rz(-0.55202019) q[2];
sx q[2];
rz(-2.7087351) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48185086) q[1];
sx q[1];
rz(-2.5999477) q[1];
sx q[1];
rz(-2.2996344) q[1];
x q[2];
rz(3.1154409) q[3];
sx q[3];
rz(-0.75687486) q[3];
sx q[3];
rz(2.7084086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40314254) q[2];
sx q[2];
rz(-0.53692997) q[2];
sx q[2];
rz(-0.37330791) q[2];
rz(2.8849854) q[3];
sx q[3];
rz(-1.4213057) q[3];
sx q[3];
rz(0.074450113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2557209) q[0];
sx q[0];
rz(-1.5626361) q[0];
sx q[0];
rz(1.7241021) q[0];
rz(1.7120842) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(1.5803554) q[2];
sx q[2];
rz(-2.0295967) q[2];
sx q[2];
rz(1.5461736) q[2];
rz(-0.96961602) q[3];
sx q[3];
rz(-1.5388699) q[3];
sx q[3];
rz(-0.13997302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
