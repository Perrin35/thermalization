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
rz(-1.6428525) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(-0.3436389) q[1];
sx q[1];
rz(2.763881) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41357562) q[0];
sx q[0];
rz(-0.33057645) q[0];
sx q[0];
rz(-2.3982993) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1324563) q[2];
sx q[2];
rz(-1.1852263) q[2];
sx q[2];
rz(1.350268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39241782) q[1];
sx q[1];
rz(-2.2274744) q[1];
sx q[1];
rz(-1.6139469) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1123766) q[3];
sx q[3];
rz(-1.8354641) q[3];
sx q[3];
rz(0.1614557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7517884) q[2];
sx q[2];
rz(-1.7916388) q[2];
sx q[2];
rz(0.26503116) q[2];
rz(2.7122279) q[3];
sx q[3];
rz(-1.2484442) q[3];
sx q[3];
rz(-0.00092367729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.86004177) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81102449) q[0];
sx q[0];
rz(-2.1386792) q[0];
sx q[0];
rz(0.83356838) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3798736) q[2];
sx q[2];
rz(-0.8702308) q[2];
sx q[2];
rz(-3.0143486) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.94446429) q[1];
sx q[1];
rz(-2.1310867) q[1];
sx q[1];
rz(0.48546882) q[1];
rz(-2.9577291) q[3];
sx q[3];
rz(-1.2753092) q[3];
sx q[3];
rz(1.1671403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5106875) q[2];
sx q[2];
rz(-2.5248933) q[2];
sx q[2];
rz(-0.56646937) q[2];
rz(0.72930068) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(0.20310371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.94839621) q[0];
sx q[0];
rz(-2.0553135) q[0];
sx q[0];
rz(2.7753944) q[0];
rz(-1.5858448) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(-0.050447024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9589466) q[0];
sx q[0];
rz(-0.0894657) q[0];
sx q[0];
rz(-2.0561809) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0773727) q[2];
sx q[2];
rz(-2.1972547) q[2];
sx q[2];
rz(-0.74365091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.09035) q[1];
sx q[1];
rz(-1.7900677) q[1];
sx q[1];
rz(1.7419113) q[1];
rz(1.9956224) q[3];
sx q[3];
rz(-1.7033938) q[3];
sx q[3];
rz(1.4736922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.067165) q[2];
sx q[2];
rz(-1.0600435) q[2];
sx q[2];
rz(1.3107497) q[2];
rz(-1.358076) q[3];
sx q[3];
rz(-0.57099968) q[3];
sx q[3];
rz(1.8324435) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27291372) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(-1.7568461) q[0];
rz(-1.9028496) q[1];
sx q[1];
rz(-1.9107995) q[1];
sx q[1];
rz(1.967427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096701972) q[0];
sx q[0];
rz(-0.47098038) q[0];
sx q[0];
rz(-1.1440008) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8205244) q[2];
sx q[2];
rz(-1.839723) q[2];
sx q[2];
rz(-3.0422701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2432855) q[1];
sx q[1];
rz(-1.8970058) q[1];
sx q[1];
rz(2.1065358) q[1];
rz(-pi) q[2];
rz(-2.1010059) q[3];
sx q[3];
rz(-2.3755666) q[3];
sx q[3];
rz(1.211973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43526402) q[2];
sx q[2];
rz(-1.8212943) q[2];
sx q[2];
rz(3.1026133) q[2];
rz(-0.6428166) q[3];
sx q[3];
rz(-2.1233852) q[3];
sx q[3];
rz(-2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103127) q[0];
sx q[0];
rz(-1.0142925) q[0];
sx q[0];
rz(-3.1211299) q[0];
rz(-0.19275716) q[1];
sx q[1];
rz(-2.5242476) q[1];
sx q[1];
rz(1.9761168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11949018) q[0];
sx q[0];
rz(-2.179232) q[0];
sx q[0];
rz(-1.0761989) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54372676) q[2];
sx q[2];
rz(-0.75060487) q[2];
sx q[2];
rz(0.0074530938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9200774) q[1];
sx q[1];
rz(-1.0672004) q[1];
sx q[1];
rz(-3.1184993) q[1];
rz(-pi) q[2];
rz(-1.7609414) q[3];
sx q[3];
rz(-1.4961637) q[3];
sx q[3];
rz(2.7281705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8274902) q[2];
sx q[2];
rz(-2.1496488) q[2];
sx q[2];
rz(0.52725434) q[2];
rz(0.54316795) q[3];
sx q[3];
rz(-0.68734622) q[3];
sx q[3];
rz(0.34272042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0813893) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(-2.7348837) q[0];
rz(-1.1609062) q[1];
sx q[1];
rz(-2.2229767) q[1];
sx q[1];
rz(1.0103753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.564931) q[0];
sx q[0];
rz(-2.1157678) q[0];
sx q[0];
rz(-0.14133639) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40596227) q[2];
sx q[2];
rz(-2.7926707) q[2];
sx q[2];
rz(2.5418848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0370673) q[1];
sx q[1];
rz(-2.7114093) q[1];
sx q[1];
rz(2.6135315) q[1];
rz(-1.3090735) q[3];
sx q[3];
rz(-0.78262586) q[3];
sx q[3];
rz(-1.4306376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0747718) q[2];
sx q[2];
rz(-1.2819042) q[2];
sx q[2];
rz(-2.1066966) q[2];
rz(1.3567989) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(-2.518173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1739625) q[0];
sx q[0];
rz(-1.6169463) q[0];
sx q[0];
rz(0.3279283) q[0];
rz(3.0294042) q[1];
sx q[1];
rz(-1.2022377) q[1];
sx q[1];
rz(-0.97253886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93720651) q[0];
sx q[0];
rz(-0.68336672) q[0];
sx q[0];
rz(2.1478189) q[0];
rz(2.8424524) q[2];
sx q[2];
rz(-1.9950331) q[2];
sx q[2];
rz(-1.5628634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2518721) q[1];
sx q[1];
rz(-2.1600284) q[1];
sx q[1];
rz(1.5934029) q[1];
x q[2];
rz(-2.6554606) q[3];
sx q[3];
rz(-0.68607578) q[3];
sx q[3];
rz(-3.0239575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5358676) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(2.4884339) q[2];
rz(0.43290916) q[3];
sx q[3];
rz(-2.631729) q[3];
sx q[3];
rz(-2.9292817) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2046278) q[0];
sx q[0];
rz(-2.300394) q[0];
sx q[0];
rz(-0.2247819) q[0];
rz(-0.79554355) q[1];
sx q[1];
rz(-1.3139775) q[1];
sx q[1];
rz(-2.0733817) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58086181) q[0];
sx q[0];
rz(-0.95968548) q[0];
sx q[0];
rz(2.1944502) q[0];
rz(0.76935591) q[2];
sx q[2];
rz(-2.0343668) q[2];
sx q[2];
rz(-0.052841436) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.873899) q[1];
sx q[1];
rz(-1.5452478) q[1];
sx q[1];
rz(-0.68768312) q[1];
rz(3.0826236) q[3];
sx q[3];
rz(-1.7921721) q[3];
sx q[3];
rz(0.27316545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3855359) q[2];
sx q[2];
rz(-1.3827366) q[2];
sx q[2];
rz(-2.5059911) q[2];
rz(-0.33291891) q[3];
sx q[3];
rz(-2.1300485) q[3];
sx q[3];
rz(-2.6788768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3391089) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(3.0185757) q[0];
rz(-1.3975551) q[1];
sx q[1];
rz(-1.3879644) q[1];
sx q[1];
rz(2.3618598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132039) q[0];
sx q[0];
rz(-1.8411978) q[0];
sx q[0];
rz(0.89288519) q[0];
rz(-0.17895584) q[2];
sx q[2];
rz(-2.4479986) q[2];
sx q[2];
rz(-3.0998942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9537264) q[1];
sx q[1];
rz(-0.426891) q[1];
sx q[1];
rz(-2.0679451) q[1];
x q[2];
rz(2.2717495) q[3];
sx q[3];
rz(-0.16582684) q[3];
sx q[3];
rz(-1.0505249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4745549) q[2];
sx q[2];
rz(-1.2660618) q[2];
sx q[2];
rz(3.0511268) q[2];
rz(2.7247834) q[3];
sx q[3];
rz(-2.3495245) q[3];
sx q[3];
rz(0.99378234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4823293) q[0];
sx q[0];
rz(-1.3220795) q[0];
sx q[0];
rz(3.0521159) q[0];
rz(-0.80884519) q[1];
sx q[1];
rz(-2.6409812) q[1];
sx q[1];
rz(2.083875) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0322005) q[0];
sx q[0];
rz(-1.9475749) q[0];
sx q[0];
rz(-0.44372875) q[0];
x q[1];
rz(1.0187947) q[2];
sx q[2];
rz(-1.5755781) q[2];
sx q[2];
rz(1.1301745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48185086) q[1];
sx q[1];
rz(-2.5999477) q[1];
sx q[1];
rz(0.84195824) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1154409) q[3];
sx q[3];
rz(-0.75687486) q[3];
sx q[3];
rz(0.43318403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40314254) q[2];
sx q[2];
rz(-2.6046627) q[2];
sx q[2];
rz(-2.7682847) q[2];
rz(0.25660723) q[3];
sx q[3];
rz(-1.720287) q[3];
sx q[3];
rz(-3.0671425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858717) q[0];
sx q[0];
rz(-1.5789565) q[0];
sx q[0];
rz(-1.4174905) q[0];
rz(1.7120842) q[1];
sx q[1];
rz(-2.7919339) q[1];
sx q[1];
rz(-2.3333593) q[1];
rz(-3.1222432) q[2];
sx q[2];
rz(-0.45889284) q[2];
sx q[2];
rz(1.5677551) q[2];
rz(1.5143916) q[3];
sx q[3];
rz(-2.5396697) q[3];
sx q[3];
rz(-1.664262) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
