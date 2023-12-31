OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(3.9217477) q[1];
sx q[1];
rz(5.2182066) q[1];
sx q[1];
rz(10.301104) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0759461) q[0];
sx q[0];
rz(-1.3407205) q[0];
sx q[0];
rz(0.40212698) q[0];
rz(-pi) q[1];
rz(0.19199065) q[2];
sx q[2];
rz(-0.99388323) q[2];
sx q[2];
rz(-0.38919762) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6688924) q[1];
sx q[1];
rz(-1.9007705) q[1];
sx q[1];
rz(-0.015923576) q[1];
rz(-0.40684367) q[3];
sx q[3];
rz(-1.4300656) q[3];
sx q[3];
rz(-1.7498457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9709388) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(-0.7286287) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(2.8858378) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(-2.2671525) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6073608) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(2.1952654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3358243) q[2];
sx q[2];
rz(-0.58832303) q[2];
sx q[2];
rz(-0.36662835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.892131) q[1];
sx q[1];
rz(-0.79626894) q[1];
sx q[1];
rz(-0.65278058) q[1];
rz(-pi) q[2];
rz(-3.087895) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(-2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(0.43593105) q[2];
rz(-2.46051) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(-0.40288231) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(-0.39594617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229189) q[0];
sx q[0];
rz(-2.4791105) q[0];
sx q[0];
rz(0.89612095) q[0];
rz(-pi) q[1];
rz(3.1275438) q[2];
sx q[2];
rz(-1.2767681) q[2];
sx q[2];
rz(-0.77169466) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0242651) q[1];
sx q[1];
rz(-1.5803442) q[1];
sx q[1];
rz(0.68318232) q[1];
rz(-pi) q[2];
rz(-3.0523236) q[3];
sx q[3];
rz(-2.8850728) q[3];
sx q[3];
rz(1.741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(0.23920693) q[2];
rz(-3.0662597) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(-0.81992942) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(0.23342361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48196402) q[0];
sx q[0];
rz(-3.0273962) q[0];
sx q[0];
rz(1.0114848) q[0];
x q[1];
rz(0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-0.13194612) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64775002) q[1];
sx q[1];
rz(-2.2550681) q[1];
sx q[1];
rz(-1.1286331) q[1];
x q[2];
rz(-2.3573973) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-2.3941669) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6915879) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(0.94379395) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(2.3805526) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2634537) q[0];
sx q[0];
rz(-1.8315151) q[0];
sx q[0];
rz(-3.1104452) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20247395) q[2];
sx q[2];
rz(-1.379181) q[2];
sx q[2];
rz(1.7970049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6619751) q[1];
sx q[1];
rz(-1.2503337) q[1];
sx q[1];
rz(-1.596405) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8499591) q[3];
sx q[3];
rz(-1.8431292) q[3];
sx q[3];
rz(-0.84701049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3349907) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(-2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(0.87310711) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(2.81566) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8317141) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(-0.563234) q[0];
x q[1];
rz(-1.2665777) q[2];
sx q[2];
rz(-1.1433257) q[2];
sx q[2];
rz(1.3104591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.95820108) q[1];
sx q[1];
rz(-2.8124053) q[1];
sx q[1];
rz(-0.42894657) q[1];
x q[2];
rz(1.5876706) q[3];
sx q[3];
rz(-1.6718739) q[3];
sx q[3];
rz(-2.8942787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(2.4874172) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.1681898) q[3];
sx q[3];
rz(-2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.6690212) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(-3.022335) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8598547) q[0];
sx q[0];
rz(-0.87346948) q[0];
sx q[0];
rz(-1.2233234) q[0];
rz(-0.10995933) q[2];
sx q[2];
rz(-1.409515) q[2];
sx q[2];
rz(1.9732628) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1489361) q[1];
sx q[1];
rz(-1.438237) q[1];
sx q[1];
rz(-2.0723144) q[1];
rz(-pi) q[2];
rz(-2.4339606) q[3];
sx q[3];
rz(-1.8579357) q[3];
sx q[3];
rz(1.5054782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(1.3593486) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(-2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(2.0781562) q[0];
rz(2.8670782) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622334) q[0];
sx q[0];
rz(-1.0842807) q[0];
sx q[0];
rz(1.3079206) q[0];
rz(-pi) q[1];
rz(-2.1503259) q[2];
sx q[2];
rz(-1.5838924) q[2];
sx q[2];
rz(0.0038557204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0091128) q[1];
sx q[1];
rz(-1.9630034) q[1];
sx q[1];
rz(1.023804) q[1];
rz(-pi) q[2];
rz(-0.76836821) q[3];
sx q[3];
rz(-2.1185015) q[3];
sx q[3];
rz(-0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(2.0098861) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(2.0057604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5539615) q[0];
sx q[0];
rz(-2.580664) q[0];
sx q[0];
rz(-1.8857303) q[0];
rz(-pi) q[1];
rz(-1.5485974) q[2];
sx q[2];
rz(-2.0054842) q[2];
sx q[2];
rz(-1.7322025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.218704) q[1];
sx q[1];
rz(-1.522038) q[1];
sx q[1];
rz(2.7152275) q[1];
x q[2];
rz(-2.2215861) q[3];
sx q[3];
rz(-1.2244867) q[3];
sx q[3];
rz(-1.8370093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-0.87289587) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(1.2783485) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(1.1970253) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7322757) q[0];
sx q[0];
rz(-2.0079552) q[0];
sx q[0];
rz(-2.2073295) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0212101) q[2];
sx q[2];
rz(-1.3345846) q[2];
sx q[2];
rz(-1.2506968) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3665109) q[1];
sx q[1];
rz(-1.0771891) q[1];
sx q[1];
rz(-2.2713186) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5246546) q[3];
sx q[3];
rz(-2.5411798) q[3];
sx q[3];
rz(2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(-1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(1.21571) q[2];
sx q[2];
rz(-1.9234895) q[2];
sx q[2];
rz(-1.1870155) q[2];
rz(-0.79896169) q[3];
sx q[3];
rz(-0.81294717) q[3];
sx q[3];
rz(-3.0019928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
