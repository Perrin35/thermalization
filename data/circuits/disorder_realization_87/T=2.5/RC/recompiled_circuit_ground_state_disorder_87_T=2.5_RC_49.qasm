OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3992231) q[0];
sx q[0];
rz(-2.7288781) q[0];
sx q[0];
rz(-0.8362008) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(1.0083415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5429687) q[0];
sx q[0];
rz(-2.7592139) q[0];
sx q[0];
rz(1.3578869) q[0];
rz(-pi) q[1];
rz(2.600619) q[2];
sx q[2];
rz(-2.7056575) q[2];
sx q[2];
rz(-1.774314) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8562682) q[1];
sx q[1];
rz(-2.6372077) q[1];
sx q[1];
rz(0.32907246) q[1];
rz(-1.4654245) q[3];
sx q[3];
rz(-2.947315) q[3];
sx q[3];
rz(2.7716178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3081554) q[2];
sx q[2];
rz(-2.1980632) q[2];
sx q[2];
rz(0.65307871) q[2];
rz(0.11450442) q[3];
sx q[3];
rz(-1.7479892) q[3];
sx q[3];
rz(-2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36650518) q[0];
sx q[0];
rz(-2.1108284) q[0];
sx q[0];
rz(-1.7979701) q[0];
rz(1.1603181) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(-2.5210099) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40784971) q[0];
sx q[0];
rz(-1.7622927) q[0];
sx q[0];
rz(1.2194077) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2389675) q[2];
sx q[2];
rz(-1.2764954) q[2];
sx q[2];
rz(-0.49605344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.008604) q[1];
sx q[1];
rz(-0.78215608) q[1];
sx q[1];
rz(-2.6166366) q[1];
rz(-pi) q[2];
rz(-2.1866131) q[3];
sx q[3];
rz(-1.338025) q[3];
sx q[3];
rz(1.7580838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6985942) q[2];
sx q[2];
rz(-1.5195165) q[2];
sx q[2];
rz(0.2300187) q[2];
rz(-1.5864141) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(2.8659099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25552937) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(-2.0607167) q[0];
rz(-2.4983662) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(3.0562775) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9832343) q[0];
sx q[0];
rz(-2.2081349) q[0];
sx q[0];
rz(2.4054224) q[0];
rz(-0.32344748) q[2];
sx q[2];
rz(-0.5964884) q[2];
sx q[2];
rz(1.310854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.18135611) q[1];
sx q[1];
rz(-1.4048368) q[1];
sx q[1];
rz(1.7039429) q[1];
x q[2];
rz(-1.5253228) q[3];
sx q[3];
rz(-2.1915428) q[3];
sx q[3];
rz(2.5705702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2049415) q[2];
sx q[2];
rz(-0.0095657883) q[2];
sx q[2];
rz(-0.19634518) q[2];
rz(-0.28228545) q[3];
sx q[3];
rz(-2.0245656) q[3];
sx q[3];
rz(-2.096368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-2.008217) q[0];
sx q[0];
rz(-2.5947925) q[0];
sx q[0];
rz(0.38544449) q[0];
rz(-1.4133981) q[1];
sx q[1];
rz(-1.9368659) q[1];
sx q[1];
rz(0.24994303) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2521533) q[0];
sx q[0];
rz(-1.982054) q[0];
sx q[0];
rz(2.2338339) q[0];
x q[1];
rz(1.5105627) q[2];
sx q[2];
rz(-1.6585729) q[2];
sx q[2];
rz(-0.46030948) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7467667) q[1];
sx q[1];
rz(-1.4020846) q[1];
sx q[1];
rz(-2.6454703) q[1];
rz(-pi) q[2];
rz(1.3718666) q[3];
sx q[3];
rz(-0.78708157) q[3];
sx q[3];
rz(-2.6347292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6803153) q[2];
sx q[2];
rz(-1.8467434) q[2];
sx q[2];
rz(-2.9052022) q[2];
rz(1.2605028) q[3];
sx q[3];
rz(-2.8915296) q[3];
sx q[3];
rz(0.67729706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33609718) q[0];
sx q[0];
rz(-1.5376872) q[0];
sx q[0];
rz(-2.6564823) q[0];
rz(1.9612034) q[1];
sx q[1];
rz(-1.5943269) q[1];
sx q[1];
rz(0.83650437) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9125875) q[0];
sx q[0];
rz(-1.8857755) q[0];
sx q[0];
rz(-2.7664685) q[0];
x q[1];
rz(2.8233158) q[2];
sx q[2];
rz(-1.1851382) q[2];
sx q[2];
rz(1.9463271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5647392) q[1];
sx q[1];
rz(-2.4702063) q[1];
sx q[1];
rz(1.6480084) q[1];
rz(0.14797609) q[3];
sx q[3];
rz(-2.736909) q[3];
sx q[3];
rz(-1.2860822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89852077) q[2];
sx q[2];
rz(-2.409755) q[2];
sx q[2];
rz(-0.76954976) q[2];
rz(0.92929333) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(0.23269674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.1997851) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(-0.76882452) q[0];
rz(-1.3517514) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(2.2519055) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79737299) q[0];
sx q[0];
rz(-1.3614206) q[0];
sx q[0];
rz(-1.9599509) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6151516) q[2];
sx q[2];
rz(-1.5694478) q[2];
sx q[2];
rz(1.1464455) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8823679) q[1];
sx q[1];
rz(-2.9881757) q[1];
sx q[1];
rz(1.2620366) q[1];
rz(-pi) q[2];
rz(2.9075648) q[3];
sx q[3];
rz(-0.46189538) q[3];
sx q[3];
rz(-1.5398825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49030226) q[2];
sx q[2];
rz(-1.4375261) q[2];
sx q[2];
rz(-1.5865405) q[2];
rz(1.8709315) q[3];
sx q[3];
rz(-0.41897604) q[3];
sx q[3];
rz(-3.0095625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1726058) q[0];
sx q[0];
rz(-2.0033328) q[0];
sx q[0];
rz(-1.1440811) q[0];
rz(-0.83089685) q[1];
sx q[1];
rz(-1.3583207) q[1];
sx q[1];
rz(-2.3544618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2427776) q[0];
sx q[0];
rz(-1.7022331) q[0];
sx q[0];
rz(-1.7064894) q[0];
rz(-pi) q[1];
rz(0.88655858) q[2];
sx q[2];
rz(-2.0066924) q[2];
sx q[2];
rz(-2.2615711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5097376) q[1];
sx q[1];
rz(-0.69586772) q[1];
sx q[1];
rz(2.3515755) q[1];
x q[2];
rz(-1.1803341) q[3];
sx q[3];
rz(-2.0956846) q[3];
sx q[3];
rz(2.386664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1181011) q[2];
sx q[2];
rz(-2.1716437) q[2];
sx q[2];
rz(-3.1084295) q[2];
rz(-1.5432594) q[3];
sx q[3];
rz(-0.71591806) q[3];
sx q[3];
rz(-1.6395114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24932662) q[0];
sx q[0];
rz(-1.5858269) q[0];
sx q[0];
rz(-1.7484885) q[0];
rz(2.6847367) q[1];
sx q[1];
rz(-1.9970048) q[1];
sx q[1];
rz(2.0410062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8580061) q[0];
sx q[0];
rz(-1.4947944) q[0];
sx q[0];
rz(1.572524) q[0];
x q[1];
rz(0.69533657) q[2];
sx q[2];
rz(-1.0461263) q[2];
sx q[2];
rz(-0.52548835) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9683084) q[1];
sx q[1];
rz(-1.6047162) q[1];
sx q[1];
rz(0.90356253) q[1];
rz(-pi) q[2];
rz(-1.4310775) q[3];
sx q[3];
rz(-1.1828239) q[3];
sx q[3];
rz(0.037484976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.042772375) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5808606) q[2];
rz(-0.76357311) q[3];
sx q[3];
rz(-1.5127425) q[3];
sx q[3];
rz(2.1055351) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7324657) q[0];
sx q[0];
rz(-2.3865073) q[0];
sx q[0];
rz(0.29148802) q[0];
rz(-1.905722) q[1];
sx q[1];
rz(-1.196685) q[1];
sx q[1];
rz(3.018697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6345469) q[0];
sx q[0];
rz(-1.6222008) q[0];
sx q[0];
rz(-2.794751) q[0];
rz(2.9192186) q[2];
sx q[2];
rz(-1.8463328) q[2];
sx q[2];
rz(-2.0982519) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1157411) q[1];
sx q[1];
rz(-0.23499712) q[1];
sx q[1];
rz(2.0629289) q[1];
x q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-1.5648309) q[3];
sx q[3];
rz(-0.98814135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.62873944) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(1.0635618) q[2];
rz(-2.9641446) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(-1.0183081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7040831) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(-1.7105239) q[0];
rz(-0.60072947) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(-2.5240555) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4040632) q[0];
sx q[0];
rz(-1.4466219) q[0];
sx q[0];
rz(1.3535208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1328636) q[2];
sx q[2];
rz(-0.8265087) q[2];
sx q[2];
rz(-1.1240608) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.70676196) q[1];
sx q[1];
rz(-1.037957) q[1];
sx q[1];
rz(0.024350655) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57281877) q[3];
sx q[3];
rz(-1.8316934) q[3];
sx q[3];
rz(0.12810055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90126976) q[2];
sx q[2];
rz(-1.0003041) q[2];
sx q[2];
rz(-0.97736248) q[2];
rz(1.0824925) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(-3.0860743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33901535) q[0];
sx q[0];
rz(-2.3727198) q[0];
sx q[0];
rz(2.4045237) q[0];
rz(0.045724178) q[1];
sx q[1];
rz(-0.8174236) q[1];
sx q[1];
rz(1.5541706) q[1];
rz(1.0411714) q[2];
sx q[2];
rz(-2.0176135) q[2];
sx q[2];
rz(1.8420646) q[2];
rz(-1.7405199) q[3];
sx q[3];
rz(-2.2206497) q[3];
sx q[3];
rz(2.7322265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
