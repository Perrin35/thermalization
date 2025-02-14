OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.709197) q[0];
sx q[0];
rz(-2.8388935) q[0];
sx q[0];
rz(2.1085289) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(1.4900102) q[1];
sx q[1];
rz(9.23041) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5759566) q[0];
sx q[0];
rz(-2.4969184) q[0];
sx q[0];
rz(-0.25523941) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4448574) q[2];
sx q[2];
rz(-1.8039304) q[2];
sx q[2];
rz(1.3145043) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3157881) q[1];
sx q[1];
rz(-1.5362843) q[1];
sx q[1];
rz(-1.6170349) q[1];
x q[2];
rz(0.62567775) q[3];
sx q[3];
rz(-1.2114269) q[3];
sx q[3];
rz(2.7070482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4091461) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(-2.5771602) q[2];
rz(-1.268528) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0821575) q[0];
sx q[0];
rz(-0.0075639021) q[0];
sx q[0];
rz(-0.91307688) q[0];
rz(-0.68436855) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(-2.2025462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6093151) q[0];
sx q[0];
rz(-1.0160722) q[0];
sx q[0];
rz(1.1251262) q[0];
rz(-pi) q[1];
rz(1.5576511) q[2];
sx q[2];
rz(-2.0347383) q[2];
sx q[2];
rz(-0.23255238) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0239183) q[1];
sx q[1];
rz(-0.54421762) q[1];
sx q[1];
rz(2.2571889) q[1];
rz(-1.9775852) q[3];
sx q[3];
rz(-1.6322246) q[3];
sx q[3];
rz(-0.2964501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3188476) q[2];
sx q[2];
rz(-0.007195909) q[2];
sx q[2];
rz(2.2491271) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(-1.8306556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398294) q[0];
sx q[0];
rz(-2.6391397) q[0];
sx q[0];
rz(-2.7318562) q[0];
rz(3.1306664) q[1];
sx q[1];
rz(-0.22053638) q[1];
sx q[1];
rz(-1.3440557) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7342012) q[0];
sx q[0];
rz(-1.0599066) q[0];
sx q[0];
rz(3.1180918) q[0];
rz(1.686269) q[2];
sx q[2];
rz(-1.6342548) q[2];
sx q[2];
rz(0.34380128) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.037087203) q[1];
sx q[1];
rz(-1.5059789) q[1];
sx q[1];
rz(-1.8291891) q[1];
rz(-pi) q[2];
rz(1.4555172) q[3];
sx q[3];
rz(-1.5805827) q[3];
sx q[3];
rz(-2.9392795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4695796) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(0.043188485) q[2];
rz(1.1399266) q[3];
sx q[3];
rz(-0.15634263) q[3];
sx q[3];
rz(1.202701) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4129055) q[0];
sx q[0];
rz(-2.7396956) q[0];
sx q[0];
rz(1.3380949) q[0];
rz(2.4463704) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-2.7241838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0510088) q[0];
sx q[0];
rz(-2.2390597) q[0];
sx q[0];
rz(-2.9221852) q[0];
x q[1];
rz(2.2171793) q[2];
sx q[2];
rz(-1.7862537) q[2];
sx q[2];
rz(0.59981031) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1811235) q[1];
sx q[1];
rz(-2.4669018) q[1];
sx q[1];
rz(1.8974278) q[1];
x q[2];
rz(-1.2715075) q[3];
sx q[3];
rz(-0.45041725) q[3];
sx q[3];
rz(1.3189486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7030316) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(-1.38928) q[2];
rz(-2.321068) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(0.70782053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97838068) q[0];
sx q[0];
rz(-0.037540171) q[0];
sx q[0];
rz(-2.1782844) q[0];
rz(0.18927255) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(0.16545573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8906358) q[0];
sx q[0];
rz(-1.5444618) q[0];
sx q[0];
rz(-0.043553003) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8420429) q[2];
sx q[2];
rz(-2.2431157) q[2];
sx q[2];
rz(1.7393665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8973863) q[1];
sx q[1];
rz(-1.4546977) q[1];
sx q[1];
rz(0.80516025) q[1];
rz(-0.31807301) q[3];
sx q[3];
rz(-1.1568562) q[3];
sx q[3];
rz(1.3596168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.307622) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(-1.0597672) q[2];
rz(-0.69093949) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(1.2714269) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5237913) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(1.6049438) q[0];
rz(-0.45283428) q[1];
sx q[1];
rz(-0.0080527877) q[1];
sx q[1];
rz(1.7013928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0217845) q[0];
sx q[0];
rz(-2.9077466) q[0];
sx q[0];
rz(2.3709557) q[0];
rz(-pi) q[1];
rz(-1.790449) q[2];
sx q[2];
rz(-1.3543324) q[2];
sx q[2];
rz(-0.6714657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78919856) q[1];
sx q[1];
rz(-1.7254313) q[1];
sx q[1];
rz(-1.5633538) q[1];
rz(0.099248107) q[3];
sx q[3];
rz(-1.4008879) q[3];
sx q[3];
rz(-2.3836294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.77233934) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(3.0625647) q[2];
rz(0.8655656) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(1.0237833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2257776) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(0.22501568) q[0];
rz(-0.30613884) q[1];
sx q[1];
rz(-0.016751079) q[1];
sx q[1];
rz(-2.2711145) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2149379) q[0];
sx q[0];
rz(-1.4995725) q[0];
sx q[0];
rz(1.5230973) q[0];
x q[1];
rz(0.6323496) q[2];
sx q[2];
rz(-0.82201695) q[2];
sx q[2];
rz(2.1123304) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5600109) q[1];
sx q[1];
rz(-2.3696179) q[1];
sx q[1];
rz(3.0549269) q[1];
rz(0.41272687) q[3];
sx q[3];
rz(-2.0214635) q[3];
sx q[3];
rz(-3.009575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47591448) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(2.1062984) q[2];
rz(-0.78694844) q[3];
sx q[3];
rz(-3.0988155) q[3];
sx q[3];
rz(0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1110474) q[0];
sx q[0];
rz(-3.1304517) q[0];
sx q[0];
rz(0.025644843) q[0];
rz(-1.2643087) q[1];
sx q[1];
rz(-0.023921078) q[1];
sx q[1];
rz(-0.69902507) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(-1.7932685) q[0];
sx q[0];
rz(-0.27357863) q[0];
rz(-2.8139171) q[2];
sx q[2];
rz(-2.6505396) q[2];
sx q[2];
rz(1.2289405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5311218) q[1];
sx q[1];
rz(-1.4969331) q[1];
sx q[1];
rz(2.7570037) q[1];
rz(-2.5949536) q[3];
sx q[3];
rz(-2.1501503) q[3];
sx q[3];
rz(2.1648615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25253025) q[2];
sx q[2];
rz(-0.74873304) q[2];
sx q[2];
rz(-0.32386455) q[2];
rz(-2.9845386) q[3];
sx q[3];
rz(-1.2374977) q[3];
sx q[3];
rz(-2.9878555) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57482982) q[0];
sx q[0];
rz(-0.014641849) q[0];
sx q[0];
rz(-2.5515442) q[0];
rz(2.4216962) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(0.86729008) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9218719) q[0];
sx q[0];
rz(-0.94014535) q[0];
sx q[0];
rz(-1.4407071) q[0];
x q[1];
rz(-1.9407104) q[2];
sx q[2];
rz(-0.6354161) q[2];
sx q[2];
rz(0.70178343) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54360753) q[1];
sx q[1];
rz(-1.5026918) q[1];
sx q[1];
rz(-1.4582514) q[1];
rz(-pi) q[2];
rz(0.18901029) q[3];
sx q[3];
rz(-2.2750686) q[3];
sx q[3];
rz(2.6420322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1934293) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(3.0286922) q[2];
rz(2.0924977) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(-0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.076544) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(-1.6381868) q[0];
rz(0.63649559) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(-1.5719302) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8419452) q[0];
sx q[0];
rz(-1.491786) q[0];
sx q[0];
rz(1.6539198) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9365065) q[2];
sx q[2];
rz(-0.0031999667) q[2];
sx q[2];
rz(0.56006685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0565785) q[1];
sx q[1];
rz(-3.1386668) q[1];
sx q[1];
rz(-3.0614733) q[1];
x q[2];
rz(2.8590389) q[3];
sx q[3];
rz(-2.6090949) q[3];
sx q[3];
rz(3.0498216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87127176) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(2.9809269) q[2];
rz(1.9130982) q[3];
sx q[3];
rz(-3.1077423) q[3];
sx q[3];
rz(0.20139774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092030839) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(-1.6231712) q[1];
sx q[1];
rz(-2.7802614) q[1];
sx q[1];
rz(-2.8526715) q[1];
rz(-3.1186947) q[2];
sx q[2];
rz(-1.2780634) q[2];
sx q[2];
rz(-2.7414049) q[2];
rz(-1.8801943) q[3];
sx q[3];
rz(-2.6831476) q[3];
sx q[3];
rz(0.89098709) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
