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
rz(-1.4323956) q[0];
sx q[0];
rz(-0.30269912) q[0];
sx q[0];
rz(-2.1085289) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(-1.6515825) q[1];
sx q[1];
rz(0.19436793) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8914999) q[0];
sx q[0];
rz(-0.95026267) q[0];
sx q[0];
rz(-1.7583855) q[0];
rz(1.871045) q[2];
sx q[2];
rz(-0.89648834) q[2];
sx q[2];
rz(-2.6943494) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3157881) q[1];
sx q[1];
rz(-1.5362843) q[1];
sx q[1];
rz(1.6170349) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5712272) q[3];
sx q[3];
rz(-0.70934848) q[3];
sx q[3];
rz(-1.589243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4091461) q[2];
sx q[2];
rz(-1.1017841) q[2];
sx q[2];
rz(-0.5644325) q[2];
rz(-1.8730646) q[3];
sx q[3];
rz(-3.1334435) q[3];
sx q[3];
rz(0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0821575) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(-0.91307688) q[0];
rz(2.4572241) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(-2.2025462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371999) q[0];
sx q[0];
rz(-0.6966204) q[0];
sx q[0];
rz(-2.533769) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5839416) q[2];
sx q[2];
rz(-2.0347383) q[2];
sx q[2];
rz(-0.23255238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1176743) q[1];
sx q[1];
rz(-2.597375) q[1];
sx q[1];
rz(2.2571889) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9775852) q[3];
sx q[3];
rz(-1.6322246) q[3];
sx q[3];
rz(0.2964501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3188476) q[2];
sx q[2];
rz(-0.007195909) q[2];
sx q[2];
rz(-0.89246559) q[2];
rz(-1.578963) q[3];
sx q[3];
rz(-3.1171418) q[3];
sx q[3];
rz(1.310937) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7398294) q[0];
sx q[0];
rz(-2.6391397) q[0];
sx q[0];
rz(-2.7318562) q[0];
rz(-3.1306664) q[1];
sx q[1];
rz(-2.9210563) q[1];
sx q[1];
rz(-1.3440557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9896802) q[0];
sx q[0];
rz(-1.5912959) q[0];
sx q[0];
rz(-2.0818039) q[0];
rz(-pi) q[1];
rz(-0.063882752) q[2];
sx q[2];
rz(-1.6860355) q[2];
sx q[2];
rz(1.921953) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.037087203) q[1];
sx q[1];
rz(-1.5059789) q[1];
sx q[1];
rz(-1.8291891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6860754) q[3];
sx q[3];
rz(-1.56101) q[3];
sx q[3];
rz(0.20231314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4695796) q[2];
sx q[2];
rz(-1.6894222) q[2];
sx q[2];
rz(-3.0984042) q[2];
rz(-2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(-1.202701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7286872) q[0];
sx q[0];
rz(-0.40189704) q[0];
sx q[0];
rz(-1.3380949) q[0];
rz(2.4463704) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-2.7241838) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5240751) q[0];
sx q[0];
rz(-1.39912) q[0];
sx q[0];
rz(0.89069946) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2171793) q[2];
sx q[2];
rz(-1.355339) q[2];
sx q[2];
rz(-2.5417823) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96046916) q[1];
sx q[1];
rz(-2.4669018) q[1];
sx q[1];
rz(1.8974278) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2715075) q[3];
sx q[3];
rz(-2.6911754) q[3];
sx q[3];
rz(-1.3189486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43856105) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(1.38928) q[2];
rz(2.321068) q[3];
sx q[3];
rz(-1.5912143) q[3];
sx q[3];
rz(-2.4337721) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97838068) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(0.96330825) q[0];
rz(-0.18927255) q[1];
sx q[1];
rz(-0.015733868) q[1];
sx q[1];
rz(0.16545573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2509569) q[0];
sx q[0];
rz(-1.5971308) q[0];
sx q[0];
rz(-3.0980397) q[0];
x q[1];
rz(-0.87617434) q[2];
sx q[2];
rz(-1.3378222) q[2];
sx q[2];
rz(-2.7829952) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.2442064) q[1];
sx q[1];
rz(-1.6868949) q[1];
sx q[1];
rz(-0.80516025) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0040163) q[3];
sx q[3];
rz(-1.8611844) q[3];
sx q[3];
rz(-0.079513915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.307622) q[2];
sx q[2];
rz(-2.622719) q[2];
sx q[2];
rz(2.0818254) q[2];
rz(-0.69093949) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(1.2714269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6178013) q[0];
sx q[0];
rz(-3.0484564) q[0];
sx q[0];
rz(1.6049438) q[0];
rz(2.6887584) q[1];
sx q[1];
rz(-0.0080527877) q[1];
sx q[1];
rz(1.7013928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.207841) q[0];
sx q[0];
rz(-1.4086723) q[0];
sx q[0];
rz(0.16926814) q[0];
rz(-pi) q[1];
rz(1.3511436) q[2];
sx q[2];
rz(-1.7872602) q[2];
sx q[2];
rz(0.6714657) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78919856) q[1];
sx q[1];
rz(-1.4161613) q[1];
sx q[1];
rz(1.5782389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0470508) q[3];
sx q[3];
rz(-0.19652995) q[3];
sx q[3];
rz(0.2257589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3692533) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(3.0625647) q[2];
rz(-2.2760271) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(1.0237833) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2257776) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(-2.916577) q[0];
rz(2.8354538) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(2.2711145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3358094) q[0];
sx q[0];
rz(-0.085698232) q[0];
sx q[0];
rz(2.5524469) q[0];
rz(-2.4267459) q[2];
sx q[2];
rz(-2.0185592) q[2];
sx q[2];
rz(3.062742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5600109) q[1];
sx q[1];
rz(-2.3696179) q[1];
sx q[1];
rz(-3.0549269) q[1];
x q[2];
rz(2.7288658) q[3];
sx q[3];
rz(-1.1201292) q[3];
sx q[3];
rz(0.13201763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47591448) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(-1.0352943) q[2];
rz(2.3546442) q[3];
sx q[3];
rz(-0.042777177) q[3];
sx q[3];
rz(2.8662445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03054522) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(3.1159478) q[0];
rz(1.2643087) q[1];
sx q[1];
rz(-3.1176716) q[1];
sx q[1];
rz(2.4425676) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5584939) q[0];
sx q[0];
rz(-1.8374658) q[0];
sx q[0];
rz(-1.8015652) q[0];
rz(2.672926) q[2];
sx q[2];
rz(-1.7231517) q[2];
sx q[2];
rz(0.050616905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6104709) q[1];
sx q[1];
rz(-1.6446596) q[1];
sx q[1];
rz(-0.38458891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5949536) q[3];
sx q[3];
rz(-0.99144236) q[3];
sx q[3];
rz(-0.97673113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8890624) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(-2.8177281) q[2];
rz(-2.9845386) q[3];
sx q[3];
rz(-1.2374977) q[3];
sx q[3];
rz(0.15373716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.5667628) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(0.59004849) q[0];
rz(2.4216962) q[1];
sx q[1];
rz(-0.059559278) q[1];
sx q[1];
rz(2.2743026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27408263) q[0];
sx q[0];
rz(-1.6757586) q[0];
sx q[0];
rz(2.5068953) q[0];
x q[1];
rz(-2.8810416) q[2];
sx q[2];
rz(-2.1571967) q[2];
sx q[2];
rz(-1.150765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5979851) q[1];
sx q[1];
rz(-1.5026918) q[1];
sx q[1];
rz(-1.4582514) q[1];
rz(-0.18901029) q[3];
sx q[3];
rz(-0.86652404) q[3];
sx q[3];
rz(2.6420322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94816339) q[2];
sx q[2];
rz(-0.89274222) q[2];
sx q[2];
rz(3.0286922) q[2];
rz(2.0924977) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(2.404864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0650487) q[0];
sx q[0];
rz(-1.5815409) q[0];
sx q[0];
rz(-1.6381868) q[0];
rz(0.63649559) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(-1.5719302) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51278567) q[0];
sx q[0];
rz(-0.11462002) q[0];
sx q[0];
rz(2.3324469) q[0];
rz(-3.1384597) q[2];
sx q[2];
rz(-1.571448) q[2];
sx q[2];
rz(2.3359483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.43409882) q[1];
sx q[1];
rz(-1.5705622) q[1];
sx q[1];
rz(-3.1386761) q[1];
rz(-pi) q[2];
rz(0.51497634) q[3];
sx q[3];
rz(-1.712821) q[3];
sx q[3];
rz(1.9076626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87127176) q[2];
sx q[2];
rz(-1.6297636) q[2];
sx q[2];
rz(-2.9809269) q[2];
rz(1.2284944) q[3];
sx q[3];
rz(-0.03385032) q[3];
sx q[3];
rz(0.20139774) q[3];
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
rz(pi/2) q[2];
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
rz(-1.5184215) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(-3.1186947) q[2];
sx q[2];
rz(-1.2780634) q[2];
sx q[2];
rz(-2.7414049) q[2];
rz(2.992441) q[3];
sx q[3];
rz(-1.1356529) q[3];
sx q[3];
rz(0.54855772) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
