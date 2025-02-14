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
rz(-1.6515825) q[1];
sx q[1];
rz(0.19436793) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5759566) q[0];
sx q[0];
rz(-0.64467421) q[0];
sx q[0];
rz(2.8863532) q[0];
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
sx q[0];
rz(pi/2) q[0];
rz(0.82580457) q[1];
sx q[1];
rz(-1.5362843) q[1];
sx q[1];
rz(-1.6170349) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5159149) q[3];
sx q[3];
rz(-1.9301658) q[3];
sx q[3];
rz(0.43454447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73244652) q[2];
sx q[2];
rz(-2.0398085) q[2];
sx q[2];
rz(-0.5644325) q[2];
rz(1.268528) q[3];
sx q[3];
rz(-0.0081491834) q[3];
sx q[3];
rz(-0.90659365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059435189) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(2.2285158) q[0];
rz(2.4572241) q[1];
sx q[1];
rz(-3.1412536) q[1];
sx q[1];
rz(2.2025462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6093151) q[0];
sx q[0];
rz(-2.1255204) q[0];
sx q[0];
rz(1.1251262) q[0];
rz(-pi) q[1];
x q[1];
rz(0.026264391) q[2];
sx q[2];
rz(-0.46411465) q[2];
sx q[2];
rz(0.26192203) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0239183) q[1];
sx q[1];
rz(-2.597375) q[1];
sx q[1];
rz(-0.88440374) q[1];
rz(0.066870832) q[3];
sx q[3];
rz(-1.9767728) q[3];
sx q[3];
rz(-1.8936881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3188476) q[2];
sx q[2];
rz(-0.007195909) q[2];
sx q[2];
rz(-2.2491271) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-3.1171418) q[3];
sx q[3];
rz(-1.310937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40176323) q[0];
sx q[0];
rz(-0.50245291) q[0];
sx q[0];
rz(0.40973642) q[0];
rz(0.01092625) q[1];
sx q[1];
rz(-2.9210563) q[1];
sx q[1];
rz(1.797537) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15191244) q[0];
sx q[0];
rz(-1.5502968) q[0];
sx q[0];
rz(2.0818039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4553237) q[2];
sx q[2];
rz(-1.5073379) q[2];
sx q[2];
rz(2.7977914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.037087203) q[1];
sx q[1];
rz(-1.6356138) q[1];
sx q[1];
rz(-1.3124036) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6860754) q[3];
sx q[3];
rz(-1.5805827) q[3];
sx q[3];
rz(2.9392795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.672013) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(3.0984042) q[2];
rz(2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(-1.9388916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7286872) q[0];
sx q[0];
rz(-2.7396956) q[0];
sx q[0];
rz(1.3380949) q[0];
rz(2.4463704) q[1];
sx q[1];
rz(-0.1278563) q[1];
sx q[1];
rz(2.7241838) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090583853) q[0];
sx q[0];
rz(-2.2390597) q[0];
sx q[0];
rz(2.9221852) q[0];
rz(-2.874006) q[2];
sx q[2];
rz(-0.94174615) q[2];
sx q[2];
rz(-0.81105328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.771811) q[1];
sx q[1];
rz(-0.93771781) q[1];
sx q[1];
rz(-2.8903583) q[1];
x q[2];
rz(-1.2715075) q[3];
sx q[3];
rz(-0.45041725) q[3];
sx q[3];
rz(-1.822644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7030316) q[2];
sx q[2];
rz(-2.265354) q[2];
sx q[2];
rz(-1.38928) q[2];
rz(0.82052463) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(0.70782053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-2.163212) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(-0.96330825) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(-2.9761369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4625797) q[0];
sx q[0];
rz(-1.5272584) q[0];
sx q[0];
rz(1.5444368) q[0];
rz(-pi) q[1];
rz(0.29954977) q[2];
sx q[2];
rz(-2.2431157) q[2];
sx q[2];
rz(1.7393665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4374809) q[1];
sx q[1];
rz(-0.81161122) q[1];
sx q[1];
rz(2.981217) q[1];
x q[2];
rz(-0.31807301) q[3];
sx q[3];
rz(-1.1568562) q[3];
sx q[3];
rz(1.3596168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8339707) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(1.0597672) q[2];
rz(2.4506532) q[3];
sx q[3];
rz(-0.86425224) q[3];
sx q[3];
rz(-1.2714269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6178013) q[0];
sx q[0];
rz(-0.093136223) q[0];
sx q[0];
rz(1.5366489) q[0];
rz(-0.45283428) q[1];
sx q[1];
rz(-3.1335399) q[1];
sx q[1];
rz(1.4401999) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0217845) q[0];
sx q[0];
rz(-0.23384604) q[0];
sx q[0];
rz(0.77063693) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9199706) q[2];
sx q[2];
rz(-1.3563507) q[2];
sx q[2];
rz(-2.1943486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3611412) q[1];
sx q[1];
rz(-1.5781501) q[1];
sx q[1];
rz(2.9869534) q[1];
rz(-2.0945419) q[3];
sx q[3];
rz(-0.19652995) q[3];
sx q[3];
rz(0.2257589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77233934) q[2];
sx q[2];
rz(-0.99268308) q[2];
sx q[2];
rz(0.079027979) q[2];
rz(0.8655656) q[3];
sx q[3];
rz(-0.95439684) q[3];
sx q[3];
rz(-1.0237833) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2257776) q[0];
sx q[0];
rz(-0.0053891698) q[0];
sx q[0];
rz(-0.22501568) q[0];
rz(-2.8354538) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(-2.2711145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8057833) q[0];
sx q[0];
rz(-3.0558944) q[0];
sx q[0];
rz(-2.5524469) q[0];
rz(-1.00433) q[2];
sx q[2];
rz(-2.2029467) q[2];
sx q[2];
rz(1.2901778) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0902033) q[1];
sx q[1];
rz(-1.5103814) q[1];
sx q[1];
rz(-0.77009542) q[1];
rz(-pi) q[2];
rz(0.41272687) q[3];
sx q[3];
rz(-1.1201292) q[3];
sx q[3];
rz(-0.13201763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6656782) q[2];
sx q[2];
rz(-0.12683882) q[2];
sx q[2];
rz(2.1062984) q[2];
rz(0.78694844) q[3];
sx q[3];
rz(-0.042777177) q[3];
sx q[3];
rz(-2.8662445) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1110474) q[0];
sx q[0];
rz(-0.011140911) q[0];
sx q[0];
rz(-0.025644843) q[0];
rz(1.2643087) q[1];
sx q[1];
rz(-3.1176716) q[1];
sx q[1];
rz(-0.69902507) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8549439) q[0];
sx q[0];
rz(-0.35084769) q[0];
sx q[0];
rz(2.4445266) q[0];
x q[1];
rz(-1.7412296) q[2];
sx q[2];
rz(-2.0336069) q[2];
sx q[2];
rz(-1.544726) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.1406724) q[1];
sx q[1];
rz(-2.7503221) q[1];
sx q[1];
rz(2.9468582) q[1];
rz(-0.54663901) q[3];
sx q[3];
rz(-0.99144236) q[3];
sx q[3];
rz(2.1648615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.25253025) q[2];
sx q[2];
rz(-2.3928596) q[2];
sx q[2];
rz(-2.8177281) q[2];
rz(0.1570541) q[3];
sx q[3];
rz(-1.9040949) q[3];
sx q[3];
rz(-0.15373716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57482982) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(-2.5515442) q[0];
rz(2.4216962) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(-2.2743026) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4380505) q[0];
sx q[0];
rz(-2.4994586) q[0];
sx q[0];
rz(0.17583986) q[0];
rz(0.9684674) q[2];
sx q[2];
rz(-1.3545389) q[2];
sx q[2];
rz(2.575084) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.572435) q[1];
sx q[1];
rz(-0.13147239) q[1];
sx q[1];
rz(1.0249895) q[1];
rz(-0.18901029) q[3];
sx q[3];
rz(-0.86652404) q[3];
sx q[3];
rz(-0.49956043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.94816339) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(3.0286922) q[2];
rz(-1.0490949) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(2.404864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.0650487) q[0];
sx q[0];
rz(-1.5600518) q[0];
sx q[0];
rz(1.5034058) q[0];
rz(-2.5050971) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(1.5696625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51278567) q[0];
sx q[0];
rz(-0.11462002) q[0];
sx q[0];
rz(-0.8091457) q[0];
rz(-1.571448) q[2];
sx q[2];
rz(-1.5739292) q[2];
sx q[2];
rz(-2.3764386) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7074938) q[1];
sx q[1];
rz(-1.5705622) q[1];
sx q[1];
rz(-3.1386761) q[1];
x q[2];
rz(2.8590389) q[3];
sx q[3];
rz(-2.6090949) q[3];
sx q[3];
rz(3.0498216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87127176) q[2];
sx q[2];
rz(-1.6297636) q[2];
sx q[2];
rz(0.16066571) q[2];
rz(-1.2284944) q[3];
sx q[3];
rz(-3.1077423) q[3];
sx q[3];
rz(0.20139774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(1.6231712) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(1.2779909) q[2];
sx q[2];
rz(-1.5488727) q[2];
sx q[2];
rz(1.9643754) q[2];
rz(2.0102228) q[3];
sx q[3];
rz(-1.4356339) q[3];
sx q[3];
rz(2.1826134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
