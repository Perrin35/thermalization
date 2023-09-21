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
rz(2.9449129) q[0];
sx q[0];
rz(11.37698) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(2.300188) q[1];
sx q[1];
rz(9.0471164) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.022764) q[0];
sx q[0];
rz(-1.6435992) q[0];
sx q[0];
rz(-1.4490436) q[0];
x q[1];
rz(-1.2871735) q[2];
sx q[2];
rz(-1.5986773) q[2];
sx q[2];
rz(-1.4717799) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79757534) q[1];
sx q[1];
rz(-1.1317778) q[1];
sx q[1];
rz(-2.3637191) q[1];
x q[2];
rz(0.30794413) q[3];
sx q[3];
rz(-0.43711284) q[3];
sx q[3];
rz(1.8969769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.010014023) q[2];
sx q[2];
rz(-0.49396124) q[2];
sx q[2];
rz(2.4689891) q[2];
rz(2.9721695) q[3];
sx q[3];
rz(-0.38893458) q[3];
sx q[3];
rz(-1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(-0.22856523) q[0];
rz(-2.9810492) q[1];
sx q[1];
rz(-1.7030145) q[1];
sx q[1];
rz(0.28796089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9790968) q[0];
sx q[0];
rz(-1.7485577) q[0];
sx q[0];
rz(1.3773247) q[0];
x q[1];
rz(-2.1177883) q[2];
sx q[2];
rz(-2.4734801) q[2];
sx q[2];
rz(-2.748865) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2634695) q[1];
sx q[1];
rz(-2.1681004) q[1];
sx q[1];
rz(-1.0695446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35034758) q[3];
sx q[3];
rz(-1.5347893) q[3];
sx q[3];
rz(-0.85067526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59445375) q[2];
sx q[2];
rz(-1.889166) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32524747) q[0];
sx q[0];
rz(-2.6419817) q[0];
sx q[0];
rz(-0.26741272) q[0];
rz(-1.422241) q[1];
sx q[1];
rz(-2.0098675) q[1];
sx q[1];
rz(0.95169383) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0591473) q[0];
sx q[0];
rz(-1.7623616) q[0];
sx q[0];
rz(-3.0491769) q[0];
rz(-1.1859602) q[2];
sx q[2];
rz(-0.54883146) q[2];
sx q[2];
rz(1.2518815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.82745951) q[1];
sx q[1];
rz(-0.810312) q[1];
sx q[1];
rz(2.00287) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6056616) q[3];
sx q[3];
rz(-1.0169573) q[3];
sx q[3];
rz(-1.2371847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.6458076) q[2];
sx q[2];
rz(-2.7974131) q[2];
rz(0.88642818) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-0.56604958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.3543509) q[1];
sx q[1];
rz(0.37277645) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5107721) q[0];
sx q[0];
rz(-1.6850867) q[0];
sx q[0];
rz(-2.3833016) q[0];
rz(-pi) q[1];
rz(0.4332306) q[2];
sx q[2];
rz(-2.3150812) q[2];
sx q[2];
rz(-2.7162958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8744295) q[1];
sx q[1];
rz(-2.0544555) q[1];
sx q[1];
rz(-2.8664385) q[1];
rz(2.8126206) q[3];
sx q[3];
rz(-2.0250468) q[3];
sx q[3];
rz(-0.3262375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5059775) q[2];
sx q[2];
rz(-1.6515235) q[2];
sx q[2];
rz(0.90488952) q[2];
rz(2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-0.99159616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6624517) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(-1.8537846) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.9961424) q[1];
sx q[1];
rz(-0.53422654) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9984765) q[0];
sx q[0];
rz(-1.8792218) q[0];
sx q[0];
rz(0.31750676) q[0];
rz(-pi) q[1];
rz(0.34110951) q[2];
sx q[2];
rz(-1.1030536) q[2];
sx q[2];
rz(-1.9175921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91671645) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(-3.0776943) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71877919) q[3];
sx q[3];
rz(-1.4712417) q[3];
sx q[3];
rz(1.0574785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6405032) q[2];
sx q[2];
rz(-1.2155632) q[2];
sx q[2];
rz(-2.9980998) q[2];
rz(1.7701373) q[3];
sx q[3];
rz(-1.8618795) q[3];
sx q[3];
rz(2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4390398) q[0];
sx q[0];
rz(-0.87600231) q[0];
sx q[0];
rz(2.904536) q[0];
rz(1.9006231) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(-1.7664849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28778827) q[0];
sx q[0];
rz(-1.3082936) q[0];
sx q[0];
rz(0.35065035) q[0];
x q[1];
rz(-1.4985956) q[2];
sx q[2];
rz(-2.318396) q[2];
sx q[2];
rz(-2.2526134) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4543537) q[1];
sx q[1];
rz(-0.41077405) q[1];
sx q[1];
rz(-0.06128581) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70011219) q[3];
sx q[3];
rz(-2.5317319) q[3];
sx q[3];
rz(1.3175347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(-3.1249629) q[3];
sx q[3];
rz(-0.36706585) q[3];
sx q[3];
rz(0.19255157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49611133) q[0];
sx q[0];
rz(-0.2922903) q[0];
sx q[0];
rz(-2.2684229) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-2.0261804) q[1];
sx q[1];
rz(-0.08392863) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8221995) q[0];
sx q[0];
rz(-0.60395741) q[0];
sx q[0];
rz(1.970406) q[0];
x q[1];
rz(-1.7882118) q[2];
sx q[2];
rz(-1.916421) q[2];
sx q[2];
rz(0.84920041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2578656) q[1];
sx q[1];
rz(-0.95974937) q[1];
sx q[1];
rz(-0.74531598) q[1];
rz(0.84079068) q[3];
sx q[3];
rz(-1.6222686) q[3];
sx q[3];
rz(1.7712902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7579047) q[2];
sx q[2];
rz(-0.44054511) q[2];
sx q[2];
rz(2.596358) q[2];
rz(2.7111354) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(-1.9301201) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-0.5935697) q[1];
sx q[1];
rz(2.1957695) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4771172) q[0];
sx q[0];
rz(-0.23453377) q[0];
sx q[0];
rz(-0.3420591) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2875167) q[2];
sx q[2];
rz(-2.1207361) q[2];
sx q[2];
rz(2.9727109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2266027) q[1];
sx q[1];
rz(-1.9493002) q[1];
sx q[1];
rz(1.539991) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66611992) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.7249858) q[2];
sx q[2];
rz(-2.8611709) q[2];
rz(0.60493207) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(0.98208565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9873001) q[0];
sx q[0];
rz(-0.91006088) q[0];
sx q[0];
rz(-2.5262685) q[0];
rz(-0.92957169) q[1];
sx q[1];
rz(-0.85837448) q[1];
sx q[1];
rz(-2.5659134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79918039) q[0];
sx q[0];
rz(-0.77501446) q[0];
sx q[0];
rz(-2.4393625) q[0];
x q[1];
rz(-3.028308) q[2];
sx q[2];
rz(-1.7656529) q[2];
sx q[2];
rz(1.726113) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5513735) q[1];
sx q[1];
rz(-3.0312523) q[1];
sx q[1];
rz(-1.8054086) q[1];
rz(2.9541295) q[3];
sx q[3];
rz(-2.1037357) q[3];
sx q[3];
rz(0.77123469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82970396) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(-3.1196307) q[2];
rz(-2.9711376) q[3];
sx q[3];
rz(-1.0237834) q[3];
sx q[3];
rz(2.8767265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72572529) q[0];
sx q[0];
rz(-2.2150345) q[0];
sx q[0];
rz(2.5073994) q[0];
rz(3.1126853) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(2.172487) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6432188) q[0];
sx q[0];
rz(-1.507483) q[0];
sx q[0];
rz(-1.4988585) q[0];
rz(-pi) q[1];
rz(-1.3766039) q[2];
sx q[2];
rz(-2.3079434) q[2];
sx q[2];
rz(3.0498691) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9173968) q[1];
sx q[1];
rz(-2.2317413) q[1];
sx q[1];
rz(-1.7978653) q[1];
rz(-pi) q[2];
rz(0.5523596) q[3];
sx q[3];
rz(-2.6354355) q[3];
sx q[3];
rz(-1.5137223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4518296) q[2];
sx q[2];
rz(-1.657106) q[2];
sx q[2];
rz(2.7815212) q[2];
rz(1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(2.578919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2071028) q[0];
sx q[0];
rz(-1.5710545) q[0];
sx q[0];
rz(1.5221773) q[0];
rz(0.044152505) q[1];
sx q[1];
rz(-1.6828729) q[1];
sx q[1];
rz(2.0353459) q[1];
rz(2.2864441) q[2];
sx q[2];
rz(-0.48968857) q[2];
sx q[2];
rz(-0.80077632) q[2];
rz(-0.45331656) q[3];
sx q[3];
rz(-1.8554056) q[3];
sx q[3];
rz(1.5649395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
