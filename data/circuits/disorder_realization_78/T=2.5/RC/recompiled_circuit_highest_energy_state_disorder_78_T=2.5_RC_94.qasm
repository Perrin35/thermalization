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
rz(-3.0566293) q[0];
sx q[0];
rz(-0.30240107) q[0];
sx q[0];
rz(0.032057134) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(-0.75262466) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40850484) q[0];
sx q[0];
rz(-1.8903539) q[0];
sx q[0];
rz(-1.7825141) q[0];
rz(-2.7787894) q[2];
sx q[2];
rz(-1.2191091) q[2];
sx q[2];
rz(-1.2376518) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51576534) q[1];
sx q[1];
rz(-1.2316362) q[1];
sx q[1];
rz(-1.8372507) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9273622) q[3];
sx q[3];
rz(-1.9211244) q[3];
sx q[3];
rz(-0.080834724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(2.3215129) q[2];
rz(0.44102937) q[3];
sx q[3];
rz(-2.0377772) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(0.41020694) q[0];
rz(-2.7606616) q[1];
sx q[1];
rz(-1.5313287) q[1];
sx q[1];
rz(1.7832696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.298734) q[0];
sx q[0];
rz(-0.72745391) q[0];
sx q[0];
rz(0.95746118) q[0];
rz(-pi) q[1];
rz(-0.83186457) q[2];
sx q[2];
rz(-0.68507776) q[2];
sx q[2];
rz(0.58704228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5839148) q[1];
sx q[1];
rz(-1.4951981) q[1];
sx q[1];
rz(1.7390523) q[1];
rz(-pi) q[2];
x q[2];
rz(2.210564) q[3];
sx q[3];
rz(-1.6954386) q[3];
sx q[3];
rz(-2.3384936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7831948) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(-2.6386293) q[2];
rz(1.8424235) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(-2.9097596) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4101039) q[0];
sx q[0];
rz(-0.51173156) q[0];
sx q[0];
rz(-0.9758392) q[0];
rz(2.2052235) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(-0.13793129) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16743037) q[0];
sx q[0];
rz(-2.5241521) q[0];
sx q[0];
rz(0.90888826) q[0];
rz(-pi) q[1];
rz(-0.27020578) q[2];
sx q[2];
rz(-0.41105726) q[2];
sx q[2];
rz(-0.11693987) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5974849) q[1];
sx q[1];
rz(-2.4307837) q[1];
sx q[1];
rz(-1.3686468) q[1];
rz(-0.46584399) q[3];
sx q[3];
rz(-1.3261822) q[3];
sx q[3];
rz(0.96058339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5522573) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(-1.662558) q[2];
rz(0.79834437) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(-3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.289157) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(-1.0668466) q[0];
rz(1.5929219) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(2.7222395) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28946149) q[0];
sx q[0];
rz(-1.3416107) q[0];
sx q[0];
rz(-2.9782692) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4801976) q[2];
sx q[2];
rz(-1.0444469) q[2];
sx q[2];
rz(1.787993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49653445) q[1];
sx q[1];
rz(-0.36318159) q[1];
sx q[1];
rz(1.3892322) q[1];
x q[2];
rz(2.1712915) q[3];
sx q[3];
rz(-0.56823778) q[3];
sx q[3];
rz(2.0214391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7833917) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(-1.5979213) q[2];
rz(1.915043) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(-1.290087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2655547) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(-0.75620404) q[0];
rz(-1.8887695) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(1.7255712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6373581) q[0];
sx q[0];
rz(-1.3535445) q[0];
sx q[0];
rz(-0.603032) q[0];
rz(2.8019287) q[2];
sx q[2];
rz(-2.77861) q[2];
sx q[2];
rz(-1.5756922) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6452613) q[1];
sx q[1];
rz(-1.0320391) q[1];
sx q[1];
rz(-1.5625728) q[1];
rz(-pi) q[2];
rz(-0.46964271) q[3];
sx q[3];
rz(-0.40294632) q[3];
sx q[3];
rz(-0.5136036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.34117928) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(-0.18079147) q[2];
rz(1.6527269) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4329231) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(-0.80023009) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-1.0601284) q[1];
sx q[1];
rz(-3.0175993) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5381966) q[0];
sx q[0];
rz(-1.5992224) q[0];
sx q[0];
rz(-1.70065) q[0];
rz(-1.8154549) q[2];
sx q[2];
rz(-0.24194939) q[2];
sx q[2];
rz(-1.9242147) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2468894) q[1];
sx q[1];
rz(-1.2661625) q[1];
sx q[1];
rz(1.6564293) q[1];
rz(-pi) q[2];
x q[2];
rz(2.769737) q[3];
sx q[3];
rz(-2.1346492) q[3];
sx q[3];
rz(-1.9052037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(3.0653811) q[2];
rz(-1.8501836) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(-0.61856234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718219) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(-2.3845657) q[0];
rz(-2.3454759) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-2.1597791) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7863203) q[0];
sx q[0];
rz(-0.85714802) q[0];
sx q[0];
rz(2.7631604) q[0];
x q[1];
rz(-3.0715838) q[2];
sx q[2];
rz(-1.677779) q[2];
sx q[2];
rz(-1.9221523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.031739) q[1];
sx q[1];
rz(-0.64424911) q[1];
sx q[1];
rz(2.2883313) q[1];
rz(-pi) q[2];
rz(-2.997274) q[3];
sx q[3];
rz(-1.2929521) q[3];
sx q[3];
rz(-1.2793737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41165274) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(-2.5977503) q[2];
rz(-0.020921556) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(0.81834832) q[3];
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
rz(2.2127317) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(0.62335706) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(-0.26434937) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0118222) q[0];
sx q[0];
rz(-1.3589824) q[0];
sx q[0];
rz(2.2925633) q[0];
rz(-1.9064205) q[2];
sx q[2];
rz(-1.5315637) q[2];
sx q[2];
rz(2.7754663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5874065) q[1];
sx q[1];
rz(-2.3465112) q[1];
sx q[1];
rz(-1.1752179) q[1];
x q[2];
rz(1.8885625) q[3];
sx q[3];
rz(-1.6971641) q[3];
sx q[3];
rz(1.4082091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.136772) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(0.47478673) q[2];
rz(-2.7759806) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356165) q[0];
sx q[0];
rz(-2.7662179) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(-1.0659418) q[1];
sx q[1];
rz(-1.424788) q[1];
sx q[1];
rz(2.8071383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1565103) q[0];
sx q[0];
rz(-1.2691536) q[0];
sx q[0];
rz(1.7443875) q[0];
rz(-0.29269258) q[2];
sx q[2];
rz(-2.6862157) q[2];
sx q[2];
rz(-2.6259881) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6600646) q[1];
sx q[1];
rz(-1.1236407) q[1];
sx q[1];
rz(1.4051953) q[1];
x q[2];
rz(-2.6869557) q[3];
sx q[3];
rz(-0.84970039) q[3];
sx q[3];
rz(-2.7434512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36755422) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(0.96662194) q[2];
rz(-1.4878081) q[3];
sx q[3];
rz(-2.8130468) q[3];
sx q[3];
rz(-3.0706792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3286572) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(-2.8712414) q[0];
rz(1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(2.5152452) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0369025) q[0];
sx q[0];
rz(-2.2459163) q[0];
sx q[0];
rz(-0.40963197) q[0];
rz(-pi) q[1];
rz(0.93339351) q[2];
sx q[2];
rz(-1.9248171) q[2];
sx q[2];
rz(2.6824981) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8737265) q[1];
sx q[1];
rz(-1.6609816) q[1];
sx q[1];
rz(1.5298046) q[1];
rz(-pi) q[2];
rz(-2.7553431) q[3];
sx q[3];
rz(-1.874141) q[3];
sx q[3];
rz(-1.4331499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7865929) q[2];
sx q[2];
rz(-0.91675106) q[2];
sx q[2];
rz(-1.5691441) q[2];
rz(-2.3908424) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(2.2641838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5759721) q[0];
sx q[0];
rz(-1.8436057) q[0];
sx q[0];
rz(2.5634503) q[0];
rz(1.3427973) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-3.1411896) q[2];
sx q[2];
rz(-0.12231356) q[2];
sx q[2];
rz(-0.077153645) q[2];
rz(-0.20082898) q[3];
sx q[3];
rz(-0.26293892) q[3];
sx q[3];
rz(0.50504167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
