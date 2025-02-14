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
rz(0.084963381) q[0];
sx q[0];
rz(3.4439937) q[0];
sx q[0];
rz(6.2511282) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(-0.75262466) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948715) q[0];
sx q[0];
rz(-1.7716452) q[0];
sx q[0];
rz(-0.3263536) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9448124) q[2];
sx q[2];
rz(-1.2311282) q[2];
sx q[2];
rz(-2.938478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.17129414) q[1];
sx q[1];
rz(-2.7134905) q[1];
sx q[1];
rz(0.64117214) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0438552) q[3];
sx q[3];
rz(-0.40832061) q[3];
sx q[3];
rz(-2.4957617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1746615) q[2];
sx q[2];
rz(-1.175468) q[2];
sx q[2];
rz(-0.82007972) q[2];
rz(2.7005633) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(2.5681514) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(2.7313857) q[0];
rz(2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.358323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5429496) q[0];
sx q[0];
rz(-0.99587593) q[0];
sx q[0];
rz(-0.47358124) q[0];
rz(-2.63851) q[2];
sx q[2];
rz(-1.0840992) q[2];
sx q[2];
rz(-2.8622735) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4315003) q[1];
sx q[1];
rz(-2.9572801) q[1];
sx q[1];
rz(-1.1460365) q[1];
rz(0.15492691) q[3];
sx q[3];
rz(-2.2048031) q[3];
sx q[3];
rz(0.67544322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.35839781) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(0.50296339) q[2];
rz(-1.8424235) q[3];
sx q[3];
rz(-0.75853577) q[3];
sx q[3];
rz(-2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7314887) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(-2.1657535) q[0];
rz(0.93636912) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(3.0036614) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5463211) q[0];
sx q[0];
rz(-2.0450651) q[0];
sx q[0];
rz(-0.41149891) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8713869) q[2];
sx q[2];
rz(-2.7305354) q[2];
sx q[2];
rz(-0.11693987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5441078) q[1];
sx q[1];
rz(-0.710809) q[1];
sx q[1];
rz(-1.3686468) q[1];
x q[2];
rz(0.46584399) q[3];
sx q[3];
rz(-1.8154105) q[3];
sx q[3];
rz(0.96058339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58933538) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(1.662558) q[2];
rz(-2.3432483) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(0.020309694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.289157) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(2.074746) q[0];
rz(-1.5929219) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(-2.7222395) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8035078) q[0];
sx q[0];
rz(-2.8610021) q[0];
sx q[0];
rz(0.96208939) q[0];
x q[1];
rz(2.9871171) q[2];
sx q[2];
rz(-0.53336582) q[2];
sx q[2];
rz(1.6090924) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49653445) q[1];
sx q[1];
rz(-2.7784111) q[1];
sx q[1];
rz(1.3892322) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97030117) q[3];
sx q[3];
rz(-0.56823778) q[3];
sx q[3];
rz(2.0214391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7833917) q[2];
sx q[2];
rz(-0.50938598) q[2];
sx q[2];
rz(-1.5979213) q[2];
rz(-1.2265497) q[3];
sx q[3];
rz(-1.2164601) q[3];
sx q[3];
rz(1.290087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.876038) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(2.3853886) q[0];
rz(1.8887695) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(1.7255712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21391103) q[0];
sx q[0];
rz(-0.98386231) q[0];
sx q[0];
rz(-1.3089433) q[0];
rz(-0.33966392) q[2];
sx q[2];
rz(-2.77861) q[2];
sx q[2];
rz(-1.5756922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.071347) q[1];
sx q[1];
rz(-1.5778549) q[1];
sx q[1];
rz(0.53877212) q[1];
rz(-2.6719499) q[3];
sx q[3];
rz(-2.7386463) q[3];
sx q[3];
rz(2.6279891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8004134) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(-0.18079147) q[2];
rz(-1.4888658) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7086696) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(-2.3413626) q[0];
rz(2.8569787) q[1];
sx q[1];
rz(-2.0814643) q[1];
sx q[1];
rz(-0.12399331) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127045) q[0];
sx q[0];
rz(-1.4409954) q[0];
sx q[0];
rz(3.1129254) q[0];
rz(-1.8154549) q[2];
sx q[2];
rz(-2.8996433) q[2];
sx q[2];
rz(1.9242147) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2468894) q[1];
sx q[1];
rz(-1.2661625) q[1];
sx q[1];
rz(-1.4851634) q[1];
rz(-pi) q[2];
rz(0.97449755) q[3];
sx q[3];
rz(-1.8829405) q[3];
sx q[3];
rz(-2.6016935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2064712) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(0.07621152) q[2];
rz(-1.291409) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(-0.61856234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718219) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(-0.75702697) q[0];
rz(2.3454759) q[1];
sx q[1];
rz(-1.7346104) q[1];
sx q[1];
rz(0.98181358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7863203) q[0];
sx q[0];
rz(-0.85714802) q[0];
sx q[0];
rz(-2.7631604) q[0];
x q[1];
rz(1.4635529) q[2];
sx q[2];
rz(-1.6404043) q[2];
sx q[2];
rz(-0.34386841) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.031739) q[1];
sx q[1];
rz(-2.4973435) q[1];
sx q[1];
rz(-2.2883313) q[1];
x q[2];
rz(1.8514093) q[3];
sx q[3];
rz(-1.4320489) q[3];
sx q[3];
rz(-0.25158238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7299399) q[2];
sx q[2];
rz(-0.81092683) q[2];
sx q[2];
rz(2.5977503) q[2];
rz(-3.1206711) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(-0.81834832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92886096) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(0.62335706) q[0];
rz(0.32132545) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(0.26434937) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0118222) q[0];
sx q[0];
rz(-1.7826102) q[0];
sx q[0];
rz(2.2925633) q[0];
rz(-pi) q[1];
rz(1.4521763) q[2];
sx q[2];
rz(-0.33782321) q[2];
sx q[2];
rz(-2.0488957) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26786727) q[1];
sx q[1];
rz(-1.2920989) q[1];
sx q[1];
rz(2.3257117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13295138) q[3];
sx q[3];
rz(-1.2556517) q[3];
sx q[3];
rz(0.12115762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.136772) q[2];
sx q[2];
rz(-0.68237582) q[2];
sx q[2];
rz(-2.6668059) q[2];
rz(-0.36561203) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30597618) q[0];
sx q[0];
rz(-2.7662179) q[0];
sx q[0];
rz(2.6557652) q[0];
rz(-2.0756508) q[1];
sx q[1];
rz(-1.424788) q[1];
sx q[1];
rz(0.33445439) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850824) q[0];
sx q[0];
rz(-1.2691536) q[0];
sx q[0];
rz(1.3972052) q[0];
rz(-pi) q[1];
rz(-2.7031519) q[2];
sx q[2];
rz(-1.4435569) q[2];
sx q[2];
rz(-1.3194989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.124467) q[1];
sx q[1];
rz(-1.7199868) q[1];
sx q[1];
rz(0.45254947) q[1];
rz(-2.0341088) q[3];
sx q[3];
rz(-0.83011887) q[3];
sx q[3];
rz(2.9024359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7740384) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(0.96662194) q[2];
rz(1.6537846) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(-0.070913471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8129355) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(-2.8712414) q[0];
rz(1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26890818) q[0];
sx q[0];
rz(-1.2546526) q[0];
sx q[0];
rz(0.85319467) q[0];
x q[1];
rz(0.93339351) q[2];
sx q[2];
rz(-1.9248171) q[2];
sx q[2];
rz(-0.4590946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2992363) q[1];
sx q[1];
rz(-1.5299712) q[1];
sx q[1];
rz(-3.051332) q[1];
rz(-0.38624951) q[3];
sx q[3];
rz(-1.2674517) q[3];
sx q[3];
rz(-1.4331499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7865929) q[2];
sx q[2];
rz(-0.91675106) q[2];
sx q[2];
rz(-1.5691441) q[2];
rz(-2.3908424) q[3];
sx q[3];
rz(-2.7995977) q[3];
sx q[3];
rz(-2.2641838) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56562051) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(-1.7987953) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-1.5708459) q[2];
sx q[2];
rz(-1.4484828) q[2];
sx q[2];
rz(3.0640329) q[2];
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
