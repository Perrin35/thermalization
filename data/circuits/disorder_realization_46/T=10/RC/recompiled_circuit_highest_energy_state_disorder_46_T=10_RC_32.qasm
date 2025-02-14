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
rz(2.1815648) q[0];
sx q[0];
rz(-0.61909827) q[0];
sx q[0];
rz(1.4033138) q[0];
rz(-0.69544855) q[1];
sx q[1];
rz(-2.5814576) q[1];
sx q[1];
rz(-0.6518031) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2346828) q[0];
sx q[0];
rz(-1.9056589) q[0];
sx q[0];
rz(-0.85483179) q[0];
rz(2.8601134) q[2];
sx q[2];
rz(-1.5764766) q[2];
sx q[2];
rz(1.8112294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.044120248) q[1];
sx q[1];
rz(-1.5481045) q[1];
sx q[1];
rz(0.23182643) q[1];
rz(-0.27068287) q[3];
sx q[3];
rz(-0.83867618) q[3];
sx q[3];
rz(1.0974778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92242509) q[2];
sx q[2];
rz(-1.3895915) q[2];
sx q[2];
rz(0.31828848) q[2];
rz(-0.90218941) q[3];
sx q[3];
rz(-0.92034322) q[3];
sx q[3];
rz(0.86377803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1282463) q[0];
sx q[0];
rz(-2.8885169) q[0];
sx q[0];
rz(-2.6213562) q[0];
rz(-3.0546313) q[1];
sx q[1];
rz(-2.088701) q[1];
sx q[1];
rz(2.3435977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32736054) q[0];
sx q[0];
rz(-1.260574) q[0];
sx q[0];
rz(-1.7807175) q[0];
x q[1];
rz(-1.410631) q[2];
sx q[2];
rz(-2.3078683) q[2];
sx q[2];
rz(-3.0708416) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4072306) q[1];
sx q[1];
rz(-1.4301908) q[1];
sx q[1];
rz(-2.749252) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4142031) q[3];
sx q[3];
rz(-0.57021457) q[3];
sx q[3];
rz(-0.19661759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.308455) q[2];
sx q[2];
rz(-1.1546346) q[2];
sx q[2];
rz(0.40033611) q[2];
rz(-0.18276246) q[3];
sx q[3];
rz(-0.57923135) q[3];
sx q[3];
rz(1.8496752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12277814) q[0];
sx q[0];
rz(-3.0570539) q[0];
sx q[0];
rz(1.2070745) q[0];
rz(1.6385551) q[1];
sx q[1];
rz(-1.8820347) q[1];
sx q[1];
rz(0.17301339) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6956222) q[0];
sx q[0];
rz(-0.39485952) q[0];
sx q[0];
rz(-1.7941712) q[0];
x q[1];
rz(3.0315057) q[2];
sx q[2];
rz(-2.1497257) q[2];
sx q[2];
rz(0.26012173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8743778) q[1];
sx q[1];
rz(-0.80433955) q[1];
sx q[1];
rz(2.4206603) q[1];
x q[2];
rz(1.7812114) q[3];
sx q[3];
rz(-2.3074779) q[3];
sx q[3];
rz(2.2115416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5485237) q[2];
sx q[2];
rz(-2.5544281) q[2];
sx q[2];
rz(2.1068088) q[2];
rz(2.2506524) q[3];
sx q[3];
rz(-0.8689298) q[3];
sx q[3];
rz(-1.9152036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064875038) q[0];
sx q[0];
rz(-2.1893976) q[0];
sx q[0];
rz(-0.12884831) q[0];
rz(0.42875641) q[1];
sx q[1];
rz(-1.4219475) q[1];
sx q[1];
rz(-1.4770003) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77177202) q[0];
sx q[0];
rz(-1.6390529) q[0];
sx q[0];
rz(-2.7959156) q[0];
rz(-2.3643983) q[2];
sx q[2];
rz(-2.5487196) q[2];
sx q[2];
rz(2.2490918) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13969888) q[1];
sx q[1];
rz(-0.47283462) q[1];
sx q[1];
rz(-1.0919149) q[1];
rz(-1.5691532) q[3];
sx q[3];
rz(-2.5503359) q[3];
sx q[3];
rz(-1.7776684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9513272) q[2];
sx q[2];
rz(-1.724388) q[2];
sx q[2];
rz(-1.8180656) q[2];
rz(-1.6773112) q[3];
sx q[3];
rz(-2.5530294) q[3];
sx q[3];
rz(2.467449) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4071963) q[0];
sx q[0];
rz(-0.36235991) q[0];
sx q[0];
rz(1.4417484) q[0];
rz(-2.7184519) q[1];
sx q[1];
rz(-0.95760456) q[1];
sx q[1];
rz(2.8513681) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8347802) q[0];
sx q[0];
rz(-2.4270822) q[0];
sx q[0];
rz(-0.45551703) q[0];
rz(2.1802203) q[2];
sx q[2];
rz(-2.0573685) q[2];
sx q[2];
rz(-2.3455462) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4729974) q[1];
sx q[1];
rz(-1.6194222) q[1];
sx q[1];
rz(-3.0134788) q[1];
rz(-1.832015) q[3];
sx q[3];
rz(-1.272837) q[3];
sx q[3];
rz(-3.1233149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71771249) q[2];
sx q[2];
rz(-2.3545357) q[2];
sx q[2];
rz(0.94669739) q[2];
rz(-0.26149073) q[3];
sx q[3];
rz(-1.3638834) q[3];
sx q[3];
rz(0.3869032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.43712) q[0];
sx q[0];
rz(-1.8081212) q[0];
sx q[0];
rz(-2.3719846) q[0];
rz(-0.17598027) q[1];
sx q[1];
rz(-1.3409921) q[1];
sx q[1];
rz(-1.5716858) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0971477) q[0];
sx q[0];
rz(-2.189069) q[0];
sx q[0];
rz(-2.3522105) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85167517) q[2];
sx q[2];
rz(-1.5181418) q[2];
sx q[2];
rz(-0.55916121) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9573201) q[1];
sx q[1];
rz(-1.9369003) q[1];
sx q[1];
rz(-1.373532) q[1];
x q[2];
rz(-0.50521781) q[3];
sx q[3];
rz(-1.1105624) q[3];
sx q[3];
rz(1.203618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1311329) q[2];
sx q[2];
rz(-0.95644462) q[2];
sx q[2];
rz(-1.825911) q[2];
rz(-1.2223505) q[3];
sx q[3];
rz(-2.0171916) q[3];
sx q[3];
rz(-0.30092064) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2216126) q[0];
sx q[0];
rz(-1.6918809) q[0];
sx q[0];
rz(1.8020887) q[0];
rz(1.9446531) q[1];
sx q[1];
rz(-2.0868128) q[1];
sx q[1];
rz(0.96492499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018835807) q[0];
sx q[0];
rz(-1.0623451) q[0];
sx q[0];
rz(-1.4623789) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2193421) q[2];
sx q[2];
rz(-1.2932475) q[2];
sx q[2];
rz(0.48257839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33789148) q[1];
sx q[1];
rz(-1.6543904) q[1];
sx q[1];
rz(0.22050942) q[1];
rz(-0.32437848) q[3];
sx q[3];
rz(-2.6357963) q[3];
sx q[3];
rz(-1.5816816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7397466) q[2];
sx q[2];
rz(-1.2082938) q[2];
sx q[2];
rz(1.9604663) q[2];
rz(3.127023) q[3];
sx q[3];
rz(-2.7404692) q[3];
sx q[3];
rz(1.9692263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75083098) q[0];
sx q[0];
rz(-1.9347235) q[0];
sx q[0];
rz(2.3857351) q[0];
rz(-2.4769056) q[1];
sx q[1];
rz(-1.942778) q[1];
sx q[1];
rz(-1.8519148) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8553365) q[0];
sx q[0];
rz(-0.77833114) q[0];
sx q[0];
rz(-1.6617213) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6014156) q[2];
sx q[2];
rz(-1.3103169) q[2];
sx q[2];
rz(1.8249701) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7482704) q[1];
sx q[1];
rz(-0.38667187) q[1];
sx q[1];
rz(2.3063117) q[1];
rz(-pi) q[2];
rz(-0.30925242) q[3];
sx q[3];
rz(-1.6669295) q[3];
sx q[3];
rz(-2.1949196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1332625) q[2];
sx q[2];
rz(-1.3349168) q[2];
sx q[2];
rz(1.5745715) q[2];
rz(-1.3957006) q[3];
sx q[3];
rz(-1.7392266) q[3];
sx q[3];
rz(0.76656669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181327) q[0];
sx q[0];
rz(-0.27579951) q[0];
sx q[0];
rz(2.7511399) q[0];
rz(-2.0402724) q[1];
sx q[1];
rz(-1.7951671) q[1];
sx q[1];
rz(0.93213814) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7926246) q[0];
sx q[0];
rz(-1.5431964) q[0];
sx q[0];
rz(-1.6076002) q[0];
rz(3.0603067) q[2];
sx q[2];
rz(-0.54078201) q[2];
sx q[2];
rz(3.0791435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.336672) q[1];
sx q[1];
rz(-1.3778566) q[1];
sx q[1];
rz(2.8590917) q[1];
rz(-pi) q[2];
x q[2];
rz(1.325439) q[3];
sx q[3];
rz(-1.1096769) q[3];
sx q[3];
rz(-2.2803015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0497047) q[2];
sx q[2];
rz(-2.4447417) q[2];
sx q[2];
rz(1.4616802) q[2];
rz(1.0852496) q[3];
sx q[3];
rz(-1.1321675) q[3];
sx q[3];
rz(0.6404883) q[3];
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
rz(-2.853448) q[0];
sx q[0];
rz(-2.4818821) q[0];
sx q[0];
rz(1.1261384) q[0];
rz(2.0938734) q[1];
sx q[1];
rz(-0.62647096) q[1];
sx q[1];
rz(-1.610021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1739222) q[0];
sx q[0];
rz(-1.5940658) q[0];
sx q[0];
rz(-1.5694005) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.80349) q[2];
sx q[2];
rz(-2.0183211) q[2];
sx q[2];
rz(1.1648503) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6781552) q[1];
sx q[1];
rz(-0.10129866) q[1];
sx q[1];
rz(-0.8763635) q[1];
rz(-pi) q[2];
rz(-3.0017742) q[3];
sx q[3];
rz(-0.56666683) q[3];
sx q[3];
rz(1.293928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.62022007) q[2];
sx q[2];
rz(-2.3282101) q[2];
sx q[2];
rz(2.9353471) q[2];
rz(-0.5591875) q[3];
sx q[3];
rz(-1.6181479) q[3];
sx q[3];
rz(2.0172113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54022057) q[0];
sx q[0];
rz(-1.843353) q[0];
sx q[0];
rz(-2.699615) q[0];
rz(1.1454918) q[1];
sx q[1];
rz(-1.3064697) q[1];
sx q[1];
rz(-1.3189955) q[1];
rz(-2.3445208) q[2];
sx q[2];
rz(-2.7596724) q[2];
sx q[2];
rz(1.8308664) q[2];
rz(2.4174684) q[3];
sx q[3];
rz(-1.4403314) q[3];
sx q[3];
rz(-2.452313) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
