OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1448016) q[0];
sx q[0];
rz(0.15455833) q[0];
sx q[0];
rz(6.9757087) q[0];
rz(1.9321631) q[1];
sx q[1];
rz(-1.2485319) q[1];
sx q[1];
rz(-1.385153) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6801075) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(2.0018342) q[0];
rz(-pi) q[1];
rz(-0.42983774) q[2];
sx q[2];
rz(-0.59525437) q[2];
sx q[2];
rz(1.0560448) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3424211) q[1];
sx q[1];
rz(-0.28563269) q[1];
sx q[1];
rz(2.530453) q[1];
rz(-pi) q[2];
rz(-0.65269835) q[3];
sx q[3];
rz(-1.9737118) q[3];
sx q[3];
rz(0.17890113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8866855) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(-0.77130476) q[3];
sx q[3];
rz(-0.78273928) q[3];
sx q[3];
rz(-2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7339864) q[0];
sx q[0];
rz(-2.3953231) q[0];
sx q[0];
rz(0.45390391) q[0];
rz(2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(1.9143547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95943806) q[0];
sx q[0];
rz(-1.6998788) q[0];
sx q[0];
rz(3.0653619) q[0];
rz(-pi) q[1];
rz(-2.0298376) q[2];
sx q[2];
rz(-1.0978205) q[2];
sx q[2];
rz(0.74726653) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7397346) q[1];
sx q[1];
rz(-1.8622073) q[1];
sx q[1];
rz(-1.3009562) q[1];
x q[2];
rz(0.20083986) q[3];
sx q[3];
rz(-1.4174995) q[3];
sx q[3];
rz(0.64985819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1002905) q[2];
sx q[2];
rz(-1.1854478) q[2];
sx q[2];
rz(-0.56742898) q[2];
rz(2.7764017) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(-0.96810961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(-0.89865249) q[0];
rz(0.99575106) q[1];
sx q[1];
rz(-1.5581222) q[1];
sx q[1];
rz(-2.8083037) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3086739) q[0];
sx q[0];
rz(-1.1090288) q[0];
sx q[0];
rz(-0.69899107) q[0];
rz(0.097924175) q[2];
sx q[2];
rz(-1.3946748) q[2];
sx q[2];
rz(-0.91913659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4132727) q[1];
sx q[1];
rz(-1.7227168) q[1];
sx q[1];
rz(-2.1876213) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0201449) q[3];
sx q[3];
rz(-1.2216976) q[3];
sx q[3];
rz(0.89494866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4553392) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(1.1509482) q[2];
rz(-0.84093705) q[3];
sx q[3];
rz(-2.0261814) q[3];
sx q[3];
rz(-1.1545198) q[3];
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
rz(-0.44160098) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(2.4568795) q[0];
rz(2.1060064) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(1.9365786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3048153) q[0];
sx q[0];
rz(-0.69201058) q[0];
sx q[0];
rz(-1.6230323) q[0];
x q[1];
rz(1.3022468) q[2];
sx q[2];
rz(-2.4797202) q[2];
sx q[2];
rz(0.72999398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11338621) q[1];
sx q[1];
rz(-0.75062597) q[1];
sx q[1];
rz(1.1651462) q[1];
rz(-2.2756696) q[3];
sx q[3];
rz(-1.5161627) q[3];
sx q[3];
rz(1.6161402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24923199) q[2];
sx q[2];
rz(-1.7148596) q[2];
sx q[2];
rz(-0.37115804) q[2];
rz(1.7403729) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(-1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.086833) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(-3.0084685) q[0];
rz(0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(2.5865119) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18300444) q[0];
sx q[0];
rz(-2.8259813) q[0];
sx q[0];
rz(-1.6118227) q[0];
x q[1];
rz(-0.06818469) q[2];
sx q[2];
rz(-1.9848595) q[2];
sx q[2];
rz(-2.806864) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75025573) q[1];
sx q[1];
rz(-1.5899982) q[1];
sx q[1];
rz(-2.0160497) q[1];
rz(-pi) q[2];
rz(2.5693232) q[3];
sx q[3];
rz(-0.096147691) q[3];
sx q[3];
rz(0.26086807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30620265) q[2];
sx q[2];
rz(-1.0058879) q[2];
sx q[2];
rz(-0.13892697) q[2];
rz(2.1991918) q[3];
sx q[3];
rz(-1.5706294) q[3];
sx q[3];
rz(-2.5901103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54365629) q[0];
sx q[0];
rz(-0.56977001) q[0];
sx q[0];
rz(-2.561835) q[0];
rz(3.014091) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(-1.5396083) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67028763) q[0];
sx q[0];
rz(-1.7738713) q[0];
sx q[0];
rz(-2.2905486) q[0];
rz(-pi) q[1];
rz(-0.13055735) q[2];
sx q[2];
rz(-1.5255543) q[2];
sx q[2];
rz(-3.0322078) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0582038) q[1];
sx q[1];
rz(-1.997588) q[1];
sx q[1];
rz(0.31174387) q[1];
x q[2];
rz(0.29069889) q[3];
sx q[3];
rz(-0.906956) q[3];
sx q[3];
rz(1.3568527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55398983) q[2];
sx q[2];
rz(-0.25045276) q[2];
sx q[2];
rz(2.8721151) q[2];
rz(0.23412165) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6948029) q[0];
sx q[0];
rz(-0.61820522) q[0];
sx q[0];
rz(-0.68429464) q[0];
rz(0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(-2.6228242) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2001901) q[0];
sx q[0];
rz(-1.4369643) q[0];
sx q[0];
rz(2.3076513) q[0];
x q[1];
rz(0.77063009) q[2];
sx q[2];
rz(-1.4118886) q[2];
sx q[2];
rz(0.7030013) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5468532) q[1];
sx q[1];
rz(-0.93187983) q[1];
sx q[1];
rz(-0.98313318) q[1];
rz(1.5090844) q[3];
sx q[3];
rz(-2.7402174) q[3];
sx q[3];
rz(-0.20460953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2542904) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(-0.051579483) q[3];
sx q[3];
rz(-1.1392461) q[3];
sx q[3];
rz(2.1896867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1812487) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(1.3990336) q[0];
rz(-0.78701204) q[1];
sx q[1];
rz(-2.0128638) q[1];
sx q[1];
rz(0.74434892) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23119152) q[0];
sx q[0];
rz(-1.7185128) q[0];
sx q[0];
rz(1.5157248) q[0];
rz(-2.9888399) q[2];
sx q[2];
rz(-1.061201) q[2];
sx q[2];
rz(-2.5324412) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0155322) q[1];
sx q[1];
rz(-1.0711526) q[1];
sx q[1];
rz(-1.8263032) q[1];
rz(-3.0934422) q[3];
sx q[3];
rz(-2.6865494) q[3];
sx q[3];
rz(0.20850785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4259592) q[2];
sx q[2];
rz(-1.3954433) q[2];
sx q[2];
rz(0.60950935) q[2];
rz(0.65731796) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(-2.8801584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(-1.5040065) q[0];
rz(1.2414744) q[1];
sx q[1];
rz(-1.1499317) q[1];
sx q[1];
rz(0.77493587) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0938213) q[0];
sx q[0];
rz(-1.1822961) q[0];
sx q[0];
rz(-0.85431487) q[0];
x q[1];
rz(2.9853285) q[2];
sx q[2];
rz(-2.0093577) q[2];
sx q[2];
rz(1.6608479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96156582) q[1];
sx q[1];
rz(-2.0225836) q[1];
sx q[1];
rz(1.7801442) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9990066) q[3];
sx q[3];
rz(-1.7877868) q[3];
sx q[3];
rz(-2.2246974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.187414) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(1.212451) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(1.2780179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-1.071844) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(-0.21433314) q[0];
rz(-0.65746039) q[1];
sx q[1];
rz(-0.22413707) q[1];
sx q[1];
rz(1.0459895) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25046529) q[0];
sx q[0];
rz(-1.2123101) q[0];
sx q[0];
rz(-1.4573218) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5997535) q[2];
sx q[2];
rz(-2.1615897) q[2];
sx q[2];
rz(-2.7383885) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6194832) q[1];
sx q[1];
rz(-1.932718) q[1];
sx q[1];
rz(-0.98720179) q[1];
rz(-pi) q[2];
rz(-0.40542116) q[3];
sx q[3];
rz(-1.5463721) q[3];
sx q[3];
rz(-3.0693808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7252698) q[2];
sx q[2];
rz(-0.060083397) q[2];
sx q[2];
rz(-2.0521169) q[2];
rz(1.5661092) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(-0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2789223) q[0];
sx q[0];
rz(-2.537732) q[0];
sx q[0];
rz(-2.296007) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(-2.5074742) q[2];
sx q[2];
rz(-0.49370439) q[2];
sx q[2];
rz(1.6903071) q[2];
rz(-2.5915626) q[3];
sx q[3];
rz(-1.5169654) q[3];
sx q[3];
rz(-2.9548002) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];