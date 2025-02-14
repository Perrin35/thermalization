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
rz(1.3697019) q[0];
sx q[0];
rz(-0.43039027) q[0];
sx q[0];
rz(-0.16464591) q[0];
rz(1.4511664) q[1];
sx q[1];
rz(3.5199447) q[1];
sx q[1];
rz(8.7149109) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9116316) q[0];
sx q[0];
rz(-2.8459186) q[0];
sx q[0];
rz(2.0333181) q[0];
rz(-pi) q[1];
rz(-2.4455382) q[2];
sx q[2];
rz(-0.6185607) q[2];
sx q[2];
rz(-1.5355009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0069847) q[1];
sx q[1];
rz(-0.35672327) q[1];
sx q[1];
rz(-1.7118042) q[1];
x q[2];
rz(-0.48810256) q[3];
sx q[3];
rz(-2.8781366) q[3];
sx q[3];
rz(1.7149705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6115438) q[2];
sx q[2];
rz(-2.7908466) q[2];
sx q[2];
rz(1.2190399) q[2];
rz(-0.48398316) q[3];
sx q[3];
rz(-2.2958906) q[3];
sx q[3];
rz(1.774196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4265863) q[0];
sx q[0];
rz(-0.33892092) q[0];
sx q[0];
rz(1.6981079) q[0];
rz(1.8679856) q[1];
sx q[1];
rz(-2.1111635) q[1];
sx q[1];
rz(2.4748306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979157) q[0];
sx q[0];
rz(-0.34404342) q[0];
sx q[0];
rz(2.6187569) q[0];
rz(-0.86033852) q[2];
sx q[2];
rz(-0.14709148) q[2];
sx q[2];
rz(1.8837613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.264111) q[1];
sx q[1];
rz(-2.5248233) q[1];
sx q[1];
rz(-2.3186604) q[1];
rz(-pi) q[2];
rz(0.81979378) q[3];
sx q[3];
rz(-1.3043376) q[3];
sx q[3];
rz(2.418892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4296809) q[2];
sx q[2];
rz(-1.9671755) q[2];
sx q[2];
rz(2.5453117) q[2];
rz(1.0239673) q[3];
sx q[3];
rz(-1.3919316) q[3];
sx q[3];
rz(-2.021324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100818) q[0];
sx q[0];
rz(-2.5937268) q[0];
sx q[0];
rz(-1.4672853) q[0];
rz(3.1027377) q[1];
sx q[1];
rz(-1.424574) q[1];
sx q[1];
rz(2.255596) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65424991) q[0];
sx q[0];
rz(-2.4375101) q[0];
sx q[0];
rz(-2.447633) q[0];
rz(-pi) q[1];
rz(-1.5019234) q[2];
sx q[2];
rz(-0.67126545) q[2];
sx q[2];
rz(-1.6948989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3275324) q[1];
sx q[1];
rz(-2.7814354) q[1];
sx q[1];
rz(0.11231695) q[1];
x q[2];
rz(-3.0795512) q[3];
sx q[3];
rz(-1.4915371) q[3];
sx q[3];
rz(-2.6798927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92750612) q[2];
sx q[2];
rz(-1.6897886) q[2];
sx q[2];
rz(-0.56373325) q[2];
rz(0.8756513) q[3];
sx q[3];
rz(-2.0893658) q[3];
sx q[3];
rz(1.0912857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5617274) q[0];
sx q[0];
rz(-2.2932597) q[0];
sx q[0];
rz(-1.0680098) q[0];
rz(0.39455286) q[1];
sx q[1];
rz(-0.5296455) q[1];
sx q[1];
rz(-1.669917) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36763299) q[0];
sx q[0];
rz(-1.9955705) q[0];
sx q[0];
rz(-1.1921117) q[0];
x q[1];
rz(1.9392233) q[2];
sx q[2];
rz(-2.8369008) q[2];
sx q[2];
rz(-0.42542377) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.010943451) q[1];
sx q[1];
rz(-0.92045451) q[1];
sx q[1];
rz(0.54231142) q[1];
rz(-0.39880347) q[3];
sx q[3];
rz(-1.670518) q[3];
sx q[3];
rz(2.698363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7431405) q[2];
sx q[2];
rz(-1.3116216) q[2];
sx q[2];
rz(0.12942806) q[2];
rz(-0.5433003) q[3];
sx q[3];
rz(-1.0933417) q[3];
sx q[3];
rz(1.7530542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1035136) q[0];
sx q[0];
rz(-1.0032126) q[0];
sx q[0];
rz(2.6858618) q[0];
rz(2.575846) q[1];
sx q[1];
rz(-2.5028298) q[1];
sx q[1];
rz(-0.21557132) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8889897) q[0];
sx q[0];
rz(-0.71160331) q[0];
sx q[0];
rz(1.8322741) q[0];
rz(-pi) q[1];
rz(2.2224765) q[2];
sx q[2];
rz(-1.0813776) q[2];
sx q[2];
rz(0.31351837) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4734772) q[1];
sx q[1];
rz(-1.6398506) q[1];
sx q[1];
rz(1.3051239) q[1];
x q[2];
rz(3.0318063) q[3];
sx q[3];
rz(-1.2874787) q[3];
sx q[3];
rz(2.7840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4244708) q[2];
sx q[2];
rz(-0.30693808) q[2];
sx q[2];
rz(1.4324987) q[2];
rz(-1.1411544) q[3];
sx q[3];
rz(-1.9211831) q[3];
sx q[3];
rz(0.46877638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6345217) q[0];
sx q[0];
rz(-0.85352007) q[0];
sx q[0];
rz(0.61990196) q[0];
rz(1.2076521) q[1];
sx q[1];
rz(-2.4687605) q[1];
sx q[1];
rz(-2.8501453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0196592) q[0];
sx q[0];
rz(-2.1284416) q[0];
sx q[0];
rz(-2.9479821) q[0];
rz(-0.3235422) q[2];
sx q[2];
rz(-1.6969883) q[2];
sx q[2];
rz(-3.0033811) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48780729) q[1];
sx q[1];
rz(-1.331662) q[1];
sx q[1];
rz(0.59822786) q[1];
x q[2];
rz(3.1278947) q[3];
sx q[3];
rz(-0.65675101) q[3];
sx q[3];
rz(-0.17431549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.42346272) q[2];
sx q[2];
rz(-0.66897696) q[2];
sx q[2];
rz(-0.97406975) q[2];
rz(-1.918321) q[3];
sx q[3];
rz(-1.7134824) q[3];
sx q[3];
rz(2.1035002) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85422) q[0];
sx q[0];
rz(-1.9725476) q[0];
sx q[0];
rz(0.30782345) q[0];
rz(-2.8616915) q[1];
sx q[1];
rz(-0.24903909) q[1];
sx q[1];
rz(-1.261796) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7578106) q[0];
sx q[0];
rz(-1.5368665) q[0];
sx q[0];
rz(-1.7050123) q[0];
rz(-1.1818188) q[2];
sx q[2];
rz(-0.87724596) q[2];
sx q[2];
rz(0.77893585) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9692291) q[1];
sx q[1];
rz(-1.5820272) q[1];
sx q[1];
rz(-1.8658691) q[1];
rz(-pi) q[2];
rz(2.0359751) q[3];
sx q[3];
rz(-1.6334459) q[3];
sx q[3];
rz(-2.986547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5032924) q[2];
sx q[2];
rz(-1.952848) q[2];
sx q[2];
rz(-2.3390181) q[2];
rz(1.5447626) q[3];
sx q[3];
rz(-2.7464726) q[3];
sx q[3];
rz(-1.6667268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.06642) q[0];
sx q[0];
rz(-1.5936699) q[0];
sx q[0];
rz(-2.5286034) q[0];
rz(-2.0892443) q[1];
sx q[1];
rz(-2.3520825) q[1];
sx q[1];
rz(1.2552415) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6832812) q[0];
sx q[0];
rz(-1.1363239) q[0];
sx q[0];
rz(2.6660454) q[0];
rz(-pi) q[1];
rz(-1.2118039) q[2];
sx q[2];
rz(-2.5796522) q[2];
sx q[2];
rz(-0.46899624) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.44574983) q[1];
sx q[1];
rz(-0.90409213) q[1];
sx q[1];
rz(-0.95335828) q[1];
rz(-pi) q[2];
rz(-1.4358712) q[3];
sx q[3];
rz(-1.3107015) q[3];
sx q[3];
rz(-1.2430199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5455948) q[2];
sx q[2];
rz(-2.2277446) q[2];
sx q[2];
rz(2.9929898) q[2];
rz(0.654486) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(1.436208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.8504836) q[0];
sx q[0];
rz(-1.5330667) q[0];
sx q[0];
rz(0.0012375687) q[0];
rz(1.2507863) q[1];
sx q[1];
rz(-2.2092399) q[1];
sx q[1];
rz(-2.1900182) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93336521) q[0];
sx q[0];
rz(-1.7386645) q[0];
sx q[0];
rz(-1.6919447) q[0];
rz(-pi) q[1];
rz(0.95979752) q[2];
sx q[2];
rz(-1.194954) q[2];
sx q[2];
rz(-1.4025337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7965895) q[1];
sx q[1];
rz(-3.0751347) q[1];
sx q[1];
rz(2.4065325) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4755842) q[3];
sx q[3];
rz(-1.9692752) q[3];
sx q[3];
rz(-2.1284895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9963659) q[2];
sx q[2];
rz(-1.9731015) q[2];
sx q[2];
rz(-2.1053947) q[2];
rz(-2.2717617) q[3];
sx q[3];
rz(-1.6706322) q[3];
sx q[3];
rz(-2.3979392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1096126) q[0];
sx q[0];
rz(-0.59015048) q[0];
sx q[0];
rz(2.0330644) q[0];
rz(0.96500665) q[1];
sx q[1];
rz(-1.3771649) q[1];
sx q[1];
rz(-2.5722497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8764502) q[0];
sx q[0];
rz(-2.5776349) q[0];
sx q[0];
rz(-0.41570681) q[0];
rz(-2.0351719) q[2];
sx q[2];
rz(-1.9074252) q[2];
sx q[2];
rz(0.48392933) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9534193) q[1];
sx q[1];
rz(-1.0480289) q[1];
sx q[1];
rz(-0.0095349113) q[1];
rz(-pi) q[2];
rz(0.17579349) q[3];
sx q[3];
rz(-0.7657649) q[3];
sx q[3];
rz(2.898223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1278594) q[2];
sx q[2];
rz(-0.49488417) q[2];
sx q[2];
rz(-0.046517046) q[2];
rz(0.30139309) q[3];
sx q[3];
rz(-1.9207585) q[3];
sx q[3];
rz(-2.9197599) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90060577) q[0];
sx q[0];
rz(-1.931668) q[0];
sx q[0];
rz(1.3264309) q[0];
rz(-1.2251414) q[1];
sx q[1];
rz(-2.6381208) q[1];
sx q[1];
rz(-2.0874964) q[1];
rz(1.3780807) q[2];
sx q[2];
rz(-0.96619923) q[2];
sx q[2];
rz(-0.33351225) q[2];
rz(0.71618373) q[3];
sx q[3];
rz(-2.2436705) q[3];
sx q[3];
rz(2.4498037) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
