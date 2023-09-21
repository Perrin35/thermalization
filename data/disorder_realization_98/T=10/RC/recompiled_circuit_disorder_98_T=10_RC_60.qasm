OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5877085) q[0];
sx q[0];
rz(-2.2061078) q[0];
sx q[0];
rz(2.5250838) q[0];
rz(-pi) q[1];
rz(0.82586536) q[2];
sx q[2];
rz(-2.4561433) q[2];
sx q[2];
rz(-0.65537383) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.95853165) q[1];
sx q[1];
rz(-1.5325938) q[1];
sx q[1];
rz(1.7512291) q[1];
x q[2];
rz(-2.1211795) q[3];
sx q[3];
rz(-1.7341511) q[3];
sx q[3];
rz(2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-0.70409888) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(-0.63252226) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(0.65223637) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-3.1104381) q[0];
sx q[0];
rz(-0.40701436) q[0];
rz(-pi) q[1];
rz(-2.9035283) q[2];
sx q[2];
rz(-1.2859584) q[2];
sx q[2];
rz(0.1711947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3398509) q[1];
sx q[1];
rz(-0.36777126) q[1];
sx q[1];
rz(-2.4778609) q[1];
rz(-pi) q[2];
rz(2.2964301) q[3];
sx q[3];
rz(-0.80544986) q[3];
sx q[3];
rz(2.1243387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(-0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.7920378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(-1.8750989) q[0];
rz(3.1042728) q[2];
sx q[2];
rz(-1.5048426) q[2];
sx q[2];
rz(-1.0730336) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9393443) q[1];
sx q[1];
rz(-1.0679686) q[1];
sx q[1];
rz(1.3015675) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1317741) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(0.94235086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(2.2375977) q[2];
rz(2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.3999456) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(1.7310463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.3316863) q[1];
sx q[1];
rz(3.1052123) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2934389) q[0];
sx q[0];
rz(-1.358195) q[0];
sx q[0];
rz(0.92377499) q[0];
rz(-pi) q[1];
rz(-1.8965917) q[2];
sx q[2];
rz(-1.0137644) q[2];
sx q[2];
rz(1.4404802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7155647) q[1];
sx q[1];
rz(-1.4392816) q[1];
sx q[1];
rz(2.2711666) q[1];
x q[2];
rz(1.6691469) q[3];
sx q[3];
rz(-1.3373168) q[3];
sx q[3];
rz(-2.1360306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(-2.4711117) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7017512) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-2.9751076) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-0.23434815) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7974632) q[0];
sx q[0];
rz(-1.2327694) q[0];
sx q[0];
rz(0.56030886) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7469823) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(-1.8793775) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6926596) q[1];
sx q[1];
rz(-1.7198791) q[1];
sx q[1];
rz(-2.3922608) q[1];
rz(1.2469532) q[3];
sx q[3];
rz(-1.6456592) q[3];
sx q[3];
rz(2.4002241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-0.22053545) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-0.76550686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(0.81714001) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(1.9979427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4139347) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(1.6137705) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8311062) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(-1.3442163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5864582) q[1];
sx q[1];
rz(-1.1291593) q[1];
sx q[1];
rz(-3.0416136) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2915217) q[3];
sx q[3];
rz(-2.0719299) q[3];
sx q[3];
rz(0.53370332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75366655) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-2.6521519) q[2];
rz(2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.2868767) q[0];
rz(0.6634179) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(1.9082665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7466) q[0];
sx q[0];
rz(-1.9548423) q[0];
sx q[0];
rz(0.011944255) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2912824) q[2];
sx q[2];
rz(-2.0680973) q[2];
sx q[2];
rz(2.8720299) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7066321) q[1];
sx q[1];
rz(-0.11206493) q[1];
sx q[1];
rz(-3.0155165) q[1];
rz(-pi) q[2];
rz(-1.2193905) q[3];
sx q[3];
rz(-2.0316342) q[3];
sx q[3];
rz(1.2724862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(0.070090381) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
rz(3.0163517) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(-0.21487543) q[1];
sx q[1];
rz(-0.75526777) q[1];
sx q[1];
rz(-1.8833556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94851516) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(-1.8817188) q[0];
rz(-pi) q[1];
rz(-0.9332946) q[2];
sx q[2];
rz(-2.1142695) q[2];
sx q[2];
rz(0.18804929) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80553493) q[1];
sx q[1];
rz(-1.6748669) q[1];
sx q[1];
rz(-2.59471) q[1];
x q[2];
rz(-0.97790896) q[3];
sx q[3];
rz(-2.3330354) q[3];
sx q[3];
rz(0.18329328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.68226472) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(2.5788467) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-2.0195885) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(-2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(3.1138611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77596091) q[0];
sx q[0];
rz(-0.86254518) q[0];
sx q[0];
rz(-2.2328949) q[0];
rz(-pi) q[1];
rz(-2.8743923) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(-2.6097678) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2003277) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(-1.4501249) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5042217) q[3];
sx q[3];
rz(-1.596631) q[3];
sx q[3];
rz(2.0500101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-2.2843602) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89850241) q[0];
sx q[0];
rz(-1.9542964) q[0];
sx q[0];
rz(-2.6570508) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72980482) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(-2.6953816) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32138667) q[1];
sx q[1];
rz(-0.69637978) q[1];
sx q[1];
rz(0.022298261) q[1];
rz(2.8994843) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(-0.20413354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-2.5349687) q[2];
rz(2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-0.84038466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3474779) q[0];
sx q[0];
rz(-1.62962) q[0];
sx q[0];
rz(-0.99123065) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(2.6735641) q[2];
sx q[2];
rz(-1.9149018) q[2];
sx q[2];
rz(0.55185774) q[2];
rz(-2.8787981) q[3];
sx q[3];
rz(-2.0398718) q[3];
sx q[3];
rz(-0.3077988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
