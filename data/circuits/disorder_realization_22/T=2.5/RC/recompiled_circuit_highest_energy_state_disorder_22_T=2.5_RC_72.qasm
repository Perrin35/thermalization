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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(1.2494614) q[0];
rz(3.0013822) q[1];
sx q[1];
rz(-1.4445211) q[1];
sx q[1];
rz(3.0797449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1309146) q[0];
sx q[0];
rz(-2.317551) q[0];
sx q[0];
rz(2.2758621) q[0];
x q[1];
rz(-3.1377162) q[2];
sx q[2];
rz(-3.0511694) q[2];
sx q[2];
rz(0.10744444) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.34644352) q[1];
sx q[1];
rz(-2.4137874) q[1];
sx q[1];
rz(-1.7176571) q[1];
x q[2];
rz(-0.40932406) q[3];
sx q[3];
rz(-2.217389) q[3];
sx q[3];
rz(-2.4149946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7141815) q[2];
sx q[2];
rz(-2.2990871) q[2];
sx q[2];
rz(-2.8851435) q[2];
rz(0.17476684) q[3];
sx q[3];
rz(-1.1656961) q[3];
sx q[3];
rz(2.3037361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3278219) q[0];
sx q[0];
rz(-2.796266) q[0];
sx q[0];
rz(-0.12369618) q[0];
rz(2.9028614) q[1];
sx q[1];
rz(-0.78014603) q[1];
sx q[1];
rz(3.0366268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5803842) q[0];
sx q[0];
rz(-2.4260902) q[0];
sx q[0];
rz(-1.8010048) q[0];
rz(-pi) q[1];
rz(2.6506422) q[2];
sx q[2];
rz(-2.6770795) q[2];
sx q[2];
rz(-2.2238942) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2626896) q[1];
sx q[1];
rz(-2.5715552) q[1];
sx q[1];
rz(-1.0703342) q[1];
x q[2];
rz(0.70521783) q[3];
sx q[3];
rz(-1.8466443) q[3];
sx q[3];
rz(2.5666756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3671942) q[2];
sx q[2];
rz(-1.6509193) q[2];
sx q[2];
rz(1.4614089) q[2];
rz(-1.8558308) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(-2.8694966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937781) q[0];
sx q[0];
rz(-2.9312134) q[0];
sx q[0];
rz(2.126597) q[0];
rz(2.2442832) q[1];
sx q[1];
rz(-1.4941447) q[1];
sx q[1];
rz(-2.4651333) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7583026) q[0];
sx q[0];
rz(-1.5967909) q[0];
sx q[0];
rz(-2.1878408) q[0];
x q[1];
rz(0.49336596) q[2];
sx q[2];
rz(-2.1378433) q[2];
sx q[2];
rz(0.21304785) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3520917) q[1];
sx q[1];
rz(-1.3695696) q[1];
sx q[1];
rz(1.9580864) q[1];
rz(-pi) q[2];
rz(-1.5514938) q[3];
sx q[3];
rz(-0.75452828) q[3];
sx q[3];
rz(1.4135264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2066388) q[2];
sx q[2];
rz(-0.97524869) q[2];
sx q[2];
rz(-2.3812531) q[2];
rz(1.9350516) q[3];
sx q[3];
rz(-2.3307266) q[3];
sx q[3];
rz(-2.2981203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43561414) q[0];
sx q[0];
rz(-0.62788457) q[0];
sx q[0];
rz(-1.5950369) q[0];
rz(-2.4226923) q[1];
sx q[1];
rz(-1.3201821) q[1];
sx q[1];
rz(2.9100606) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0615094) q[0];
sx q[0];
rz(-0.69893796) q[0];
sx q[0];
rz(-2.438758) q[0];
x q[1];
rz(0.79265742) q[2];
sx q[2];
rz(-1.5680299) q[2];
sx q[2];
rz(0.94209988) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.097024767) q[1];
sx q[1];
rz(-0.49355727) q[1];
sx q[1];
rz(-0.19792168) q[1];
rz(-pi) q[2];
rz(1.9636964) q[3];
sx q[3];
rz(-2.1390954) q[3];
sx q[3];
rz(-1.5400122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.71451688) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(0.60849774) q[2];
rz(-0.19615873) q[3];
sx q[3];
rz(-1.720263) q[3];
sx q[3];
rz(2.1430446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2732368) q[0];
sx q[0];
rz(-2.4899857) q[0];
sx q[0];
rz(-2.3625145) q[0];
rz(-0.92998663) q[1];
sx q[1];
rz(-1.168074) q[1];
sx q[1];
rz(1.9680061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54544696) q[0];
sx q[0];
rz(-1.9088424) q[0];
sx q[0];
rz(3.1374187) q[0];
rz(-pi) q[1];
rz(-1.7601682) q[2];
sx q[2];
rz(-2.5177285) q[2];
sx q[2];
rz(0.80846918) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.25949892) q[1];
sx q[1];
rz(-0.32392392) q[1];
sx q[1];
rz(-1.367635) q[1];
rz(-0.42901943) q[3];
sx q[3];
rz(-1.7439505) q[3];
sx q[3];
rz(-0.28399434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8194627) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(2.4665311) q[2];
rz(2.2760462) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(1.9338098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.753767) q[0];
sx q[0];
rz(-0.017711552) q[0];
sx q[0];
rz(-0.87131635) q[0];
rz(0.21698347) q[1];
sx q[1];
rz(-1.5419518) q[1];
sx q[1];
rz(1.3380231) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3335452) q[0];
sx q[0];
rz(-1.4879972) q[0];
sx q[0];
rz(-0.6525349) q[0];
rz(-2.2964444) q[2];
sx q[2];
rz(-1.4652227) q[2];
sx q[2];
rz(-2.0279864) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8585526) q[1];
sx q[1];
rz(-1.5225019) q[1];
sx q[1];
rz(-1.7142536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7609672) q[3];
sx q[3];
rz(-2.0479879) q[3];
sx q[3];
rz(-0.3482477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1285105) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(-0.80005542) q[2];
rz(-0.99177805) q[3];
sx q[3];
rz(-0.9477152) q[3];
sx q[3];
rz(0.94314027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25471383) q[0];
sx q[0];
rz(-0.44097057) q[0];
sx q[0];
rz(2.9898341) q[0];
rz(1.4166547) q[1];
sx q[1];
rz(-0.85985008) q[1];
sx q[1];
rz(-2.3053665) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72513103) q[0];
sx q[0];
rz(-0.90354474) q[0];
sx q[0];
rz(-0.50531549) q[0];
rz(-pi) q[1];
rz(1.5791527) q[2];
sx q[2];
rz(-2.1474693) q[2];
sx q[2];
rz(1.6461262) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6265833) q[1];
sx q[1];
rz(-1.8191511) q[1];
sx q[1];
rz(-0.088661389) q[1];
rz(-pi) q[2];
rz(1.7015905) q[3];
sx q[3];
rz(-2.0328641) q[3];
sx q[3];
rz(0.82926428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6553216) q[2];
sx q[2];
rz(-3.0781015) q[2];
sx q[2];
rz(3.0625694) q[2];
rz(0.71845734) q[3];
sx q[3];
rz(-1.6301165) q[3];
sx q[3];
rz(-2.2169936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543168) q[0];
sx q[0];
rz(-2.5434255) q[0];
sx q[0];
rz(-0.86791903) q[0];
rz(2.4527841) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(-2.205663) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3132738) q[0];
sx q[0];
rz(-0.5567282) q[0];
sx q[0];
rz(-2.7931045) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8008119) q[2];
sx q[2];
rz(-1.9592957) q[2];
sx q[2];
rz(2.4430371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0088499) q[1];
sx q[1];
rz(-2.8121083) q[1];
sx q[1];
rz(1.4047755) q[1];
x q[2];
rz(-0.74347382) q[3];
sx q[3];
rz(-1.5718565) q[3];
sx q[3];
rz(-0.27976945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9800637) q[2];
sx q[2];
rz(-0.65902013) q[2];
sx q[2];
rz(2.8847983) q[2];
rz(-1.0181381) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(-2.170678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37860206) q[0];
sx q[0];
rz(-0.55736962) q[0];
sx q[0];
rz(-2.2879404) q[0];
rz(1.254982) q[1];
sx q[1];
rz(-1.9957142) q[1];
sx q[1];
rz(-0.34057158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7434692) q[0];
sx q[0];
rz(-0.2532244) q[0];
sx q[0];
rz(-2.2608093) q[0];
rz(-1.8179632) q[2];
sx q[2];
rz(-2.2296612) q[2];
sx q[2];
rz(0.62550046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3601968) q[1];
sx q[1];
rz(-0.92222795) q[1];
sx q[1];
rz(0.99332033) q[1];
rz(-pi) q[2];
rz(2.0331618) q[3];
sx q[3];
rz(-1.3232854) q[3];
sx q[3];
rz(2.4570297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.42369947) q[2];
sx q[2];
rz(-1.7607949) q[2];
sx q[2];
rz(-0.37810668) q[2];
rz(-2.9041491) q[3];
sx q[3];
rz(-1.0708555) q[3];
sx q[3];
rz(0.28089359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0963999) q[0];
sx q[0];
rz(-1.0843596) q[0];
sx q[0];
rz(2.4943446) q[0];
rz(1.1680394) q[1];
sx q[1];
rz(-0.85473514) q[1];
sx q[1];
rz(-0.38356575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952571) q[0];
sx q[0];
rz(-1.5492166) q[0];
sx q[0];
rz(-3.1109592) q[0];
rz(-1.9012544) q[2];
sx q[2];
rz(-1.0631732) q[2];
sx q[2];
rz(-0.79703813) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3473222) q[1];
sx q[1];
rz(-1.8350661) q[1];
sx q[1];
rz(-1.7842954) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0553352) q[3];
sx q[3];
rz(-1.664242) q[3];
sx q[3];
rz(-2.6938113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8770404) q[2];
sx q[2];
rz(-2.2525747) q[2];
sx q[2];
rz(-1.5846579) q[2];
rz(-1.6282188) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(0.42678601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7814816) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(-1.9602641) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(2.0170586) q[2];
sx q[2];
rz(-2.6206803) q[2];
sx q[2];
rz(-0.16823106) q[2];
rz(-2.6570126) q[3];
sx q[3];
rz(-1.5427586) q[3];
sx q[3];
rz(1.2324738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
