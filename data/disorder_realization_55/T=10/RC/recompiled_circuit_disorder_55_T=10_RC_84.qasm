OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5947333) q[0];
sx q[0];
rz(-1.5164627) q[0];
sx q[0];
rz(-2.8773142) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7289294) q[0];
sx q[0];
rz(-2.465415) q[0];
sx q[0];
rz(2.9039608) q[0];
rz(-pi) q[1];
rz(-2.0779607) q[2];
sx q[2];
rz(-2.1858366) q[2];
sx q[2];
rz(2.087649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4525675) q[1];
sx q[1];
rz(-1.2847632) q[1];
sx q[1];
rz(-0.16098117) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5988889) q[3];
sx q[3];
rz(-2.1510604) q[3];
sx q[3];
rz(-2.7147646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4720817) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(-1.2708698) q[3];
sx q[3];
rz(-1.9138252) q[3];
sx q[3];
rz(-0.27403533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1141777) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(2.6787058) q[1];
sx q[1];
rz(-1.0375689) q[1];
sx q[1];
rz(-2.8754821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5824495) q[0];
sx q[0];
rz(-0.83280116) q[0];
sx q[0];
rz(-2.0931431) q[0];
rz(-pi) q[1];
rz(1.3219464) q[2];
sx q[2];
rz(-2.5517533) q[2];
sx q[2];
rz(0.8283386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7539566) q[1];
sx q[1];
rz(-1.0322744) q[1];
sx q[1];
rz(2.3073879) q[1];
x q[2];
rz(-0.18272419) q[3];
sx q[3];
rz(-2.1696739) q[3];
sx q[3];
rz(-0.11174186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(0.51149386) q[2];
rz(0.809508) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(2.8539343) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70616102) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(0.92873746) q[0];
rz(-1.7354895) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(-1.3471289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0179694) q[0];
sx q[0];
rz(-2.6051913) q[0];
sx q[0];
rz(-0.53039741) q[0];
rz(-2.9834619) q[2];
sx q[2];
rz(-1.2889382) q[2];
sx q[2];
rz(-0.49027298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.128708) q[1];
sx q[1];
rz(-1.8599659) q[1];
sx q[1];
rz(-0.94990001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15123983) q[3];
sx q[3];
rz(-2.2721707) q[3];
sx q[3];
rz(2.410694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1469664) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(0.49595293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3494444) q[0];
sx q[0];
rz(-1.9932207) q[0];
sx q[0];
rz(0.96570063) q[0];
rz(2.4194338) q[1];
sx q[1];
rz(-1.5042217) q[1];
sx q[1];
rz(0.55975634) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20630079) q[0];
sx q[0];
rz(-1.9959404) q[0];
sx q[0];
rz(1.6598808) q[0];
rz(0.81360929) q[2];
sx q[2];
rz(-2.1289729) q[2];
sx q[2];
rz(1.3249601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2019129) q[1];
sx q[1];
rz(-0.61720467) q[1];
sx q[1];
rz(-1.4918785) q[1];
x q[2];
rz(2.9066554) q[3];
sx q[3];
rz(-0.80312356) q[3];
sx q[3];
rz(2.0149751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42797783) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.4245865) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150862) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(-0.36636233) q[0];
rz(-1.5461961) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(-0.34367925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3643091) q[0];
sx q[0];
rz(-1.0316348) q[0];
sx q[0];
rz(-0.76215141) q[0];
x q[1];
rz(0.91243773) q[2];
sx q[2];
rz(-2.3235333) q[2];
sx q[2];
rz(-0.91615265) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74663631) q[1];
sx q[1];
rz(-2.6379105) q[1];
sx q[1];
rz(-1.4175182) q[1];
rz(-pi) q[2];
x q[2];
rz(3.120317) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(-1.4847886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3917824) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(3.0878477) q[2];
rz(1.404445) q[3];
sx q[3];
rz(-2.6440547) q[3];
sx q[3];
rz(-0.29156175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6546201) q[0];
sx q[0];
rz(-2.4916861) q[0];
sx q[0];
rz(0.75575954) q[0];
rz(3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(-0.25973928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1007337) q[0];
sx q[0];
rz(-1.9319527) q[0];
sx q[0];
rz(0.33613899) q[0];
rz(-pi) q[1];
rz(1.2345215) q[2];
sx q[2];
rz(-1.5521126) q[2];
sx q[2];
rz(2.3980354) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7163135) q[1];
sx q[1];
rz(-1.7967766) q[1];
sx q[1];
rz(1.9081566) q[1];
rz(-pi) q[2];
rz(1.1214244) q[3];
sx q[3];
rz(-1.5869889) q[3];
sx q[3];
rz(-1.3390954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(0.10061131) q[2];
rz(-2.9595024) q[3];
sx q[3];
rz(-2.263335) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9942193) q[0];
sx q[0];
rz(-1.4202776) q[0];
sx q[0];
rz(2.8314262) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(2.535634) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2180442) q[0];
sx q[0];
rz(-0.96203066) q[0];
sx q[0];
rz(-1.2952842) q[0];
rz(-pi) q[1];
rz(1.7714959) q[2];
sx q[2];
rz(-1.2846652) q[2];
sx q[2];
rz(1.2933033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76304945) q[1];
sx q[1];
rz(-2.0705283) q[1];
sx q[1];
rz(0.39079697) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0244272) q[3];
sx q[3];
rz(-2.0351279) q[3];
sx q[3];
rz(-1.6950316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(-0.85285464) q[2];
rz(1.3700221) q[3];
sx q[3];
rz(-1.6829237) q[3];
sx q[3];
rz(-2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9119499) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(-1.461331) q[0];
rz(-1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(-3.0775552) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56211573) q[0];
sx q[0];
rz(-2.1770283) q[0];
sx q[0];
rz(-0.22867815) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0649101) q[2];
sx q[2];
rz(-1.3008899) q[2];
sx q[2];
rz(3.0907061) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.774051) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(-2.7864085) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.039489432) q[3];
sx q[3];
rz(-1.2590623) q[3];
sx q[3];
rz(-1.5453651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0503851) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(1.5554265) q[2];
rz(-0.88820109) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(-0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9973307) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-2.0137537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9040684) q[0];
sx q[0];
rz(-1.5547817) q[0];
sx q[0];
rz(2.9928656) q[0];
rz(-2.2531829) q[2];
sx q[2];
rz(-0.75196224) q[2];
sx q[2];
rz(-2.106452) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7833264) q[1];
sx q[1];
rz(-1.3291385) q[1];
sx q[1];
rz(-3.1296455) q[1];
rz(2.1209675) q[3];
sx q[3];
rz(-2.5084825) q[3];
sx q[3];
rz(1.8819295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-2.2884559) q[2];
sx q[2];
rz(-0.17364994) q[2];
rz(-2.8052143) q[3];
sx q[3];
rz(-1.9189546) q[3];
sx q[3];
rz(0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4353452) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(-0.82410518) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-0.5724268) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7596282) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(2.1419924) q[0];
rz(-1.1119214) q[2];
sx q[2];
rz(-1.3238812) q[2];
sx q[2];
rz(-2.0725476) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6483557) q[1];
sx q[1];
rz(-1.8093458) q[1];
sx q[1];
rz(1.3484713) q[1];
rz(-0.70865788) q[3];
sx q[3];
rz(-1.7950247) q[3];
sx q[3];
rz(1.538016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77999014) q[2];
sx q[2];
rz(-0.47452351) q[2];
sx q[2];
rz(-0.53722107) q[2];
rz(2.0843263) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2873516) q[0];
sx q[0];
rz(-1.9762522) q[0];
sx q[0];
rz(1.5594788) q[0];
rz(2.6782425) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(0.097975227) q[2];
sx q[2];
rz(-2.2231979) q[2];
sx q[2];
rz(-2.4616432) q[2];
rz(-1.6871917) q[3];
sx q[3];
rz(-2.6506861) q[3];
sx q[3];
rz(-0.41674137) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
