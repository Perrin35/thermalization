OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.54685932) q[0];
sx q[0];
rz(-1.62513) q[0];
sx q[0];
rz(-0.2642785) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(-1.2101313) q[1];
sx q[1];
rz(0.73524737) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0300232) q[0];
sx q[0];
rz(-0.91696793) q[0];
sx q[0];
rz(-1.3841188) q[0];
rz(-2.0779607) q[2];
sx q[2];
rz(-2.1858366) q[2];
sx q[2];
rz(2.087649) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83598614) q[1];
sx q[1];
rz(-1.7251833) q[1];
sx q[1];
rz(1.2812213) q[1];
rz(-pi) q[2];
rz(1.5427038) q[3];
sx q[3];
rz(-2.1510604) q[3];
sx q[3];
rz(0.42682808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.66951093) q[2];
sx q[2];
rz(-1.3005723) q[2];
sx q[2];
rz(2.0377339) q[2];
rz(1.8707229) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(-2.8675573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0274149) q[0];
sx q[0];
rz(-2.6129621) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(0.46288681) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(0.26611051) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55914315) q[0];
sx q[0];
rz(-2.3087915) q[0];
sx q[0];
rz(2.0931431) q[0];
x q[1];
rz(2.9782148) q[2];
sx q[2];
rz(-2.1401569) q[2];
sx q[2];
rz(2.6100104) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7539566) q[1];
sx q[1];
rz(-2.1093183) q[1];
sx q[1];
rz(0.83420475) q[1];
rz(-1.3105884) q[3];
sx q[3];
rz(-2.5187413) q[3];
sx q[3];
rz(-0.42850307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(0.51149386) q[2];
rz(-2.3320847) q[3];
sx q[3];
rz(-1.5313238) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4354316) q[0];
sx q[0];
rz(-1.5478739) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(1.7354895) q[1];
sx q[1];
rz(-2.4415253) q[1];
sx q[1];
rz(1.7944638) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72235332) q[0];
sx q[0];
rz(-2.0273211) q[0];
sx q[0];
rz(-1.8629575) q[0];
rz(-pi) q[1];
rz(-1.8560266) q[2];
sx q[2];
rz(-1.7226379) q[2];
sx q[2];
rz(-1.0361995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9007064) q[1];
sx q[1];
rz(-0.97929231) q[1];
sx q[1];
rz(-0.35066168) q[1];
rz(-pi) q[2];
rz(1.3942765) q[3];
sx q[3];
rz(-0.71478292) q[3];
sx q[3];
rz(2.178758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(-0.26432031) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(-2.6456397) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(0.96570063) q[0];
rz(2.4194338) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(-0.55975634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7402732) q[0];
sx q[0];
rz(-1.4896605) q[0];
sx q[0];
rz(2.7149537) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70978769) q[2];
sx q[2];
rz(-2.1927532) q[2];
sx q[2];
rz(-0.70993916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2019129) q[1];
sx q[1];
rz(-0.61720467) q[1];
sx q[1];
rz(1.6497142) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9066554) q[3];
sx q[3];
rz(-2.3384691) q[3];
sx q[3];
rz(2.0149751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7136148) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(-1.7170061) q[2];
rz(-0.26040855) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99073064) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(-2.7752303) q[0];
rz(1.5461961) q[1];
sx q[1];
rz(-0.55391824) q[1];
sx q[1];
rz(-0.34367925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2491964) q[0];
sx q[0];
rz(-0.93660347) q[0];
sx q[0];
rz(-2.2618494) q[0];
rz(-0.91243773) q[2];
sx q[2];
rz(-0.81805938) q[2];
sx q[2];
rz(-0.91615265) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5695565) q[1];
sx q[1];
rz(-2.0680288) q[1];
sx q[1];
rz(-3.0576502) q[1];
x q[2];
rz(-3.120317) q[3];
sx q[3];
rz(-2.2145503) q[3];
sx q[3];
rz(-1.6568041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3917824) q[2];
sx q[2];
rz(-0.85931531) q[2];
sx q[2];
rz(0.053744944) q[2];
rz(1.7371477) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(-0.29156175) q[3];
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
rz(-pi) q[3];
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
rz(2.6546201) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(0.75575954) q[0];
rz(-0.02515633) q[1];
sx q[1];
rz(-2.2143366) q[1];
sx q[1];
rz(0.25973928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.320967) q[0];
sx q[0];
rz(-0.48829406) q[0];
sx q[0];
rz(-2.2886306) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1218006) q[2];
sx q[2];
rz(-1.2345825) q[2];
sx q[2];
rz(2.3078231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7163135) q[1];
sx q[1];
rz(-1.7967766) q[1];
sx q[1];
rz(-1.9081566) q[1];
rz(-pi) q[2];
rz(2.0201683) q[3];
sx q[3];
rz(-1.5869889) q[3];
sx q[3];
rz(1.3390954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48866895) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-3.0409813) q[2];
rz(-0.18209022) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.7939059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1473734) q[0];
sx q[0];
rz(-1.7213151) q[0];
sx q[0];
rz(-0.31016645) q[0];
rz(-0.50225964) q[1];
sx q[1];
rz(-2.6175833) q[1];
sx q[1];
rz(-0.60595864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46442407) q[0];
sx q[0];
rz(-0.66093984) q[0];
sx q[0];
rz(2.7695157) q[0];
x q[1];
rz(1.7714959) q[2];
sx q[2];
rz(-1.2846652) q[2];
sx q[2];
rz(1.2933033) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1389321) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(2.1041811) q[1];
rz(-pi) q[2];
rz(-2.3378387) q[3];
sx q[3];
rz(-2.440212) q[3];
sx q[3];
rz(-2.3826117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24017748) q[2];
sx q[2];
rz(-1.1837974) q[2];
sx q[2];
rz(2.288738) q[2];
rz(1.7715706) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(-2.9366233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9119499) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(1.6802616) q[0];
rz(-1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(0.064037474) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56211573) q[0];
sx q[0];
rz(-0.96456438) q[0];
sx q[0];
rz(-2.9129145) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0649101) q[2];
sx q[2];
rz(-1.3008899) q[2];
sx q[2];
rz(-3.0907061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.774051) q[1];
sx q[1];
rz(-1.5343752) q[1];
sx q[1];
rz(0.35518412) q[1];
x q[2];
rz(1.6927034) q[3];
sx q[3];
rz(-2.8274483) q[3];
sx q[3];
rz(1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.091207592) q[2];
sx q[2];
rz(-2.51077) q[2];
sx q[2];
rz(-1.5554265) q[2];
rz(-2.2533916) q[3];
sx q[3];
rz(-1.9017838) q[3];
sx q[3];
rz(-0.92938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9973307) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.7096747) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(-1.127839) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8059212) q[0];
sx q[0];
rz(-1.4220884) q[0];
sx q[0];
rz(-1.554603) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1986507) q[2];
sx q[2];
rz(-1.1254416) q[2];
sx q[2];
rz(3.1415591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3582663) q[1];
sx q[1];
rz(-1.3291385) q[1];
sx q[1];
rz(-3.1296455) q[1];
rz(-pi) q[2];
rz(-2.1298218) q[3];
sx q[3];
rz(-1.8852919) q[3];
sx q[3];
rz(-0.77034706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.613712) q[2];
sx q[2];
rz(-0.85313672) q[2];
sx q[2];
rz(2.9679427) q[2];
rz(-0.33637834) q[3];
sx q[3];
rz(-1.222638) q[3];
sx q[3];
rz(0.39150795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-2.5543537) q[0];
sx q[0];
rz(-1.6760814) q[0];
rz(2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-0.5724268) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7596282) q[0];
sx q[0];
rz(-2.1158764) q[0];
sx q[0];
rz(-0.99960021) q[0];
rz(-pi) q[1];
rz(-2.8675251) q[2];
sx q[2];
rz(-1.1268508) q[2];
sx q[2];
rz(0.62192164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.010667) q[1];
sx q[1];
rz(-1.7867242) q[1];
sx q[1];
rz(2.8972577) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4329348) q[3];
sx q[3];
rz(-1.346568) q[3];
sx q[3];
rz(-1.6035767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(-2.6043716) q[2];
rz(2.0843263) q[3];
sx q[3];
rz(-0.89151645) q[3];
sx q[3];
rz(-2.4479772) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.854241) q[0];
sx q[0];
rz(-1.1653405) q[0];
sx q[0];
rz(-1.5821138) q[0];
rz(-0.46335012) q[1];
sx q[1];
rz(-0.87711038) q[1];
sx q[1];
rz(-1.6323485) q[1];
rz(-2.2255185) q[2];
sx q[2];
rz(-1.6486042) q[2];
sx q[2];
rz(-0.83124607) q[2];
rz(-0.061999576) q[3];
sx q[3];
rz(-1.0835032) q[3];
sx q[3];
rz(-0.54855357) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];