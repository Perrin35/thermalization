OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6150317) q[0];
sx q[0];
rz(-0.57305133) q[0];
sx q[0];
rz(-2.2990062) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18544491) q[0];
sx q[0];
rz(-1.7831793) q[0];
sx q[0];
rz(0.74622112) q[0];
x q[1];
rz(-0.027878472) q[2];
sx q[2];
rz(-2.5275702) q[2];
sx q[2];
rz(0.30464722) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7588501) q[1];
sx q[1];
rz(-1.1945801) q[1];
sx q[1];
rz(1.43169) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4487212) q[3];
sx q[3];
rz(-2.1616518) q[3];
sx q[3];
rz(-1.5199682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52790102) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.2256631) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(0.82204449) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(-0.32546145) q[0];
rz(-1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.1546086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75630674) q[0];
sx q[0];
rz(-0.025408832) q[0];
sx q[0];
rz(0.73861648) q[0];
rz(-1.0500533) q[2];
sx q[2];
rz(-2.4467231) q[2];
sx q[2];
rz(0.75887242) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7696015) q[1];
sx q[1];
rz(-2.3725315) q[1];
sx q[1];
rz(0.12188697) q[1];
x q[2];
rz(-1.8335908) q[3];
sx q[3];
rz(-1.3920708) q[3];
sx q[3];
rz(2.5457515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(0.78770351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(-1.6261684) q[0];
rz(2.5405163) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(-2.0498958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9804304) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(0.79361332) q[0];
rz(-1.0513564) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(-1.6348334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25635438) q[1];
sx q[1];
rz(-1.2043722) q[1];
sx q[1];
rz(-0.79856915) q[1];
rz(-pi) q[2];
rz(1.6352429) q[3];
sx q[3];
rz(-1.4430874) q[3];
sx q[3];
rz(2.8116022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(1.7679924) q[3];
sx q[3];
rz(-1.6146086) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3110733) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(-0.4367035) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-0.31035796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46582857) q[0];
sx q[0];
rz(-0.40256631) q[0];
sx q[0];
rz(2.7990544) q[0];
rz(-1.6967625) q[2];
sx q[2];
rz(-2.2519886) q[2];
sx q[2];
rz(-1.3931526) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1102317) q[1];
sx q[1];
rz(-0.35846113) q[1];
sx q[1];
rz(-1.5178174) q[1];
rz(-1.2918606) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13005304) q[2];
sx q[2];
rz(-2.4235642) q[2];
sx q[2];
rz(2.0641573) q[2];
rz(-3.0854026) q[3];
sx q[3];
rz(-0.63779938) q[3];
sx q[3];
rz(-1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524277) q[0];
sx q[0];
rz(-2.0963033) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0349883) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(-0.58017054) q[0];
rz(2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(-0.86415926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.00759) q[1];
sx q[1];
rz(-2.3569199) q[1];
sx q[1];
rz(-2.6919634) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16964511) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(-0.58562216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(-2.8379748) q[3];
sx q[3];
rz(-1.2250591) q[3];
sx q[3];
rz(1.3195066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-0.93739167) q[0];
sx q[0];
rz(2.8836024) q[0];
rz(2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(1.4917096) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314635) q[0];
sx q[0];
rz(-1.9995814) q[0];
sx q[0];
rz(1.6786806) q[0];
rz(-0.67955534) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(2.213775) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3867492) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(1.8317354) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(2.6747243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(3.0498665) q[2];
rz(2.2979459) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3180852) q[0];
sx q[0];
rz(-1.2473236) q[0];
sx q[0];
rz(-0.41123018) q[0];
rz(-2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-3.1076028) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5916409) q[0];
sx q[0];
rz(-0.80054987) q[0];
sx q[0];
rz(-2.9151025) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2529644) q[2];
sx q[2];
rz(-2.3805328) q[2];
sx q[2];
rz(1.6737446) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36044869) q[1];
sx q[1];
rz(-0.89372674) q[1];
sx q[1];
rz(0.064868049) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61492413) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(-1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(0.87654385) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2440764) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(-2.774033) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(-1.6961018) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.53795) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(0.005714697) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58480279) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(2.1540097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8135012) q[1];
sx q[1];
rz(-0.5760759) q[1];
sx q[1];
rz(0.2920132) q[1];
rz(-0.94533841) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8470856) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-2.7344446) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.8898213) q[0];
rz(2.4720526) q[1];
sx q[1];
rz(-1.957683) q[1];
sx q[1];
rz(0.30977419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1447434) q[0];
sx q[0];
rz(-1.0486756) q[0];
sx q[0];
rz(-1.4247308) q[0];
rz(-pi) q[1];
rz(3.1321399) q[2];
sx q[2];
rz(-1.7477034) q[2];
sx q[2];
rz(-2.0131907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.0032776912) q[1];
sx q[1];
rz(-1.7764846) q[1];
sx q[1];
rz(3.0823207) q[1];
rz(-pi) q[2];
rz(1.314332) q[3];
sx q[3];
rz(-0.98955742) q[3];
sx q[3];
rz(2.6687711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70242515) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(-1.9343728) q[2];
rz(-2.1045945) q[3];
sx q[3];
rz(-1.2456649) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(-0.58569113) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(-1.4996128) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6584872) q[0];
sx q[0];
rz(-0.30880901) q[0];
sx q[0];
rz(-1.2312066) q[0];
rz(-0.76403107) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(1.3438091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26350281) q[1];
sx q[1];
rz(-0.51551688) q[1];
sx q[1];
rz(2.0104736) q[1];
x q[2];
rz(1.3979982) q[3];
sx q[3];
rz(-1.5683335) q[3];
sx q[3];
rz(0.60009225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(-0.45599109) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35836999) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(2.6772066) q[2];
sx q[2];
rz(-2.0225564) q[2];
sx q[2];
rz(0.49973942) q[2];
rz(1.0740888) q[3];
sx q[3];
rz(-2.2278193) q[3];
sx q[3];
rz(-0.57886119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
