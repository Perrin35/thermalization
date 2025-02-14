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
rz(2.7895522) q[0];
sx q[0];
rz(-2.3166603) q[0];
sx q[0];
rz(2.6068249) q[0];
rz(-2.1575902) q[1];
sx q[1];
rz(-0.50552955) q[1];
sx q[1];
rz(-1.2619789) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5231294) q[0];
sx q[0];
rz(-0.97181706) q[0];
sx q[0];
rz(1.8904786) q[0];
rz(2.2625669) q[2];
sx q[2];
rz(-2.3874385) q[2];
sx q[2];
rz(1.9944432) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.64691862) q[1];
sx q[1];
rz(-0.54975408) q[1];
sx q[1];
rz(-0.099262909) q[1];
rz(0.38930157) q[3];
sx q[3];
rz(-1.3817781) q[3];
sx q[3];
rz(2.1019328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8925605) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(-1.4830291) q[2];
rz(0.28462166) q[3];
sx q[3];
rz(-2.5183545) q[3];
sx q[3];
rz(2.0774138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738709) q[0];
sx q[0];
rz(-2.4148648) q[0];
sx q[0];
rz(-3.100585) q[0];
rz(1.2076591) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(1.6809195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6448633) q[0];
sx q[0];
rz(-1.4251727) q[0];
sx q[0];
rz(-1.54099) q[0];
x q[1];
rz(-1.0233489) q[2];
sx q[2];
rz(-2.259575) q[2];
sx q[2];
rz(0.42519409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1665312) q[1];
sx q[1];
rz(-1.6094065) q[1];
sx q[1];
rz(-1.8273942) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19537851) q[3];
sx q[3];
rz(-0.88802216) q[3];
sx q[3];
rz(-1.186779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4760806) q[2];
sx q[2];
rz(-1.4613232) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(-0.0811854) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(1.7056874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3944893) q[0];
sx q[0];
rz(-0.012270027) q[0];
sx q[0];
rz(2.5315206) q[0];
rz(-1.3129781) q[1];
sx q[1];
rz(-2.4459631) q[1];
sx q[1];
rz(2.5909766) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2064821) q[0];
sx q[0];
rz(-1.2121823) q[0];
sx q[0];
rz(-2.4367843) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32260334) q[2];
sx q[2];
rz(-1.8292973) q[2];
sx q[2];
rz(-2.5532364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0455835) q[1];
sx q[1];
rz(-2.7852045) q[1];
sx q[1];
rz(-2.136904) q[1];
rz(-pi) q[2];
rz(-1.7254225) q[3];
sx q[3];
rz(-2.3904413) q[3];
sx q[3];
rz(0.42645833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6637471) q[2];
sx q[2];
rz(-2.9691594) q[2];
sx q[2];
rz(1.0663859) q[2];
rz(-1.1586698) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(0.042044736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.31024194) q[0];
sx q[0];
rz(-0.61315918) q[0];
sx q[0];
rz(-0.9170652) q[0];
rz(-2.4861368) q[1];
sx q[1];
rz(-2.2323445) q[1];
sx q[1];
rz(-0.38958946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3026019) q[0];
sx q[0];
rz(-0.20379681) q[0];
sx q[0];
rz(-1.1558644) q[0];
rz(-2.5823103) q[2];
sx q[2];
rz(-0.96508316) q[2];
sx q[2];
rz(0.16597151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.37871088) q[1];
sx q[1];
rz(-0.40232752) q[1];
sx q[1];
rz(-1.6454003) q[1];
x q[2];
rz(-2.5960931) q[3];
sx q[3];
rz(-2.0974468) q[3];
sx q[3];
rz(-0.57034501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83170825) q[2];
sx q[2];
rz(-2.1990621) q[2];
sx q[2];
rz(1.3523098) q[2];
rz(2.3934707) q[3];
sx q[3];
rz(-1.8375405) q[3];
sx q[3];
rz(1.6149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-2.5717773) q[0];
sx q[0];
rz(-1.8001528) q[0];
rz(-0.98126423) q[1];
sx q[1];
rz(-1.2813247) q[1];
sx q[1];
rz(0.70995465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7959952) q[0];
sx q[0];
rz(-1.5516722) q[0];
sx q[0];
rz(2.8157793) q[0];
rz(-pi) q[1];
rz(1.9290438) q[2];
sx q[2];
rz(-0.60645559) q[2];
sx q[2];
rz(-2.0403634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4057096) q[1];
sx q[1];
rz(-0.27389474) q[1];
sx q[1];
rz(2.2475858) q[1];
x q[2];
rz(-3.0026765) q[3];
sx q[3];
rz(-0.60231042) q[3];
sx q[3];
rz(-1.7642782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.183341) q[2];
sx q[2];
rz(-2.3812582) q[2];
sx q[2];
rz(2.23488) q[2];
rz(2.9592311) q[3];
sx q[3];
rz(-1.4401108) q[3];
sx q[3];
rz(0.85842925) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7957423) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(0.20172754) q[0];
rz(-1.3564823) q[1];
sx q[1];
rz(-2.4701665) q[1];
sx q[1];
rz(-1.1511525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7089823) q[0];
sx q[0];
rz(-1.1469736) q[0];
sx q[0];
rz(0.2291542) q[0];
rz(-0.62080748) q[2];
sx q[2];
rz(-0.70638958) q[2];
sx q[2];
rz(-2.8548129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21143895) q[1];
sx q[1];
rz(-1.0027998) q[1];
sx q[1];
rz(2.0536973) q[1];
rz(2.6309507) q[3];
sx q[3];
rz(-1.770442) q[3];
sx q[3];
rz(-2.8996403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1543033) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(0.34745535) q[2];
rz(-0.82184982) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(1.52389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3407985) q[0];
sx q[0];
rz(-0.64685416) q[0];
sx q[0];
rz(1.1232173) q[0];
rz(-2.7128291) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(-2.9926328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65229416) q[0];
sx q[0];
rz(-2.6241488) q[0];
sx q[0];
rz(1.5259471) q[0];
x q[1];
rz(1.7725905) q[2];
sx q[2];
rz(-1.861146) q[2];
sx q[2];
rz(-1.5648016) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.625714) q[1];
sx q[1];
rz(-0.82059723) q[1];
sx q[1];
rz(-1.3681812) q[1];
x q[2];
rz(-1.4895053) q[3];
sx q[3];
rz(-2.7346161) q[3];
sx q[3];
rz(-1.5291884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.011977) q[2];
sx q[2];
rz(-2.7836383) q[2];
sx q[2];
rz(0.32148662) q[2];
rz(-2.0251515) q[3];
sx q[3];
rz(-2.2827086) q[3];
sx q[3];
rz(1.4139676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9880992) q[0];
sx q[0];
rz(-1.7970002) q[0];
sx q[0];
rz(3.1186812) q[0];
rz(1.5921536) q[1];
sx q[1];
rz(-1.0433082) q[1];
sx q[1];
rz(1.9291838) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1003636) q[0];
sx q[0];
rz(-1.0335021) q[0];
sx q[0];
rz(2.8365718) q[0];
x q[1];
rz(-2.7801082) q[2];
sx q[2];
rz(-2.2687842) q[2];
sx q[2];
rz(2.7113376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7914932) q[1];
sx q[1];
rz(-1.9684569) q[1];
sx q[1];
rz(1.8997764) q[1];
rz(-1.4597662) q[3];
sx q[3];
rz(-2.4250406) q[3];
sx q[3];
rz(2.6727303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2108078) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(-0.99736324) q[2];
rz(-1.8831683) q[3];
sx q[3];
rz(-1.1910028) q[3];
sx q[3];
rz(1.2303801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9835984) q[0];
sx q[0];
rz(-0.65584922) q[0];
sx q[0];
rz(0.47437814) q[0];
rz(2.5899218) q[1];
sx q[1];
rz(-0.76849476) q[1];
sx q[1];
rz(-2.4840568) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0484587) q[0];
sx q[0];
rz(-3.0952929) q[0];
sx q[0];
rz(-2.6048624) q[0];
x q[1];
rz(-0.16993292) q[2];
sx q[2];
rz(-1.9714154) q[2];
sx q[2];
rz(-2.0037494) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.20831693) q[1];
sx q[1];
rz(-0.36872702) q[1];
sx q[1];
rz(1.5255724) q[1];
rz(-0.10513427) q[3];
sx q[3];
rz(-1.7145559) q[3];
sx q[3];
rz(-1.2105699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23003301) q[2];
sx q[2];
rz(-0.40731373) q[2];
sx q[2];
rz(-1.8302906) q[2];
rz(0.71410549) q[3];
sx q[3];
rz(-1.8592535) q[3];
sx q[3];
rz(0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0062200935) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(-2.5332992) q[0];
rz(0.81781203) q[1];
sx q[1];
rz(-1.9655971) q[1];
sx q[1];
rz(-0.83293319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9271256) q[0];
sx q[0];
rz(-1.9916849) q[0];
sx q[0];
rz(1.5141928) q[0];
rz(-pi) q[1];
rz(-0.3141381) q[2];
sx q[2];
rz(-1.3262265) q[2];
sx q[2];
rz(2.7121674) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91580931) q[1];
sx q[1];
rz(-1.7034917) q[1];
sx q[1];
rz(-0.58376329) q[1];
rz(-1.6021483) q[3];
sx q[3];
rz(-2.1636181) q[3];
sx q[3];
rz(1.5167936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9587162) q[2];
sx q[2];
rz(-0.35280886) q[2];
sx q[2];
rz(-2.555441) q[2];
rz(-0.90306774) q[3];
sx q[3];
rz(-0.62246263) q[3];
sx q[3];
rz(-3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027973) q[0];
sx q[0];
rz(-1.1302523) q[0];
sx q[0];
rz(-0.94373066) q[0];
rz(-2.0441652) q[1];
sx q[1];
rz(-0.25249093) q[1];
sx q[1];
rz(2.9846334) q[1];
rz(1.1949933) q[2];
sx q[2];
rz(-2.7766418) q[2];
sx q[2];
rz(3.1212213) q[2];
rz(-2.2678866) q[3];
sx q[3];
rz(-1.9601964) q[3];
sx q[3];
rz(-1.8662069) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
