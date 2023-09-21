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
rz(4.6580553) q[0];
sx q[0];
rz(9.1604995) q[0];
rz(-0.9737941) q[1];
sx q[1];
rz(5.073054) q[1];
sx q[1];
rz(10.160025) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0300232) q[0];
sx q[0];
rz(-2.2246247) q[0];
sx q[0];
rz(1.3841188) q[0];
rz(-pi) q[1];
rz(-2.4618857) q[2];
sx q[2];
rz(-1.1628816) q[2];
sx q[2];
rz(2.9349875) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.211261) q[1];
sx q[1];
rz(-2.8144565) q[1];
sx q[1];
rz(-1.0717908) q[1];
rz(-pi) q[2];
rz(0.04282184) q[3];
sx q[3];
rz(-2.560727) q[3];
sx q[3];
rz(0.37561852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66951093) q[2];
sx q[2];
rz(-1.8410204) q[2];
sx q[2];
rz(-1.1038587) q[2];
rz(-1.2708698) q[3];
sx q[3];
rz(-1.2277675) q[3];
sx q[3];
rz(0.27403533) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0274149) q[0];
sx q[0];
rz(-0.52863055) q[0];
sx q[0];
rz(-0.43637481) q[0];
rz(-0.46288681) q[1];
sx q[1];
rz(-2.1040237) q[1];
sx q[1];
rz(2.8754821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9929745) q[0];
sx q[0];
rz(-2.2668112) q[0];
sx q[0];
rz(0.50177411) q[0];
x q[1];
rz(1.3219464) q[2];
sx q[2];
rz(-0.58983931) q[2];
sx q[2];
rz(2.313254) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6978584) q[1];
sx q[1];
rz(-0.88159544) q[1];
sx q[1];
rz(-2.2976848) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3105884) q[3];
sx q[3];
rz(-2.5187413) q[3];
sx q[3];
rz(-2.7130896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6767072) q[2];
sx q[2];
rz(-1.8647944) q[2];
sx q[2];
rz(-2.6300988) q[2];
rz(2.3320847) q[3];
sx q[3];
rz(-1.6102689) q[3];
sx q[3];
rz(2.8539343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4354316) q[0];
sx q[0];
rz(-1.5937188) q[0];
sx q[0];
rz(-0.92873746) q[0];
rz(1.4061032) q[1];
sx q[1];
rz(-0.7000674) q[1];
sx q[1];
rz(-1.3471289) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1613306) q[0];
sx q[0];
rz(-1.8322924) q[0];
sx q[0];
rz(-2.6677092) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0727097) q[2];
sx q[2];
rz(-0.32215298) q[2];
sx q[2];
rz(-1.0108394) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9007064) q[1];
sx q[1];
rz(-2.1623003) q[1];
sx q[1];
rz(-2.790931) q[1];
x q[2];
rz(0.8637572) q[3];
sx q[3];
rz(-1.45544) q[3];
sx q[3];
rz(0.74187169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9946263) q[2];
sx q[2];
rz(-1.7549843) q[2];
sx q[2];
rz(1.6195126) q[2];
rz(2.8772723) q[3];
sx q[3];
rz(-2.1051354) q[3];
sx q[3];
rz(-2.6456397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79214823) q[0];
sx q[0];
rz(-1.148372) q[0];
sx q[0];
rz(-2.175892) q[0];
rz(-0.72215885) q[1];
sx q[1];
rz(-1.637371) q[1];
sx q[1];
rz(2.5818363) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352919) q[0];
sx q[0];
rz(-1.9959404) q[0];
sx q[0];
rz(-1.6598808) q[0];
rz(-pi) q[1];
rz(-0.70978769) q[2];
sx q[2];
rz(-2.1927532) q[2];
sx q[2];
rz(0.70993916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8430082) q[1];
sx q[1];
rz(-0.95579879) q[1];
sx q[1];
rz(3.0857012) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78919952) q[3];
sx q[3];
rz(-1.7390828) q[3];
sx q[3];
rz(-2.5326953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7136148) q[2];
sx q[2];
rz(-1.5087936) q[2];
sx q[2];
rz(1.7170061) q[2];
rz(2.8811841) q[3];
sx q[3];
rz(-1.3924761) q[3];
sx q[3];
rz(2.7105455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.150862) q[0];
sx q[0];
rz(-1.4980415) q[0];
sx q[0];
rz(2.7752303) q[0];
rz(-1.5953966) q[1];
sx q[1];
rz(-2.5876744) q[1];
sx q[1];
rz(0.34367925) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2491964) q[0];
sx q[0];
rz(-2.2049892) q[0];
sx q[0];
rz(-0.87974324) q[0];
rz(-pi) q[1];
rz(2.2291549) q[2];
sx q[2];
rz(-0.81805938) q[2];
sx q[2];
rz(-0.91615265) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74663631) q[1];
sx q[1];
rz(-0.50368217) q[1];
sx q[1];
rz(-1.7240745) q[1];
rz(-pi) q[2];
rz(-1.5424535) q[3];
sx q[3];
rz(-0.64405555) q[3];
sx q[3];
rz(-1.6213662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7498103) q[2];
sx q[2];
rz(-2.2822773) q[2];
sx q[2];
rz(-3.0878477) q[2];
rz(-1.404445) q[3];
sx q[3];
rz(-0.497538) q[3];
sx q[3];
rz(2.8500309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4869726) q[0];
sx q[0];
rz(-0.64990652) q[0];
sx q[0];
rz(-2.3858331) q[0];
rz(-3.1164363) q[1];
sx q[1];
rz(-0.92725602) q[1];
sx q[1];
rz(0.25973928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4070968) q[0];
sx q[0];
rz(-1.2571113) q[0];
sx q[0];
rz(-1.9514027) q[0];
rz(0.01979205) q[2];
sx q[2];
rz(-1.2345825) q[2];
sx q[2];
rz(-2.3078231) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0670843) q[1];
sx q[1];
rz(-1.2423406) q[1];
sx q[1];
rz(2.902608) q[1];
rz(-pi) q[2];
x q[2];
rz(1.608058) q[3];
sx q[3];
rz(-0.44964368) q[3];
sx q[3];
rz(0.26526181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6529237) q[2];
sx q[2];
rz(-2.3355464) q[2];
sx q[2];
rz(-0.10061131) q[2];
rz(-2.9595024) q[3];
sx q[3];
rz(-0.87825769) q[3];
sx q[3];
rz(-1.3476868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9942193) q[0];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46442407) q[0];
sx q[0];
rz(-0.66093984) q[0];
sx q[0];
rz(-0.37207694) q[0];
rz(2.546054) q[2];
sx q[2];
rz(-0.34791246) q[2];
sx q[2];
rz(-2.4728342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0026605) q[1];
sx q[1];
rz(-1.9117038) q[1];
sx q[1];
rz(-2.1041811) q[1];
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
rz(-2.9014152) q[2];
sx q[2];
rz(-1.9577953) q[2];
sx q[2];
rz(0.85285464) q[2];
rz(-1.7715706) q[3];
sx q[3];
rz(-1.4586689) q[3];
sx q[3];
rz(-0.20496932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2296427) q[0];
sx q[0];
rz(-0.59597534) q[0];
sx q[0];
rz(-1.6802616) q[0];
rz(1.7386859) q[1];
sx q[1];
rz(-2.1673514) q[1];
sx q[1];
rz(3.0775552) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1916699) q[0];
sx q[0];
rz(-2.4987698) q[0];
sx q[0];
rz(-1.2547917) q[0];
x q[1];
rz(1.0649101) q[2];
sx q[2];
rz(-1.3008899) q[2];
sx q[2];
rz(-0.050886521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0361955) q[1];
sx q[1];
rz(-2.7846249) q[1];
sx q[1];
rz(3.0371975) q[1];
rz(1.6927034) q[3];
sx q[3];
rz(-0.31414437) q[3];
sx q[3];
rz(-1.4172518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.091207592) q[2];
sx q[2];
rz(-0.63082266) q[2];
sx q[2];
rz(1.5861661) q[2];
rz(2.2533916) q[3];
sx q[3];
rz(-1.2398088) q[3];
sx q[3];
rz(2.2122038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14426194) q[0];
sx q[0];
rz(-2.0781131) q[0];
sx q[0];
rz(-1.4319179) q[0];
rz(-0.56888467) q[1];
sx q[1];
rz(-2.6061997) q[1];
sx q[1];
rz(2.0137537) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2375243) q[0];
sx q[0];
rz(-1.5547817) q[0];
sx q[0];
rz(-2.9928656) q[0];
rz(-pi) q[1];
rz(2.1986507) q[2];
sx q[2];
rz(-2.0161511) q[2];
sx q[2];
rz(-3.359059e-05) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7334426) q[1];
sx q[1];
rz(-0.24194716) q[1];
sx q[1];
rz(-1.5223632) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7752152) q[3];
sx q[3];
rz(-1.0421703) q[3];
sx q[3];
rz(-2.5322994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.52788064) q[2];
sx q[2];
rz(-2.2884559) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7062475) q[0];
sx q[0];
rz(-0.58723891) q[0];
sx q[0];
rz(1.6760814) q[0];
rz(2.3174875) q[1];
sx q[1];
rz(-1.6128287) q[1];
sx q[1];
rz(-0.5724268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7596282) q[0];
sx q[0];
rz(-1.0257162) q[0];
sx q[0];
rz(-0.99960021) q[0];
rz(-2.8675251) q[2];
sx q[2];
rz(-1.1268508) q[2];
sx q[2];
rz(0.62192164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4116284) q[1];
sx q[1];
rz(-2.8169544) q[1];
sx q[1];
rz(2.4050729) q[1];
rz(1.8626067) q[3];
sx q[3];
rz(-0.88340532) q[3];
sx q[3];
rz(0.22112267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3616025) q[2];
sx q[2];
rz(-2.6670691) q[2];
sx q[2];
rz(-2.6043716) q[2];
rz(-2.0843263) q[3];
sx q[3];
rz(-2.2500762) q[3];
sx q[3];
rz(0.69361544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2873516) q[0];
sx q[0];
rz(-1.9762522) q[0];
sx q[0];
rz(1.5594788) q[0];
rz(0.46335012) q[1];
sx q[1];
rz(-2.2644823) q[1];
sx q[1];
rz(1.5092441) q[1];
rz(1.6981381) q[2];
sx q[2];
rz(-0.6586532) q[2];
sx q[2];
rz(-2.3011343) q[2];
rz(2.0588856) q[3];
sx q[3];
rz(-1.516021) q[3];
sx q[3];
rz(1.0513023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];