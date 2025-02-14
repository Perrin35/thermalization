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
rz(-2.4432776) q[0];
sx q[0];
rz(-2.7628216) q[0];
sx q[0];
rz(1.9842499) q[0];
rz(10.209822) q[1];
sx q[1];
rz(0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3279545) q[0];
sx q[0];
rz(-2.1305973) q[0];
sx q[0];
rz(0.35388057) q[0];
rz(-pi) q[1];
rz(-0.73424642) q[2];
sx q[2];
rz(-0.99319211) q[2];
sx q[2];
rz(2.2478888) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68702811) q[1];
sx q[1];
rz(-0.44378456) q[1];
sx q[1];
rz(0.084974809) q[1];
rz(-1.6185226) q[3];
sx q[3];
rz(-2.5807947) q[3];
sx q[3];
rz(-0.58855295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36246768) q[2];
sx q[2];
rz(-2.1802826) q[2];
sx q[2];
rz(-0.15164068) q[2];
rz(-2.9996297) q[3];
sx q[3];
rz(-1.6715334) q[3];
sx q[3];
rz(-1.1375455) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66459429) q[0];
sx q[0];
rz(-2.4095896) q[0];
sx q[0];
rz(0.81749302) q[0];
rz(-0.11257653) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(-0.89961019) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7821723) q[0];
sx q[0];
rz(-0.36672584) q[0];
sx q[0];
rz(0.72215778) q[0];
rz(-2.4075872) q[2];
sx q[2];
rz(-0.81799928) q[2];
sx q[2];
rz(0.98246511) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4790597) q[1];
sx q[1];
rz(-2.3182437) q[1];
sx q[1];
rz(1.4940628) q[1];
x q[2];
rz(-0.894015) q[3];
sx q[3];
rz(-1.1543659) q[3];
sx q[3];
rz(2.3125966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.67109913) q[2];
sx q[2];
rz(-1.6464536) q[2];
sx q[2];
rz(-1.0279083) q[2];
rz(-0.73728621) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(-3.1173053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0055493) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.9480202) q[0];
rz(-1.8434803) q[1];
sx q[1];
rz(-1.4793414) q[1];
sx q[1];
rz(-1.4847635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5201621) q[0];
sx q[0];
rz(-1.545816) q[0];
sx q[0];
rz(-1.6170184) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0151691) q[2];
sx q[2];
rz(-1.3609386) q[2];
sx q[2];
rz(-2.1239547) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.130467) q[1];
sx q[1];
rz(-0.36917205) q[1];
sx q[1];
rz(-2.7586731) q[1];
x q[2];
rz(-1.3562702) q[3];
sx q[3];
rz(-2.7709922) q[3];
sx q[3];
rz(1.4166946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.16126157) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(-0.70183357) q[2];
rz(3.0724604) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(0.70772901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0858916) q[0];
sx q[0];
rz(-2.0022855) q[0];
sx q[0];
rz(-2.5162146) q[0];
rz(-2.0471795) q[1];
sx q[1];
rz(-1.6182599) q[1];
sx q[1];
rz(1.8249003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8024361) q[0];
sx q[0];
rz(-0.56471741) q[0];
sx q[0];
rz(1.8403649) q[0];
rz(-pi) q[1];
rz(2.8595692) q[2];
sx q[2];
rz(-2.8057753) q[2];
sx q[2];
rz(1.5329602) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91176468) q[1];
sx q[1];
rz(-2.266464) q[1];
sx q[1];
rz(-1.8063481) q[1];
rz(-pi) q[2];
rz(-2.5685124) q[3];
sx q[3];
rz(-2.2551558) q[3];
sx q[3];
rz(-2.7671368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40654287) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(3.1114846) q[2];
rz(0.41664577) q[3];
sx q[3];
rz(-1.4285587) q[3];
sx q[3];
rz(3.0863975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3629214) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(1.2177421) q[0];
rz(-0.22449224) q[1];
sx q[1];
rz(-1.6480564) q[1];
sx q[1];
rz(2.7405558) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5767956) q[0];
sx q[0];
rz(-0.57379299) q[0];
sx q[0];
rz(-2.3032715) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46868344) q[2];
sx q[2];
rz(-2.4309553) q[2];
sx q[2];
rz(0.24220322) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6575091) q[1];
sx q[1];
rz(-0.24166688) q[1];
sx q[1];
rz(2.9422791) q[1];
rz(0.16245654) q[3];
sx q[3];
rz(-1.6704644) q[3];
sx q[3];
rz(-1.398744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2991221) q[2];
sx q[2];
rz(-1.2229908) q[2];
sx q[2];
rz(0.50745884) q[2];
rz(-2.2211645) q[3];
sx q[3];
rz(-2.7046552) q[3];
sx q[3];
rz(-0.18032716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4949263) q[0];
sx q[0];
rz(-0.05519069) q[0];
sx q[0];
rz(-2.1221509) q[0];
rz(-1.9693718) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(1.9196462) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9393765) q[0];
sx q[0];
rz(-0.76986137) q[0];
sx q[0];
rz(-1.059993) q[0];
x q[1];
rz(-3.1062791) q[2];
sx q[2];
rz(-1.0562293) q[2];
sx q[2];
rz(-3.0701612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1091369) q[1];
sx q[1];
rz(-2.1137521) q[1];
sx q[1];
rz(1.892966) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99626793) q[3];
sx q[3];
rz(-1.822851) q[3];
sx q[3];
rz(-1.1561405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6601861) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(2.0484203) q[2];
rz(2.6045351) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(0.66948906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2375803) q[0];
sx q[0];
rz(-0.31929382) q[0];
sx q[0];
rz(-1.4599266) q[0];
rz(1.6319252) q[1];
sx q[1];
rz(-0.62848148) q[1];
sx q[1];
rz(-0.07930886) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58188841) q[0];
sx q[0];
rz(-1.6515344) q[0];
sx q[0];
rz(-1.8084007) q[0];
rz(0.20360501) q[2];
sx q[2];
rz(-1.4362659) q[2];
sx q[2];
rz(-2.8367538) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7056337) q[1];
sx q[1];
rz(-2.6323032) q[1];
sx q[1];
rz(-0.58756709) q[1];
rz(-2.9856624) q[3];
sx q[3];
rz(-2.4178388) q[3];
sx q[3];
rz(2.2035905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1377533) q[2];
sx q[2];
rz(-1.1649818) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(-1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(2.7220461) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15692784) q[0];
sx q[0];
rz(-0.65021896) q[0];
sx q[0];
rz(-2.8562163) q[0];
rz(-0.05263075) q[1];
sx q[1];
rz(-1.7335408) q[1];
sx q[1];
rz(0.1836798) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4405555) q[0];
sx q[0];
rz(-1.6822878) q[0];
sx q[0];
rz(-0.035758408) q[0];
x q[1];
rz(-0.48641522) q[2];
sx q[2];
rz(-0.32054311) q[2];
sx q[2];
rz(2.965791) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3425828) q[1];
sx q[1];
rz(-1.6470121) q[1];
sx q[1];
rz(-2.5460495) q[1];
rz(1.8623137) q[3];
sx q[3];
rz(-2.1983578) q[3];
sx q[3];
rz(-1.4087536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8949184) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(-1.0651917) q[2];
rz(-1.6861606) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(-0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10373779) q[0];
sx q[0];
rz(-0.81473628) q[0];
sx q[0];
rz(1.6023585) q[0];
rz(-2.1922951) q[1];
sx q[1];
rz(-1.0485336) q[1];
sx q[1];
rz(-1.7112214) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024139013) q[0];
sx q[0];
rz(-1.1012337) q[0];
sx q[0];
rz(0.59654327) q[0];
x q[1];
rz(-0.075242356) q[2];
sx q[2];
rz(-2.2981567) q[2];
sx q[2];
rz(2.5204349) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51586823) q[1];
sx q[1];
rz(-1.6004381) q[1];
sx q[1];
rz(2.8457706) q[1];
rz(-pi) q[2];
rz(-1.999302) q[3];
sx q[3];
rz(-1.8290038) q[3];
sx q[3];
rz(2.9203897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1324233) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(0.12651786) q[2];
rz(-0.33454076) q[3];
sx q[3];
rz(-0.25944513) q[3];
sx q[3];
rz(0.87219605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3009406) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(-0.51710039) q[0];
rz(-0.05489796) q[1];
sx q[1];
rz(-1.6322735) q[1];
sx q[1];
rz(-0.089769207) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15799284) q[0];
sx q[0];
rz(-1.9741892) q[0];
sx q[0];
rz(2.1703815) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1316809) q[2];
sx q[2];
rz(-1.5176306) q[2];
sx q[2];
rz(1.9694984) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5489784) q[1];
sx q[1];
rz(-1.7294809) q[1];
sx q[1];
rz(2.8864546) q[1];
rz(0.2281727) q[3];
sx q[3];
rz(-1.2163278) q[3];
sx q[3];
rz(1.6443524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.53139293) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(-0.72719491) q[2];
rz(-0.7555035) q[3];
sx q[3];
rz(-0.38755363) q[3];
sx q[3];
rz(0.21272794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1995734) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(1.3923116) q[1];
sx q[1];
rz(-1.6284457) q[1];
sx q[1];
rz(1.5096691) q[1];
rz(0.86633273) q[2];
sx q[2];
rz(-1.6523747) q[2];
sx q[2];
rz(-2.0518641) q[2];
rz(2.200442) q[3];
sx q[3];
rz(-2.0654021) q[3];
sx q[3];
rz(0.13765814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
