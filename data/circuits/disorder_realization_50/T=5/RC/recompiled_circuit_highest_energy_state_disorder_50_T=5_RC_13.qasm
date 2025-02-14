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
rz(0.69831508) q[0];
sx q[0];
rz(-0.3787711) q[0];
sx q[0];
rz(1.1573428) q[0];
rz(0.78504374) q[1];
sx q[1];
rz(-0.68612376) q[1];
sx q[1];
rz(2.6148028) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9357932) q[0];
sx q[0];
rz(-0.65202921) q[0];
sx q[0];
rz(1.0656641) q[0];
rz(-pi) q[1];
rz(0.85028591) q[2];
sx q[2];
rz(-0.97480259) q[2];
sx q[2];
rz(-1.13499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59299027) q[1];
sx q[1];
rz(-1.1287264) q[1];
sx q[1];
rz(-1.5304687) q[1];
x q[2];
rz(0.029954682) q[3];
sx q[3];
rz(-2.1308793) q[3];
sx q[3];
rz(-0.64489472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.779125) q[2];
sx q[2];
rz(-0.96131009) q[2];
sx q[2];
rz(0.15164068) q[2];
rz(-2.9996297) q[3];
sx q[3];
rz(-1.4700593) q[3];
sx q[3];
rz(-2.0040472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66459429) q[0];
sx q[0];
rz(-0.73200309) q[0];
sx q[0];
rz(-0.81749302) q[0];
rz(0.11257653) q[1];
sx q[1];
rz(-1.6855626) q[1];
sx q[1];
rz(0.89961019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647659) q[0];
sx q[0];
rz(-1.810084) q[0];
sx q[0];
rz(-0.28061687) q[0];
x q[1];
rz(0.67018761) q[2];
sx q[2];
rz(-2.0815583) q[2];
sx q[2];
rz(-2.0005884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1020722) q[1];
sx q[1];
rz(-1.6270492) q[1];
sx q[1];
rz(-0.7489167) q[1];
rz(-0.51606222) q[3];
sx q[3];
rz(-2.1806697) q[3];
sx q[3];
rz(-0.42760951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4704935) q[2];
sx q[2];
rz(-1.6464536) q[2];
sx q[2];
rz(-2.1136843) q[2];
rz(-2.4043064) q[3];
sx q[3];
rz(-1.6643915) q[3];
sx q[3];
rz(-0.024287311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0055493) q[0];
sx q[0];
rz(-2.5856954) q[0];
sx q[0];
rz(-1.1935724) q[0];
rz(1.8434803) q[1];
sx q[1];
rz(-1.6622512) q[1];
sx q[1];
rz(1.6568291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62143059) q[0];
sx q[0];
rz(-1.545816) q[0];
sx q[0];
rz(1.6170184) q[0];
x q[1];
rz(1.7822927) q[2];
sx q[2];
rz(-1.4471608) q[2];
sx q[2];
rz(0.57963003) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3419108) q[1];
sx q[1];
rz(-1.7060301) q[1];
sx q[1];
rz(-2.7970201) q[1];
x q[2];
rz(1.3562702) q[3];
sx q[3];
rz(-0.37060043) q[3];
sx q[3];
rz(-1.724898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16126157) q[2];
sx q[2];
rz(-1.9183466) q[2];
sx q[2];
rz(0.70183357) q[2];
rz(-0.069132239) q[3];
sx q[3];
rz(-2.2528503) q[3];
sx q[3];
rz(-2.4338636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0858916) q[0];
sx q[0];
rz(-1.1393071) q[0];
sx q[0];
rz(0.62537801) q[0];
rz(1.0944132) q[1];
sx q[1];
rz(-1.5233327) q[1];
sx q[1];
rz(1.3166924) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4863281) q[0];
sx q[0];
rz(-1.0287971) q[0];
sx q[0];
rz(0.16714759) q[0];
rz(-pi) q[1];
rz(1.473963) q[2];
sx q[2];
rz(-1.8928573) q[2];
sx q[2];
rz(-1.2351954) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.229828) q[1];
sx q[1];
rz(-2.266464) q[1];
sx q[1];
rz(-1.8063481) q[1];
rz(-pi) q[2];
rz(-2.1573477) q[3];
sx q[3];
rz(-2.2797425) q[3];
sx q[3];
rz(-1.9714485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40654287) q[2];
sx q[2];
rz(-0.64261618) q[2];
sx q[2];
rz(0.030108062) q[2];
rz(-2.7249469) q[3];
sx q[3];
rz(-1.7130339) q[3];
sx q[3];
rz(0.055195181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77867126) q[0];
sx q[0];
rz(-1.8632977) q[0];
sx q[0];
rz(-1.9238506) q[0];
rz(0.22449224) q[1];
sx q[1];
rz(-1.4935363) q[1];
sx q[1];
rz(2.7405558) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7571791) q[0];
sx q[0];
rz(-1.1553488) q[0];
sx q[0];
rz(-2.7336304) q[0];
rz(-pi) q[1];
rz(-1.941576) q[2];
sx q[2];
rz(-0.94964281) q[2];
sx q[2];
rz(-0.3467243) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48408357) q[1];
sx q[1];
rz(-0.24166688) q[1];
sx q[1];
rz(-2.9422791) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55373904) q[3];
sx q[3];
rz(-0.19036346) q[3];
sx q[3];
rz(-2.7680264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.84247056) q[2];
sx q[2];
rz(-1.2229908) q[2];
sx q[2];
rz(-2.6341338) q[2];
rz(-0.92042813) q[3];
sx q[3];
rz(-0.43693742) q[3];
sx q[3];
rz(2.9612655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64666635) q[0];
sx q[0];
rz(-3.086402) q[0];
sx q[0];
rz(-2.1221509) q[0];
rz(1.1722209) q[1];
sx q[1];
rz(-1.1727138) q[1];
sx q[1];
rz(-1.2219465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2766314) q[0];
sx q[0];
rz(-0.91827938) q[0];
sx q[0];
rz(2.6990273) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5084294) q[2];
sx q[2];
rz(-0.51566873) q[2];
sx q[2];
rz(0.0002278681) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.032455803) q[1];
sx q[1];
rz(-1.0278406) q[1];
sx q[1];
rz(-1.2486267) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0133408) q[3];
sx q[3];
rz(-0.62164111) q[3];
sx q[3];
rz(0.78237247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4814066) q[2];
sx q[2];
rz(-2.47611) q[2];
sx q[2];
rz(1.0931724) q[2];
rz(-0.53705755) q[3];
sx q[3];
rz(-1.6780746) q[3];
sx q[3];
rz(-2.4721036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2375803) q[0];
sx q[0];
rz(-2.8222988) q[0];
sx q[0];
rz(-1.681666) q[0];
rz(1.5096674) q[1];
sx q[1];
rz(-2.5131112) q[1];
sx q[1];
rz(3.0622838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597042) q[0];
sx q[0];
rz(-1.4900582) q[0];
sx q[0];
rz(-1.8084007) q[0];
rz(0.58987381) q[2];
sx q[2];
rz(-0.24352077) q[2];
sx q[2];
rz(2.4517454) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39187688) q[1];
sx q[1];
rz(-1.2971216) q[1];
sx q[1];
rz(0.43507149) q[1];
rz(-pi) q[2];
rz(1.4344057) q[3];
sx q[3];
rz(-0.85771424) q[3];
sx q[3];
rz(0.73120414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0038393) q[2];
sx q[2];
rz(-1.9766108) q[2];
sx q[2];
rz(-0.4129146) q[2];
rz(1.2952992) q[3];
sx q[3];
rz(-2.2173939) q[3];
sx q[3];
rz(0.41954654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15692784) q[0];
sx q[0];
rz(-0.65021896) q[0];
sx q[0];
rz(-0.28537634) q[0];
rz(-3.0889619) q[1];
sx q[1];
rz(-1.4080518) q[1];
sx q[1];
rz(-2.9579128) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.015332) q[0];
sx q[0];
rz(-1.53526) q[0];
sx q[0];
rz(1.6823586) q[0];
x q[1];
rz(0.48641522) q[2];
sx q[2];
rz(-0.32054311) q[2];
sx q[2];
rz(-2.965791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.86183) q[1];
sx q[1];
rz(-2.1643736) q[1];
sx q[1];
rz(-1.4788126) q[1];
x q[2];
rz(1.2792789) q[3];
sx q[3];
rz(-0.94323483) q[3];
sx q[3];
rz(1.732839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24667428) q[2];
sx q[2];
rz(-1.7343212) q[2];
sx q[2];
rz(2.0764009) q[2];
rz(1.4554321) q[3];
sx q[3];
rz(-1.287241) q[3];
sx q[3];
rz(-0.58568946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10373779) q[0];
sx q[0];
rz(-2.3268564) q[0];
sx q[0];
rz(1.6023585) q[0];
rz(2.1922951) q[1];
sx q[1];
rz(-2.0930591) q[1];
sx q[1];
rz(-1.7112214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8447782) q[0];
sx q[0];
rz(-2.0955293) q[0];
sx q[0];
rz(2.120976) q[0];
rz(-0.075242356) q[2];
sx q[2];
rz(-2.2981567) q[2];
sx q[2];
rz(2.5204349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0956968) q[1];
sx q[1];
rz(-1.2751082) q[1];
sx q[1];
rz(-1.5398094) q[1];
rz(-pi) q[2];
rz(1.999302) q[3];
sx q[3];
rz(-1.3125889) q[3];
sx q[3];
rz(2.9203897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1324233) q[2];
sx q[2];
rz(-1.3045661) q[2];
sx q[2];
rz(-3.0150748) q[2];
rz(-0.33454076) q[3];
sx q[3];
rz(-2.8821475) q[3];
sx q[3];
rz(2.2693966) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8406521) q[0];
sx q[0];
rz(-2.2569188) q[0];
sx q[0];
rz(2.6244923) q[0];
rz(0.05489796) q[1];
sx q[1];
rz(-1.5093191) q[1];
sx q[1];
rz(3.0518234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9835998) q[0];
sx q[0];
rz(-1.9741892) q[0];
sx q[0];
rz(-0.97121111) q[0];
rz(-pi) q[1];
rz(0.062762063) q[2];
sx q[2];
rz(-1.0107991) q[2];
sx q[2];
rz(-2.7095209) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.43347142) q[1];
sx q[1];
rz(-2.8420487) q[1];
sx q[1];
rz(0.56510651) q[1];
rz(-1.9338984) q[3];
sx q[3];
rz(-1.3570367) q[3];
sx q[3];
rz(-0.0068706415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6101997) q[2];
sx q[2];
rz(-1.1343845) q[2];
sx q[2];
rz(-0.72719491) q[2];
rz(2.3860892) q[3];
sx q[3];
rz(-0.38755363) q[3];
sx q[3];
rz(-2.9288647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1995734) q[0];
sx q[0];
rz(-1.5409536) q[0];
sx q[0];
rz(1.090747) q[0];
rz(1.749281) q[1];
sx q[1];
rz(-1.5131469) q[1];
sx q[1];
rz(-1.6319235) q[1];
rz(2.2752599) q[2];
sx q[2];
rz(-1.489218) q[2];
sx q[2];
rz(1.0897286) q[2];
rz(2.5531664) q[3];
sx q[3];
rz(-2.1156433) q[3];
sx q[3];
rz(2.0414203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
