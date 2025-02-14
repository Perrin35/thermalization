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
rz(-3.0814085) q[0];
sx q[0];
rz(-1.0589851) q[0];
sx q[0];
rz(1.0101779) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(-0.29616907) q[1];
sx q[1];
rz(0.30997601) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.298807) q[0];
sx q[0];
rz(-1.9901384) q[0];
sx q[0];
rz(-1.4489277) q[0];
rz(-pi) q[1];
rz(1.5484137) q[2];
sx q[2];
rz(-1.3655919) q[2];
sx q[2];
rz(-0.6485282) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0204326) q[1];
sx q[1];
rz(-1.766664) q[1];
sx q[1];
rz(-1.2365667) q[1];
rz(-pi) q[2];
rz(-0.83155379) q[3];
sx q[3];
rz(-1.7208916) q[3];
sx q[3];
rz(0.41285601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4787204) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(-3.0459246) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-2.2662558) q[3];
sx q[3];
rz(-1.7101425) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2410759) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(2.526793) q[0];
rz(2.2531033) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(2.6420171) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7730656) q[0];
sx q[0];
rz(-1.8963666) q[0];
sx q[0];
rz(-0.21296176) q[0];
rz(-pi) q[1];
rz(-1.4128628) q[2];
sx q[2];
rz(-2.0896032) q[2];
sx q[2];
rz(0.11473303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37812472) q[1];
sx q[1];
rz(-0.92510447) q[1];
sx q[1];
rz(0.21685361) q[1];
rz(0.98827117) q[3];
sx q[3];
rz(-1.2856021) q[3];
sx q[3];
rz(-2.5799283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8932314) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(-0.020817967) q[2];
rz(-1.9013532) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(-1.9564691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5027387) q[0];
sx q[0];
rz(-0.26832142) q[0];
sx q[0];
rz(0.33367208) q[0];
rz(-2.4036713) q[1];
sx q[1];
rz(-1.2260022) q[1];
sx q[1];
rz(3.0121682) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20724587) q[0];
sx q[0];
rz(-0.6967623) q[0];
sx q[0];
rz(1.8201226) q[0];
rz(-pi) q[1];
rz(2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(-2.3455623) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63505748) q[1];
sx q[1];
rz(-1.4806857) q[1];
sx q[1];
rz(1.4485301) q[1];
x q[2];
rz(0.78227104) q[3];
sx q[3];
rz(-1.2201628) q[3];
sx q[3];
rz(0.80814894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(1.370149) q[2];
rz(-0.74639368) q[3];
sx q[3];
rz(-1.3179702) q[3];
sx q[3];
rz(-3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.2097825) q[0];
sx q[0];
rz(0.56513894) q[0];
rz(2.7124229) q[1];
sx q[1];
rz(-0.64650911) q[1];
sx q[1];
rz(-0.82829222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1278616) q[0];
sx q[0];
rz(-0.7452226) q[0];
sx q[0];
rz(3.0237314) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28675927) q[2];
sx q[2];
rz(-1.1017403) q[2];
sx q[2];
rz(-3.0931851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2678623) q[1];
sx q[1];
rz(-0.56247382) q[1];
sx q[1];
rz(-2.6270694) q[1];
x q[2];
rz(-0.24803646) q[3];
sx q[3];
rz(-1.0674849) q[3];
sx q[3];
rz(0.65500427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4565178) q[2];
sx q[2];
rz(-2.0802616) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(-0.00069869839) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(-1.8121207) q[0];
rz(1.9895408) q[1];
sx q[1];
rz(-1.0877129) q[1];
sx q[1];
rz(0.00024814127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084293289) q[0];
sx q[0];
rz(-1.3868185) q[0];
sx q[0];
rz(1.3217682) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0543112) q[2];
sx q[2];
rz(-1.0382663) q[2];
sx q[2];
rz(-2.0244348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8793638) q[1];
sx q[1];
rz(-1.2254224) q[1];
sx q[1];
rz(1.0443347) q[1];
x q[2];
rz(3.0609691) q[3];
sx q[3];
rz(-0.73314694) q[3];
sx q[3];
rz(2.6997363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0033215) q[2];
sx q[2];
rz(-1.891581) q[2];
sx q[2];
rz(3.1357583) q[2];
rz(-0.88998574) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(1.1267004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(-1.4877315) q[0];
rz(-0.030390175) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(0.94246513) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0069601428) q[0];
sx q[0];
rz(-2.83036) q[0];
sx q[0];
rz(2.4853021) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29877383) q[2];
sx q[2];
rz(-1.818383) q[2];
sx q[2];
rz(0.58591671) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9760175) q[1];
sx q[1];
rz(-1.2846795) q[1];
sx q[1];
rz(2.0598434) q[1];
rz(2.8811137) q[3];
sx q[3];
rz(-1.1090058) q[3];
sx q[3];
rz(-0.33212979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5421062) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(1.3055118) q[2];
rz(0.54715884) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(0.80593306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18043537) q[0];
sx q[0];
rz(-1.0338217) q[0];
sx q[0];
rz(-3.0410774) q[0];
rz(-0.48249498) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(-2.2844792) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6109638) q[0];
sx q[0];
rz(-0.93994323) q[0];
sx q[0];
rz(1.1683589) q[0];
rz(-1.6202507) q[2];
sx q[2];
rz(-2.2170545) q[2];
sx q[2];
rz(-0.44909278) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3639314) q[1];
sx q[1];
rz(-0.45501935) q[1];
sx q[1];
rz(-2.6520686) q[1];
rz(-1.1825652) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(-1.5608112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40019217) q[2];
sx q[2];
rz(-2.1705748) q[2];
sx q[2];
rz(2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-1.0023578) q[3];
sx q[3];
rz(-1.8625331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.361146) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(-3.1106023) q[0];
rz(-0.57394761) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.8642289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4784067) q[0];
sx q[0];
rz(-0.48106964) q[0];
sx q[0];
rz(2.3548954) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9979111) q[2];
sx q[2];
rz(-0.53887212) q[2];
sx q[2];
rz(-0.32921916) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0884235) q[1];
sx q[1];
rz(-2.649247) q[1];
sx q[1];
rz(1.4150934) q[1];
rz(1.8559009) q[3];
sx q[3];
rz(-1.5114956) q[3];
sx q[3];
rz(-1.3467195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0370827) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(0.48745298) q[2];
rz(2.4449352) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(-0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071851991) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(2.790614) q[0];
rz(-2.9674496) q[1];
sx q[1];
rz(-1.5085647) q[1];
sx q[1];
rz(-1.7399656) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202209) q[0];
sx q[0];
rz(-1.3231771) q[0];
sx q[0];
rz(-1.4445452) q[0];
rz(-pi) q[1];
rz(0.75058896) q[2];
sx q[2];
rz(-1.9074512) q[2];
sx q[2];
rz(2.2845993) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40608866) q[1];
sx q[1];
rz(-1.5417409) q[1];
sx q[1];
rz(-1.3872223) q[1];
rz(-pi) q[2];
rz(2.4003567) q[3];
sx q[3];
rz(-1.110807) q[3];
sx q[3];
rz(2.7335087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.031124) q[2];
sx q[2];
rz(-0.68507552) q[2];
sx q[2];
rz(-2.5049211) q[2];
rz(1.1104256) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(3.1029347) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-2.0775332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.80262) q[0];
sx q[0];
rz(-1.5674599) q[0];
sx q[0];
rz(-0.0044857684) q[0];
x q[1];
rz(0.27711192) q[2];
sx q[2];
rz(-1.2925576) q[2];
sx q[2];
rz(-0.81467512) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0088996) q[1];
sx q[1];
rz(-1.6311967) q[1];
sx q[1];
rz(-1.6113043) q[1];
rz(-pi) q[2];
rz(0.070399447) q[3];
sx q[3];
rz(-2.681085) q[3];
sx q[3];
rz(-1.2706336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(0.074020298) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114125) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(-0.61182712) q[1];
sx q[1];
rz(-0.71129967) q[1];
sx q[1];
rz(1.5536972) q[1];
rz(1.189497) q[2];
sx q[2];
rz(-1.0972037) q[2];
sx q[2];
rz(0.9158132) q[2];
rz(-1.0985804) q[3];
sx q[3];
rz(-2.2618812) q[3];
sx q[3];
rz(-1.19899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
