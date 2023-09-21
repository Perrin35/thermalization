OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5207198) q[0];
sx q[0];
rz(-1.7680661) q[0];
sx q[0];
rz(-1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(2.642282) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141657) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(2.8660197) q[0];
x q[1];
rz(2.9688641) q[2];
sx q[2];
rz(-2.2842801) q[2];
sx q[2];
rz(1.2781065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5159113) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(0.73966571) q[1];
rz(-pi) q[2];
rz(1.4673759) q[3];
sx q[3];
rz(-1.5973583) q[3];
sx q[3];
rz(-0.51277044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(2.1271558) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4734128) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53045814) q[0];
sx q[0];
rz(-2.3728752) q[0];
sx q[0];
rz(2.4918633) q[0];
x q[1];
rz(2.9687256) q[2];
sx q[2];
rz(-1.7555408) q[2];
sx q[2];
rz(2.6181521) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4304632) q[1];
sx q[1];
rz(-1.4158447) q[1];
sx q[1];
rz(2.1610545) q[1];
rz(-2.5428883) q[3];
sx q[3];
rz(-1.8667392) q[3];
sx q[3];
rz(0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(-3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012264) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(0.40584392) q[0];
x q[1];
rz(-1.1728889) q[2];
sx q[2];
rz(-1.298561) q[2];
sx q[2];
rz(2.9583601) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1556342) q[1];
sx q[1];
rz(-1.8556719) q[1];
sx q[1];
rz(-0.023521544) q[1];
rz(-pi) q[2];
rz(-1.6885353) q[3];
sx q[3];
rz(-2.2265834) q[3];
sx q[3];
rz(-1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3016004) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(-2.7123614) q[2];
rz(2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4317959) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(-2.5298932) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(2.591419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94565369) q[0];
sx q[0];
rz(-1.494207) q[0];
sx q[0];
rz(-2.4038195) q[0];
rz(2.8470464) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(-1.7812658) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78754567) q[1];
sx q[1];
rz(-2.5552632) q[1];
sx q[1];
rz(1.7802618) q[1];
rz(-1.6226193) q[3];
sx q[3];
rz(-1.341815) q[3];
sx q[3];
rz(2.1932972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-2.6553314) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(3.0299305) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(0.87669796) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(0.91526389) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6633776) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(0.94232725) q[0];
rz(1.5259597) q[2];
sx q[2];
rz(-1.0460639) q[2];
sx q[2];
rz(0.63477883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76449672) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(1.3678958) q[1];
rz(-pi) q[2];
rz(2.4004585) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(2.0625431) q[0];
rz(2.5462529) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(-0.39658305) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3298033) q[0];
sx q[0];
rz(-1.8237231) q[0];
sx q[0];
rz(1.1008218) q[0];
rz(-pi) q[1];
rz(-2.6655212) q[2];
sx q[2];
rz(-1.2928315) q[2];
sx q[2];
rz(1.6518041) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.086416883) q[1];
sx q[1];
rz(-2.2284818) q[1];
sx q[1];
rz(-3.0228826) q[1];
rz(-pi) q[2];
rz(-1.3271689) q[3];
sx q[3];
rz(-2.1414087) q[3];
sx q[3];
rz(-2.7057196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2369909) q[0];
sx q[0];
rz(-2.8303574) q[0];
sx q[0];
rz(-0.29445946) q[0];
rz(-pi) q[1];
rz(2.0790714) q[2];
sx q[2];
rz(-1.695343) q[2];
sx q[2];
rz(1.5073656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.78325242) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.3621484) q[1];
rz(0.095110006) q[3];
sx q[3];
rz(-1.2246119) q[3];
sx q[3];
rz(0.1971052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37041935) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2643124) q[0];
sx q[0];
rz(-2.7240629) q[0];
sx q[0];
rz(2.2420922) q[0];
x q[1];
rz(-2.6997386) q[2];
sx q[2];
rz(-1.0519069) q[2];
sx q[2];
rz(-1.2043124) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36383648) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(-1.9883518) q[1];
rz(1.1912187) q[3];
sx q[3];
rz(-1.4879585) q[3];
sx q[3];
rz(2.6759202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1774566) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(1.1671676) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.099982925) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(-1.6149678) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-2.4885978) q[1];
sx q[1];
rz(2.656235) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.582726) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(-0.38082122) q[0];
rz(-pi) q[1];
rz(-2.543407) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(-0.26192947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(1.1361213) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5076809) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(2.4661635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(0.1594485) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-3.1260417) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.7074701) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-2.1914296) q[1];
sx q[1];
rz(-0.5982582) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
rz(-pi) q[1];
rz(0.60600772) q[2];
sx q[2];
rz(-0.25642828) q[2];
sx q[2];
rz(2.3026349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4926589) q[1];
sx q[1];
rz(-1.6655386) q[1];
sx q[1];
rz(1.2050864) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9833097) q[3];
sx q[3];
rz(-1.4083107) q[3];
sx q[3];
rz(-0.76009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3020637) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-0.17157208) q[2];
sx q[2];
rz(-0.62660672) q[2];
sx q[2];
rz(-1.3852711) q[2];
rz(2.0486352) q[3];
sx q[3];
rz(-0.7328877) q[3];
sx q[3];
rz(2.4536798) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];