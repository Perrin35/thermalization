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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141657) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(-0.27557296) q[0];
x q[1];
rz(-0.84988611) q[2];
sx q[2];
rz(-1.7011142) q[2];
sx q[2];
rz(2.7352114) q[2];
sx q[3];
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
rz(-1.4673759) q[3];
sx q[3];
rz(-1.5442344) q[3];
sx q[3];
rz(2.6288222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(-1.0144368) q[2];
rz(2.9075918) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(1.0043253) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(2.1388334) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5402055) q[0];
sx q[0];
rz(-2.0048855) q[0];
sx q[0];
rz(-0.65625221) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3833212) q[2];
sx q[2];
rz(-1.7406929) q[2];
sx q[2];
rz(2.1263009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36675378) q[1];
sx q[1];
rz(-2.5336821) q[1];
sx q[1];
rz(1.2971836) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59870436) q[3];
sx q[3];
rz(-1.2748534) q[3];
sx q[3];
rz(0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(-2.9906452) q[2];
rz(2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(3.0533561) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(2.3625968) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.2483968) q[1];
sx q[1];
rz(2.2580106) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2471836) q[0];
sx q[0];
rz(-2.7054225) q[0];
sx q[0];
rz(2.7437074) q[0];
rz(-1.9687037) q[2];
sx q[2];
rz(-1.8430317) q[2];
sx q[2];
rz(-0.18323252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1556342) q[1];
sx q[1];
rz(-1.8556719) q[1];
sx q[1];
rz(0.023521544) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15150841) q[3];
sx q[3];
rz(-0.66473367) q[3];
sx q[3];
rz(0.96707771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3016004) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(-0.99003506) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(0.5501737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54116762) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(-0.11359544) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8407211) q[2];
sx q[2];
rz(-2.8340325) q[2];
sx q[2];
rz(-0.077066271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6039227) q[1];
sx q[1];
rz(-2.1426755) q[1];
sx q[1];
rz(0.13725431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21869603) q[3];
sx q[3];
rz(-0.23467206) q[3];
sx q[3];
rz(-2.4179539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4438641) q[0];
sx q[0];
rz(-1.0794909) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1137431) q[0];
sx q[0];
rz(-2.5115974) q[0];
sx q[0];
rz(-1.6508474) q[0];
x q[1];
rz(0.5251685) q[2];
sx q[2];
rz(-1.6095973) q[2];
sx q[2];
rz(-2.2280488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3770959) q[1];
sx q[1];
rz(-2.9829574) q[1];
sx q[1];
rz(-1.7736969) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4004585) q[3];
sx q[3];
rz(-2.1949265) q[3];
sx q[3];
rz(-2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(-0.56435895) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(-2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-1.0790496) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-2.3880312) q[1];
sx q[1];
rz(-2.7450096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2561589) q[0];
sx q[0];
rz(-2.0246756) q[0];
sx q[0];
rz(2.8594349) q[0];
rz(1.8814538) q[2];
sx q[2];
rz(-1.1144181) q[2];
sx q[2];
rz(0.059547101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5844333) q[1];
sx q[1];
rz(-1.4769308) q[1];
sx q[1];
rz(2.2319016) q[1];
rz(-pi) q[2];
rz(2.5572705) q[3];
sx q[3];
rz(-1.366368) q[3];
sx q[3];
rz(1.0014597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(1.4554626) q[3];
sx q[3];
rz(-2.5374135) q[3];
sx q[3];
rz(-0.82715183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7431188) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(-1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267245) q[0];
sx q[0];
rz(-1.4818026) q[0];
sx q[0];
rz(0.29863775) q[0];
x q[1];
rz(1.8225841) q[2];
sx q[2];
rz(-0.52201027) q[2];
sx q[2];
rz(-2.8587647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3583402) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.3621484) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3133615) q[3];
sx q[3];
rz(-0.35850393) q[3];
sx q[3];
rz(3.0646216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.209098) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(0.15360019) q[2];
rz(2.8364654) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(-0.11288189) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(-3.1088366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(2.2420922) q[0];
x q[1];
rz(0.92808) q[2];
sx q[2];
rz(-0.66814458) q[2];
sx q[2];
rz(-1.966114) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7777562) q[1];
sx q[1];
rz(-1.5118074) q[1];
sx q[1];
rz(-1.1532409) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.350358) q[3];
sx q[3];
rz(-2.7535097) q[3];
sx q[3];
rz(-2.2409852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1774566) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-0.39045236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(-1.6149678) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-2.656235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55886666) q[0];
sx q[0];
rz(-2.9968046) q[0];
sx q[0];
rz(2.7607714) q[0];
rz(2.543407) q[2];
sx q[2];
rz(-1.7028114) q[2];
sx q[2];
rz(-0.26192947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.607) q[1];
sx q[1];
rz(-1.2619164) q[1];
sx q[1];
rz(-2.328518) q[1];
x q[2];
rz(-0.58166196) q[3];
sx q[3];
rz(-1.9546486) q[3];
sx q[3];
rz(-1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(-1.738328) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(-0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6195246) q[0];
sx q[0];
rz(-0.63729006) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78865096) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(1.3593332) q[0];
x q[1];
rz(1.7190476) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(1.6809747) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88565247) q[1];
sx q[1];
rz(-1.9347895) q[1];
sx q[1];
rz(3.0401858) q[1];
x q[2];
rz(2.3365799) q[3];
sx q[3];
rz(-0.22634889) q[3];
sx q[3];
rz(-1.602802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.082211994) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-0.65199488) q[2];
rz(0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(-1.1317071) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(1.9051753) q[1];
sx q[1];
rz(-0.9691144) q[1];
sx q[1];
rz(-1.2126927) q[1];
rz(-1.6937704) q[2];
sx q[2];
rz(-2.1868144) q[2];
sx q[2];
rz(-1.1745324) q[2];
rz(-0.39246172) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
