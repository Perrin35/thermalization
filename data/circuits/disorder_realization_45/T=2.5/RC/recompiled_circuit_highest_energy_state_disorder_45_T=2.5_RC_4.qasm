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
rz(-0.51195872) q[0];
sx q[0];
rz(-2.7380138) q[0];
sx q[0];
rz(-0.0041178218) q[0];
rz(-2.4601958) q[1];
sx q[1];
rz(-1.8713142) q[1];
sx q[1];
rz(-1.4407925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82872921) q[0];
sx q[0];
rz(-1.2932737) q[0];
sx q[0];
rz(-0.1779495) q[0];
rz(0.75136331) q[2];
sx q[2];
rz(-1.3168502) q[2];
sx q[2];
rz(-0.51515173) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.6328675) q[1];
sx q[1];
rz(-2.3687216) q[1];
sx q[1];
rz(2.6231324) q[1];
x q[2];
rz(0.012083455) q[3];
sx q[3];
rz(-1.7444495) q[3];
sx q[3];
rz(2.1326667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19930856) q[2];
sx q[2];
rz(-1.0465304) q[2];
sx q[2];
rz(0.47041565) q[2];
rz(2.2453902) q[3];
sx q[3];
rz(-1.4678518) q[3];
sx q[3];
rz(-1.5508274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.6386221) q[0];
sx q[0];
rz(-0.48415411) q[0];
sx q[0];
rz(0.36889398) q[0];
rz(2.6781354) q[1];
sx q[1];
rz(-1.8577441) q[1];
sx q[1];
rz(-0.2643815) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48369147) q[0];
sx q[0];
rz(-0.92934791) q[0];
sx q[0];
rz(0.19199706) q[0];
rz(2.8910341) q[2];
sx q[2];
rz(-1.7224285) q[2];
sx q[2];
rz(0.79094584) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2575327) q[1];
sx q[1];
rz(-2.7142378) q[1];
sx q[1];
rz(-0.76622574) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78231298) q[3];
sx q[3];
rz(-1.6564329) q[3];
sx q[3];
rz(2.8748363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5276864) q[2];
sx q[2];
rz(-0.1656342) q[2];
sx q[2];
rz(-1.3045093) q[2];
rz(0.3197318) q[3];
sx q[3];
rz(-2.3535862) q[3];
sx q[3];
rz(0.50362292) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9673135) q[0];
sx q[0];
rz(-0.19198424) q[0];
sx q[0];
rz(-2.018003) q[0];
rz(2.7536821) q[1];
sx q[1];
rz(-1.265641) q[1];
sx q[1];
rz(2.8111615) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0200121) q[0];
sx q[0];
rz(-2.8224484) q[0];
sx q[0];
rz(-0.97346755) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60744384) q[2];
sx q[2];
rz(-1.4943549) q[2];
sx q[2];
rz(0.7668524) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36281113) q[1];
sx q[1];
rz(-1.8900202) q[1];
sx q[1];
rz(-0.37324639) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73883812) q[3];
sx q[3];
rz(-0.78966245) q[3];
sx q[3];
rz(1.9146321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.419751) q[2];
sx q[2];
rz(-1.9998877) q[2];
sx q[2];
rz(1.9640131) q[2];
rz(1.2490595) q[3];
sx q[3];
rz(-1.306059) q[3];
sx q[3];
rz(-3.1031928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053442001) q[0];
sx q[0];
rz(-1.7120687) q[0];
sx q[0];
rz(-3.0528659) q[0];
rz(2.7170722) q[1];
sx q[1];
rz(-0.71966925) q[1];
sx q[1];
rz(-1.7321865) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78414375) q[0];
sx q[0];
rz(-1.3848834) q[0];
sx q[0];
rz(-1.9108514) q[0];
rz(-pi) q[1];
rz(2.652651) q[2];
sx q[2];
rz(-1.0283111) q[2];
sx q[2];
rz(-0.20636339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3864435) q[1];
sx q[1];
rz(-2.1514822) q[1];
sx q[1];
rz(-2.0232852) q[1];
rz(-2.0921246) q[3];
sx q[3];
rz(-1.8531257) q[3];
sx q[3];
rz(1.2897593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38922629) q[2];
sx q[2];
rz(-1.0800635) q[2];
sx q[2];
rz(-2.7092095) q[2];
rz(0.45807517) q[3];
sx q[3];
rz(-1.658344) q[3];
sx q[3];
rz(1.0079916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3573989) q[0];
sx q[0];
rz(-2.9965239) q[0];
sx q[0];
rz(-0.31132895) q[0];
rz(1.9642824) q[1];
sx q[1];
rz(-1.177634) q[1];
sx q[1];
rz(-2.6572773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84660778) q[0];
sx q[0];
rz(-1.7292882) q[0];
sx q[0];
rz(-0.096317795) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2782477) q[2];
sx q[2];
rz(-1.6451837) q[2];
sx q[2];
rz(2.0015581) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53953417) q[1];
sx q[1];
rz(-1.1975529) q[1];
sx q[1];
rz(2.9887524) q[1];
x q[2];
rz(-1.442904) q[3];
sx q[3];
rz(-1.4223243) q[3];
sx q[3];
rz(-3.0612891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3207265) q[2];
sx q[2];
rz(-1.4297239) q[2];
sx q[2];
rz(2.7166264) q[2];
rz(-3.037437) q[3];
sx q[3];
rz(-2.7724373) q[3];
sx q[3];
rz(2.4764376) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3686309) q[0];
sx q[0];
rz(-2.5017128) q[0];
sx q[0];
rz(0.1567008) q[0];
rz(-2.8796097) q[1];
sx q[1];
rz(-1.1527088) q[1];
sx q[1];
rz(1.6798457) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91297594) q[0];
sx q[0];
rz(-1.3973087) q[0];
sx q[0];
rz(0.96488379) q[0];
x q[1];
rz(-0.22905519) q[2];
sx q[2];
rz(-0.88729268) q[2];
sx q[2];
rz(-2.7102269) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9380381) q[1];
sx q[1];
rz(-2.9239836) q[1];
sx q[1];
rz(-0.59534351) q[1];
rz(-0.42243345) q[3];
sx q[3];
rz(-2.3452873) q[3];
sx q[3];
rz(2.5811206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4713952) q[2];
sx q[2];
rz(-1.6452226) q[2];
sx q[2];
rz(-0.71693286) q[2];
rz(-0.21643058) q[3];
sx q[3];
rz(-1.8559034) q[3];
sx q[3];
rz(1.9560248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4099429) q[0];
sx q[0];
rz(-2.8493632) q[0];
sx q[0];
rz(-1.2662079) q[0];
rz(-0.75824291) q[1];
sx q[1];
rz(-1.1781324) q[1];
sx q[1];
rz(2.0936802) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3558725) q[0];
sx q[0];
rz(-1.2730755) q[0];
sx q[0];
rz(2.1773318) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1282461) q[2];
sx q[2];
rz(-1.1801022) q[2];
sx q[2];
rz(-1.7683795) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3057237) q[1];
sx q[1];
rz(-0.36986884) q[1];
sx q[1];
rz(-1.3772411) q[1];
rz(-pi) q[2];
rz(1.326845) q[3];
sx q[3];
rz(-1.1317963) q[3];
sx q[3];
rz(1.148996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3332112) q[2];
sx q[2];
rz(-1.0071249) q[2];
sx q[2];
rz(0.0086616596) q[2];
rz(-1.0041142) q[3];
sx q[3];
rz(-2.6486371) q[3];
sx q[3];
rz(-2.139411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36326161) q[0];
sx q[0];
rz(-0.14334981) q[0];
sx q[0];
rz(2.1891201) q[0];
rz(-1.6078423) q[1];
sx q[1];
rz(-1.8550823) q[1];
sx q[1];
rz(1.0362157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8995953) q[0];
sx q[0];
rz(-1.138666) q[0];
sx q[0];
rz(2.5431387) q[0];
rz(-pi) q[1];
rz(1.9373193) q[2];
sx q[2];
rz(-2.1749305) q[2];
sx q[2];
rz(2.0458536) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0091102) q[1];
sx q[1];
rz(-0.81071893) q[1];
sx q[1];
rz(-0.85217635) q[1];
rz(-1.1187068) q[3];
sx q[3];
rz(-0.89632672) q[3];
sx q[3];
rz(1.3687641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2967534) q[2];
sx q[2];
rz(-2.4256458) q[2];
sx q[2];
rz(1.0003662) q[2];
rz(-1.865546) q[3];
sx q[3];
rz(-1.4714656) q[3];
sx q[3];
rz(-2.0744417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0412096) q[0];
sx q[0];
rz(-0.90827933) q[0];
sx q[0];
rz(2.8874183) q[0];
rz(-1.8036448) q[1];
sx q[1];
rz(-2.281052) q[1];
sx q[1];
rz(2.3635704) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4648424) q[0];
sx q[0];
rz(-1.8942231) q[0];
sx q[0];
rz(2.4535562) q[0];
x q[1];
rz(-2.4619589) q[2];
sx q[2];
rz(-0.25654116) q[2];
sx q[2];
rz(0.40479615) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3624448) q[1];
sx q[1];
rz(-2.1333284) q[1];
sx q[1];
rz(2.0570807) q[1];
rz(1.3555525) q[3];
sx q[3];
rz(-1.6674893) q[3];
sx q[3];
rz(-1.7413063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.079387) q[2];
sx q[2];
rz(-1.6239245) q[2];
sx q[2];
rz(-3.0089231) q[2];
rz(-2.0343871) q[3];
sx q[3];
rz(-2.5532494) q[3];
sx q[3];
rz(1.6975105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4001813) q[0];
sx q[0];
rz(-1.9453229) q[0];
sx q[0];
rz(2.3924526) q[0];
rz(-2.6465042) q[1];
sx q[1];
rz(-1.9063213) q[1];
sx q[1];
rz(1.3912158) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.063569) q[0];
sx q[0];
rz(-1.1207523) q[0];
sx q[0];
rz(0.70722945) q[0];
rz(1.130213) q[2];
sx q[2];
rz(-1.9398707) q[2];
sx q[2];
rz(2.1598494) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5500945) q[1];
sx q[1];
rz(-1.6143394) q[1];
sx q[1];
rz(-0.40710475) q[1];
rz(-pi) q[2];
rz(-0.81328525) q[3];
sx q[3];
rz(-0.67768598) q[3];
sx q[3];
rz(2.4285674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5007925) q[2];
sx q[2];
rz(-1.4747138) q[2];
sx q[2];
rz(-0.87727171) q[2];
rz(0.52214617) q[3];
sx q[3];
rz(-0.79972655) q[3];
sx q[3];
rz(1.6845901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7397639) q[0];
sx q[0];
rz(-0.54720989) q[0];
sx q[0];
rz(-1.3450958) q[0];
rz(-1.2783891) q[1];
sx q[1];
rz(-0.31247333) q[1];
sx q[1];
rz(-0.023275274) q[1];
rz(-1.5398434) q[2];
sx q[2];
rz(-2.2014115) q[2];
sx q[2];
rz(1.9104107) q[2];
rz(2.4458281) q[3];
sx q[3];
rz(-0.56752612) q[3];
sx q[3];
rz(-0.023434536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
