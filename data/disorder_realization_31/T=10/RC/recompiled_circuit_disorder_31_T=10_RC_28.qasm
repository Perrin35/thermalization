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
rz(0.84258643) q[0];
rz(2.1057582) q[1];
sx q[1];
rz(8.3254568) q[1];
sx q[1];
rz(7.96666) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1928756) q[0];
sx q[0];
rz(-2.2964381) q[0];
sx q[0];
rz(1.2851508) q[0];
x q[1];
rz(-3.1137142) q[2];
sx q[2];
rz(-0.61402245) q[2];
sx q[2];
rz(0.30464722) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.018651389) q[1];
sx q[1];
rz(-2.7416347) q[1];
sx q[1];
rz(2.8040228) q[1];
rz(-2.9620142) q[3];
sx q[3];
rz(-0.60185963) q[3];
sx q[3];
rz(1.7367401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6136916) q[2];
sx q[2];
rz(-1.0062904) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(-1.9159296) q[3];
sx q[3];
rz(-1.3464728) q[3];
sx q[3];
rz(2.3195482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74801385) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(0.32546145) q[0];
rz(1.7851967) q[1];
sx q[1];
rz(-1.0486832) q[1];
sx q[1];
rz(-1.9869841) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75630674) q[0];
sx q[0];
rz(-3.1161838) q[0];
sx q[0];
rz(2.4029762) q[0];
rz(-pi) q[1];
rz(-2.7484659) q[2];
sx q[2];
rz(-2.1596585) q[2];
sx q[2];
rz(1.7413505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.11101152) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(-2.3762523) q[1];
rz(0.96315893) q[3];
sx q[3];
rz(-0.31664407) q[3];
sx q[3];
rz(-2.7505927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4521728) q[2];
sx q[2];
rz(-1.8916811) q[2];
sx q[2];
rz(2.2581805) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-2.3538891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.5939012) q[1];
sx q[1];
rz(-1.0916969) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1611623) q[0];
sx q[0];
rz(-0.83320252) q[0];
sx q[0];
rz(2.3479793) q[0];
x q[1];
rz(0.91471471) q[2];
sx q[2];
rz(-1.9064184) q[2];
sx q[2];
rz(2.6732973) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25635438) q[1];
sx q[1];
rz(-1.9372205) q[1];
sx q[1];
rz(-2.3430235) q[1];
rz(-pi) q[2];
rz(3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(-1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8213356) q[2];
sx q[2];
rz(-2.6358423) q[2];
sx q[2];
rz(2.2606405) q[2];
rz(1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-1.0176456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83051935) q[0];
sx q[0];
rz(-1.7493462) q[0];
sx q[0];
rz(-2.7048892) q[0];
rz(2.9084335) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(-2.8312347) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46582857) q[0];
sx q[0];
rz(-2.7390263) q[0];
sx q[0];
rz(0.34253828) q[0];
rz(1.6967625) q[2];
sx q[2];
rz(-2.2519886) q[2];
sx q[2];
rz(-1.7484401) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0879285) q[1];
sx q[1];
rz(-1.9287319) q[1];
sx q[1];
rz(0.019836516) q[1];
rz(-pi) q[2];
rz(1.849732) q[3];
sx q[3];
rz(-0.12860563) q[3];
sx q[3];
rz(2.0898553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13005304) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(-1.0774353) q[2];
rz(0.056190101) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(1.5475387) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.189165) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(-0.24965832) q[0];
rz(-1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(0.87019428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0349883) q[0];
sx q[0];
rz(-2.7634794) q[0];
sx q[0];
rz(2.5614221) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9179847) q[2];
sx q[2];
rz(-0.70906559) q[2];
sx q[2];
rz(-2.2774334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2493077) q[1];
sx q[1];
rz(-1.2586437) q[1];
sx q[1];
rz(-2.409163) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16964511) q[3];
sx q[3];
rz(-1.0153474) q[3];
sx q[3];
rz(-2.5559705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.3262649) q[2];
sx q[2];
rz(0.67374054) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(1.3195066) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(-2.8836024) q[0];
rz(0.42516431) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(1.649883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314635) q[0];
sx q[0];
rz(-1.9995814) q[0];
sx q[0];
rz(-1.462912) q[0];
rz(-pi) q[1];
rz(-0.50478023) q[2];
sx q[2];
rz(-2.3961888) q[2];
sx q[2];
rz(2.8842852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7548435) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(-pi) q[2];
rz(-0.10156472) q[3];
sx q[3];
rz(-1.2079117) q[3];
sx q[3];
rz(0.18728072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1288746) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(0.84364676) q[3];
sx q[3];
rz(-0.97976145) q[3];
sx q[3];
rz(0.89404026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(-0.86589083) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(3.1076028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2798529) q[0];
sx q[0];
rz(-1.4089157) q[0];
sx q[0];
rz(0.78761657) q[0];
x q[1];
rz(-0.54079536) q[2];
sx q[2];
rz(-1.0058837) q[2];
sx q[2];
rz(0.62513798) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36044869) q[1];
sx q[1];
rz(-2.2478659) q[1];
sx q[1];
rz(3.0767246) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5266685) q[3];
sx q[3];
rz(-1.691754) q[3];
sx q[3];
rz(-1.3161591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-0.59843439) q[2];
sx q[2];
rz(-0.87654385) q[2];
rz(-0.34902469) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8975163) q[0];
sx q[0];
rz(-1.4325457) q[0];
sx q[0];
rz(0.39392719) q[0];
rz(2.774033) q[1];
sx q[1];
rz(-1.7575248) q[1];
sx q[1];
rz(-1.6961018) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1720393) q[0];
sx q[0];
rz(-1.5656099) q[0];
sx q[0];
rz(-2.0041549) q[0];
rz(-pi) q[1];
rz(2.1922853) q[2];
sx q[2];
rz(-1.0770123) q[2];
sx q[2];
rz(0.91044237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.98390025) q[1];
sx q[1];
rz(-1.0219814) q[1];
sx q[1];
rz(-1.3859315) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94533841) q[3];
sx q[3];
rz(-2.1350386) q[3];
sx q[3];
rz(-0.64627796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2945071) q[2];
sx q[2];
rz(-2.2512348) q[2];
sx q[2];
rz(2.7344446) q[2];
rz(-1.5173222) q[3];
sx q[3];
rz(-1.1573236) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80609926) q[0];
sx q[0];
rz(-0.51500106) q[0];
sx q[0];
rz(-1.2517713) q[0];
rz(2.4720526) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(2.8318185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71000242) q[0];
sx q[0];
rz(-2.6012523) q[0];
sx q[0];
rz(2.8938328) q[0];
x q[1];
rz(1.5179713) q[2];
sx q[2];
rz(-0.17715684) q[2];
sx q[2];
rz(1.0747386) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.138315) q[1];
sx q[1];
rz(-1.365108) q[1];
sx q[1];
rz(0.059271952) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5449949) q[3];
sx q[3];
rz(-1.784424) q[3];
sx q[3];
rz(-1.240977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4391675) q[2];
sx q[2];
rz(-2.4283786) q[2];
sx q[2];
rz(1.2072198) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.050215125) q[0];
sx q[0];
rz(-1.3239048) q[0];
sx q[0];
rz(1.9357095) q[0];
rz(2.5559015) q[1];
sx q[1];
rz(-2.060545) q[1];
sx q[1];
rz(1.6419798) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6584872) q[0];
sx q[0];
rz(-2.8327836) q[0];
sx q[0];
rz(1.2312066) q[0];
rz(2.3775616) q[2];
sx q[2];
rz(-1.521763) q[2];
sx q[2];
rz(1.3438091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2227576) q[1];
sx q[1];
rz(-1.3593874) q[1];
sx q[1];
rz(-1.097015) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5851192) q[3];
sx q[3];
rz(-2.9687772) q[3];
sx q[3];
rz(2.1567791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88400921) q[2];
sx q[2];
rz(-1.7929701) q[2];
sx q[2];
rz(-0.94669) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.5654516) q[3];
sx q[3];
rz(0.45599109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7832227) q[0];
sx q[0];
rz(-1.1483648) q[0];
sx q[0];
rz(-0.4233465) q[0];
rz(-0.070925698) q[1];
sx q[1];
rz(-1.4535041) q[1];
sx q[1];
rz(2.8765875) q[1];
rz(-1.0735738) q[2];
sx q[2];
rz(-1.9855269) q[2];
sx q[2];
rz(-1.2863458) q[2];
rz(-2.4214217) q[3];
sx q[3];
rz(-1.9577033) q[3];
sx q[3];
rz(-2.4693558) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
