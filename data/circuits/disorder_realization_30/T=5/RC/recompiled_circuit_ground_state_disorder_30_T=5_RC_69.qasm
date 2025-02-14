OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(2.6240786) q[0];
rz(-1.8139047) q[1];
sx q[1];
rz(-0.89193901) q[1];
sx q[1];
rz(2.5995624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21236496) q[0];
sx q[0];
rz(-2.3327418) q[0];
sx q[0];
rz(1.8083284) q[0];
rz(1.8250685) q[2];
sx q[2];
rz(-0.82077998) q[2];
sx q[2];
rz(-0.17707846) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65329018) q[1];
sx q[1];
rz(-0.79217109) q[1];
sx q[1];
rz(-1.8503055) q[1];
rz(-pi) q[2];
rz(-2.9085701) q[3];
sx q[3];
rz(-1.7400422) q[3];
sx q[3];
rz(0.47692933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0077670495) q[2];
sx q[2];
rz(-2.0962891) q[2];
sx q[2];
rz(0.93057752) q[2];
rz(-0.13115701) q[3];
sx q[3];
rz(-2.7635062) q[3];
sx q[3];
rz(-1.4906496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5432878) q[0];
sx q[0];
rz(-0.85619339) q[0];
sx q[0];
rz(-1.244586) q[0];
rz(2.2159131) q[1];
sx q[1];
rz(-1.8857748) q[1];
sx q[1];
rz(1.6474887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5410575) q[0];
sx q[0];
rz(-0.78406683) q[0];
sx q[0];
rz(-2.0406141) q[0];
rz(-pi) q[1];
rz(-1.235512) q[2];
sx q[2];
rz(-0.92520324) q[2];
sx q[2];
rz(-0.52926999) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0665613) q[1];
sx q[1];
rz(-1.2258471) q[1];
sx q[1];
rz(-2.6607772) q[1];
rz(-pi) q[2];
rz(0.027493942) q[3];
sx q[3];
rz(-2.8373233) q[3];
sx q[3];
rz(-2.9472443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.99574789) q[2];
sx q[2];
rz(-0.44439849) q[2];
sx q[2];
rz(2.2399529) q[2];
rz(-1.3580648) q[3];
sx q[3];
rz(-1.3172904) q[3];
sx q[3];
rz(2.5942514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667592) q[0];
sx q[0];
rz(-1.1792553) q[0];
sx q[0];
rz(-1.1460079) q[0];
rz(-2.215812) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(-2.944223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5615047) q[0];
sx q[0];
rz(-0.60778032) q[0];
sx q[0];
rz(1.32919) q[0];
rz(-0.80855753) q[2];
sx q[2];
rz(-0.41038358) q[2];
sx q[2];
rz(1.9661599) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6655953) q[1];
sx q[1];
rz(-1.531812) q[1];
sx q[1];
rz(-1.798756) q[1];
rz(-pi) q[2];
rz(1.3696543) q[3];
sx q[3];
rz(-1.0363058) q[3];
sx q[3];
rz(-2.8902658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4027412) q[2];
sx q[2];
rz(-0.98029843) q[2];
sx q[2];
rz(2.951238) q[2];
rz(3.0800152) q[3];
sx q[3];
rz(-2.172251) q[3];
sx q[3];
rz(1.0524582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9848118) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(0.62830997) q[0];
rz(1.620694) q[1];
sx q[1];
rz(-1.2077121) q[1];
sx q[1];
rz(-0.77888387) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.223004) q[0];
sx q[0];
rz(-2.1445334) q[0];
sx q[0];
rz(-0.45044915) q[0];
rz(-0.65108786) q[2];
sx q[2];
rz(-2.6991803) q[2];
sx q[2];
rz(-2.173732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0592093) q[1];
sx q[1];
rz(-1.6184855) q[1];
sx q[1];
rz(-1.5573182) q[1];
x q[2];
rz(-0.86872673) q[3];
sx q[3];
rz(-1.0208289) q[3];
sx q[3];
rz(-1.3438674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4817619) q[2];
sx q[2];
rz(-2.063282) q[2];
sx q[2];
rz(-3.1415494) q[2];
rz(-1.0846042) q[3];
sx q[3];
rz(-1.9856039) q[3];
sx q[3];
rz(0.46561766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67149177) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(-1.2985562) q[0];
rz(-0.57688722) q[1];
sx q[1];
rz(-1.8678317) q[1];
sx q[1];
rz(2.6947122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90335411) q[0];
sx q[0];
rz(-2.2132259) q[0];
sx q[0];
rz(1.7967671) q[0];
rz(1.8700897) q[2];
sx q[2];
rz(-1.4993877) q[2];
sx q[2];
rz(-1.4324607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2858823) q[1];
sx q[1];
rz(-1.1709164) q[1];
sx q[1];
rz(2.9215496) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6303667) q[3];
sx q[3];
rz(-1.7559663) q[3];
sx q[3];
rz(-0.51051729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.4345066) q[2];
sx q[2];
rz(-2.862317) q[2];
sx q[2];
rz(0.2198098) q[2];
rz(-1.48742) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(2.4421104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(0.18911067) q[0];
sx q[0];
rz(-0.43287745) q[0];
sx q[0];
rz(0.89749807) q[0];
rz(1.3709566) q[1];
sx q[1];
rz(-2.1015344) q[1];
sx q[1];
rz(0.8955566) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5286551) q[0];
sx q[0];
rz(-2.2870758) q[0];
sx q[0];
rz(-1.3837293) q[0];
rz(-0.61195749) q[2];
sx q[2];
rz(-1.0494119) q[2];
sx q[2];
rz(-2.251791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5472153) q[1];
sx q[1];
rz(-1.2378232) q[1];
sx q[1];
rz(-0.39526387) q[1];
rz(-2.6704253) q[3];
sx q[3];
rz(-2.0063324) q[3];
sx q[3];
rz(2.7903231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36876496) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(-0.15711288) q[2];
rz(0.67240063) q[3];
sx q[3];
rz(-1.3998569) q[3];
sx q[3];
rz(-0.11252832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86033487) q[0];
sx q[0];
rz(-2.4761138) q[0];
sx q[0];
rz(-1.6819287) q[0];
rz(0.36059391) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(1.8416539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7304218) q[0];
sx q[0];
rz(-1.9864559) q[0];
sx q[0];
rz(-1.704373) q[0];
rz(-pi) q[1];
rz(-1.1088013) q[2];
sx q[2];
rz(-2.1401775) q[2];
sx q[2];
rz(1.5608567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2092525) q[1];
sx q[1];
rz(-2.8289218) q[1];
sx q[1];
rz(-2.8220842) q[1];
x q[2];
rz(0.54706562) q[3];
sx q[3];
rz(-2.4506042) q[3];
sx q[3];
rz(0.83741659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.931687) q[2];
sx q[2];
rz(-2.1114712) q[2];
sx q[2];
rz(0.21200655) q[2];
rz(1.0509521) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(3.0534548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6050922) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(2.9686046) q[0];
rz(-1.1146924) q[1];
sx q[1];
rz(-1.0429017) q[1];
sx q[1];
rz(-3.0317422) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3152471) q[0];
sx q[0];
rz(-1.7638806) q[0];
sx q[0];
rz(-0.56435926) q[0];
rz(-pi) q[1];
rz(-0.82288701) q[2];
sx q[2];
rz(-1.7639065) q[2];
sx q[2];
rz(-0.88179526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.122095) q[1];
sx q[1];
rz(-1.8260584) q[1];
sx q[1];
rz(0.17664102) q[1];
rz(-1.4509216) q[3];
sx q[3];
rz(-1.0195512) q[3];
sx q[3];
rz(-2.3012183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7383808) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(2.476725) q[2];
rz(-2.6730149) q[3];
sx q[3];
rz(-1.5237619) q[3];
sx q[3];
rz(-0.81336462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68798962) q[0];
sx q[0];
rz(-2.5316694) q[0];
sx q[0];
rz(0.74158057) q[0];
rz(1.1874416) q[1];
sx q[1];
rz(-1.5267742) q[1];
sx q[1];
rz(2.5724519) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5461972) q[0];
sx q[0];
rz(-2.6962207) q[0];
sx q[0];
rz(0.73308848) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9958565) q[2];
sx q[2];
rz(-0.45782858) q[2];
sx q[2];
rz(-2.482058) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5829825) q[1];
sx q[1];
rz(-2.4931968) q[1];
sx q[1];
rz(2.3984539) q[1];
x q[2];
rz(-2.3145202) q[3];
sx q[3];
rz(-2.2571392) q[3];
sx q[3];
rz(-2.4857869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49447122) q[2];
sx q[2];
rz(-2.1529866) q[2];
sx q[2];
rz(-0.76623255) q[2];
rz(-1.1726441) q[3];
sx q[3];
rz(-1.5394883) q[3];
sx q[3];
rz(-2.9197781) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2072993) q[0];
sx q[0];
rz(-0.3485637) q[0];
sx q[0];
rz(-2.9497414) q[0];
rz(0.81709298) q[1];
sx q[1];
rz(-0.25661883) q[1];
sx q[1];
rz(-2.4339035) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548508) q[0];
sx q[0];
rz(-1.3258717) q[0];
sx q[0];
rz(-1.3913416) q[0];
rz(-pi) q[1];
rz(2.0470263) q[2];
sx q[2];
rz(-0.64733887) q[2];
sx q[2];
rz(0.29905427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7308049) q[1];
sx q[1];
rz(-1.9390701) q[1];
sx q[1];
rz(-2.1965697) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97174208) q[3];
sx q[3];
rz(-1.069456) q[3];
sx q[3];
rz(-0.80162345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91844687) q[2];
sx q[2];
rz(-0.94474363) q[2];
sx q[2];
rz(-2.4165912) q[2];
rz(0.21507344) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(2.422629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216777) q[0];
sx q[0];
rz(-0.81659962) q[0];
sx q[0];
rz(0.30246977) q[0];
rz(2.0401781) q[1];
sx q[1];
rz(-1.2322203) q[1];
sx q[1];
rz(0.89422918) q[1];
rz(-0.11753043) q[2];
sx q[2];
rz(-1.0385658) q[2];
sx q[2];
rz(-0.081296878) q[2];
rz(-0.054417944) q[3];
sx q[3];
rz(-2.6734753) q[3];
sx q[3];
rz(-2.5772167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
