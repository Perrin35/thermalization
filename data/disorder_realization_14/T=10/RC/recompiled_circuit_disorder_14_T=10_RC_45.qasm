OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(4.598773) q[1];
sx q[1];
rz(7.4814319) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449013) q[0];
sx q[0];
rz(-0.20151073) q[0];
sx q[0];
rz(-2.6411396) q[0];
rz(-2.3660907) q[2];
sx q[2];
rz(-1.3368133) q[2];
sx q[2];
rz(-2.6393059) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6557505) q[1];
sx q[1];
rz(-1.9006923) q[1];
sx q[1];
rz(-0.54539036) q[1];
rz(-3.0953193) q[3];
sx q[3];
rz(-1.7414021) q[3];
sx q[3];
rz(1.5266649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9464232) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-0.18134376) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-0.75631022) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392035) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(-2.7547577) q[0];
rz(0.50239262) q[1];
sx q[1];
rz(-0.97351176) q[1];
sx q[1];
rz(-1.5418672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4138448) q[0];
sx q[0];
rz(-1.9808597) q[0];
sx q[0];
rz(-2.4983665) q[0];
x q[1];
rz(2.0273655) q[2];
sx q[2];
rz(-1.7184988) q[2];
sx q[2];
rz(0.6870803) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6501573) q[1];
sx q[1];
rz(-0.037089247) q[1];
sx q[1];
rz(2.8703719) q[1];
x q[2];
rz(-2.0412444) q[3];
sx q[3];
rz(-1.1332266) q[3];
sx q[3];
rz(0.13089422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5388422) q[2];
sx q[2];
rz(-0.87783146) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(2.291262) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.26281115) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(-0.61022726) q[0];
rz(1.2894851) q[1];
sx q[1];
rz(-2.162343) q[1];
sx q[1];
rz(-2.1496225) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1410071) q[0];
sx q[0];
rz(-1.7551384) q[0];
sx q[0];
rz(2.1328451) q[0];
rz(-1.7730373) q[2];
sx q[2];
rz(-2.2291406) q[2];
sx q[2];
rz(0.973268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5132644) q[1];
sx q[1];
rz(-2.6918415) q[1];
sx q[1];
rz(2.9343534) q[1];
x q[2];
rz(-2.6842327) q[3];
sx q[3];
rz(-1.3501581) q[3];
sx q[3];
rz(1.8822576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.758574) q[2];
sx q[2];
rz(-0.16123161) q[2];
sx q[2];
rz(-2.8733011) q[2];
rz(-0.39408436) q[3];
sx q[3];
rz(-1.9106617) q[3];
sx q[3];
rz(-0.18850732) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2878993) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(0.40507856) q[0];
rz(0.69008094) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(2.4437723) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0426548) q[0];
sx q[0];
rz(-1.5942759) q[0];
sx q[0];
rz(1.4809181) q[0];
rz(-pi) q[1];
rz(-1.9525098) q[2];
sx q[2];
rz(-2.082798) q[2];
sx q[2];
rz(-0.94101671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.3623558) q[1];
sx q[1];
rz(-0.52473611) q[1];
sx q[1];
rz(1.1974105) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.97045578) q[3];
sx q[3];
rz(-2.6451023) q[3];
sx q[3];
rz(-2.0127726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71965376) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.3556577) q[0];
sx q[0];
rz(2.1160545) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(0.62932032) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9176661) q[0];
sx q[0];
rz(-1.9804945) q[0];
sx q[0];
rz(-1.2439338) q[0];
rz(-1.2567369) q[2];
sx q[2];
rz(-2.4410508) q[2];
sx q[2];
rz(-0.26740057) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3931261) q[1];
sx q[1];
rz(-1.6854291) q[1];
sx q[1];
rz(-0.51490358) q[1];
x q[2];
rz(0.96962813) q[3];
sx q[3];
rz(-1.789635) q[3];
sx q[3];
rz(0.18274433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59262529) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(-2.7992451) q[2];
rz(-1.7957548) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(0.19601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65258566) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(-0.18519369) q[0];
rz(-1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(1.3669744) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6706086) q[0];
sx q[0];
rz(-0.99262041) q[0];
sx q[0];
rz(-1.3060119) q[0];
rz(-pi) q[1];
rz(0.30250678) q[2];
sx q[2];
rz(-2.2777646) q[2];
sx q[2];
rz(-2.6940341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1320912) q[1];
sx q[1];
rz(-2.9638634) q[1];
sx q[1];
rz(-1.0264261) q[1];
rz(-pi) q[2];
rz(-0.6188789) q[3];
sx q[3];
rz(-0.66580171) q[3];
sx q[3];
rz(-0.57951365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(0.883376) q[2];
rz(-1.4128489) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(0.81378716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8900523) q[0];
sx q[0];
rz(-3.0391356) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(-3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(-1.3758804) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517075) q[0];
sx q[0];
rz(-2.2635248) q[0];
sx q[0];
rz(-0.38498199) q[0];
rz(-pi) q[1];
rz(1.6872348) q[2];
sx q[2];
rz(-0.47633007) q[2];
sx q[2];
rz(1.6778698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77586506) q[1];
sx q[1];
rz(-2.9850246) q[1];
sx q[1];
rz(-0.74128976) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7405628) q[3];
sx q[3];
rz(-2.8068636) q[3];
sx q[3];
rz(1.3917805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(-0.73105556) q[2];
rz(-0.11080065) q[3];
sx q[3];
rz(-1.585107) q[3];
sx q[3];
rz(-0.69537648) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-2.2667363) q[0];
sx q[0];
rz(2.4080283) q[0];
rz(-0.60797524) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(-0.2342934) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5406571) q[0];
sx q[0];
rz(-2.0458851) q[0];
sx q[0];
rz(-0.29151543) q[0];
x q[1];
rz(0.90227622) q[2];
sx q[2];
rz(-2.1736645) q[2];
sx q[2];
rz(-0.55842802) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.81208166) q[1];
sx q[1];
rz(-0.78204621) q[1];
sx q[1];
rz(-1.300632) q[1];
x q[2];
rz(1.0844564) q[3];
sx q[3];
rz(-0.64389766) q[3];
sx q[3];
rz(-0.38734303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.12604788) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(1.7162494) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(1.19338) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(-2.1967922) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1244922) q[0];
sx q[0];
rz(-1.955535) q[0];
sx q[0];
rz(2.4186633) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2009732) q[2];
sx q[2];
rz(-0.47626469) q[2];
sx q[2];
rz(-0.95048743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8069357) q[1];
sx q[1];
rz(-1.8493223) q[1];
sx q[1];
rz(1.3929699) q[1];
rz(-pi) q[2];
rz(1.2148083) q[3];
sx q[3];
rz(-1.755852) q[3];
sx q[3];
rz(-1.2026279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9251359) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(0.28820583) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-1.0498485) q[3];
sx q[3];
rz(1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(2.2286041) q[0];
rz(-2.7669725) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(0.8909117) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64542949) q[0];
sx q[0];
rz(-1.2362288) q[0];
sx q[0];
rz(-1.8339001) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.183379) q[2];
sx q[2];
rz(-1.805086) q[2];
sx q[2];
rz(-2.2189552) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18894698) q[1];
sx q[1];
rz(-1.9083605) q[1];
sx q[1];
rz(1.360421) q[1];
rz(-pi) q[2];
rz(1.6463514) q[3];
sx q[3];
rz(-1.75845) q[3];
sx q[3];
rz(-0.56896602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23641071) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(-0.50160828) q[2];
rz(1.8524528) q[3];
sx q[3];
rz(-1.6882378) q[3];
sx q[3];
rz(-2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63507737) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(-2.0093909) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-0.51433993) q[2];
sx q[2];
rz(-0.30403501) q[2];
sx q[2];
rz(-2.4652849) q[2];
rz(2.1174259) q[3];
sx q[3];
rz(-1.6129941) q[3];
sx q[3];
rz(2.2254406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
