OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(3.6462311) q[0];
sx q[0];
rz(12.991821) q[0];
rz(0.57234859) q[1];
sx q[1];
rz(-2.0444137) q[1];
sx q[1];
rz(-1.9414577) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.67681) q[0];
sx q[0];
rz(-1.3199782) q[0];
sx q[0];
rz(-1.5808015) q[0];
rz(0.91325735) q[2];
sx q[2];
rz(-2.7139671) q[2];
sx q[2];
rz(-1.0853801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7537743) q[1];
sx q[1];
rz(-1.942831) q[1];
sx q[1];
rz(1.8005936) q[1];
rz(-pi) q[2];
rz(0.67114222) q[3];
sx q[3];
rz(-1.6832388) q[3];
sx q[3];
rz(0.46026106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9878865) q[2];
sx q[2];
rz(-1.1517297) q[2];
sx q[2];
rz(-0.64725867) q[2];
rz(-0.17793947) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(-1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4030289) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(-0.52363288) q[0];
rz(-1.2298443) q[1];
sx q[1];
rz(-0.74408999) q[1];
sx q[1];
rz(-2.8935166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9427633) q[0];
sx q[0];
rz(-0.83434767) q[0];
sx q[0];
rz(-1.426986) q[0];
rz(-2.7047728) q[2];
sx q[2];
rz(-1.9736276) q[2];
sx q[2];
rz(0.4536597) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4011127) q[1];
sx q[1];
rz(-1.4251627) q[1];
sx q[1];
rz(0.039197368) q[1];
rz(0.6389736) q[3];
sx q[3];
rz(-2.3666214) q[3];
sx q[3];
rz(-2.4955622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6123885) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(2.6017792) q[2];
rz(-2.9811033) q[3];
sx q[3];
rz(-1.5925946) q[3];
sx q[3];
rz(-0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983343) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(-1.5191822) q[1];
sx q[1];
rz(-1.2956053) q[1];
sx q[1];
rz(0.79230961) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80817079) q[0];
sx q[0];
rz(-1.9229888) q[0];
sx q[0];
rz(-0.19547021) q[0];
rz(0.14459212) q[2];
sx q[2];
rz(-1.7011614) q[2];
sx q[2];
rz(3.0693288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6349259) q[1];
sx q[1];
rz(-2.185195) q[1];
sx q[1];
rz(2.3022815) q[1];
rz(-pi) q[2];
rz(0.00040690502) q[3];
sx q[3];
rz(-0.30012977) q[3];
sx q[3];
rz(-2.2426734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85492674) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(2.3112042) q[2];
rz(1.7928127) q[3];
sx q[3];
rz(-2.1068137) q[3];
sx q[3];
rz(2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26375672) q[0];
sx q[0];
rz(-2.1649375) q[0];
sx q[0];
rz(-2.7184955) q[0];
rz(1.1571723) q[1];
sx q[1];
rz(-1.6559699) q[1];
sx q[1];
rz(2.5194936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4870479) q[0];
sx q[0];
rz(-1.9584987) q[0];
sx q[0];
rz(0.88519208) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2086925) q[2];
sx q[2];
rz(-1.9819753) q[2];
sx q[2];
rz(1.5087939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3286666) q[1];
sx q[1];
rz(-1.5541967) q[1];
sx q[1];
rz(-1.9783201) q[1];
rz(-pi) q[2];
rz(1.6016209) q[3];
sx q[3];
rz(-0.72576952) q[3];
sx q[3];
rz(0.14528615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28079924) q[2];
sx q[2];
rz(-1.4192899) q[2];
sx q[2];
rz(-0.61817509) q[2];
rz(-1.6587967) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16267714) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(2.3090114) q[0];
rz(0.32514462) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(-0.064373374) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58080855) q[0];
sx q[0];
rz(-2.7940895) q[0];
sx q[0];
rz(-2.7827713) q[0];
rz(-pi) q[1];
rz(-0.60466296) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(-0.99550216) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87582237) q[1];
sx q[1];
rz(-2.0811305) q[1];
sx q[1];
rz(2.2052664) q[1];
rz(-pi) q[2];
rz(0.35803087) q[3];
sx q[3];
rz(-0.97409407) q[3];
sx q[3];
rz(-0.95355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9084106) q[2];
sx q[2];
rz(-0.57047129) q[2];
sx q[2];
rz(1.2916267) q[2];
rz(0.49606797) q[3];
sx q[3];
rz(-2.3521164) q[3];
sx q[3];
rz(-1.1424278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9555776) q[0];
sx q[0];
rz(-0.29009524) q[0];
sx q[0];
rz(2.7225323) q[0];
rz(-2.8298607) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(2.9980803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9547894) q[0];
sx q[0];
rz(-1.1214897) q[0];
sx q[0];
rz(1.1727929) q[0];
rz(1.3355876) q[2];
sx q[2];
rz(-1.303788) q[2];
sx q[2];
rz(0.63797229) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0106502) q[1];
sx q[1];
rz(-0.82320628) q[1];
sx q[1];
rz(-2.7014616) q[1];
x q[2];
rz(2.8815202) q[3];
sx q[3];
rz(-1.4315769) q[3];
sx q[3];
rz(1.8782488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.36025563) q[2];
sx q[2];
rz(-0.90136734) q[2];
sx q[2];
rz(-2.0085013) q[2];
rz(-1.1059149) q[3];
sx q[3];
rz(-1.7224885) q[3];
sx q[3];
rz(0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7853506) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(1.6957977) q[0];
rz(-0.71714199) q[1];
sx q[1];
rz(-1.8627867) q[1];
sx q[1];
rz(0.76146567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324156) q[0];
sx q[0];
rz(-1.7178365) q[0];
sx q[0];
rz(2.5086286) q[0];
rz(-pi) q[1];
rz(-2.285073) q[2];
sx q[2];
rz(-1.1978784) q[2];
sx q[2];
rz(0.97717092) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95132212) q[1];
sx q[1];
rz(-2.1441064) q[1];
sx q[1];
rz(-0.058752937) q[1];
rz(-pi) q[2];
rz(-2.2753115) q[3];
sx q[3];
rz(-0.24476642) q[3];
sx q[3];
rz(-2.1921128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7722499) q[2];
sx q[2];
rz(-0.49392924) q[2];
sx q[2];
rz(0.93144766) q[2];
rz(-0.61819589) q[3];
sx q[3];
rz(-1.6341011) q[3];
sx q[3];
rz(-0.96562323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6922927) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(1.3344673) q[0];
rz(2.6284699) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(1.2293053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216825) q[0];
sx q[0];
rz(-1.1510885) q[0];
sx q[0];
rz(-0.76292636) q[0];
rz(-3.0487829) q[2];
sx q[2];
rz(-1.3121737) q[2];
sx q[2];
rz(-1.5802204) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8905764) q[1];
sx q[1];
rz(-2.9843247) q[1];
sx q[1];
rz(3.0096439) q[1];
x q[2];
rz(-1.5531814) q[3];
sx q[3];
rz(-2.2329139) q[3];
sx q[3];
rz(0.78126844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3076155) q[2];
sx q[2];
rz(-2.6145085) q[2];
sx q[2];
rz(2.9311467) q[2];
rz(3.0757507) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(0.79894799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54302067) q[0];
sx q[0];
rz(-2.4871171) q[0];
sx q[0];
rz(1.8684813) q[0];
rz(-1.8674564) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(-2.2575016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9719111) q[0];
sx q[0];
rz(-0.95392862) q[0];
sx q[0];
rz(0.00029460987) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8294677) q[2];
sx q[2];
rz(-2.0805367) q[2];
sx q[2];
rz(2.3296251) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6940004) q[1];
sx q[1];
rz(-1.2939556) q[1];
sx q[1];
rz(2.0029699) q[1];
x q[2];
rz(-2.4482449) q[3];
sx q[3];
rz(-0.39284387) q[3];
sx q[3];
rz(-1.0885914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1511496) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(1.3746064) q[2];
rz(0.19045842) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16354887) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(-1.8769886) q[0];
rz(0.073607445) q[1];
sx q[1];
rz(-1.0818447) q[1];
sx q[1];
rz(-0.61990613) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8148851) q[0];
sx q[0];
rz(-2.305079) q[0];
sx q[0];
rz(2.919038) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5307199) q[2];
sx q[2];
rz(-1.9985285) q[2];
sx q[2];
rz(1.7868702) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1312841) q[1];
sx q[1];
rz(-0.29228739) q[1];
sx q[1];
rz(0.75914219) q[1];
rz(-pi) q[2];
rz(-0.4510433) q[3];
sx q[3];
rz(-2.6313734) q[3];
sx q[3];
rz(3.1376145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2747043) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(2.9956024) q[2];
rz(-2.919096) q[3];
sx q[3];
rz(-2.9166418) q[3];
sx q[3];
rz(2.4116624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9324026) q[0];
sx q[0];
rz(-1.9575735) q[0];
sx q[0];
rz(0.14458543) q[0];
rz(-1.2296386) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(-0.59496224) q[2];
sx q[2];
rz(-1.3887726) q[2];
sx q[2];
rz(-0.61979823) q[2];
rz(-1.360523) q[3];
sx q[3];
rz(-1.539264) q[3];
sx q[3];
rz(2.8217583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
