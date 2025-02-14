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
rz(0.49864545) q[0];
sx q[0];
rz(3.6874229) q[0];
sx q[0];
rz(9.3064718) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(4.4448648) q[1];
sx q[1];
rz(8.2104609) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34284232) q[0];
sx q[0];
rz(-1.3590706) q[0];
sx q[0];
rz(-2.0131074) q[0];
x q[1];
rz(0.6251752) q[2];
sx q[2];
rz(-2.0707651) q[2];
sx q[2];
rz(-1.878405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2976297) q[1];
sx q[1];
rz(-1.7073892) q[1];
sx q[1];
rz(-0.94874391) q[1];
x q[2];
rz(-3.1200527) q[3];
sx q[3];
rz(-1.2976754) q[3];
sx q[3];
rz(2.5813951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.168657) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(1.0966148) q[2];
rz(0.24122572) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(2.296804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0519401) q[0];
sx q[0];
rz(-1.0900499) q[0];
sx q[0];
rz(0.3845149) q[0];
rz(-0.3081201) q[1];
sx q[1];
rz(-0.48079753) q[1];
sx q[1];
rz(0.62517977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70474941) q[0];
sx q[0];
rz(-1.539121) q[0];
sx q[0];
rz(-0.28625536) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4433631) q[2];
sx q[2];
rz(-2.2837167) q[2];
sx q[2];
rz(-2.8486203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58142606) q[1];
sx q[1];
rz(-2.4871965) q[1];
sx q[1];
rz(0.1398211) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2143638) q[3];
sx q[3];
rz(-0.72690287) q[3];
sx q[3];
rz(-0.34360269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.008931) q[2];
sx q[2];
rz(-2.2098358) q[2];
sx q[2];
rz(2.8399732) q[2];
rz(2.610142) q[3];
sx q[3];
rz(-0.88574946) q[3];
sx q[3];
rz(-1.1227932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5136435) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(1.1847786) q[0];
rz(-2.9966677) q[1];
sx q[1];
rz(-1.8791608) q[1];
sx q[1];
rz(-1.5473993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14159053) q[0];
sx q[0];
rz(-0.98391082) q[0];
sx q[0];
rz(-0.97372719) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0058262) q[2];
sx q[2];
rz(-2.7270881) q[2];
sx q[2];
rz(2.5786825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6559534) q[1];
sx q[1];
rz(-1.9536634) q[1];
sx q[1];
rz(-0.98321557) q[1];
rz(-0.66158847) q[3];
sx q[3];
rz(-2.6433655) q[3];
sx q[3];
rz(0.029904043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80921119) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(1.3261718) q[2];
rz(-2.2281846) q[3];
sx q[3];
rz(-1.2181506) q[3];
sx q[3];
rz(-2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42295414) q[0];
sx q[0];
rz(-2.664743) q[0];
sx q[0];
rz(-1.3612716) q[0];
rz(-0.71296972) q[1];
sx q[1];
rz(-1.1402592) q[1];
sx q[1];
rz(2.4580809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26091247) q[0];
sx q[0];
rz(-0.99755462) q[0];
sx q[0];
rz(-0.011996847) q[0];
x q[1];
rz(0.81869469) q[2];
sx q[2];
rz(-2.8398189) q[2];
sx q[2];
rz(1.8723328) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5665413) q[1];
sx q[1];
rz(-2.0861826) q[1];
sx q[1];
rz(-2.383286) q[1];
rz(-2.7587851) q[3];
sx q[3];
rz(-0.26255739) q[3];
sx q[3];
rz(-2.5298339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6354562) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(-0.48640856) q[2];
rz(-2.6247315) q[3];
sx q[3];
rz(-1.1812527) q[3];
sx q[3];
rz(2.1080871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0559693) q[0];
sx q[0];
rz(-1.6317246) q[0];
sx q[0];
rz(-1.3729209) q[0];
rz(1.4068475) q[1];
sx q[1];
rz(-0.56766784) q[1];
sx q[1];
rz(2.6645606) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45390391) q[0];
sx q[0];
rz(-2.6756508) q[0];
sx q[0];
rz(2.0201319) q[0];
rz(2.1978756) q[2];
sx q[2];
rz(-0.6155799) q[2];
sx q[2];
rz(-1.0301539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7030695) q[1];
sx q[1];
rz(-1.0756936) q[1];
sx q[1];
rz(-2.7994521) q[1];
rz(2.2955016) q[3];
sx q[3];
rz(-0.93069221) q[3];
sx q[3];
rz(2.7654331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(-0.46727115) q[2];
rz(-2.9723736) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(0.30538487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071103) q[0];
sx q[0];
rz(-2.2619673) q[0];
sx q[0];
rz(-2.95209) q[0];
rz(-0.27319187) q[1];
sx q[1];
rz(-2.0676282) q[1];
sx q[1];
rz(2.3409519) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439745) q[0];
sx q[0];
rz(-1.6139784) q[0];
sx q[0];
rz(1.8042685) q[0];
x q[1];
rz(2.8932299) q[2];
sx q[2];
rz(-1.4830605) q[2];
sx q[2];
rz(1.9875634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.86294293) q[1];
sx q[1];
rz(-2.9065296) q[1];
sx q[1];
rz(-2.1549015) q[1];
rz(2.6012035) q[3];
sx q[3];
rz(-2.0107993) q[3];
sx q[3];
rz(0.42127452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2022986) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(-2.6046216) q[2];
rz(1.8592853) q[3];
sx q[3];
rz(-1.8351646) q[3];
sx q[3];
rz(-1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.110431) q[0];
sx q[0];
rz(-0.81542504) q[0];
sx q[0];
rz(1.331331) q[0];
rz(1.4415461) q[1];
sx q[1];
rz(-1.4794289) q[1];
sx q[1];
rz(1.2225245) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6922369) q[0];
sx q[0];
rz(-0.48654443) q[0];
sx q[0];
rz(-1.0445717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4961081) q[2];
sx q[2];
rz(-1.0202346) q[2];
sx q[2];
rz(2.2979743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0646831) q[1];
sx q[1];
rz(-1.7864702) q[1];
sx q[1];
rz(-2.2806858) q[1];
x q[2];
rz(-3.0223373) q[3];
sx q[3];
rz(-1.1290765) q[3];
sx q[3];
rz(-2.2252803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20087251) q[2];
sx q[2];
rz(-2.6437289) q[2];
sx q[2];
rz(-0.3328003) q[2];
rz(0.29916549) q[3];
sx q[3];
rz(-1.2406415) q[3];
sx q[3];
rz(0.021473264) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0126295) q[0];
sx q[0];
rz(-2.3733932) q[0];
sx q[0];
rz(2.8666038) q[0];
rz(-0.081427447) q[1];
sx q[1];
rz(-0.80137253) q[1];
sx q[1];
rz(-1.3224695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0596432) q[0];
sx q[0];
rz(-0.23898757) q[0];
sx q[0];
rz(-1.8247402) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36825387) q[2];
sx q[2];
rz(-2.2241631) q[2];
sx q[2];
rz(-0.39619941) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53912932) q[1];
sx q[1];
rz(-0.88383288) q[1];
sx q[1];
rz(-1.0844356) q[1];
rz(-pi) q[2];
rz(-2.4422936) q[3];
sx q[3];
rz(-0.85996503) q[3];
sx q[3];
rz(-3.1131203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.72121173) q[2];
sx q[2];
rz(-0.12004852) q[2];
sx q[2];
rz(1.4018641) q[2];
rz(-1.2841355) q[3];
sx q[3];
rz(-1.9522342) q[3];
sx q[3];
rz(3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.5818802) q[0];
sx q[0];
rz(-1.548883) q[0];
sx q[0];
rz(2.6522719) q[0];
rz(-3.1210461) q[1];
sx q[1];
rz(-1.9053883) q[1];
sx q[1];
rz(-2.4403341) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72557455) q[0];
sx q[0];
rz(-1.3301992) q[0];
sx q[0];
rz(0.98576905) q[0];
x q[1];
rz(2.2743938) q[2];
sx q[2];
rz(-1.6072306) q[2];
sx q[2];
rz(1.4796204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8119231) q[1];
sx q[1];
rz(-2.3282302) q[1];
sx q[1];
rz(-1.8977036) q[1];
x q[2];
rz(0.55072983) q[3];
sx q[3];
rz(-2.3756873) q[3];
sx q[3];
rz(-2.285241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7265085) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(0.48039594) q[2];
rz(1.7324309) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(0.33717808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3163863) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(-2.7301042) q[0];
rz(-1.2443789) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(2.8172475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8702983) q[0];
sx q[0];
rz(-2.76519) q[0];
sx q[0];
rz(1.7145304) q[0];
rz(1.2191992) q[2];
sx q[2];
rz(-1.2484387) q[2];
sx q[2];
rz(0.11478648) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7601519) q[1];
sx q[1];
rz(-2.5709929) q[1];
sx q[1];
rz(1.4165611) q[1];
rz(-0.12121157) q[3];
sx q[3];
rz(-1.5577661) q[3];
sx q[3];
rz(-0.37076326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1856498) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-0.11697098) q[2];
rz(-2.1864435) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(-0.54168934) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87880001) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(2.0211438) q[1];
sx q[1];
rz(-2.3255377) q[1];
sx q[1];
rz(-1.9768523) q[1];
rz(-1.2016313) q[2];
sx q[2];
rz(-2.6045447) q[2];
sx q[2];
rz(2.4640026) q[2];
rz(1.3001623) q[3];
sx q[3];
rz(-1.5951372) q[3];
sx q[3];
rz(-1.5230877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
