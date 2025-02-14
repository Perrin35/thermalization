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
rz(2.565413) q[0];
sx q[0];
rz(-1.6771069) q[0];
sx q[0];
rz(0.65281868) q[0];
rz(2.6989812) q[1];
sx q[1];
rz(-1.1224597) q[1];
sx q[1];
rz(0.62243661) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1718423) q[0];
sx q[0];
rz(-1.6873708) q[0];
sx q[0];
rz(1.8844834) q[0];
rz(0.60127778) q[2];
sx q[2];
rz(-1.8493422) q[2];
sx q[2];
rz(0.33180607) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6316526) q[1];
sx q[1];
rz(-1.7544244) q[1];
sx q[1];
rz(-2.0765523) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7951598) q[3];
sx q[3];
rz(-2.5208559) q[3];
sx q[3];
rz(-0.94514314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5114674) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(2.5959065) q[2];
rz(2.4557579) q[3];
sx q[3];
rz(-2.4617709) q[3];
sx q[3];
rz(2.1849476) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086455258) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(-1.0954683) q[0];
rz(2.144004) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(-1.197061) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67162469) q[0];
sx q[0];
rz(-1.3177682) q[0];
sx q[0];
rz(2.1718911) q[0];
x q[1];
rz(0.62792553) q[2];
sx q[2];
rz(-0.88191477) q[2];
sx q[2];
rz(0.21180001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6857288) q[1];
sx q[1];
rz(-1.2867754) q[1];
sx q[1];
rz(-1.0953254) q[1];
rz(1.3838816) q[3];
sx q[3];
rz(-1.8946365) q[3];
sx q[3];
rz(-1.8139386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7257488) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(1.3843298) q[2];
rz(-1.7978801) q[3];
sx q[3];
rz(-1.5680771) q[3];
sx q[3];
rz(1.869092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62190732) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(0.4271048) q[0];
rz(-1.9097795) q[1];
sx q[1];
rz(-1.4991263) q[1];
sx q[1];
rz(-2.5274091) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2285143) q[0];
sx q[0];
rz(-0.43134637) q[0];
sx q[0];
rz(2.449663) q[0];
rz(-pi) q[1];
rz(0.74684192) q[2];
sx q[2];
rz(-1.421442) q[2];
sx q[2];
rz(0.17490444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78313561) q[1];
sx q[1];
rz(-1.6880717) q[1];
sx q[1];
rz(-3.0774119) q[1];
rz(-0.50602405) q[3];
sx q[3];
rz(-1.3576686) q[3];
sx q[3];
rz(1.0552849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4517639) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(0.27099398) q[2];
rz(-2.6594035) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(-1.2838001) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2778306) q[0];
sx q[0];
rz(-1.8067124) q[0];
sx q[0];
rz(0.41896391) q[0];
rz(-0.22467443) q[1];
sx q[1];
rz(-1.2089665) q[1];
sx q[1];
rz(-3.0640501) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6200752) q[0];
sx q[0];
rz(-1.2780315) q[0];
sx q[0];
rz(-1.0241246) q[0];
rz(-pi) q[1];
rz(-2.3691142) q[2];
sx q[2];
rz(-1.8260406) q[2];
sx q[2];
rz(0.30568631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78752199) q[1];
sx q[1];
rz(-1.4284627) q[1];
sx q[1];
rz(-2.700129) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7794533) q[3];
sx q[3];
rz(-1.28445) q[3];
sx q[3];
rz(1.3680259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0858687) q[2];
sx q[2];
rz(-2.8102977) q[2];
sx q[2];
rz(-1.2223988) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.2374249) q[3];
sx q[3];
rz(0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772575) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(-1.3539535) q[0];
rz(2.1063781) q[1];
sx q[1];
rz(-1.2151006) q[1];
sx q[1];
rz(-1.530102) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28883067) q[0];
sx q[0];
rz(-1.0438598) q[0];
sx q[0];
rz(2.2599392) q[0];
x q[1];
rz(-2.8132854) q[2];
sx q[2];
rz(-1.4151148) q[2];
sx q[2];
rz(2.1552483) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93188796) q[1];
sx q[1];
rz(-1.5173638) q[1];
sx q[1];
rz(2.7393952) q[1];
rz(1.6282998) q[3];
sx q[3];
rz(-1.981145) q[3];
sx q[3];
rz(1.8156605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3349541) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(-0.063892603) q[2];
rz(-2.4334) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(-1.3341058) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.719249) q[0];
sx q[0];
rz(-0.023119211) q[0];
sx q[0];
rz(2.0656021) q[0];
rz(3.0883582) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(-2.7630973) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94476766) q[0];
sx q[0];
rz(-1.7009957) q[0];
sx q[0];
rz(-2.8848335) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6588062) q[2];
sx q[2];
rz(-1.529379) q[2];
sx q[2];
rz(-0.27829042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82847257) q[1];
sx q[1];
rz(-2.8250029) q[1];
sx q[1];
rz(-0.60917882) q[1];
x q[2];
rz(0.76734425) q[3];
sx q[3];
rz(-1.9986614) q[3];
sx q[3];
rz(1.8024753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1899015) q[2];
sx q[2];
rz(-1.6804164) q[2];
sx q[2];
rz(-2.107673) q[2];
rz(0.90384358) q[3];
sx q[3];
rz(-2.240182) q[3];
sx q[3];
rz(-0.044053642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0251625) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(-0.12475363) q[1];
sx q[1];
rz(-1.1589103) q[1];
sx q[1];
rz(2.6483026) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83324821) q[0];
sx q[0];
rz(-1.253729) q[0];
sx q[0];
rz(2.8494016) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9134995) q[2];
sx q[2];
rz(-1.5950632) q[2];
sx q[2];
rz(-0.14957854) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1715162) q[1];
sx q[1];
rz(-1.0708717) q[1];
sx q[1];
rz(-2.058567) q[1];
rz(-2.7760963) q[3];
sx q[3];
rz(-2.8474244) q[3];
sx q[3];
rz(-2.1902167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.3844246) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(-1.1473131) q[2];
rz(2.3186963) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(2.4790922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61854521) q[0];
sx q[0];
rz(-1.2724027) q[0];
sx q[0];
rz(-0.4185032) q[0];
rz(3.0283527) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(-1.07771) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0200054) q[0];
sx q[0];
rz(-1.5527871) q[0];
sx q[0];
rz(-1.2822582) q[0];
x q[1];
rz(-0.79506008) q[2];
sx q[2];
rz(-1.803298) q[2];
sx q[2];
rz(0.40337053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5716963) q[1];
sx q[1];
rz(-1.2087617) q[1];
sx q[1];
rz(-3.0891332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4850278) q[3];
sx q[3];
rz(-1.5235008) q[3];
sx q[3];
rz(-1.8790203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8354127) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(-1.4554679) q[2];
rz(2.9742187) q[3];
sx q[3];
rz(-2.6264329) q[3];
sx q[3];
rz(-2.0696056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48285943) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(0.40153781) q[0];
rz(-0.43831476) q[1];
sx q[1];
rz(-0.79151789) q[1];
sx q[1];
rz(1.4567136) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7244723) q[0];
sx q[0];
rz(-3.1062685) q[0];
sx q[0];
rz(1.2958584) q[0];
rz(-pi) q[1];
rz(2.8598665) q[2];
sx q[2];
rz(-2.4664219) q[2];
sx q[2];
rz(3.0670241) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6746725) q[1];
sx q[1];
rz(-2.2385257) q[1];
sx q[1];
rz(-2.2064462) q[1];
rz(2.4026186) q[3];
sx q[3];
rz(-2.4904076) q[3];
sx q[3];
rz(2.2923451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7598286) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(0.71167243) q[2];
rz(-0.034218637) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(1.155352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782368) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(-0.79886287) q[0];
rz(2.0955739) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(-0.23652133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72457658) q[0];
sx q[0];
rz(-1.5931221) q[0];
sx q[0];
rz(1.6315559) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78153083) q[2];
sx q[2];
rz(-0.81630125) q[2];
sx q[2];
rz(0.81354173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7812209) q[1];
sx q[1];
rz(-2.0016333) q[1];
sx q[1];
rz(1.9442417) q[1];
x q[2];
rz(1.4528689) q[3];
sx q[3];
rz(-0.8632568) q[3];
sx q[3];
rz(-2.4411069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.11655) q[2];
sx q[2];
rz(-0.5566842) q[2];
sx q[2];
rz(0.64209783) q[2];
rz(-2.5721278) q[3];
sx q[3];
rz(-0.48782188) q[3];
sx q[3];
rz(-0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4702598) q[0];
sx q[0];
rz(-1.8076121) q[0];
sx q[0];
rz(-2.1847834) q[0];
rz(1.9500465) q[1];
sx q[1];
rz(-1.5622495) q[1];
sx q[1];
rz(-1.539485) q[1];
rz(0.12267648) q[2];
sx q[2];
rz(-0.79570607) q[2];
sx q[2];
rz(-0.88655587) q[2];
rz(0.33470086) q[3];
sx q[3];
rz(-1.1664433) q[3];
sx q[3];
rz(0.56952624) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
