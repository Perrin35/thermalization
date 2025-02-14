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
rz(-2.488774) q[0];
rz(-0.44261143) q[1];
sx q[1];
rz(-2.0191329) q[1];
sx q[1];
rz(2.519156) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25678062) q[0];
sx q[0];
rz(-2.8076165) q[0];
sx q[0];
rz(-1.9335174) q[0];
rz(0.4680674) q[2];
sx q[2];
rz(-2.4862346) q[2];
sx q[2];
rz(1.5214024) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5099401) q[1];
sx q[1];
rz(-1.7544244) q[1];
sx q[1];
rz(-1.0650403) q[1];
rz(-pi) q[2];
x q[2];
rz(0.962019) q[3];
sx q[3];
rz(-1.441027) q[3];
sx q[3];
rz(-0.4421086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6301253) q[2];
sx q[2];
rz(-1.3518535) q[2];
sx q[2];
rz(2.5959065) q[2];
rz(-2.4557579) q[3];
sx q[3];
rz(-0.67982173) q[3];
sx q[3];
rz(2.1849476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0551374) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(2.0461244) q[0];
rz(0.99758863) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(-1.9445317) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2491978) q[0];
sx q[0];
rz(-0.64606842) q[0];
sx q[0];
rz(-1.9996253) q[0];
x q[1];
rz(-2.1904693) q[2];
sx q[2];
rz(-0.89604267) q[2];
sx q[2];
rz(2.5017966) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.45586387) q[1];
sx q[1];
rz(-1.8548172) q[1];
sx q[1];
rz(2.0462672) q[1];
x q[2];
rz(1.7577111) q[3];
sx q[3];
rz(-1.2469562) q[3];
sx q[3];
rz(-1.8139386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7257488) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(-1.3843298) q[2];
rz(-1.7978801) q[3];
sx q[3];
rz(-1.5680771) q[3];
sx q[3];
rz(-1.2725007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62190732) q[0];
sx q[0];
rz(-0.78751957) q[0];
sx q[0];
rz(-2.7144879) q[0];
rz(1.2318132) q[1];
sx q[1];
rz(-1.4991263) q[1];
sx q[1];
rz(-2.5274091) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2285143) q[0];
sx q[0];
rz(-2.7102463) q[0];
sx q[0];
rz(-2.449663) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7730447) q[2];
sx q[2];
rz(-2.3073811) q[2];
sx q[2];
rz(1.6088161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.358457) q[1];
sx q[1];
rz(-1.453521) q[1];
sx q[1];
rz(-3.0774119) q[1];
x q[2];
rz(0.41992979) q[3];
sx q[3];
rz(-2.5961317) q[3];
sx q[3];
rz(-0.88014102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4517639) q[2];
sx q[2];
rz(-2.557705) q[2];
sx q[2];
rz(-0.27099398) q[2];
rz(-2.6594035) q[3];
sx q[3];
rz(-0.8420344) q[3];
sx q[3];
rz(1.8577925) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2778306) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(-2.7226287) q[0];
rz(2.9169182) q[1];
sx q[1];
rz(-1.2089665) q[1];
sx q[1];
rz(0.077542543) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6200752) q[0];
sx q[0];
rz(-1.2780315) q[0];
sx q[0];
rz(2.1174681) q[0];
rz(0.77247844) q[2];
sx q[2];
rz(-1.8260406) q[2];
sx q[2];
rz(-2.8359063) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.78752199) q[1];
sx q[1];
rz(-1.4284627) q[1];
sx q[1];
rz(0.44146367) q[1];
rz(-pi) q[2];
rz(2.5285012) q[3];
sx q[3];
rz(-0.35260751) q[3];
sx q[3];
rz(-2.4168454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0858687) q[2];
sx q[2];
rz(-2.8102977) q[2];
sx q[2];
rz(-1.9191939) q[2];
rz(2.870765) q[3];
sx q[3];
rz(-1.2374249) q[3];
sx q[3];
rz(0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772575) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(1.7876392) q[0];
rz(-1.0352146) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.6114906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2524231) q[0];
sx q[0];
rz(-0.98888654) q[0];
sx q[0];
rz(-2.4956368) q[0];
rz(-0.45299977) q[2];
sx q[2];
rz(-0.36213798) q[2];
sx q[2];
rz(2.1299794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2097047) q[1];
sx q[1];
rz(-1.5173638) q[1];
sx q[1];
rz(0.40219743) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7306385) q[3];
sx q[3];
rz(-1.5180713) q[3];
sx q[3];
rz(2.8737674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80663854) q[2];
sx q[2];
rz(-0.97794509) q[2];
sx q[2];
rz(0.063892603) q[2];
rz(-0.70819267) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(-1.8074869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42234364) q[0];
sx q[0];
rz(-3.1184734) q[0];
sx q[0];
rz(1.0759906) q[0];
rz(-3.0883582) q[1];
sx q[1];
rz(-2.2884171) q[1];
sx q[1];
rz(2.7630973) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814893) q[0];
sx q[0];
rz(-1.8253339) q[0];
sx q[0];
rz(-1.7053563) q[0];
rz(-0.041578023) q[2];
sx q[2];
rz(-1.6587306) q[2];
sx q[2];
rz(1.2961594) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9847451) q[1];
sx q[1];
rz(-1.3917006) q[1];
sx q[1];
rz(-2.8791134) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5603797) q[3];
sx q[3];
rz(-0.85678116) q[3];
sx q[3];
rz(-2.9670144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9516912) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(1.0339197) q[2];
rz(-0.90384358) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(-0.044053642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164301) q[0];
sx q[0];
rz(-1.4465541) q[0];
sx q[0];
rz(0.59301162) q[0];
rz(-0.12475363) q[1];
sx q[1];
rz(-1.9826823) q[1];
sx q[1];
rz(-2.6483026) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3083444) q[0];
sx q[0];
rz(-1.8878636) q[0];
sx q[0];
rz(-2.8494016) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9134995) q[2];
sx q[2];
rz(-1.5950632) q[2];
sx q[2];
rz(0.14957854) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.806099) q[1];
sx q[1];
rz(-0.68365288) q[1];
sx q[1];
rz(0.70913507) q[1];
x q[2];
rz(-0.36549632) q[3];
sx q[3];
rz(-0.29416829) q[3];
sx q[3];
rz(0.95137596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7571681) q[2];
sx q[2];
rz(-1.1606471) q[2];
sx q[2];
rz(-1.1473131) q[2];
rz(-0.82289639) q[3];
sx q[3];
rz(-1.9673248) q[3];
sx q[3];
rz(-0.66250044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61854521) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(2.7230895) q[0];
rz(3.0283527) q[1];
sx q[1];
rz(-1.4694045) q[1];
sx q[1];
rz(1.07771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5961469) q[0];
sx q[0];
rz(-1.8592863) q[0];
sx q[0];
rz(-0.01878564) q[0];
rz(0.79506008) q[2];
sx q[2];
rz(-1.3382947) q[2];
sx q[2];
rz(0.40337053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.019494836) q[1];
sx q[1];
rz(-1.6198525) q[1];
sx q[1];
rz(8*pi/13) q[1];
x q[2];
rz(1.5111132) q[3];
sx q[3];
rz(-2.2264997) q[3];
sx q[3];
rz(0.34464097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.30617994) q[2];
sx q[2];
rz(-2.4330008) q[2];
sx q[2];
rz(1.4554679) q[2];
rz(-2.9742187) q[3];
sx q[3];
rz(-2.6264329) q[3];
sx q[3];
rz(2.0696056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587332) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(2.7400548) q[0];
rz(-2.7032779) q[1];
sx q[1];
rz(-0.79151789) q[1];
sx q[1];
rz(1.6848791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4171203) q[0];
sx q[0];
rz(-0.03532413) q[0];
sx q[0];
rz(-1.2958584) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4859547) q[2];
sx q[2];
rz(-1.7454503) q[2];
sx q[2];
rz(1.2740335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6091793) q[1];
sx q[1];
rz(-2.0557771) q[1];
sx q[1];
rz(0.77528016) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4026186) q[3];
sx q[3];
rz(-2.4904076) q[3];
sx q[3];
rz(0.84924752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38176408) q[2];
sx q[2];
rz(-0.20918736) q[2];
sx q[2];
rz(-2.4299202) q[2];
rz(-3.107374) q[3];
sx q[3];
rz(-0.70845571) q[3];
sx q[3];
rz(1.155352) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.782368) q[0];
sx q[0];
rz(-1.490626) q[0];
sx q[0];
rz(-0.79886287) q[0];
rz(-1.0460188) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(-0.23652133) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4170161) q[0];
sx q[0];
rz(-1.5931221) q[0];
sx q[0];
rz(1.6315559) q[0];
rz(-pi) q[1];
rz(-0.92774074) q[2];
sx q[2];
rz(-2.1143713) q[2];
sx q[2];
rz(1.3613995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0276841) q[1];
sx q[1];
rz(-2.5792173) q[1];
sx q[1];
rz(2.4706868) q[1];
x q[2];
rz(-3.0048851) q[3];
sx q[3];
rz(-2.4259704) q[3];
sx q[3];
rz(2.2608044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.11655) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(2.4994948) q[2];
rz(-0.56946483) q[3];
sx q[3];
rz(-0.48782188) q[3];
sx q[3];
rz(-3.0552982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4702598) q[0];
sx q[0];
rz(-1.8076121) q[0];
sx q[0];
rz(-2.1847834) q[0];
rz(-1.9500465) q[1];
sx q[1];
rz(-1.5793431) q[1];
sx q[1];
rz(1.6021077) q[1];
rz(1.446522) q[2];
sx q[2];
rz(-0.7827324) q[2];
sx q[2];
rz(2.0806349) q[2];
rz(0.91611923) q[3];
sx q[3];
rz(-0.51894938) q[3];
sx q[3];
rz(2.9872158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
