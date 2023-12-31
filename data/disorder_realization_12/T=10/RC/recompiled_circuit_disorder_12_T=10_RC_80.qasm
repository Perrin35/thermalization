OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(1.4847423) q[0];
rz(4.3447189) q[1];
sx q[1];
rz(0.523518) q[1];
sx q[1];
rz(8.5365774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523525) q[0];
sx q[0];
rz(-1.7416735) q[0];
sx q[0];
rz(-1.6160374) q[0];
x q[1];
rz(1.0339123) q[2];
sx q[2];
rz(-1.7134588) q[2];
sx q[2];
rz(-2.2906274) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8468195) q[1];
sx q[1];
rz(-1.5702489) q[1];
sx q[1];
rz(0.16695395) q[1];
rz(-pi) q[2];
rz(0.039460823) q[3];
sx q[3];
rz(-2.4426115) q[3];
sx q[3];
rz(-2.3208502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730826) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(-0.65482393) q[0];
rz(-1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-3.0156946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7315642) q[0];
sx q[0];
rz(-2.0429987) q[0];
sx q[0];
rz(2.7239951) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89181487) q[2];
sx q[2];
rz(-2.0795155) q[2];
sx q[2];
rz(1.5430792) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1344578) q[1];
sx q[1];
rz(-1.4440698) q[1];
sx q[1];
rz(-3.111372) q[1];
rz(2.294299) q[3];
sx q[3];
rz(-1.5627268) q[3];
sx q[3];
rz(1.9891667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(2.753567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5858784) q[0];
sx q[0];
rz(-0.782574) q[0];
sx q[0];
rz(0.079332381) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(-1.1598587) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9281611) q[0];
sx q[0];
rz(-2.0800989) q[0];
sx q[0];
rz(0.094497735) q[0];
rz(-0.50767501) q[2];
sx q[2];
rz(-2.1360364) q[2];
sx q[2];
rz(0.6914247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68418903) q[1];
sx q[1];
rz(-2.0261129) q[1];
sx q[1];
rz(-1.1475569) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0852774) q[3];
sx q[3];
rz(-1.0259797) q[3];
sx q[3];
rz(-2.44851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7362061) q[2];
sx q[2];
rz(-1.8182886) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(-2.4978499) q[3];
sx q[3];
rz(-1.0780004) q[3];
sx q[3];
rz(-0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9765587) q[0];
sx q[0];
rz(-1.7465916) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(2.4404793) q[1];
sx q[1];
rz(-1.8891524) q[1];
sx q[1];
rz(-0.28465095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92589256) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(-1.0924933) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20547159) q[2];
sx q[2];
rz(-1.2887508) q[2];
sx q[2];
rz(-0.13946433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2440566) q[1];
sx q[1];
rz(-1.4923864) q[1];
sx q[1];
rz(1.4439911) q[1];
x q[2];
rz(-2.1553667) q[3];
sx q[3];
rz(-0.75891906) q[3];
sx q[3];
rz(-2.330247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9099137) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(-0.56751928) q[2];
rz(-2.7275758) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48859566) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(-1.3265142) q[0];
rz(-1.9891706) q[1];
sx q[1];
rz(-1.3782586) q[1];
sx q[1];
rz(-0.93793905) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2269103) q[0];
sx q[0];
rz(-1.662928) q[0];
sx q[0];
rz(3.1196306) q[0];
rz(-pi) q[1];
rz(-2.6229834) q[2];
sx q[2];
rz(-1.0728288) q[2];
sx q[2];
rz(-1.2155611) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.282498) q[1];
sx q[1];
rz(-1.5364093) q[1];
sx q[1];
rz(0.026245898) q[1];
rz(-pi) q[2];
rz(1.7004847) q[3];
sx q[3];
rz(-1.218759) q[3];
sx q[3];
rz(-0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-1.1338736) q[3];
sx q[3];
rz(-1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(-2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(2.5851137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0628478) q[0];
sx q[0];
rz(-1.1947462) q[0];
sx q[0];
rz(1.7929121) q[0];
x q[1];
rz(0.17369341) q[2];
sx q[2];
rz(-0.86280338) q[2];
sx q[2];
rz(-3.0922572) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1408128) q[1];
sx q[1];
rz(-1.8218578) q[1];
sx q[1];
rz(-2.954133) q[1];
rz(-1.6218833) q[3];
sx q[3];
rz(-2.3458977) q[3];
sx q[3];
rz(1.5339799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0625793) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(3.0986837) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(2.6409805) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6734877) q[0];
sx q[0];
rz(-1.7485432) q[0];
sx q[0];
rz(0.034981473) q[0];
rz(-pi) q[1];
rz(0.31502864) q[2];
sx q[2];
rz(-2.0909967) q[2];
sx q[2];
rz(0.5859642) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4088879) q[1];
sx q[1];
rz(-1.9703456) q[1];
sx q[1];
rz(-1.7016181) q[1];
rz(-pi) q[2];
x q[2];
rz(2.385419) q[3];
sx q[3];
rz(-1.2740967) q[3];
sx q[3];
rz(2.365436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(0.67561692) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(2.0741529) q[0];
rz(-0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(1.4656461) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059793652) q[0];
sx q[0];
rz(-1.1104715) q[0];
sx q[0];
rz(2.9914809) q[0];
x q[1];
rz(-2.9550214) q[2];
sx q[2];
rz(-1.4636883) q[2];
sx q[2];
rz(-2.4466116) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9626179) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(0.38260539) q[1];
rz(-pi) q[2];
rz(2.1727209) q[3];
sx q[3];
rz(-2.1224009) q[3];
sx q[3];
rz(-1.5495891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33891588) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.476293) q[2];
rz(2.2680797) q[3];
sx q[3];
rz(-0.87681186) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2575689) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-0.73721686) q[0];
rz(0.018741477) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(-2.2163056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5697445) q[0];
sx q[0];
rz(-2.0664555) q[0];
sx q[0];
rz(-2.3832688) q[0];
x q[1];
rz(-2.4783752) q[2];
sx q[2];
rz(-1.466201) q[2];
sx q[2];
rz(1.6246206) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9463897) q[1];
sx q[1];
rz(-0.22559799) q[1];
sx q[1];
rz(-0.619508) q[1];
x q[2];
rz(-0.85261811) q[3];
sx q[3];
rz(-1.8947621) q[3];
sx q[3];
rz(2.1924803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(0.090099661) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1019679) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(-1.8819303) q[0];
rz(-2.8758077) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(-2.3840747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.789061) q[0];
sx q[0];
rz(-1.7081982) q[0];
sx q[0];
rz(1.765541) q[0];
rz(-1.8025814) q[2];
sx q[2];
rz(-0.60283649) q[2];
sx q[2];
rz(-1.1834809) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.47093686) q[1];
sx q[1];
rz(-0.48479143) q[1];
sx q[1];
rz(-1.0792586) q[1];
x q[2];
rz(3.1086139) q[3];
sx q[3];
rz(-2.4313201) q[3];
sx q[3];
rz(-0.23746333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-2.5027067) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(2.1209774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(-1.5007301) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(0.43754229) q[2];
sx q[2];
rz(-1.2699288) q[2];
sx q[2];
rz(-1.5581836) q[2];
rz(-0.32140857) q[3];
sx q[3];
rz(-2.1264429) q[3];
sx q[3];
rz(0.85254729) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
