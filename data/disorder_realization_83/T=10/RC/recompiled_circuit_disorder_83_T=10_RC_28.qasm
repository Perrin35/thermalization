OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(2.3798556) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0640472) q[0];
sx q[0];
rz(-0.24928688) q[0];
sx q[0];
rz(-1.3011978) q[0];
x q[1];
rz(0.6526297) q[2];
sx q[2];
rz(-1.666781) q[2];
sx q[2];
rz(-1.3053615) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8929157) q[1];
sx q[1];
rz(-2.6495653) q[1];
sx q[1];
rz(-2.1881275) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7573962) q[3];
sx q[3];
rz(-2.7717234) q[3];
sx q[3];
rz(-0.44934011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(0.10786954) q[2];
rz(-0.14262959) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-0.67255783) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(2.9108677) q[0];
rz(1.327286) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-3.1006295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9757257) q[0];
sx q[0];
rz(-1.3250933) q[0];
sx q[0];
rz(-2.9742572) q[0];
rz(-pi) q[1];
rz(-1.2626921) q[2];
sx q[2];
rz(-1.3034046) q[2];
sx q[2];
rz(1.845713) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.47536182) q[1];
sx q[1];
rz(-2.1365709) q[1];
sx q[1];
rz(-1.0223785) q[1];
x q[2];
rz(1.1702234) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(-2.7924201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.162398) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(2.6290821) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-2.5797243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(0.67726642) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30401858) q[0];
sx q[0];
rz(-1.7390334) q[0];
sx q[0];
rz(1.6234682) q[0];
rz(2.9708423) q[2];
sx q[2];
rz(-0.0064592529) q[2];
sx q[2];
rz(-0.016575459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22563572) q[1];
sx q[1];
rz(-1.3924053) q[1];
sx q[1];
rz(-0.13326463) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1797818) q[3];
sx q[3];
rz(-2.0174332) q[3];
sx q[3];
rz(3.1004268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(0.051483367) q[0];
rz(0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.1725918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6068891) q[0];
sx q[0];
rz(-2.4692315) q[0];
sx q[0];
rz(0.21557233) q[0];
rz(-1.555221) q[2];
sx q[2];
rz(-0.96484631) q[2];
sx q[2];
rz(-1.3941927) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7395775) q[1];
sx q[1];
rz(-0.47164729) q[1];
sx q[1];
rz(-2.2775843) q[1];
rz(-1.7309534) q[3];
sx q[3];
rz(-2.7048116) q[3];
sx q[3];
rz(0.76524759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2525758) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(-1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(3.1161599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(-1.8160965) q[0];
rz(-2.5505113) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7514873) q[0];
sx q[0];
rz(-0.57919466) q[0];
sx q[0];
rz(0.67436995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3833369) q[2];
sx q[2];
rz(-1.1584632) q[2];
sx q[2];
rz(0.86134855) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7655189) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(-2.5942624) q[1];
rz(-pi) q[2];
rz(-1.4773024) q[3];
sx q[3];
rz(-1.9545385) q[3];
sx q[3];
rz(-1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6902265) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(0.29904547) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4371571) q[0];
sx q[0];
rz(-0.91059443) q[0];
sx q[0];
rz(3.0755088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5619179) q[2];
sx q[2];
rz(-1.056676) q[2];
sx q[2];
rz(-0.74619734) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1491579) q[1];
sx q[1];
rz(-0.5545534) q[1];
sx q[1];
rz(2.1402332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3214345) q[3];
sx q[3];
rz(-1.2994088) q[3];
sx q[3];
rz(0.67941487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6254639) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(0.027475474) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(0.97822899) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(-0.79089975) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1066125) q[0];
sx q[0];
rz(-0.62823409) q[0];
sx q[0];
rz(-1.7321222) q[0];
rz(-pi) q[1];
rz(-1.9060471) q[2];
sx q[2];
rz(-0.8555879) q[2];
sx q[2];
rz(-2.7170979) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0485718) q[1];
sx q[1];
rz(-1.934811) q[1];
sx q[1];
rz(-3.0572592) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7780276) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(0.082106575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-0.86501914) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-0.46736026) q[0];
rz(-2.6043747) q[1];
sx q[1];
rz(-2.1610114) q[1];
sx q[1];
rz(2.8875202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8763037) q[0];
sx q[0];
rz(-2.9348001) q[0];
sx q[0];
rz(0.32949038) q[0];
rz(-3.1363856) q[2];
sx q[2];
rz(-1.569283) q[2];
sx q[2];
rz(0.13841275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1797807) q[1];
sx q[1];
rz(-1.4439549) q[1];
sx q[1];
rz(-0.3133874) q[1];
rz(3.0244175) q[3];
sx q[3];
rz(-2.5391038) q[3];
sx q[3];
rz(0.52091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6415928) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083387233) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(2.0876419) q[0];
rz(2.9846233) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(0.98186791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69753416) q[0];
sx q[0];
rz(-2.7148348) q[0];
sx q[0];
rz(2.8696637) q[0];
rz(-1.7392731) q[2];
sx q[2];
rz(-0.57541621) q[2];
sx q[2];
rz(-1.1489431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6135785) q[1];
sx q[1];
rz(-2.947223) q[1];
sx q[1];
rz(-1.2685052) q[1];
rz(-pi) q[2];
rz(0.28717678) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2173569) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(2.1441933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-0.034328073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85703325) q[0];
sx q[0];
rz(-2.4301346) q[0];
sx q[0];
rz(-2.6249264) q[0];
rz(0.38756469) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(-0.26836747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3049406) q[0];
sx q[0];
rz(-2.1450844) q[0];
sx q[0];
rz(0.71026295) q[0];
x q[1];
rz(2.7029413) q[2];
sx q[2];
rz(-1.029656) q[2];
sx q[2];
rz(-1.3676608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8473709) q[1];
sx q[1];
rz(-0.98681824) q[1];
sx q[1];
rz(-1.8933312) q[1];
rz(-0.99276944) q[3];
sx q[3];
rz(-1.7939995) q[3];
sx q[3];
rz(2.2900667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(1.1901723) q[3];
sx q[3];
rz(-0.95919132) q[3];
sx q[3];
rz(0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47678369) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-0.51856002) q[2];
sx q[2];
rz(-1.3052365) q[2];
sx q[2];
rz(-2.4076751) q[2];
rz(-0.56223829) q[3];
sx q[3];
rz(-0.36459618) q[3];
sx q[3];
rz(0.96095745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
