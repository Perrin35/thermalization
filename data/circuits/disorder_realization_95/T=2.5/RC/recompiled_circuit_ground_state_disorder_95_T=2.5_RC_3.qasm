OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9606544) q[0];
sx q[0];
rz(-0.063194312) q[0];
sx q[0];
rz(-0.59046459) q[0];
rz(-0.2904627) q[1];
sx q[1];
rz(-1.3757179) q[1];
sx q[1];
rz(0.003493316) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0079903) q[0];
sx q[0];
rz(-0.65235814) q[0];
sx q[0];
rz(-0.86835394) q[0];
x q[1];
rz(0.89031808) q[2];
sx q[2];
rz(-1.3290452) q[2];
sx q[2];
rz(1.3554475) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.063979186) q[1];
sx q[1];
rz(-3.1099456) q[1];
sx q[1];
rz(0.74546234) q[1];
x q[2];
rz(0.9736075) q[3];
sx q[3];
rz(-1.191621) q[3];
sx q[3];
rz(1.7969064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2575839) q[2];
sx q[2];
rz(-2.282892) q[2];
sx q[2];
rz(-2.1250471) q[2];
rz(0.86333418) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(-1.6421002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6739552) q[0];
sx q[0];
rz(-1.4731982) q[0];
sx q[0];
rz(-0.13482811) q[0];
rz(0.22678953) q[1];
sx q[1];
rz(-0.062954523) q[1];
sx q[1];
rz(-2.8917868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86827474) q[0];
sx q[0];
rz(-0.76486838) q[0];
sx q[0];
rz(-2.5165783) q[0];
x q[1];
rz(-3.089014) q[2];
sx q[2];
rz(-1.0710395) q[2];
sx q[2];
rz(-1.6249958) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0909101) q[1];
sx q[1];
rz(-1.5600648) q[1];
sx q[1];
rz(-0.077535943) q[1];
rz(1.2050912) q[3];
sx q[3];
rz(-1.4451175) q[3];
sx q[3];
rz(2.2179536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9164385) q[2];
sx q[2];
rz(-3.0735922) q[2];
sx q[2];
rz(-0.98006836) q[2];
rz(0.89503908) q[3];
sx q[3];
rz(-2.3990302) q[3];
sx q[3];
rz(1.9612954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5880244) q[0];
sx q[0];
rz(-0.50694412) q[0];
sx q[0];
rz(1.6059599) q[0];
rz(1.5060679) q[1];
sx q[1];
rz(-0.82553828) q[1];
sx q[1];
rz(-0.95605409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0900537) q[0];
sx q[0];
rz(-2.9057626) q[0];
sx q[0];
rz(-2.7084742) q[0];
rz(3.0676715) q[2];
sx q[2];
rz(-0.73237458) q[2];
sx q[2];
rz(0.074490212) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3967665) q[1];
sx q[1];
rz(-1.6770207) q[1];
sx q[1];
rz(0.83858072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5713137) q[3];
sx q[3];
rz(-0.59123838) q[3];
sx q[3];
rz(1.6275004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8354336) q[2];
sx q[2];
rz(-2.4296438) q[2];
sx q[2];
rz(2.5095059) q[2];
rz(-0.88405526) q[3];
sx q[3];
rz(-0.017280936) q[3];
sx q[3];
rz(2.250905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.66568351) q[0];
sx q[0];
rz(-1.2983687) q[0];
sx q[0];
rz(-2.9744398) q[0];
rz(-1.2627603) q[1];
sx q[1];
rz(-2.1985168) q[1];
sx q[1];
rz(1.4080338) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4564877) q[0];
sx q[0];
rz(-0.52336687) q[0];
sx q[0];
rz(2.531764) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5535545) q[2];
sx q[2];
rz(-0.21176007) q[2];
sx q[2];
rz(-1.0946158) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1442529) q[1];
sx q[1];
rz(-1.3648811) q[1];
sx q[1];
rz(1.3763672) q[1];
rz(0.77045958) q[3];
sx q[3];
rz(-1.3589199) q[3];
sx q[3];
rz(2.6765064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7888646) q[2];
sx q[2];
rz(-3.126324) q[2];
sx q[2];
rz(0.35066476) q[2];
rz(-2.8777425) q[3];
sx q[3];
rz(-3.1407472) q[3];
sx q[3];
rz(-1.9381757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8769787) q[0];
sx q[0];
rz(-2.3961841) q[0];
sx q[0];
rz(-1.718148) q[0];
rz(0.28165948) q[1];
sx q[1];
rz(-1.5080844) q[1];
sx q[1];
rz(1.8219832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21319178) q[0];
sx q[0];
rz(-1.1494779) q[0];
sx q[0];
rz(-2.1405959) q[0];
rz(-pi) q[1];
x q[1];
rz(0.020242029) q[2];
sx q[2];
rz(-1.5432095) q[2];
sx q[2];
rz(2.8351912) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3710821) q[1];
sx q[1];
rz(-0.96643172) q[1];
sx q[1];
rz(1.1649067) q[1];
rz(-pi) q[2];
rz(-1.9044456) q[3];
sx q[3];
rz(-1.9556442) q[3];
sx q[3];
rz(-2.9001482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8772956) q[2];
sx q[2];
rz(-3.0946315) q[2];
sx q[2];
rz(1.4287255) q[2];
rz(-1.9127539) q[3];
sx q[3];
rz(-0.052534025) q[3];
sx q[3];
rz(-0.088168941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0838715) q[0];
sx q[0];
rz(-0.11744048) q[0];
sx q[0];
rz(0.55026662) q[0];
rz(-0.30217198) q[1];
sx q[1];
rz(-1.3928394) q[1];
sx q[1];
rz(-0.43513939) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9254375) q[0];
sx q[0];
rz(-1.7820246) q[0];
sx q[0];
rz(1.5351377) q[0];
x q[1];
rz(1.5698293) q[2];
sx q[2];
rz(-1.5633601) q[2];
sx q[2];
rz(2.815757) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3966538) q[1];
sx q[1];
rz(-2.7135773) q[1];
sx q[1];
rz(-1.6950083) q[1];
x q[2];
rz(-0.72316678) q[3];
sx q[3];
rz(-1.0042913) q[3];
sx q[3];
rz(0.97149937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90007323) q[2];
sx q[2];
rz(-3.1406431) q[2];
sx q[2];
rz(-0.97236902) q[2];
rz(2.1420245) q[3];
sx q[3];
rz(-0.0071503706) q[3];
sx q[3];
rz(1.0424559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480176) q[0];
sx q[0];
rz(-1.8997718) q[0];
sx q[0];
rz(1.0663363) q[0];
rz(1.713133) q[1];
sx q[1];
rz(-2.8722873) q[1];
sx q[1];
rz(1.2880633) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.903164) q[0];
sx q[0];
rz(-2.3580436) q[0];
sx q[0];
rz(-2.6859693) q[0];
rz(-pi) q[1];
rz(-1.1797899) q[2];
sx q[2];
rz(-0.76863747) q[2];
sx q[2];
rz(2.1360508) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3235229) q[1];
sx q[1];
rz(-1.5841636) q[1];
sx q[1];
rz(1.6105777) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24434511) q[3];
sx q[3];
rz(-1.0536031) q[3];
sx q[3];
rz(1.2822267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.11918934) q[2];
sx q[2];
rz(-0.069312118) q[2];
sx q[2];
rz(-2.8663087) q[2];
rz(-0.22661181) q[3];
sx q[3];
rz(-0.35609069) q[3];
sx q[3];
rz(-0.4376469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.3771123) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(2.5544033) q[0];
rz(-1.6181234) q[1];
sx q[1];
rz(-1.1174997) q[1];
sx q[1];
rz(1.5709741) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8507152) q[0];
sx q[0];
rz(-1.1505011) q[0];
sx q[0];
rz(-0.72776003) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4373407) q[2];
sx q[2];
rz(-0.46758662) q[2];
sx q[2];
rz(0.0011006265) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60698673) q[1];
sx q[1];
rz(-3.0879277) q[1];
sx q[1];
rz(0.95602687) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0471025) q[3];
sx q[3];
rz(-2.5137551) q[3];
sx q[3];
rz(1.353457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1945343) q[2];
sx q[2];
rz(-2.9583866) q[2];
sx q[2];
rz(1.3018695) q[2];
rz(2.7070847) q[3];
sx q[3];
rz(-0.040381581) q[3];
sx q[3];
rz(2.6935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.7942363) q[0];
sx q[0];
rz(-0.11581049) q[0];
sx q[0];
rz(1.3809123) q[0];
rz(-1.5877089) q[1];
sx q[1];
rz(-0.96276182) q[1];
sx q[1];
rz(3.0316321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1292353) q[0];
sx q[0];
rz(-1.7596272) q[0];
sx q[0];
rz(-0.83602943) q[0];
rz(-2.8822363) q[2];
sx q[2];
rz(-2.4217941) q[2];
sx q[2];
rz(2.6754745) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.546098) q[1];
sx q[1];
rz(-1.5459041) q[1];
sx q[1];
rz(-0.81374641) q[1];
x q[2];
rz(-0.31153714) q[3];
sx q[3];
rz(-1.7230125) q[3];
sx q[3];
rz(-1.0675586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3515557) q[2];
sx q[2];
rz(-1.0645126) q[2];
sx q[2];
rz(2.7891187) q[2];
rz(-2.500109) q[3];
sx q[3];
rz(-3.0951169) q[3];
sx q[3];
rz(2.7789153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27956692) q[0];
sx q[0];
rz(-2.8643705) q[0];
sx q[0];
rz(-2.5773881) q[0];
rz(-1.5203681) q[1];
sx q[1];
rz(-2.11026) q[1];
sx q[1];
rz(-3.0618111) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.736981) q[0];
sx q[0];
rz(-2.4546162) q[0];
sx q[0];
rz(1.3300598) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3529813) q[2];
sx q[2];
rz(-0.067316003) q[2];
sx q[2];
rz(-2.2612417) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3120739) q[1];
sx q[1];
rz(-0.79527277) q[1];
sx q[1];
rz(1.2088695) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8989927) q[3];
sx q[3];
rz(-1.651847) q[3];
sx q[3];
rz(-1.9417861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0148049) q[2];
sx q[2];
rz(-0.0099651907) q[2];
sx q[2];
rz(1.0753746) q[2];
rz(-0.77438313) q[3];
sx q[3];
rz(-3.1175678) q[3];
sx q[3];
rz(0.27124852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2992235) q[0];
sx q[0];
rz(-1.5689701) q[0];
sx q[0];
rz(-1.5690621) q[0];
rz(0.52539274) q[1];
sx q[1];
rz(-3.0699984) q[1];
sx q[1];
rz(2.933554) q[1];
rz(-1.2425497) q[2];
sx q[2];
rz(-0.59292159) q[2];
sx q[2];
rz(-2.4688724) q[2];
rz(-2.0490859) q[3];
sx q[3];
rz(-0.94259613) q[3];
sx q[3];
rz(-2.5362956) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
