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
rz(-2.0237145) q[0];
sx q[0];
rz(-1.0853381) q[0];
sx q[0];
rz(2.7827652) q[0];
rz(1.2598414) q[1];
sx q[1];
rz(2.0460195) q[1];
sx q[1];
rz(10.464315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0721999) q[0];
sx q[0];
rz(-0.93749267) q[0];
sx q[0];
rz(2.0149489) q[0];
rz(3.1141236) q[2];
sx q[2];
rz(-1.5284163) q[2];
sx q[2];
rz(-0.24951631) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4182836) q[1];
sx q[1];
rz(-1.9305349) q[1];
sx q[1];
rz(0.32029256) q[1];
rz(-pi) q[2];
x q[2];
rz(1.667439) q[3];
sx q[3];
rz(-2.2386239) q[3];
sx q[3];
rz(2.5092595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.871668) q[2];
sx q[2];
rz(-2.2871127) q[2];
sx q[2];
rz(-2.975115) q[2];
rz(-0.13655937) q[3];
sx q[3];
rz(-0.86524335) q[3];
sx q[3];
rz(1.6816007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.0110375) q[0];
sx q[0];
rz(-0.32259536) q[0];
sx q[0];
rz(-2.7815797) q[0];
rz(-0.75045466) q[1];
sx q[1];
rz(-0.89025918) q[1];
sx q[1];
rz(3.0976565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.982807) q[0];
sx q[0];
rz(-1.7097397) q[0];
sx q[0];
rz(-1.7278633) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6360128) q[2];
sx q[2];
rz(-0.64526099) q[2];
sx q[2];
rz(1.1087904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3219731) q[1];
sx q[1];
rz(-0.83998806) q[1];
sx q[1];
rz(1.2801484) q[1];
x q[2];
rz(-2.9437795) q[3];
sx q[3];
rz(-1.3105243) q[3];
sx q[3];
rz(-2.71432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3540078) q[2];
sx q[2];
rz(-0.44469357) q[2];
sx q[2];
rz(-2.1237874) q[2];
rz(2.9338037) q[3];
sx q[3];
rz(-0.60616797) q[3];
sx q[3];
rz(0.35473216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6924234) q[0];
sx q[0];
rz(-2.848931) q[0];
sx q[0];
rz(2.0617275) q[0];
rz(-2.1366468) q[1];
sx q[1];
rz(-0.69098538) q[1];
sx q[1];
rz(1.8963922) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8522135) q[0];
sx q[0];
rz(-2.820008) q[0];
sx q[0];
rz(-2.8711161) q[0];
x q[1];
rz(-0.4246906) q[2];
sx q[2];
rz(-1.7199368) q[2];
sx q[2];
rz(1.5678659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9745404) q[1];
sx q[1];
rz(-0.61950942) q[1];
sx q[1];
rz(1.0654679) q[1];
x q[2];
rz(-1.336634) q[3];
sx q[3];
rz(-2.5076205) q[3];
sx q[3];
rz(-2.0869107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9006742) q[2];
sx q[2];
rz(-1.8697048) q[2];
sx q[2];
rz(2.7980878) q[2];
rz(-0.014723012) q[3];
sx q[3];
rz(-0.69535178) q[3];
sx q[3];
rz(1.9007614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13046509) q[0];
sx q[0];
rz(-1.0437597) q[0];
sx q[0];
rz(2.4567132) q[0];
rz(-0.20826805) q[1];
sx q[1];
rz(-1.3083649) q[1];
sx q[1];
rz(0.79316717) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858636) q[0];
sx q[0];
rz(-1.0287713) q[0];
sx q[0];
rz(0.11491187) q[0];
rz(2.9689413) q[2];
sx q[2];
rz(-1.358629) q[2];
sx q[2];
rz(1.4868975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.316721) q[1];
sx q[1];
rz(-0.61144704) q[1];
sx q[1];
rz(1.5427521) q[1];
rz(2.9248167) q[3];
sx q[3];
rz(-1.9329251) q[3];
sx q[3];
rz(-2.0538052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3785279) q[2];
sx q[2];
rz(-1.4497764) q[2];
sx q[2];
rz(0.48279631) q[2];
rz(1.637508) q[3];
sx q[3];
rz(-0.62862325) q[3];
sx q[3];
rz(-2.8211236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046431635) q[0];
sx q[0];
rz(-0.23155364) q[0];
sx q[0];
rz(-3.0401373) q[0];
rz(-2.9265112) q[1];
sx q[1];
rz(-1.2098034) q[1];
sx q[1];
rz(-1.869841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6658507) q[0];
sx q[0];
rz(-1.2870732) q[0];
sx q[0];
rz(-1.0918777) q[0];
x q[1];
rz(-2.099718) q[2];
sx q[2];
rz(-2.1019962) q[2];
sx q[2];
rz(0.94106969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.78465652) q[1];
sx q[1];
rz(-1.9627092) q[1];
sx q[1];
rz(2.980113) q[1];
rz(-pi) q[2];
rz(1.029483) q[3];
sx q[3];
rz(-2.7446163) q[3];
sx q[3];
rz(3.027806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35470206) q[2];
sx q[2];
rz(-2.6933647) q[2];
sx q[2];
rz(2.9917713) q[2];
rz(-2.8938854) q[3];
sx q[3];
rz(-0.37023538) q[3];
sx q[3];
rz(-2.3410334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54514068) q[0];
sx q[0];
rz(-0.62262028) q[0];
sx q[0];
rz(1.2860292) q[0];
rz(-0.81895685) q[1];
sx q[1];
rz(-2.8352234) q[1];
sx q[1];
rz(0.18316306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4990627) q[0];
sx q[0];
rz(-2.1234491) q[0];
sx q[0];
rz(1.8091317) q[0];
rz(-pi) q[1];
rz(2.8122453) q[2];
sx q[2];
rz(-0.19528772) q[2];
sx q[2];
rz(-1.5415217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0350041) q[1];
sx q[1];
rz(-2.2640806) q[1];
sx q[1];
rz(0.50030746) q[1];
x q[2];
rz(2.6425903) q[3];
sx q[3];
rz(-1.2899644) q[3];
sx q[3];
rz(-0.71597354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5554819) q[2];
sx q[2];
rz(-1.4376983) q[2];
sx q[2];
rz(0.5371896) q[2];
rz(1.7389343) q[3];
sx q[3];
rz(-1.2638998) q[3];
sx q[3];
rz(-2.5130443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1527767) q[0];
sx q[0];
rz(-2.7429136) q[0];
sx q[0];
rz(-3.0231349) q[0];
rz(2.0351824) q[1];
sx q[1];
rz(-1.2245919) q[1];
sx q[1];
rz(2.4647663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40672228) q[0];
sx q[0];
rz(-3.0483735) q[0];
sx q[0];
rz(2.5073476) q[0];
rz(-2.4873892) q[2];
sx q[2];
rz(-1.6779644) q[2];
sx q[2];
rz(0.079710641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54296498) q[1];
sx q[1];
rz(-1.7522546) q[1];
sx q[1];
rz(1.7398029) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5854225) q[3];
sx q[3];
rz(-1.5717603) q[3];
sx q[3];
rz(-2.2617676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38830882) q[2];
sx q[2];
rz(-2.1948185) q[2];
sx q[2];
rz(3.0906265) q[2];
rz(0.1554337) q[3];
sx q[3];
rz(-1.4916462) q[3];
sx q[3];
rz(1.0203993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.257457) q[0];
sx q[0];
rz(-1.9866885) q[0];
sx q[0];
rz(-0.98404032) q[0];
rz(0.46977305) q[1];
sx q[1];
rz(-1.4642986) q[1];
sx q[1];
rz(-0.22107548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2673938) q[0];
sx q[0];
rz(-1.4765863) q[0];
sx q[0];
rz(-1.5896958) q[0];
rz(-pi) q[1];
rz(0.65024489) q[2];
sx q[2];
rz(-1.8762372) q[2];
sx q[2];
rz(-2.1540856) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1578678) q[1];
sx q[1];
rz(-2.7642784) q[1];
sx q[1];
rz(0.91767101) q[1];
rz(-2.2443767) q[3];
sx q[3];
rz(-0.69475896) q[3];
sx q[3];
rz(-1.0396797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2880963) q[2];
sx q[2];
rz(-2.4595538) q[2];
sx q[2];
rz(-1.3300396) q[2];
rz(2.5174777) q[3];
sx q[3];
rz(-0.95804405) q[3];
sx q[3];
rz(-0.38930711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5617111) q[0];
sx q[0];
rz(-1.6498673) q[0];
sx q[0];
rz(-2.0910182) q[0];
rz(2.1613278) q[1];
sx q[1];
rz(-1.0824208) q[1];
sx q[1];
rz(1.6260653) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1001871) q[0];
sx q[0];
rz(-0.55313191) q[0];
sx q[0];
rz(2.83908) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43584906) q[2];
sx q[2];
rz(-0.51700355) q[2];
sx q[2];
rz(-0.85384254) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2962471) q[1];
sx q[1];
rz(-1.913684) q[1];
sx q[1];
rz(-1.898349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.6330746) q[3];
sx q[3];
rz(-2.4313723) q[3];
sx q[3];
rz(-2.4274277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3773697) q[2];
sx q[2];
rz(-1.2219967) q[2];
sx q[2];
rz(-2.4019305) q[2];
rz(0.1554648) q[3];
sx q[3];
rz(-2.7630617) q[3];
sx q[3];
rz(2.9445061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71843201) q[0];
sx q[0];
rz(-0.53373706) q[0];
sx q[0];
rz(1.1370132) q[0];
rz(-2.5088189) q[1];
sx q[1];
rz(-2.4686765) q[1];
sx q[1];
rz(0.64579642) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5006913) q[0];
sx q[0];
rz(-2.7703022) q[0];
sx q[0];
rz(-0.47993203) q[0];
rz(-0.77984758) q[2];
sx q[2];
rz(-1.4044242) q[2];
sx q[2];
rz(0.8188687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2254722) q[1];
sx q[1];
rz(-1.5016306) q[1];
sx q[1];
rz(-0.81838496) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0876376) q[3];
sx q[3];
rz(-2.5541383) q[3];
sx q[3];
rz(1.7980391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48468963) q[2];
sx q[2];
rz(-0.36269665) q[2];
sx q[2];
rz(-0.038385782) q[2];
rz(0.11135993) q[3];
sx q[3];
rz(-2.7146118) q[3];
sx q[3];
rz(-0.67489433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4511694) q[0];
sx q[0];
rz(-1.7043865) q[0];
sx q[0];
rz(-1.7519328) q[0];
rz(-1.8545064) q[1];
sx q[1];
rz(-0.99936395) q[1];
sx q[1];
rz(-2.0157464) q[1];
rz(-2.1846175) q[2];
sx q[2];
rz(-1.0697094) q[2];
sx q[2];
rz(-0.30164837) q[2];
rz(1.3948364) q[3];
sx q[3];
rz(-2.3223489) q[3];
sx q[3];
rz(-1.7391218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
