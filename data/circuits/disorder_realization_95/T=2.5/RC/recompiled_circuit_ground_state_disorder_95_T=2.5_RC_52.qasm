OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18093827) q[0];
sx q[0];
rz(-3.0783983) q[0];
sx q[0];
rz(0.59046459) q[0];
rz(-0.2904627) q[1];
sx q[1];
rz(-1.3757179) q[1];
sx q[1];
rz(0.003493316) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1549415) q[0];
sx q[0];
rz(-1.9738324) q[0];
sx q[0];
rz(2.0986845) q[0];
x q[1];
rz(-2.8344049) q[2];
sx q[2];
rz(-2.2279539) q[2];
sx q[2];
rz(0.023935723) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3319015) q[1];
sx q[1];
rz(-1.5475447) q[1];
sx q[1];
rz(-1.5922668) q[1];
rz(-pi) q[2];
rz(-2.1679852) q[3];
sx q[3];
rz(-1.9499717) q[3];
sx q[3];
rz(-1.7969064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2575839) q[2];
sx q[2];
rz(-0.85870063) q[2];
sx q[2];
rz(-2.1250471) q[2];
rz(-0.86333418) q[3];
sx q[3];
rz(-3.1035247) q[3];
sx q[3];
rz(-1.6421002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6739552) q[0];
sx q[0];
rz(-1.6683945) q[0];
sx q[0];
rz(3.0067645) q[0];
rz(-0.22678953) q[1];
sx q[1];
rz(-0.062954523) q[1];
sx q[1];
rz(-0.24980587) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0587388) q[0];
sx q[0];
rz(-2.1670411) q[0];
sx q[0];
rz(1.0591177) q[0];
rz(-pi) q[1];
rz(-1.6667566) q[2];
sx q[2];
rz(-2.6393106) q[2];
sx q[2];
rz(1.7343777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6625401) q[1];
sx q[1];
rz(-1.4932649) q[1];
sx q[1];
rz(1.5600325) q[1];
rz(-pi) q[2];
rz(-1.9365015) q[3];
sx q[3];
rz(-1.4451175) q[3];
sx q[3];
rz(2.2179536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9164385) q[2];
sx q[2];
rz(-3.0735922) q[2];
sx q[2];
rz(0.98006836) q[2];
rz(-2.2465536) q[3];
sx q[3];
rz(-2.3990302) q[3];
sx q[3];
rz(-1.1802973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5535683) q[0];
sx q[0];
rz(-2.6346485) q[0];
sx q[0];
rz(-1.5356327) q[0];
rz(-1.5060679) q[1];
sx q[1];
rz(-2.3160544) q[1];
sx q[1];
rz(2.1855386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5339622) q[0];
sx q[0];
rz(-1.3570983) q[0];
sx q[0];
rz(-1.4702809) q[0];
rz(-pi) q[1];
rz(1.6371085) q[2];
sx q[2];
rz(-0.8408747) q[2];
sx q[2];
rz(0.17375565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0625831) q[1];
sx q[1];
rz(-0.843636) q[1];
sx q[1];
rz(2.9991902) q[1];
x q[2];
rz(2.5713137) q[3];
sx q[3];
rz(-0.59123838) q[3];
sx q[3];
rz(1.6275004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8354336) q[2];
sx q[2];
rz(-0.71194887) q[2];
sx q[2];
rz(-2.5095059) q[2];
rz(-0.88405526) q[3];
sx q[3];
rz(-3.1243117) q[3];
sx q[3];
rz(0.89068762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4759091) q[0];
sx q[0];
rz(-1.8432239) q[0];
sx q[0];
rz(0.16715288) q[0];
rz(-1.2627603) q[1];
sx q[1];
rz(-0.94307584) q[1];
sx q[1];
rz(1.7335588) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4830925) q[0];
sx q[0];
rz(-1.2804872) q[0];
sx q[0];
rz(0.4418375) q[0];
x q[1];
rz(-1.5880382) q[2];
sx q[2];
rz(-2.9298326) q[2];
sx q[2];
rz(2.0469768) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1442529) q[1];
sx q[1];
rz(-1.7767116) q[1];
sx q[1];
rz(-1.7652254) q[1];
x q[2];
rz(2.8420343) q[3];
sx q[3];
rz(-0.79323461) q[3];
sx q[3];
rz(1.8225074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7888646) q[2];
sx q[2];
rz(-3.126324) q[2];
sx q[2];
rz(-2.7909279) q[2];
rz(2.8777425) q[3];
sx q[3];
rz(-3.1407472) q[3];
sx q[3];
rz(1.9381757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8769787) q[0];
sx q[0];
rz(-0.74540859) q[0];
sx q[0];
rz(1.4234446) q[0];
rz(-2.8599332) q[1];
sx q[1];
rz(-1.5080844) q[1];
sx q[1];
rz(-1.3196094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78948259) q[0];
sx q[0];
rz(-2.4470959) q[0];
sx q[0];
rz(2.2640034) q[0];
rz(-pi) q[1];
rz(3.1213506) q[2];
sx q[2];
rz(-1.5983832) q[2];
sx q[2];
rz(-0.30640145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4180243) q[1];
sx q[1];
rz(-0.71360525) q[1];
sx q[1];
rz(-0.51939555) q[1];
x q[2];
rz(-1.9044456) q[3];
sx q[3];
rz(-1.1859484) q[3];
sx q[3];
rz(-0.24144444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8772956) q[2];
sx q[2];
rz(-3.0946315) q[2];
sx q[2];
rz(-1.7128672) q[2];
rz(-1.2288387) q[3];
sx q[3];
rz(-3.0890586) q[3];
sx q[3];
rz(3.0534237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0577211) q[0];
sx q[0];
rz(-3.0241522) q[0];
sx q[0];
rz(-2.591326) q[0];
rz(-2.8394207) q[1];
sx q[1];
rz(-1.3928394) q[1];
sx q[1];
rz(0.43513939) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2161552) q[0];
sx q[0];
rz(-1.7820246) q[0];
sx q[0];
rz(-1.5351377) q[0];
rz(-1.5698293) q[2];
sx q[2];
rz(-1.5633601) q[2];
sx q[2];
rz(-2.815757) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2602712) q[1];
sx q[1];
rz(-1.146293) q[1];
sx q[1];
rz(0.05646245) q[1];
rz(-pi) q[2];
rz(-0.76558785) q[3];
sx q[3];
rz(-0.88589719) q[3];
sx q[3];
rz(1.9960038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90007323) q[2];
sx q[2];
rz(-0.0009495432) q[2];
sx q[2];
rz(0.97236902) q[2];
rz(2.1420245) q[3];
sx q[3];
rz(-0.0071503706) q[3];
sx q[3];
rz(1.0424559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59357506) q[0];
sx q[0];
rz(-1.2418208) q[0];
sx q[0];
rz(1.0663363) q[0];
rz(-1.4284596) q[1];
sx q[1];
rz(-2.8722873) q[1];
sx q[1];
rz(1.2880633) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.47511) q[0];
sx q[0];
rz(-1.2550071) q[0];
sx q[0];
rz(2.411729) q[0];
rz(2.3003103) q[2];
sx q[2];
rz(-1.8389353) q[2];
sx q[2];
rz(2.8644305) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.076723531) q[1];
sx q[1];
rz(-3.0996267) q[1];
sx q[1];
rz(-1.8950592) q[1];
x q[2];
rz(-1.1686822) q[3];
sx q[3];
rz(-0.56723226) q[3];
sx q[3];
rz(-0.81517787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11918934) q[2];
sx q[2];
rz(-0.069312118) q[2];
sx q[2];
rz(-2.8663087) q[2];
rz(-2.9149808) q[3];
sx q[3];
rz(-0.35609069) q[3];
sx q[3];
rz(0.4376469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76448035) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(0.58718938) q[0];
rz(1.6181234) q[1];
sx q[1];
rz(-1.1174997) q[1];
sx q[1];
rz(-1.5709741) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0728517) q[0];
sx q[0];
rz(-2.223461) q[0];
sx q[0];
rz(1.0314343) q[0];
x q[1];
rz(-1.704252) q[2];
sx q[2];
rz(-2.674006) q[2];
sx q[2];
rz(-0.0011006265) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5636922) q[1];
sx q[1];
rz(-1.5398539) q[1];
sx q[1];
rz(-1.5269431) q[1];
x q[2];
rz(-2.8203378) q[3];
sx q[3];
rz(-2.1200088) q[3];
sx q[3];
rz(-2.3556701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1945343) q[2];
sx q[2];
rz(-2.9583866) q[2];
sx q[2];
rz(-1.3018695) q[2];
rz(2.7070847) q[3];
sx q[3];
rz(-3.1012111) q[3];
sx q[3];
rz(0.44809189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3473564) q[0];
sx q[0];
rz(-3.0257822) q[0];
sx q[0];
rz(-1.7606803) q[0];
rz(-1.5538838) q[1];
sx q[1];
rz(-0.96276182) q[1];
sx q[1];
rz(-3.0316321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5320265) q[0];
sx q[0];
rz(-2.2896165) q[0];
sx q[0];
rz(0.25208313) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3496341) q[2];
sx q[2];
rz(-0.87995565) q[2];
sx q[2];
rz(-3.0146989) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.546098) q[1];
sx q[1];
rz(-1.5956886) q[1];
sx q[1];
rz(2.3278462) q[1];
x q[2];
rz(0.31153714) q[3];
sx q[3];
rz(-1.4185801) q[3];
sx q[3];
rz(-1.0675586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79003698) q[2];
sx q[2];
rz(-2.07708) q[2];
sx q[2];
rz(2.7891187) q[2];
rz(0.6414837) q[3];
sx q[3];
rz(-0.04647579) q[3];
sx q[3];
rz(-2.7789153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8620257) q[0];
sx q[0];
rz(-0.27722219) q[0];
sx q[0];
rz(-0.56420457) q[0];
rz(1.6212246) q[1];
sx q[1];
rz(-2.11026) q[1];
sx q[1];
rz(-3.0618111) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1201982) q[0];
sx q[0];
rz(-1.4190089) q[0];
sx q[0];
rz(0.89812507) q[0];
rz(-pi) q[1];
rz(-2.3529813) q[2];
sx q[2];
rz(-0.067316003) q[2];
sx q[2];
rz(-2.2612417) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3120739) q[1];
sx q[1];
rz(-0.79527277) q[1];
sx q[1];
rz(-1.2088695) q[1];
rz(-pi) q[2];
rz(0.24259992) q[3];
sx q[3];
rz(-1.651847) q[3];
sx q[3];
rz(-1.9417861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12678777) q[2];
sx q[2];
rz(-3.1316275) q[2];
sx q[2];
rz(1.0753746) q[2];
rz(-0.77438313) q[3];
sx q[3];
rz(-0.024024809) q[3];
sx q[3];
rz(-0.27124852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(1.2992235) q[0];
sx q[0];
rz(-1.5689701) q[0];
sx q[0];
rz(-1.5690621) q[0];
rz(2.6161999) q[1];
sx q[1];
rz(-0.071594302) q[1];
sx q[1];
rz(-0.20803861) q[1];
rz(-0.21389773) q[2];
sx q[2];
rz(-1.0134309) q[2];
sx q[2];
rz(0.2830563) q[2];
rz(-0.68572509) q[3];
sx q[3];
rz(-1.9523806) q[3];
sx q[3];
rz(2.4718391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
