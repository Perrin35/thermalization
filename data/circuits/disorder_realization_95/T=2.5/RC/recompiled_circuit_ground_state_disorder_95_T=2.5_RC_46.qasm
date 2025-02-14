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
rz(2.5511281) q[0];
rz(2.85113) q[1];
sx q[1];
rz(-1.7658748) q[1];
sx q[1];
rz(3.1380993) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0079903) q[0];
sx q[0];
rz(-0.65235814) q[0];
sx q[0];
rz(0.86835394) q[0];
x q[1];
rz(1.1972869) q[2];
sx q[2];
rz(-0.71564681) q[2];
sx q[2];
rz(2.6386821) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.063979186) q[1];
sx q[1];
rz(-0.031647041) q[1];
sx q[1];
rz(-2.3961303) q[1];
rz(-pi) q[2];
rz(-0.9736075) q[3];
sx q[3];
rz(-1.191621) q[3];
sx q[3];
rz(-1.7969064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8840088) q[2];
sx q[2];
rz(-0.85870063) q[2];
sx q[2];
rz(-1.0165455) q[2];
rz(-0.86333418) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(1.6421002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.6739552) q[0];
sx q[0];
rz(-1.6683945) q[0];
sx q[0];
rz(-0.13482811) q[0];
rz(2.9148031) q[1];
sx q[1];
rz(-3.0786381) q[1];
sx q[1];
rz(0.24980587) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86827474) q[0];
sx q[0];
rz(-2.3767243) q[0];
sx q[0];
rz(2.5165783) q[0];
rz(-pi) q[1];
rz(-1.6667566) q[2];
sx q[2];
rz(-0.50228206) q[2];
sx q[2];
rz(1.407215) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.47905251) q[1];
sx q[1];
rz(-1.4932649) q[1];
sx q[1];
rz(-1.5815602) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2050912) q[3];
sx q[3];
rz(-1.4451175) q[3];
sx q[3];
rz(-0.92363908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2251542) q[2];
sx q[2];
rz(-0.068000451) q[2];
sx q[2];
rz(-0.98006836) q[2];
rz(-2.2465536) q[3];
sx q[3];
rz(-0.74256247) q[3];
sx q[3];
rz(1.1802973) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5880244) q[0];
sx q[0];
rz(-2.6346485) q[0];
sx q[0];
rz(-1.6059599) q[0];
rz(-1.5060679) q[1];
sx q[1];
rz(-0.82553828) q[1];
sx q[1];
rz(0.95605409) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0900537) q[0];
sx q[0];
rz(-2.9057626) q[0];
sx q[0];
rz(0.43311849) q[0];
rz(-1.5044841) q[2];
sx q[2];
rz(-0.8408747) q[2];
sx q[2];
rz(0.17375565) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0625831) q[1];
sx q[1];
rz(-2.2979567) q[1];
sx q[1];
rz(-0.14240245) q[1];
rz(-pi) q[2];
rz(-1.2230824) q[3];
sx q[3];
rz(-1.082431) q[3];
sx q[3];
rz(2.2852416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.306159) q[2];
sx q[2];
rz(-2.4296438) q[2];
sx q[2];
rz(0.63208675) q[2];
rz(2.2575374) q[3];
sx q[3];
rz(-0.017280936) q[3];
sx q[3];
rz(2.250905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.66568351) q[0];
sx q[0];
rz(-1.2983687) q[0];
sx q[0];
rz(0.16715288) q[0];
rz(-1.2627603) q[1];
sx q[1];
rz(-0.94307584) q[1];
sx q[1];
rz(-1.4080338) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77771132) q[0];
sx q[0];
rz(-1.9929307) q[0];
sx q[0];
rz(1.889983) q[0];
rz(-pi) q[1];
rz(-1.5535545) q[2];
sx q[2];
rz(-0.21176007) q[2];
sx q[2];
rz(2.0469768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.38630292) q[1];
sx q[1];
rz(-1.3805256) q[1];
sx q[1];
rz(-2.931837) q[1];
x q[2];
rz(-0.77045958) q[3];
sx q[3];
rz(-1.3589199) q[3];
sx q[3];
rz(0.4650863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7888646) q[2];
sx q[2];
rz(-3.126324) q[2];
sx q[2];
rz(0.35066476) q[2];
rz(2.8777425) q[3];
sx q[3];
rz(-0.00084547384) q[3];
sx q[3];
rz(-1.9381757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8769787) q[0];
sx q[0];
rz(-0.74540859) q[0];
sx q[0];
rz(-1.4234446) q[0];
rz(2.8599332) q[1];
sx q[1];
rz(-1.5080844) q[1];
sx q[1];
rz(-1.8219832) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9284009) q[0];
sx q[0];
rz(-1.1494779) q[0];
sx q[0];
rz(-2.1405959) q[0];
rz(0.93792589) q[2];
sx q[2];
rz(-3.1073776) q[2];
sx q[2];
rz(0.93955112) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4180243) q[1];
sx q[1];
rz(-0.71360525) q[1];
sx q[1];
rz(-0.51939555) q[1];
rz(-pi) q[2];
rz(1.237147) q[3];
sx q[3];
rz(-1.9556442) q[3];
sx q[3];
rz(0.24144444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8772956) q[2];
sx q[2];
rz(-3.0946315) q[2];
sx q[2];
rz(-1.4287255) q[2];
rz(-1.9127539) q[3];
sx q[3];
rz(-0.052534025) q[3];
sx q[3];
rz(-0.088168941) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-0.30217198) q[1];
sx q[1];
rz(-1.3928394) q[1];
sx q[1];
rz(-0.43513939) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36212039) q[0];
sx q[0];
rz(-1.5359306) q[0];
sx q[0];
rz(-0.21135862) q[0];
x q[1];
rz(1.5717634) q[2];
sx q[2];
rz(-1.5633601) q[2];
sx q[2];
rz(0.32583562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8813215) q[1];
sx q[1];
rz(-1.146293) q[1];
sx q[1];
rz(0.05646245) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76558785) q[3];
sx q[3];
rz(-2.2556955) q[3];
sx q[3];
rz(1.9960038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90007323) q[2];
sx q[2];
rz(-3.1406431) q[2];
sx q[2];
rz(-2.1692236) q[2];
rz(2.1420245) q[3];
sx q[3];
rz(-3.1344423) q[3];
sx q[3];
rz(-1.0424559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59357506) q[0];
sx q[0];
rz(-1.8997718) q[0];
sx q[0];
rz(1.0663363) q[0];
rz(1.4284596) q[1];
sx q[1];
rz(-2.8722873) q[1];
sx q[1];
rz(-1.2880633) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6664827) q[0];
sx q[0];
rz(-1.8865856) q[0];
sx q[0];
rz(-0.72986365) q[0];
rz(-pi) q[1];
rz(-0.35310642) q[2];
sx q[2];
rz(-2.2688008) q[2];
sx q[2];
rz(-1.6153796) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8937872) q[1];
sx q[1];
rz(-1.6105741) q[1];
sx q[1];
rz(-0.013377845) q[1];
rz(1.0405447) q[3];
sx q[3];
rz(-1.3589348) q[3];
sx q[3];
rz(2.7303641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0224033) q[2];
sx q[2];
rz(-0.069312118) q[2];
sx q[2];
rz(0.27528396) q[2];
rz(2.9149808) q[3];
sx q[3];
rz(-2.785502) q[3];
sx q[3];
rz(-2.7039458) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3771123) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(2.5544033) q[0];
rz(-1.5234692) q[1];
sx q[1];
rz(-1.1174997) q[1];
sx q[1];
rz(1.5706185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2908774) q[0];
sx q[0];
rz(-1.1505011) q[0];
sx q[0];
rz(0.72776003) q[0];
x q[1];
rz(-0.067085548) q[2];
sx q[2];
rz(-1.1076945) q[2];
sx q[2];
rz(2.9912134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0084616737) q[1];
sx q[1];
rz(-1.6146286) q[1];
sx q[1];
rz(0.030972199) q[1];
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
rz(-pi) q[1];
rz(1.9470584) q[2];
sx q[2];
rz(-2.9583866) q[2];
sx q[2];
rz(1.8397231) q[2];
rz(-0.434508) q[3];
sx q[3];
rz(-3.1012111) q[3];
sx q[3];
rz(-2.6935008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3473564) q[0];
sx q[0];
rz(-0.11581049) q[0];
sx q[0];
rz(1.7606803) q[0];
rz(-1.5538838) q[1];
sx q[1];
rz(-0.96276182) q[1];
sx q[1];
rz(0.10996058) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9048683) q[0];
sx q[0];
rz(-2.3873608) q[0];
sx q[0];
rz(-1.2931025) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70297697) q[2];
sx q[2];
rz(-1.7406782) q[2];
sx q[2];
rz(1.301606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1404071) q[1];
sx q[1];
rz(-2.3275536) q[1];
sx q[1];
rz(0.034239727) q[1];
rz(-pi) q[2];
rz(1.7305829) q[3];
sx q[3];
rz(-1.8786123) q[3];
sx q[3];
rz(0.45444835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3515557) q[2];
sx q[2];
rz(-1.0645126) q[2];
sx q[2];
rz(-2.7891187) q[2];
rz(-2.500109) q[3];
sx q[3];
rz(-3.0951169) q[3];
sx q[3];
rz(-0.36267734) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8620257) q[0];
sx q[0];
rz(-2.8643705) q[0];
sx q[0];
rz(0.56420457) q[0];
rz(1.5203681) q[1];
sx q[1];
rz(-1.0313326) q[1];
sx q[1];
rz(0.079781562) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.736981) q[0];
sx q[0];
rz(-0.68697646) q[0];
sx q[0];
rz(-1.3300598) q[0];
x q[1];
rz(2.3529813) q[2];
sx q[2];
rz(-3.0742767) q[2];
sx q[2];
rz(0.88035098) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0003766) q[1];
sx q[1];
rz(-1.8264007) q[1];
sx q[1];
rz(-0.80900286) q[1];
x q[2];
rz(-1.4873119) q[3];
sx q[3];
rz(-1.8125839) q[3];
sx q[3];
rz(0.35095645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0148049) q[2];
sx q[2];
rz(-0.0099651907) q[2];
sx q[2];
rz(-1.0753746) q[2];
rz(-2.3672095) q[3];
sx q[3];
rz(-3.1175678) q[3];
sx q[3];
rz(-0.27124852) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423691) q[0];
sx q[0];
rz(-1.5726226) q[0];
sx q[0];
rz(1.5725305) q[0];
rz(-0.52539274) q[1];
sx q[1];
rz(-0.071594302) q[1];
sx q[1];
rz(-0.20803861) q[1];
rz(1.8990429) q[2];
sx q[2];
rz(-0.59292159) q[2];
sx q[2];
rz(-2.4688724) q[2];
rz(2.0490859) q[3];
sx q[3];
rz(-2.1989965) q[3];
sx q[3];
rz(0.60529706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
