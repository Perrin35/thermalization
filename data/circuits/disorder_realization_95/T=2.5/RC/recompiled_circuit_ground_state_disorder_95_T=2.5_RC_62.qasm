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
rz(-3.1380993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19102819) q[0];
sx q[0];
rz(-1.0890245) q[0];
sx q[0];
rz(0.45848565) q[0];
x q[1];
rz(1.9443058) q[2];
sx q[2];
rz(-0.71564681) q[2];
sx q[2];
rz(0.50291059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3799882) q[1];
sx q[1];
rz(-1.5922609) q[1];
sx q[1];
rz(3.1183356) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95432561) q[3];
sx q[3];
rz(-2.4467957) q[3];
sx q[3];
rz(-0.27237005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2575839) q[2];
sx q[2];
rz(-2.282892) q[2];
sx q[2];
rz(-1.0165455) q[2];
rz(0.86333418) q[3];
sx q[3];
rz(-0.038067929) q[3];
sx q[3];
rz(1.4994924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6739552) q[0];
sx q[0];
rz(-1.4731982) q[0];
sx q[0];
rz(0.13482811) q[0];
rz(2.9148031) q[1];
sx q[1];
rz(-0.062954523) q[1];
sx q[1];
rz(-0.24980587) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1824719) q[0];
sx q[0];
rz(-1.9879436) q[0];
sx q[0];
rz(-0.66138256) q[0];
x q[1];
rz(-1.4748361) q[2];
sx q[2];
rz(-2.6393106) q[2];
sx q[2];
rz(1.407215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5244478) q[1];
sx q[1];
rz(-3.063319) q[1];
sx q[1];
rz(-0.13767482) q[1];
x q[2];
rz(1.2311801) q[3];
sx q[3];
rz(-0.38577739) q[3];
sx q[3];
rz(0.3308109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2251542) q[2];
sx q[2];
rz(-3.0735922) q[2];
sx q[2];
rz(-2.1615243) q[2];
rz(-2.2465536) q[3];
sx q[3];
rz(-0.74256247) q[3];
sx q[3];
rz(1.1802973) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5880244) q[0];
sx q[0];
rz(-2.6346485) q[0];
sx q[0];
rz(1.5356327) q[0];
rz(1.6355248) q[1];
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
rz(-1.0900537) q[0];
sx q[0];
rz(-2.9057626) q[0];
sx q[0];
rz(-0.43311849) q[0];
rz(1.5044841) q[2];
sx q[2];
rz(-0.8408747) q[2];
sx q[2];
rz(-0.17375565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29143476) q[1];
sx q[1];
rz(-2.4031284) q[1];
sx q[1];
rz(-1.4126331) q[1];
rz(-pi) q[2];
rz(1.9185103) q[3];
sx q[3];
rz(-2.0591617) q[3];
sx q[3];
rz(0.85635105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.306159) q[2];
sx q[2];
rz(-0.71194887) q[2];
sx q[2];
rz(2.5095059) q[2];
rz(-2.2575374) q[3];
sx q[3];
rz(-3.1243117) q[3];
sx q[3];
rz(-0.89068762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4759091) q[0];
sx q[0];
rz(-1.2983687) q[0];
sx q[0];
rz(-2.9744398) q[0];
rz(1.2627603) q[1];
sx q[1];
rz(-2.1985168) q[1];
sx q[1];
rz(1.7335588) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65850019) q[0];
sx q[0];
rz(-1.8611055) q[0];
sx q[0];
rz(2.6997552) q[0];
x q[1];
rz(3.1378861) q[2];
sx q[2];
rz(-1.7825244) q[2];
sx q[2];
rz(-2.0293411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9973397) q[1];
sx q[1];
rz(-1.3648811) q[1];
sx q[1];
rz(1.3763672) q[1];
rz(-pi) q[2];
rz(-2.3711331) q[3];
sx q[3];
rz(-1.7826728) q[3];
sx q[3];
rz(0.4650863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.35272804) q[2];
sx q[2];
rz(-0.015268607) q[2];
sx q[2];
rz(-0.35066476) q[2];
rz(2.8777425) q[3];
sx q[3];
rz(-3.1407472) q[3];
sx q[3];
rz(-1.2034169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8769787) q[0];
sx q[0];
rz(-0.74540859) q[0];
sx q[0];
rz(1.718148) q[0];
rz(2.8599332) q[1];
sx q[1];
rz(-1.6335082) q[1];
sx q[1];
rz(-1.3196094) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9284009) q[0];
sx q[0];
rz(-1.1494779) q[0];
sx q[0];
rz(-1.0009967) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5432038) q[2];
sx q[2];
rz(-1.5910307) q[2];
sx q[2];
rz(1.8777562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77051059) q[1];
sx q[1];
rz(-0.96643172) q[1];
sx q[1];
rz(1.9766859) q[1];
rz(0.67992391) q[3];
sx q[3];
rz(-2.6377684) q[3];
sx q[3];
rz(2.6375204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.264297) q[2];
sx q[2];
rz(-3.0946315) q[2];
sx q[2];
rz(-1.4287255) q[2];
rz(-1.9127539) q[3];
sx q[3];
rz(-3.0890586) q[3];
sx q[3];
rz(0.088168941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0838715) q[0];
sx q[0];
rz(-3.0241522) q[0];
sx q[0];
rz(-2.591326) q[0];
rz(0.30217198) q[1];
sx q[1];
rz(-1.3928394) q[1];
sx q[1];
rz(0.43513939) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9254375) q[0];
sx q[0];
rz(-1.3595681) q[0];
sx q[0];
rz(1.5351377) q[0];
x q[1];
rz(3.0122717) q[2];
sx q[2];
rz(-3.1340938) q[2];
sx q[2];
rz(-0.19651112) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8813215) q[1];
sx q[1];
rz(-1.9952996) q[1];
sx q[1];
rz(3.0851302) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3760048) q[3];
sx q[3];
rz(-0.88589719) q[3];
sx q[3];
rz(1.9960038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90007323) q[2];
sx q[2];
rz(-0.0009495432) q[2];
sx q[2];
rz(-2.1692236) q[2];
rz(-0.9995681) q[3];
sx q[3];
rz(-3.1344423) q[3];
sx q[3];
rz(-1.0424559) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59357506) q[0];
sx q[0];
rz(-1.8997718) q[0];
sx q[0];
rz(2.0752564) q[0];
rz(1.713133) q[1];
sx q[1];
rz(-0.26930535) q[1];
sx q[1];
rz(-1.2880633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36670524) q[0];
sx q[0];
rz(-2.2572491) q[0];
sx q[0];
rz(1.1576325) q[0];
x q[1];
rz(-2.3003103) q[2];
sx q[2];
rz(-1.3026574) q[2];
sx q[2];
rz(-0.27716217) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8180698) q[1];
sx q[1];
rz(-1.5574291) q[1];
sx q[1];
rz(-1.531015) q[1];
rz(-pi) q[2];
rz(0.24434511) q[3];
sx q[3];
rz(-1.0536031) q[3];
sx q[3];
rz(1.2822267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11918934) q[2];
sx q[2];
rz(-0.069312118) q[2];
sx q[2];
rz(0.27528396) q[2];
rz(-0.22661181) q[3];
sx q[3];
rz(-0.35609069) q[3];
sx q[3];
rz(2.7039458) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76448035) q[0];
sx q[0];
rz(-0.28700101) q[0];
sx q[0];
rz(2.5544033) q[0];
rz(-1.5234692) q[1];
sx q[1];
rz(-1.1174997) q[1];
sx q[1];
rz(1.5706185) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8507152) q[0];
sx q[0];
rz(-1.9910915) q[0];
sx q[0];
rz(-0.72776003) q[0];
rz(-pi) q[1];
rz(-1.4373407) q[2];
sx q[2];
rz(-0.46758662) q[2];
sx q[2];
rz(-0.0011006265) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5636922) q[1];
sx q[1];
rz(-1.6017388) q[1];
sx q[1];
rz(-1.5269431) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0471025) q[3];
sx q[3];
rz(-0.62783754) q[3];
sx q[3];
rz(1.7881356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1945343) q[2];
sx q[2];
rz(-0.18320601) q[2];
sx q[2];
rz(1.8397231) q[2];
rz(-2.7070847) q[3];
sx q[3];
rz(-3.1012111) q[3];
sx q[3];
rz(-0.44809189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7942363) q[0];
sx q[0];
rz(-0.11581049) q[0];
sx q[0];
rz(1.7606803) q[0];
rz(1.5538838) q[1];
sx q[1];
rz(-0.96276182) q[1];
sx q[1];
rz(3.0316321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048683) q[0];
sx q[0];
rz(-0.75423181) q[0];
sx q[0];
rz(-1.2931025) q[0];
rz(-0.70297697) q[2];
sx q[2];
rz(-1.7406782) q[2];
sx q[2];
rz(1.8399866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.090558) q[1];
sx q[1];
rz(-0.75737774) q[1];
sx q[1];
rz(-1.6070328) q[1];
rz(-pi) q[2];
rz(1.4110097) q[3];
sx q[3];
rz(-1.8786123) q[3];
sx q[3];
rz(-0.45444835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79003698) q[2];
sx q[2];
rz(-1.0645126) q[2];
sx q[2];
rz(-2.7891187) q[2];
rz(-0.6414837) q[3];
sx q[3];
rz(-3.0951169) q[3];
sx q[3];
rz(0.36267734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27956692) q[0];
sx q[0];
rz(-2.8643705) q[0];
sx q[0];
rz(2.5773881) q[0];
rz(-1.5203681) q[1];
sx q[1];
rz(-2.11026) q[1];
sx q[1];
rz(-3.0618111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4046116) q[0];
sx q[0];
rz(-0.68697646) q[0];
sx q[0];
rz(1.3300598) q[0];
rz(3.0941102) q[2];
sx q[2];
rz(-1.523062) q[2];
sx q[2];
rz(1.477923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0003766) q[1];
sx q[1];
rz(-1.8264007) q[1];
sx q[1];
rz(-0.80900286) q[1];
rz(-pi) q[2];
rz(-2.8989927) q[3];
sx q[3];
rz(-1.651847) q[3];
sx q[3];
rz(-1.9417861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0148049) q[2];
sx q[2];
rz(-0.0099651907) q[2];
sx q[2];
rz(1.0753746) q[2];
rz(2.3672095) q[3];
sx q[3];
rz(-0.024024809) q[3];
sx q[3];
rz(-0.27124852) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423691) q[0];
sx q[0];
rz(-1.5726226) q[0];
sx q[0];
rz(1.5725305) q[0];
rz(2.6161999) q[1];
sx q[1];
rz(-0.071594302) q[1];
sx q[1];
rz(-0.20803861) q[1];
rz(0.21389773) q[2];
sx q[2];
rz(-2.1281617) q[2];
sx q[2];
rz(-2.8585363) q[2];
rz(1.0925068) q[3];
sx q[3];
rz(-0.94259613) q[3];
sx q[3];
rz(-2.5362956) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
