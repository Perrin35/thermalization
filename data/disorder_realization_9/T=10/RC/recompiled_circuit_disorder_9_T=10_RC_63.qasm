OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(2.467632) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(2.34692) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2406143) q[0];
sx q[0];
rz(-2.6128747) q[0];
sx q[0];
rz(-0.76529495) q[0];
rz(-pi) q[1];
rz(1.4235714) q[2];
sx q[2];
rz(-1.0861673) q[2];
sx q[2];
rz(1.1046315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1311878) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(-2.5351934) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13866339) q[3];
sx q[3];
rz(-1.1274459) q[3];
sx q[3];
rz(1.5821622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-2.2241212) q[3];
sx q[3];
rz(0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(0.95970884) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(0.62746343) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9976881) q[0];
sx q[0];
rz(-1.5630957) q[0];
sx q[0];
rz(1.5639923) q[0];
x q[1];
rz(-1.5487899) q[2];
sx q[2];
rz(-1.7823879) q[2];
sx q[2];
rz(1.6197268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6154502) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(2.1089094) q[1];
x q[2];
rz(1.1153568) q[3];
sx q[3];
rz(-3.1133828) q[3];
sx q[3];
rz(1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(0.97989782) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(-1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58127789) q[0];
sx q[0];
rz(-0.53239765) q[0];
sx q[0];
rz(2.6530932) q[0];
rz(1.5930487) q[2];
sx q[2];
rz(-1.77102) q[2];
sx q[2];
rz(-2.3929838) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2587535) q[1];
sx q[1];
rz(-1.6707318) q[1];
sx q[1];
rz(0.12211166) q[1];
rz(0.090431902) q[3];
sx q[3];
rz(-1.3867497) q[3];
sx q[3];
rz(-0.66891608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-2.7518318) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-2.5770082) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72258654) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(1.0976085) q[0];
rz(-0.98323804) q[2];
sx q[2];
rz(-1.431627) q[2];
sx q[2];
rz(2.522574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.035186471) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(0.10722843) q[1];
x q[2];
rz(2.5921949) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(2.3466563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(-2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-2.1767298) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(0.66666493) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4283838) q[0];
sx q[0];
rz(-0.41233006) q[0];
sx q[0];
rz(-1.3165738) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(-1.0416043) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95496817) q[1];
sx q[1];
rz(-2.4116227) q[1];
sx q[1];
rz(2.7093922) q[1];
rz(-2.9186451) q[3];
sx q[3];
rz(-1.2811536) q[3];
sx q[3];
rz(1.2900316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.6113575) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0479193) q[0];
sx q[0];
rz(-1.1512655) q[0];
sx q[0];
rz(1.6506667) q[0];
x q[1];
rz(-0.62263454) q[2];
sx q[2];
rz(-2.6175675) q[2];
sx q[2];
rz(-0.084409075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46576071) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(-2.0783706) q[1];
rz(-pi) q[2];
rz(-3.0909782) q[3];
sx q[3];
rz(-1.317987) q[3];
sx q[3];
rz(1.8265754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.6507089) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-0.18280612) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144007) q[0];
sx q[0];
rz(-1.6806024) q[0];
sx q[0];
rz(-0.26392428) q[0];
x q[1];
rz(-0.73166087) q[2];
sx q[2];
rz(-0.96207843) q[2];
sx q[2];
rz(-0.70639709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.70799815) q[1];
sx q[1];
rz(-0.79213789) q[1];
sx q[1];
rz(-2.5670693) q[1];
rz(-pi) q[2];
rz(-1.0484344) q[3];
sx q[3];
rz(-1.8494542) q[3];
sx q[3];
rz(3.0184961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(-0.42256045) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(1.0890695) q[0];
rz(2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7437744) q[0];
sx q[0];
rz(-1.213307) q[0];
sx q[0];
rz(-2.9767569) q[0];
rz(-pi) q[1];
rz(2.0681357) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(-1.9286326) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(2.5177588) q[1];
rz(-pi) q[2];
rz(0.20181228) q[3];
sx q[3];
rz(-1.0045856) q[3];
sx q[3];
rz(-1.4334397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(-0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-0.91264549) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(2.7776921) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(-2.231853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7649987) q[0];
sx q[0];
rz(-1.8928327) q[0];
sx q[0];
rz(-1.2850719) q[0];
rz(-pi) q[1];
rz(1.7889195) q[2];
sx q[2];
rz(-1.911273) q[2];
sx q[2];
rz(-0.097629658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7418356) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(-0.21106212) q[1];
x q[2];
rz(-2.8839888) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(-0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45458083) q[0];
sx q[0];
rz(-0.63475383) q[0];
sx q[0];
rz(1.4995585) q[0];
rz(1.1012494) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(-1.7061526) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2243832) q[1];
sx q[1];
rz(-1.6341126) q[1];
sx q[1];
rz(1.2095832) q[1];
x q[2];
rz(-1.7747202) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5932896) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(0.9085761) q[2];
sx q[2];
rz(-0.77975811) q[2];
sx q[2];
rz(-2.1644885) q[2];
rz(-1.2699119) q[3];
sx q[3];
rz(-1.6172998) q[3];
sx q[3];
rz(-0.87340364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];