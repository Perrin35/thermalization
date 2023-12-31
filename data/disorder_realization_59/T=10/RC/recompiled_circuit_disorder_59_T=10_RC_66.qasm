OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(2.3595915) q[0];
sx q[0];
rz(10.695988) q[0];
rz(-2.864569) q[1];
sx q[1];
rz(-2.6695873) q[1];
sx q[1];
rz(-0.0013874887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3256677) q[0];
sx q[0];
rz(-2.7357833) q[0];
sx q[0];
rz(1.9947467) q[0];
rz(-3.0543047) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(2.0729614) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6359771) q[1];
sx q[1];
rz(-2.639289) q[1];
sx q[1];
rz(2.99519) q[1];
rz(-pi) q[2];
rz(-0.9790768) q[3];
sx q[3];
rz(-0.43833971) q[3];
sx q[3];
rz(1.0867659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9871621) q[2];
sx q[2];
rz(-2.5240832) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(2.1253712) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(2.7367676) q[3];
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
rz(-1.7063023) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(0.81545365) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1708508) q[0];
sx q[0];
rz(-2.4903957) q[0];
sx q[0];
rz(-0.09911508) q[0];
rz(0.72950659) q[2];
sx q[2];
rz(-1.7765877) q[2];
sx q[2];
rz(-1.6934998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.38824575) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(-1.6080329) q[1];
x q[2];
rz(2.8498597) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(0.093689703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6796391) q[2];
sx q[2];
rz(-1.5660428) q[2];
sx q[2];
rz(0.63278502) q[2];
rz(1.1535545) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(-0.87483037) q[0];
rz(1.3300928) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777269) q[0];
sx q[0];
rz(-0.3361055) q[0];
sx q[0];
rz(2.0396114) q[0];
rz(-pi) q[1];
rz(-0.61727662) q[2];
sx q[2];
rz(-1.0059788) q[2];
sx q[2];
rz(-1.1519943) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1332902) q[1];
sx q[1];
rz(-1.3160719) q[1];
sx q[1];
rz(-0.50526527) q[1];
x q[2];
rz(-3.0619377) q[3];
sx q[3];
rz(-2.0783391) q[3];
sx q[3];
rz(-1.5505276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76628768) q[0];
sx q[0];
rz(-1.5910609) q[0];
sx q[0];
rz(0.91039175) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(2.8667563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1823605) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(-1.6596646) q[0];
rz(-pi) q[1];
rz(-1.574013) q[2];
sx q[2];
rz(-2.6811757) q[2];
sx q[2];
rz(-0.38052961) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6064925) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(2.5938354) q[1];
x q[2];
rz(1.4196017) q[3];
sx q[3];
rz(-2.4393743) q[3];
sx q[3];
rz(0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3466907) q[2];
sx q[2];
rz(-2.0236423) q[2];
sx q[2];
rz(-1.6332731) q[2];
rz(-1.1446965) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65790025) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(-2.143798) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.654637) q[1];
sx q[1];
rz(-1.516974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8109587) q[0];
sx q[0];
rz(-1.9348382) q[0];
sx q[0];
rz(-0.90322687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60913182) q[2];
sx q[2];
rz(-2.5172148) q[2];
sx q[2];
rz(-0.49027157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5114054) q[1];
sx q[1];
rz(-1.9747707) q[1];
sx q[1];
rz(1.1955111) q[1];
rz(2.153271) q[3];
sx q[3];
rz(-2.0867996) q[3];
sx q[3];
rz(-0.9048681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(2.8175763) q[2];
rz(-1.8185395) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76535392) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(1.8776241) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.2011386) q[1];
sx q[1];
rz(-0.34067672) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7621988) q[0];
sx q[0];
rz(-1.5389991) q[0];
sx q[0];
rz(1.5048774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75603007) q[2];
sx q[2];
rz(-1.3544193) q[2];
sx q[2];
rz(1.1083958) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1813982) q[1];
sx q[1];
rz(-1.3409233) q[1];
sx q[1];
rz(2.6541436) q[1];
x q[2];
rz(2.0601222) q[3];
sx q[3];
rz(-1.9483856) q[3];
sx q[3];
rz(-1.1946354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.95057758) q[2];
sx q[2];
rz(-2.5577929) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(0.54780444) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1691549) q[0];
sx q[0];
rz(-1.4470402) q[0];
sx q[0];
rz(-0.34564885) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-0.47880104) q[1];
sx q[1];
rz(2.6766434) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3195254) q[0];
sx q[0];
rz(-1.3789346) q[0];
sx q[0];
rz(2.0454387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49352383) q[2];
sx q[2];
rz(-2.0888622) q[2];
sx q[2];
rz(-1.354419) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.17864922) q[1];
sx q[1];
rz(-1.4649676) q[1];
sx q[1];
rz(2.9220198) q[1];
x q[2];
rz(0.025860272) q[3];
sx q[3];
rz(-1.6169294) q[3];
sx q[3];
rz(-2.9010454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1376301) q[2];
sx q[2];
rz(-1.6370862) q[2];
sx q[2];
rz(0.31759343) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-2.8542744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(-0.28433329) q[0];
rz(-0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5057482) q[0];
sx q[0];
rz(-0.76857476) q[0];
sx q[0];
rz(1.4710674) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3381091) q[2];
sx q[2];
rz(-0.33499559) q[2];
sx q[2];
rz(-1.6840881) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6137177) q[1];
sx q[1];
rz(-0.67786874) q[1];
sx q[1];
rz(1.8498969) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7275229) q[3];
sx q[3];
rz(-1.6046451) q[3];
sx q[3];
rz(1.0364929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-0.58319432) q[3];
sx q[3];
rz(-1.1116894) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(-0.051368512) q[0];
rz(2.2180166) q[1];
sx q[1];
rz(-0.66134614) q[1];
sx q[1];
rz(0.87402469) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6286205) q[0];
sx q[0];
rz(-0.7681094) q[0];
sx q[0];
rz(-0.012499768) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85272312) q[2];
sx q[2];
rz(-0.52041473) q[2];
sx q[2];
rz(-2.4783217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.45174949) q[1];
sx q[1];
rz(-2.6362231) q[1];
sx q[1];
rz(-1.7187353) q[1];
rz(-pi) q[2];
rz(2.6690528) q[3];
sx q[3];
rz(-2.4774385) q[3];
sx q[3];
rz(-0.00045517552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.727227) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.254436) q[3];
sx q[3];
rz(1.363389) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1353564) q[0];
sx q[0];
rz(-2.0993435) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(1.9627409) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26319474) q[0];
sx q[0];
rz(-2.7976755) q[0];
sx q[0];
rz(-0.92349903) q[0];
rz(-pi) q[1];
rz(-0.63209052) q[2];
sx q[2];
rz(-2.9846016) q[2];
sx q[2];
rz(2.0096411) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2420173) q[1];
sx q[1];
rz(-1.3841108) q[1];
sx q[1];
rz(0.94355299) q[1];
x q[2];
rz(-0.52272777) q[3];
sx q[3];
rz(-1.3739112) q[3];
sx q[3];
rz(-2.6678391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.05802352) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(-2.2422092) q[2];
rz(0.87456885) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62109229) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-2.3836366) q[1];
sx q[1];
rz(-1.6592204) q[1];
sx q[1];
rz(1.6827676) q[1];
rz(-1.6451251) q[2];
sx q[2];
rz(-0.44114124) q[2];
sx q[2];
rz(0.91976358) q[2];
rz(0.27579565) q[3];
sx q[3];
rz(-1.4281359) q[3];
sx q[3];
rz(1.7208163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
