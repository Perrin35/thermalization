OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.5503791) q[0];
sx q[0];
rz(3.1449218) q[0];
sx q[0];
rz(9.7700906) q[0];
rz(-1.3357063) q[1];
sx q[1];
rz(-0.3392646) q[1];
sx q[1];
rz(-2.8621434) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0246668) q[0];
sx q[0];
rz(-1.7755839) q[0];
sx q[0];
rz(1.8860399) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3596256) q[2];
sx q[2];
rz(-1.8188393) q[2];
sx q[2];
rz(2.9458407) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.63424078) q[1];
sx q[1];
rz(-1.5841286) q[1];
sx q[1];
rz(-1.3807382) q[1];
x q[2];
rz(2.4996098) q[3];
sx q[3];
rz(-1.6025474) q[3];
sx q[3];
rz(-0.25105219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7001069) q[2];
sx q[2];
rz(-1.4725279) q[2];
sx q[2];
rz(2.1842365) q[2];
rz(2.8422614) q[3];
sx q[3];
rz(-0.39756164) q[3];
sx q[3];
rz(-0.41199747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1608202) q[0];
sx q[0];
rz(-0.94962025) q[0];
sx q[0];
rz(0.41369307) q[0];
rz(1.7970239) q[1];
sx q[1];
rz(-0.78318703) q[1];
sx q[1];
rz(-2.5059674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5168415) q[0];
sx q[0];
rz(-0.13953129) q[0];
sx q[0];
rz(1.8952888) q[0];
rz(-pi) q[1];
rz(1.4388678) q[2];
sx q[2];
rz(-1.6172098) q[2];
sx q[2];
rz(-0.012243587) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0483413) q[1];
sx q[1];
rz(-1.3158568) q[1];
sx q[1];
rz(-1.821372) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2222744) q[3];
sx q[3];
rz(-0.98005664) q[3];
sx q[3];
rz(2.7880993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7522493) q[2];
sx q[2];
rz(-1.2717335) q[2];
sx q[2];
rz(-2.6129369) q[2];
rz(-1.2403437) q[3];
sx q[3];
rz(-0.35959187) q[3];
sx q[3];
rz(0.24578978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56025958) q[0];
sx q[0];
rz(-1.9868877) q[0];
sx q[0];
rz(-0.2581968) q[0];
rz(-1.6351581) q[1];
sx q[1];
rz(-2.5841027) q[1];
sx q[1];
rz(2.7071276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397944) q[0];
sx q[0];
rz(-0.5895624) q[0];
sx q[0];
rz(0.88296417) q[0];
rz(-pi) q[1];
rz(-1.4591818) q[2];
sx q[2];
rz(-0.53225213) q[2];
sx q[2];
rz(0.75512952) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9745969) q[1];
sx q[1];
rz(-1.3612862) q[1];
sx q[1];
rz(-0.4619044) q[1];
rz(-2.4934019) q[3];
sx q[3];
rz(-2.18581) q[3];
sx q[3];
rz(-2.0730413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7210641) q[2];
sx q[2];
rz(-1.8177744) q[2];
sx q[2];
rz(-0.055796441) q[2];
rz(-0.93786401) q[3];
sx q[3];
rz(-2.8580229) q[3];
sx q[3];
rz(2.8985033) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043512251) q[0];
sx q[0];
rz(-2.1976017) q[0];
sx q[0];
rz(-3.0010624) q[0];
rz(-0.17164104) q[1];
sx q[1];
rz(-1.834603) q[1];
sx q[1];
rz(-2.8780639) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182501) q[0];
sx q[0];
rz(-1.468987) q[0];
sx q[0];
rz(-2.5725767) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66480555) q[2];
sx q[2];
rz(-1.9812599) q[2];
sx q[2];
rz(2.7523506) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7100704) q[1];
sx q[1];
rz(-0.78774161) q[1];
sx q[1];
rz(2.6452933) q[1];
x q[2];
rz(0.88818355) q[3];
sx q[3];
rz(-2.7215951) q[3];
sx q[3];
rz(-0.94204599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0984829) q[2];
sx q[2];
rz(-0.3442328) q[2];
sx q[2];
rz(-1.2496703) q[2];
rz(-1.194687) q[3];
sx q[3];
rz(-1.0281111) q[3];
sx q[3];
rz(2.4627114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93976218) q[0];
sx q[0];
rz(-1.8562466) q[0];
sx q[0];
rz(0.029065954) q[0];
rz(1.7395696) q[1];
sx q[1];
rz(-0.35750917) q[1];
sx q[1];
rz(2.8129541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11360725) q[0];
sx q[0];
rz(-2.0562045) q[0];
sx q[0];
rz(3.0263682) q[0];
x q[1];
rz(1.7476014) q[2];
sx q[2];
rz(-1.8362152) q[2];
sx q[2];
rz(-0.058942827) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8481816) q[1];
sx q[1];
rz(-1.1383346) q[1];
sx q[1];
rz(1.5007988) q[1];
x q[2];
rz(-2.394649) q[3];
sx q[3];
rz(-1.9455823) q[3];
sx q[3];
rz(-0.67583109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11792004) q[2];
sx q[2];
rz(-2.6866044) q[2];
sx q[2];
rz(-0.16061352) q[2];
rz(1.1462071) q[3];
sx q[3];
rz(-1.3137772) q[3];
sx q[3];
rz(-2.5207991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6495431) q[0];
sx q[0];
rz(-2.2821125) q[0];
sx q[0];
rz(2.3550526) q[0];
rz(0.37711626) q[1];
sx q[1];
rz(-0.84723324) q[1];
sx q[1];
rz(2.699111) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0241644) q[0];
sx q[0];
rz(-2.6167343) q[0];
sx q[0];
rz(2.3470122) q[0];
x q[1];
rz(-2.6857576) q[2];
sx q[2];
rz(-2.5172533) q[2];
sx q[2];
rz(2.3221743) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9128531) q[1];
sx q[1];
rz(-0.40816669) q[1];
sx q[1];
rz(-2.2250882) q[1];
rz(-1.2908859) q[3];
sx q[3];
rz(-0.32752447) q[3];
sx q[3];
rz(0.97400986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68142146) q[2];
sx q[2];
rz(-0.57839102) q[2];
sx q[2];
rz(-0.9712514) q[2];
rz(-2.4462637) q[3];
sx q[3];
rz(-2.9054502) q[3];
sx q[3];
rz(2.7142081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5044395) q[0];
sx q[0];
rz(-1.9313066) q[0];
sx q[0];
rz(2.1321645) q[0];
rz(0.55074739) q[1];
sx q[1];
rz(-1.8661205) q[1];
sx q[1];
rz(-1.1632464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7834085) q[0];
sx q[0];
rz(-1.9005214) q[0];
sx q[0];
rz(0.14915906) q[0];
rz(-pi) q[1];
rz(-2.8909056) q[2];
sx q[2];
rz(-2.2633268) q[2];
sx q[2];
rz(-2.3547949) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5442859) q[1];
sx q[1];
rz(-1.7339216) q[1];
sx q[1];
rz(-2.9189928) q[1];
x q[2];
rz(0.30927741) q[3];
sx q[3];
rz(-0.89865548) q[3];
sx q[3];
rz(-2.6291763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.358868) q[2];
sx q[2];
rz(-1.1512558) q[2];
sx q[2];
rz(-0.7981832) q[2];
rz(-2.362137) q[3];
sx q[3];
rz(-0.53648406) q[3];
sx q[3];
rz(-1.1727758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9629795) q[0];
sx q[0];
rz(-2.6586752) q[0];
sx q[0];
rz(3.0187507) q[0];
rz(-0.12610647) q[1];
sx q[1];
rz(-1.5051196) q[1];
sx q[1];
rz(-1.925148) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0053596) q[0];
sx q[0];
rz(-1.981712) q[0];
sx q[0];
rz(-0.67816011) q[0];
rz(2.8100796) q[2];
sx q[2];
rz(-0.47795948) q[2];
sx q[2];
rz(-0.31994672) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6718037) q[1];
sx q[1];
rz(-0.43699139) q[1];
sx q[1];
rz(0.48204084) q[1];
x q[2];
rz(1.3835195) q[3];
sx q[3];
rz(-2.1547199) q[3];
sx q[3];
rz(0.54959471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.2609666) q[2];
sx q[2];
rz(-2.5239021) q[2];
sx q[2];
rz(3.0984666) q[2];
rz(0.17523266) q[3];
sx q[3];
rz(-2.2641116) q[3];
sx q[3];
rz(-1.5847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6373428) q[0];
sx q[0];
rz(-3.0872587) q[0];
sx q[0];
rz(-0.20877008) q[0];
rz(-1.4978706) q[1];
sx q[1];
rz(-1.6521963) q[1];
sx q[1];
rz(1.1245022) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18692423) q[0];
sx q[0];
rz(-2.9483729) q[0];
sx q[0];
rz(0.78878553) q[0];
x q[1];
rz(-2.8378216) q[2];
sx q[2];
rz(-2.4035932) q[2];
sx q[2];
rz(-0.27429013) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8090325) q[1];
sx q[1];
rz(-0.13339116) q[1];
sx q[1];
rz(-1.5796214) q[1];
x q[2];
rz(-2.1379925) q[3];
sx q[3];
rz(-0.78046679) q[3];
sx q[3];
rz(-1.6571852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1407397) q[2];
sx q[2];
rz(-2.4177987) q[2];
sx q[2];
rz(-0.17803426) q[2];
rz(1.595165) q[3];
sx q[3];
rz(-1.833257) q[3];
sx q[3];
rz(-2.6509638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7889325) q[0];
sx q[0];
rz(-1.1505609) q[0];
sx q[0];
rz(-1.0349405) q[0];
rz(2.3433698) q[1];
sx q[1];
rz(-0.98033506) q[1];
sx q[1];
rz(-0.14990526) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7184188) q[0];
sx q[0];
rz(-2.383854) q[0];
sx q[0];
rz(0.54990479) q[0];
rz(-pi) q[1];
rz(-2.855004) q[2];
sx q[2];
rz(-1.5865241) q[2];
sx q[2];
rz(2.36731) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1356218) q[1];
sx q[1];
rz(-1.0711728) q[1];
sx q[1];
rz(-0.47132229) q[1];
x q[2];
rz(-0.67930119) q[3];
sx q[3];
rz(-1.7656109) q[3];
sx q[3];
rz(-0.65115813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7297111) q[2];
sx q[2];
rz(-1.9103266) q[2];
sx q[2];
rz(2.7330772) q[2];
rz(0.92489964) q[3];
sx q[3];
rz(-2.3291406) q[3];
sx q[3];
rz(-2.6598721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2395441) q[0];
sx q[0];
rz(-1.5477381) q[0];
sx q[0];
rz(1.5292194) q[0];
rz(-1.3760024) q[1];
sx q[1];
rz(-1.1767495) q[1];
sx q[1];
rz(-1.8935988) q[1];
rz(-2.3464936) q[2];
sx q[2];
rz(-1.3394525) q[2];
sx q[2];
rz(0.48223334) q[2];
rz(2.8077447) q[3];
sx q[3];
rz(-0.6469938) q[3];
sx q[3];
rz(-2.7713431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
