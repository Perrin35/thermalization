OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.85668808) q[0];
sx q[0];
rz(4.9520725) q[0];
sx q[0];
rz(7.0926275) q[0];
rz(-2.5419905) q[1];
sx q[1];
rz(2.0166346) q[1];
sx q[1];
rz(12.044608) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4057465) q[0];
sx q[0];
rz(-3.0015115) q[0];
sx q[0];
rz(1.5790042) q[0];
rz(-pi) q[1];
rz(-0.56057616) q[2];
sx q[2];
rz(-1.3180705) q[2];
sx q[2];
rz(1.8552903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.45568902) q[1];
sx q[1];
rz(-2.0200811) q[1];
sx q[1];
rz(-0.30251578) q[1];
x q[2];
rz(1.4843602) q[3];
sx q[3];
rz(-2.1051791) q[3];
sx q[3];
rz(-2.5680281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9894422) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(1.36261) q[2];
rz(-1.4710434) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(-0.41489261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3926369) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(2.0718527) q[0];
rz(-2.3906129) q[1];
sx q[1];
rz(-0.90194482) q[1];
sx q[1];
rz(2.353207) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8892944) q[0];
sx q[0];
rz(-1.7579334) q[0];
sx q[0];
rz(-2.5146913) q[0];
rz(-0.66547243) q[2];
sx q[2];
rz(-1.5147527) q[2];
sx q[2];
rz(1.9201345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0505817) q[1];
sx q[1];
rz(-1.3772703) q[1];
sx q[1];
rz(1.0440318) q[1];
x q[2];
rz(-1.1343003) q[3];
sx q[3];
rz(-1.1009163) q[3];
sx q[3];
rz(1.7738455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1145733) q[2];
sx q[2];
rz(-0.35832778) q[2];
sx q[2];
rz(-0.9066073) q[2];
rz(1.00057) q[3];
sx q[3];
rz(-2.0228736) q[3];
sx q[3];
rz(0.51893836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24666102) q[0];
sx q[0];
rz(-0.38235679) q[0];
sx q[0];
rz(1.3038127) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.3273393) q[1];
sx q[1];
rz(1.4132285) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50158635) q[0];
sx q[0];
rz(-1.4297856) q[0];
sx q[0];
rz(1.4835351) q[0];
x q[1];
rz(-2.5557466) q[2];
sx q[2];
rz(-2.4163462) q[2];
sx q[2];
rz(-2.9278529) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4354463) q[1];
sx q[1];
rz(-2.6455952) q[1];
sx q[1];
rz(2.7643327) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85573825) q[3];
sx q[3];
rz(-1.6659587) q[3];
sx q[3];
rz(-2.9449938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.383519) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(-3.0124532) q[2];
rz(-0.50390759) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(-0.57607877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1103766) q[0];
sx q[0];
rz(-1.2126558) q[0];
sx q[0];
rz(2.2953798) q[0];
rz(0.57202488) q[1];
sx q[1];
rz(-1.5049479) q[1];
sx q[1];
rz(-1.9073073) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0219824) q[0];
sx q[0];
rz(-1.2391157) q[0];
sx q[0];
rz(-1.673438) q[0];
rz(-pi) q[1];
rz(-0.28469601) q[2];
sx q[2];
rz(-1.1276334) q[2];
sx q[2];
rz(0.62847947) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14951359) q[1];
sx q[1];
rz(-1.7798335) q[1];
sx q[1];
rz(-1.6531709) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3448496) q[3];
sx q[3];
rz(-1.9588184) q[3];
sx q[3];
rz(-1.4772367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7531551) q[2];
sx q[2];
rz(-0.93418241) q[2];
sx q[2];
rz(1.4446806) q[2];
rz(-1.8709024) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-1.1893893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020849) q[0];
sx q[0];
rz(-1.6725699) q[0];
sx q[0];
rz(-1.0953267) q[0];
rz(-0.0630088) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(2.1953886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61990935) q[0];
sx q[0];
rz(-2.5961726) q[0];
sx q[0];
rz(-1.5167461) q[0];
rz(-2.5054362) q[2];
sx q[2];
rz(-1.2604179) q[2];
sx q[2];
rz(-2.7458422) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7530147) q[1];
sx q[1];
rz(-1.565897) q[1];
sx q[1];
rz(-2.6304662) q[1];
rz(-2.0949515) q[3];
sx q[3];
rz(-0.44455179) q[3];
sx q[3];
rz(-1.1994282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.63518628) q[2];
sx q[2];
rz(-2.0621641) q[2];
sx q[2];
rz(1.5080473) q[2];
rz(-2.9888198) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(2.2085786) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2372811) q[0];
sx q[0];
rz(-0.58552423) q[0];
sx q[0];
rz(3.0329419) q[0];
rz(3.0142504) q[1];
sx q[1];
rz(-1.8904949) q[1];
sx q[1];
rz(0.88643518) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21289794) q[0];
sx q[0];
rz(-1.587364) q[0];
sx q[0];
rz(-2.6221656) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3662874) q[2];
sx q[2];
rz(-1.5123774) q[2];
sx q[2];
rz(1.2893334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.242524) q[1];
sx q[1];
rz(-1.3521114) q[1];
sx q[1];
rz(-0.058700801) q[1];
rz(-pi) q[2];
rz(2.9720794) q[3];
sx q[3];
rz(-1.4614786) q[3];
sx q[3];
rz(1.6131356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2966557) q[2];
sx q[2];
rz(-0.79534328) q[2];
sx q[2];
rz(-0.33357683) q[2];
rz(-0.9345471) q[3];
sx q[3];
rz(-2.3881674) q[3];
sx q[3];
rz(-0.88753382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4535256) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(0.30853477) q[0];
rz(2.535179) q[1];
sx q[1];
rz(-2.2589222) q[1];
sx q[1];
rz(-2.6920614) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7768984) q[0];
sx q[0];
rz(-2.6230271) q[0];
sx q[0];
rz(-3.0125951) q[0];
rz(-2.3541127) q[2];
sx q[2];
rz(-1.228294) q[2];
sx q[2];
rz(-0.94405789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36360221) q[1];
sx q[1];
rz(-2.5153179) q[1];
sx q[1];
rz(1.7454778) q[1];
rz(-pi) q[2];
rz(0.54517643) q[3];
sx q[3];
rz(-0.94674331) q[3];
sx q[3];
rz(1.6069309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.061607925) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(1.1770581) q[2];
rz(2.4401149) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(0.85566068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6845282) q[0];
sx q[0];
rz(-2.6236911) q[0];
sx q[0];
rz(-1.6494226) q[0];
rz(2.0860784) q[1];
sx q[1];
rz(-1.5857453) q[1];
sx q[1];
rz(-2.7141056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1558091) q[0];
sx q[0];
rz(-2.0844578) q[0];
sx q[0];
rz(-2.2751121) q[0];
rz(1.4094818) q[2];
sx q[2];
rz(-1.9967134) q[2];
sx q[2];
rz(-0.49954712) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44459773) q[1];
sx q[1];
rz(-0.31365221) q[1];
sx q[1];
rz(-3.0693377) q[1];
rz(-0.70722039) q[3];
sx q[3];
rz(-0.40509352) q[3];
sx q[3];
rz(1.9735379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9588354) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(-1.2048362) q[2];
rz(-3.0086573) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(2.0601864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28250113) q[0];
sx q[0];
rz(-1.7318672) q[0];
sx q[0];
rz(-2.6891563) q[0];
rz(0.85894194) q[1];
sx q[1];
rz(-0.57933885) q[1];
sx q[1];
rz(-1.4525684) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70266028) q[0];
sx q[0];
rz(-1.0947252) q[0];
sx q[0];
rz(-2.4804546) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33272393) q[2];
sx q[2];
rz(-2.1004268) q[2];
sx q[2];
rz(0.67654787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4367797) q[1];
sx q[1];
rz(-1.5734993) q[1];
sx q[1];
rz(2.0767434) q[1];
x q[2];
rz(-0.77073578) q[3];
sx q[3];
rz(-0.49388921) q[3];
sx q[3];
rz(1.0115185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23898827) q[2];
sx q[2];
rz(-1.8391823) q[2];
sx q[2];
rz(0.39546173) q[2];
rz(1.0836733) q[3];
sx q[3];
rz(-0.62267059) q[3];
sx q[3];
rz(-1.7284988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85635066) q[0];
sx q[0];
rz(-3.0265891) q[0];
sx q[0];
rz(1.2913936) q[0];
rz(0.77766386) q[1];
sx q[1];
rz(-1.5823369) q[1];
sx q[1];
rz(1.8596328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34436953) q[0];
sx q[0];
rz(-0.31569052) q[0];
sx q[0];
rz(-2.3395721) q[0];
rz(1.3488814) q[2];
sx q[2];
rz(-2.2063428) q[2];
sx q[2];
rz(-1.2254305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.79445542) q[1];
sx q[1];
rz(-2.0973118) q[1];
sx q[1];
rz(1.305425) q[1];
rz(-pi) q[2];
rz(-1.3814209) q[3];
sx q[3];
rz(-0.52703372) q[3];
sx q[3];
rz(-0.98637146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.97734863) q[2];
sx q[2];
rz(-0.35173309) q[2];
sx q[2];
rz(-0.39883167) q[2];
rz(-0.93419689) q[3];
sx q[3];
rz(-2.4103006) q[3];
sx q[3];
rz(-0.092223316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4451404) q[0];
sx q[0];
rz(-2.2284989) q[0];
sx q[0];
rz(-3.1365119) q[0];
rz(-0.65771163) q[1];
sx q[1];
rz(-1.7291768) q[1];
sx q[1];
rz(-1.9531858) q[1];
rz(-1.2451386) q[2];
sx q[2];
rz(-0.79380582) q[2];
sx q[2];
rz(0.028887916) q[2];
rz(-1.4767968) q[3];
sx q[3];
rz(-0.58782676) q[3];
sx q[3];
rz(1.6997433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
