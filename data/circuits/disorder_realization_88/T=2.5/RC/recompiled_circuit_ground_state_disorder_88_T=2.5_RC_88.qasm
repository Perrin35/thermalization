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
rz(0.59960214) q[1];
sx q[1];
rz(-2.0166346) q[1];
sx q[1];
rz(-0.52176276) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7441352) q[0];
sx q[0];
rz(-1.7108727) q[0];
sx q[0];
rz(-0.0011573275) q[0];
rz(-2.6894301) q[2];
sx q[2];
rz(-0.60930646) q[2];
sx q[2];
rz(0.094560187) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45568902) q[1];
sx q[1];
rz(-1.1215116) q[1];
sx q[1];
rz(2.8390769) q[1];
x q[2];
rz(2.9967446) q[3];
sx q[3];
rz(-0.54065767) q[3];
sx q[3];
rz(0.40504211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9894422) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(1.7789827) q[2];
rz(1.6705492) q[3];
sx q[3];
rz(-1.5971284) q[3];
sx q[3];
rz(-2.7267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3926369) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(-1.0697399) q[0];
rz(2.3906129) q[1];
sx q[1];
rz(-2.2396478) q[1];
sx q[1];
rz(2.353207) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9570558) q[0];
sx q[0];
rz(-2.1851077) q[0];
sx q[0];
rz(-1.8004824) q[0];
rz(-0.66547243) q[2];
sx q[2];
rz(-1.5147527) q[2];
sx q[2];
rz(-1.2214582) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.091011) q[1];
sx q[1];
rz(-1.7643223) q[1];
sx q[1];
rz(-1.0440318) q[1];
rz(-pi) q[2];
rz(2.0072924) q[3];
sx q[3];
rz(-1.1009163) q[3];
sx q[3];
rz(1.7738455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.02701935) q[2];
sx q[2];
rz(-0.35832778) q[2];
sx q[2];
rz(-2.2349854) q[2];
rz(-2.1410227) q[3];
sx q[3];
rz(-1.118719) q[3];
sx q[3];
rz(-0.51893836) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24666102) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(1.83778) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(1.7283641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055186836) q[0];
sx q[0];
rz(-0.16567437) q[0];
sx q[0];
rz(0.55060668) q[0];
rz(-pi) q[1];
rz(-2.0264852) q[2];
sx q[2];
rz(-2.1564061) q[2];
sx q[2];
rz(0.93916303) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6709137) q[1];
sx q[1];
rz(-1.7470198) q[1];
sx q[1];
rz(2.6754968) q[1];
rz(-0.85573825) q[3];
sx q[3];
rz(-1.475634) q[3];
sx q[3];
rz(-0.19659886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7580737) q[2];
sx q[2];
rz(-1.588409) q[2];
sx q[2];
rz(-0.12913945) q[2];
rz(-2.6376851) q[3];
sx q[3];
rz(-1.4174856) q[3];
sx q[3];
rz(2.5655139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.031216089) q[0];
sx q[0];
rz(-1.9289368) q[0];
sx q[0];
rz(-0.84621286) q[0];
rz(2.5695678) q[1];
sx q[1];
rz(-1.6366448) q[1];
sx q[1];
rz(-1.9073073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0219824) q[0];
sx q[0];
rz(-1.2391157) q[0];
sx q[0];
rz(1.4681547) q[0];
rz(2.8568966) q[2];
sx q[2];
rz(-1.1276334) q[2];
sx q[2];
rz(0.62847947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.438414) q[1];
sx q[1];
rz(-1.6513737) q[1];
sx q[1];
rz(0.2097278) q[1];
rz(-2.6223141) q[3];
sx q[3];
rz(-2.2745345) q[3];
sx q[3];
rz(-0.26070538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7531551) q[2];
sx q[2];
rz(-0.93418241) q[2];
sx q[2];
rz(1.6969121) q[2];
rz(1.2706903) q[3];
sx q[3];
rz(-1.5710187) q[3];
sx q[3];
rz(-1.9522033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.63950771) q[0];
sx q[0];
rz(-1.4690228) q[0];
sx q[0];
rz(2.0462659) q[0];
rz(3.0785839) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(-0.94620401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5216833) q[0];
sx q[0];
rz(-0.54542002) q[0];
sx q[0];
rz(-1.6248466) q[0];
x q[1];
rz(-2.5054362) q[2];
sx q[2];
rz(-1.8811748) q[2];
sx q[2];
rz(2.7458422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19095382) q[1];
sx q[1];
rz(-0.51114782) q[1];
sx q[1];
rz(0.010015476) q[1];
x q[2];
rz(-0.23403731) q[3];
sx q[3];
rz(-1.952301) q[3];
sx q[3];
rz(-0.62988076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63518628) q[2];
sx q[2];
rz(-1.0794285) q[2];
sx q[2];
rz(-1.5080473) q[2];
rz(0.15277282) q[3];
sx q[3];
rz(-1.2341876) q[3];
sx q[3];
rz(2.2085786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9286947) q[0];
sx q[0];
rz(-1.5542287) q[0];
sx q[0];
rz(2.6221656) q[0];
rz(-pi) q[1];
rz(-0.059659307) q[2];
sx q[2];
rz(-1.7749514) q[2];
sx q[2];
rz(2.8722389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.163609) q[1];
sx q[1];
rz(-2.915288) q[1];
sx q[1];
rz(1.8288906) q[1];
rz(-2.9720794) q[3];
sx q[3];
rz(-1.4614786) q[3];
sx q[3];
rz(-1.6131356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.844937) q[2];
sx q[2];
rz(-0.79534328) q[2];
sx q[2];
rz(-2.8080158) q[2];
rz(-2.2070456) q[3];
sx q[3];
rz(-0.75342527) q[3];
sx q[3];
rz(-0.88753382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4535256) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(0.30853477) q[0];
rz(-2.535179) q[1];
sx q[1];
rz(-0.88267046) q[1];
sx q[1];
rz(-2.6920614) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9251557) q[0];
sx q[0];
rz(-1.0569658) q[0];
sx q[0];
rz(-1.4975182) q[0];
x q[1];
rz(2.3541127) q[2];
sx q[2];
rz(-1.228294) q[2];
sx q[2];
rz(0.94405789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7779904) q[1];
sx q[1];
rz(-2.5153179) q[1];
sx q[1];
rz(-1.3961149) q[1];
x q[2];
rz(0.94664871) q[3];
sx q[3];
rz(-2.3376541) q[3];
sx q[3];
rz(2.338666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.061607925) q[2];
sx q[2];
rz(-1.4915165) q[2];
sx q[2];
rz(-1.1770581) q[2];
rz(0.70147771) q[3];
sx q[3];
rz(-2.8886815) q[3];
sx q[3];
rz(-2.285932) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6845282) q[0];
sx q[0];
rz(-2.6236911) q[0];
sx q[0];
rz(1.6494226) q[0];
rz(2.0860784) q[1];
sx q[1];
rz(-1.5857453) q[1];
sx q[1];
rz(-2.7141056) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2023017) q[0];
sx q[0];
rz(-0.84502506) q[0];
sx q[0];
rz(0.85406749) q[0];
x q[1];
rz(0.43085499) q[2];
sx q[2];
rz(-1.7175894) q[2];
sx q[2];
rz(-2.1374747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0841353) q[1];
sx q[1];
rz(-1.5485207) q[1];
sx q[1];
rz(-0.31288625) q[1];
x q[2];
rz(1.8425138) q[3];
sx q[3];
rz(-1.2665367) q[3];
sx q[3];
rz(2.7226686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9588354) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(1.2048362) q[2];
rz(-3.0086573) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(-1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28250113) q[0];
sx q[0];
rz(-1.4097255) q[0];
sx q[0];
rz(2.6891563) q[0];
rz(-0.85894194) q[1];
sx q[1];
rz(-0.57933885) q[1];
sx q[1];
rz(-1.6890242) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4389324) q[0];
sx q[0];
rz(-2.0468674) q[0];
sx q[0];
rz(2.4804546) q[0];
x q[1];
rz(2.8088687) q[2];
sx q[2];
rz(-1.0411658) q[2];
sx q[2];
rz(-2.4650448) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4367797) q[1];
sx q[1];
rz(-1.5734993) q[1];
sx q[1];
rz(1.0648492) q[1];
x q[2];
rz(2.3708569) q[3];
sx q[3];
rz(-0.49388921) q[3];
sx q[3];
rz(1.0115185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23898827) q[2];
sx q[2];
rz(-1.3024104) q[2];
sx q[2];
rz(2.7461309) q[2];
rz(1.0836733) q[3];
sx q[3];
rz(-2.5189221) q[3];
sx q[3];
rz(1.7284988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.285242) q[0];
sx q[0];
rz(-3.0265891) q[0];
sx q[0];
rz(1.2913936) q[0];
rz(-2.3639288) q[1];
sx q[1];
rz(-1.5823369) q[1];
sx q[1];
rz(1.8596328) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6586475) q[0];
sx q[0];
rz(-1.3532257) q[0];
sx q[0];
rz(1.3402142) q[0];
rz(-0.6474495) q[2];
sx q[2];
rz(-1.7488591) q[2];
sx q[2];
rz(2.9293729) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91207764) q[1];
sx q[1];
rz(-1.34206) q[1];
sx q[1];
rz(-2.5994615) q[1];
rz(-1.3814209) q[3];
sx q[3];
rz(-2.6145589) q[3];
sx q[3];
rz(-2.1552212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.164244) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(-2.742761) q[2];
rz(0.93419689) q[3];
sx q[3];
rz(-2.4103006) q[3];
sx q[3];
rz(0.092223316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6964523) q[0];
sx q[0];
rz(-0.91309375) q[0];
sx q[0];
rz(0.0050807411) q[0];
rz(-2.483881) q[1];
sx q[1];
rz(-1.4124159) q[1];
sx q[1];
rz(1.1884069) q[1];
rz(-1.896454) q[2];
sx q[2];
rz(-2.3477868) q[2];
sx q[2];
rz(-3.1127047) q[2];
rz(1.6647958) q[3];
sx q[3];
rz(-0.58782676) q[3];
sx q[3];
rz(1.6997433) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
