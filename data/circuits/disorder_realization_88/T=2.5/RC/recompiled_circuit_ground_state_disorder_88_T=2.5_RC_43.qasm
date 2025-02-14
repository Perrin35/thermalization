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
rz(-1.3974575) q[0];
sx q[0];
rz(-1.43072) q[0];
sx q[0];
rz(-3.1404353) q[0];
rz(-pi) q[1];
rz(2.5810165) q[2];
sx q[2];
rz(-1.8235221) q[2];
sx q[2];
rz(1.2863024) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45568902) q[1];
sx q[1];
rz(-1.1215116) q[1];
sx q[1];
rz(2.8390769) q[1];
rz(-0.5360236) q[3];
sx q[3];
rz(-1.6451577) q[3];
sx q[3];
rz(-2.1002567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.9894422) q[2];
sx q[2];
rz(-0.70789727) q[2];
sx q[2];
rz(1.7789827) q[2];
rz(-1.6705492) q[3];
sx q[3];
rz(-1.5444642) q[3];
sx q[3];
rz(-2.7267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7489557) q[0];
sx q[0];
rz(-1.3480659) q[0];
sx q[0];
rz(1.0697399) q[0];
rz(-2.3906129) q[1];
sx q[1];
rz(-2.2396478) q[1];
sx q[1];
rz(0.78838563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1845368) q[0];
sx q[0];
rz(-2.1851077) q[0];
sx q[0];
rz(1.3411103) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66547243) q[2];
sx q[2];
rz(-1.5147527) q[2];
sx q[2];
rz(1.2214582) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0505817) q[1];
sx q[1];
rz(-1.7643223) q[1];
sx q[1];
rz(-2.0975608) q[1];
x q[2];
rz(-0.69425868) q[3];
sx q[3];
rz(-2.5115847) q[3];
sx q[3];
rz(-0.97351551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.02701935) q[2];
sx q[2];
rz(-2.7832649) q[2];
sx q[2];
rz(-2.2349854) q[2];
rz(-1.00057) q[3];
sx q[3];
rz(-2.0228736) q[3];
sx q[3];
rz(-0.51893836) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24666102) q[0];
sx q[0];
rz(-2.7592359) q[0];
sx q[0];
rz(-1.3038127) q[0];
rz(1.9815824) q[1];
sx q[1];
rz(-1.8142533) q[1];
sx q[1];
rz(-1.4132285) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055186836) q[0];
sx q[0];
rz(-2.9759183) q[0];
sx q[0];
rz(2.590986) q[0];
rz(-2.0264852) q[2];
sx q[2];
rz(-2.1564061) q[2];
sx q[2];
rz(0.93916303) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47067898) q[1];
sx q[1];
rz(-1.7470198) q[1];
sx q[1];
rz(2.6754968) q[1];
rz(-pi) q[2];
rz(0.85573825) q[3];
sx q[3];
rz(-1.6659587) q[3];
sx q[3];
rz(-0.19659886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.383519) q[2];
sx q[2];
rz(-1.5531837) q[2];
sx q[2];
rz(-3.0124532) q[2];
rz(-2.6376851) q[3];
sx q[3];
rz(-1.7241071) q[3];
sx q[3];
rz(0.57607877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.031216089) q[0];
sx q[0];
rz(-1.9289368) q[0];
sx q[0];
rz(-2.2953798) q[0];
rz(-0.57202488) q[1];
sx q[1];
rz(-1.5049479) q[1];
sx q[1];
rz(-1.2342854) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11961022) q[0];
sx q[0];
rz(-1.9024769) q[0];
sx q[0];
rz(1.673438) q[0];
x q[1];
rz(2.8568966) q[2];
sx q[2];
rz(-2.0139593) q[2];
sx q[2];
rz(-0.62847947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22913904) q[1];
sx q[1];
rz(-0.22446147) q[1];
sx q[1];
rz(-2.771586) q[1];
x q[2];
rz(2.3448496) q[3];
sx q[3];
rz(-1.9588184) q[3];
sx q[3];
rz(-1.6643559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3884376) q[2];
sx q[2];
rz(-2.2074102) q[2];
sx q[2];
rz(-1.4446806) q[2];
rz(-1.2706903) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-1.9522033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5020849) q[0];
sx q[0];
rz(-1.4690228) q[0];
sx q[0];
rz(-2.0462659) q[0];
rz(0.0630088) q[1];
sx q[1];
rz(-1.6687702) q[1];
sx q[1];
rz(0.94620401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444855) q[0];
sx q[0];
rz(-1.5988264) q[0];
sx q[0];
rz(2.1155684) q[0];
rz(0.63615648) q[2];
sx q[2];
rz(-1.2604179) q[2];
sx q[2];
rz(-2.7458422) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7530147) q[1];
sx q[1];
rz(-1.565897) q[1];
sx q[1];
rz(-0.51112643) q[1];
x q[2];
rz(2.0949515) q[3];
sx q[3];
rz(-2.6970409) q[3];
sx q[3];
rz(-1.1994282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5064064) q[2];
sx q[2];
rz(-1.0794285) q[2];
sx q[2];
rz(1.5080473) q[2];
rz(0.15277282) q[3];
sx q[3];
rz(-1.907405) q[3];
sx q[3];
rz(-2.2085786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9043115) q[0];
sx q[0];
rz(-0.58552423) q[0];
sx q[0];
rz(3.0329419) q[0];
rz(0.12734224) q[1];
sx q[1];
rz(-1.8904949) q[1];
sx q[1];
rz(-0.88643518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3868635) q[0];
sx q[0];
rz(-0.51966705) q[0];
sx q[0];
rz(-3.1082252) q[0];
x q[1];
rz(-3.0819333) q[2];
sx q[2];
rz(-1.3666412) q[2];
sx q[2];
rz(2.8722389) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.163609) q[1];
sx q[1];
rz(-2.915288) q[1];
sx q[1];
rz(-1.3127021) q[1];
rz(-pi) q[2];
rz(2.9720794) q[3];
sx q[3];
rz(-1.680114) q[3];
sx q[3];
rz(1.5284571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.844937) q[2];
sx q[2];
rz(-0.79534328) q[2];
sx q[2];
rz(0.33357683) q[2];
rz(-2.2070456) q[3];
sx q[3];
rz(-0.75342527) q[3];
sx q[3];
rz(-0.88753382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4535256) q[0];
sx q[0];
rz(-1.6418566) q[0];
sx q[0];
rz(-0.30853477) q[0];
rz(-2.535179) q[1];
sx q[1];
rz(-0.88267046) q[1];
sx q[1];
rz(0.44953129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7768984) q[0];
sx q[0];
rz(-2.6230271) q[0];
sx q[0];
rz(3.0125951) q[0];
rz(2.6753898) q[2];
sx q[2];
rz(-2.297819) q[2];
sx q[2];
rz(0.94949338) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7923741) q[1];
sx q[1];
rz(-1.6728396) q[1];
sx q[1];
rz(-0.95178638) q[1];
rz(-pi) q[2];
rz(2.1949439) q[3];
sx q[3];
rz(-2.3376541) q[3];
sx q[3];
rz(-2.338666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0799847) q[2];
sx q[2];
rz(-1.6500762) q[2];
sx q[2];
rz(-1.1770581) q[2];
rz(-0.70147771) q[3];
sx q[3];
rz(-0.25291118) q[3];
sx q[3];
rz(0.85566068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-2.6845282) q[0];
sx q[0];
rz(-2.6236911) q[0];
sx q[0];
rz(1.6494226) q[0];
rz(2.0860784) q[1];
sx q[1];
rz(-1.5558473) q[1];
sx q[1];
rz(-0.42748705) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857835) q[0];
sx q[0];
rz(-2.0844578) q[0];
sx q[0];
rz(-2.2751121) q[0];
x q[1];
rz(2.801339) q[2];
sx q[2];
rz(-0.45368735) q[2];
sx q[2];
rz(-0.87475529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.44459773) q[1];
sx q[1];
rz(-0.31365221) q[1];
sx q[1];
rz(0.072254953) q[1];
rz(-0.70722039) q[3];
sx q[3];
rz(-2.7364991) q[3];
sx q[3];
rz(1.1680548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18275729) q[2];
sx q[2];
rz(-1.3311102) q[2];
sx q[2];
rz(1.2048362) q[2];
rz(3.0086573) q[3];
sx q[3];
rz(-3.0429621) q[3];
sx q[3];
rz(1.0814063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
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
rz(-2.2826507) q[1];
sx q[1];
rz(-0.57933885) q[1];
sx q[1];
rz(1.6890242) q[1];
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
rz(0.66113801) q[0];
rz(-pi) q[1];
rz(-0.33272393) q[2];
sx q[2];
rz(-2.1004268) q[2];
sx q[2];
rz(-0.67654787) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4367797) q[1];
sx q[1];
rz(-1.5680934) q[1];
sx q[1];
rz(1.0648492) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77073578) q[3];
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
rz(-2.9026044) q[2];
sx q[2];
rz(-1.3024104) q[2];
sx q[2];
rz(-2.7461309) q[2];
rz(-2.0579193) q[3];
sx q[3];
rz(-0.62267059) q[3];
sx q[3];
rz(-1.7284988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.285242) q[0];
sx q[0];
rz(-0.1150035) q[0];
sx q[0];
rz(1.850199) q[0];
rz(0.77766386) q[1];
sx q[1];
rz(-1.5592557) q[1];
sx q[1];
rz(1.2819598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4829452) q[0];
sx q[0];
rz(-1.3532257) q[0];
sx q[0];
rz(1.3402142) q[0];
rz(-pi) q[1];
rz(2.4941432) q[2];
sx q[2];
rz(-1.7488591) q[2];
sx q[2];
rz(2.9293729) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.229515) q[1];
sx q[1];
rz(-1.7995326) q[1];
sx q[1];
rz(-2.5994615) q[1];
x q[2];
rz(2.090023) q[3];
sx q[3];
rz(-1.4759721) q[3];
sx q[3];
rz(0.42025305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97734863) q[2];
sx q[2];
rz(-2.7898596) q[2];
sx q[2];
rz(-2.742761) q[2];
rz(-0.93419689) q[3];
sx q[3];
rz(-2.4103006) q[3];
sx q[3];
rz(3.0493693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4451404) q[0];
sx q[0];
rz(-0.91309375) q[0];
sx q[0];
rz(0.0050807411) q[0];
rz(-2.483881) q[1];
sx q[1];
rz(-1.4124159) q[1];
sx q[1];
rz(1.1884069) q[1];
rz(0.31455518) q[2];
sx q[2];
rz(-2.3125074) q[2];
sx q[2];
rz(-0.41994647) q[2];
rz(-3.0791238) q[3];
sx q[3];
rz(-2.1556839) q[3];
sx q[3];
rz(1.5869303) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
