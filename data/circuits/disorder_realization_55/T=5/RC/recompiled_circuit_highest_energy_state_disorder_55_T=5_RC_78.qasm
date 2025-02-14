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
rz(0.060184181) q[0];
sx q[0];
rz(4.2005778) q[0];
sx q[0];
rz(8.4146001) q[0];
rz(2.3330359) q[1];
sx q[1];
rz(-2.8454236) q[1];
sx q[1];
rz(-0.30997601) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55059075) q[0];
sx q[0];
rz(-0.43568107) q[0];
sx q[0];
rz(0.26623078) q[0];
rz(-pi) q[1];
rz(1.5484137) q[2];
sx q[2];
rz(-1.7760008) q[2];
sx q[2];
rz(0.6485282) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0204326) q[1];
sx q[1];
rz(-1.3749287) q[1];
sx q[1];
rz(1.2365667) q[1];
rz(-pi) q[2];
rz(-2.9397291) q[3];
sx q[3];
rz(-0.84175693) q[3];
sx q[3];
rz(-1.0224467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4787204) q[2];
sx q[2];
rz(-0.24820776) q[2];
sx q[2];
rz(3.0459246) q[2];
rz(-0.11224789) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2410759) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(0.61479968) q[0];
rz(-2.2531033) q[1];
sx q[1];
rz(-1.5526086) q[1];
sx q[1];
rz(2.6420171) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.963012) q[0];
sx q[0];
rz(-0.38694978) q[0];
sx q[0];
rz(-1.011417) q[0];
x q[1];
rz(0.26878727) q[2];
sx q[2];
rz(-2.6013881) q[2];
sx q[2];
rz(-0.42551431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.412498) q[1];
sx q[1];
rz(-0.67616588) q[1];
sx q[1];
rz(-1.8489643) q[1];
x q[2];
rz(2.1533215) q[3];
sx q[3];
rz(-1.8559905) q[3];
sx q[3];
rz(-2.5799283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8932314) q[2];
sx q[2];
rz(-1.9042559) q[2];
sx q[2];
rz(3.1207747) q[2];
rz(-1.9013532) q[3];
sx q[3];
rz(-1.6720684) q[3];
sx q[3];
rz(-1.1851236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5027387) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(2.8079206) q[0];
rz(-2.4036713) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(0.12942448) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5851916) q[0];
sx q[0];
rz(-1.4117774) q[0];
sx q[0];
rz(-0.88944737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.470181) q[2];
sx q[2];
rz(-1.08403) q[2];
sx q[2];
rz(-2.0109796) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1947965) q[1];
sx q[1];
rz(-1.4490286) q[1];
sx q[1];
rz(3.0508079) q[1];
rz(-pi) q[2];
rz(-2.0468726) q[3];
sx q[3];
rz(-0.84730803) q[3];
sx q[3];
rz(2.0499961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1239329) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(-1.7714436) q[2];
rz(-0.74639368) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0081886) q[0];
sx q[0];
rz(-1.9318102) q[0];
sx q[0];
rz(-2.5764537) q[0];
rz(2.7124229) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(-2.3133004) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9681184) q[0];
sx q[0];
rz(-0.83195639) q[0];
sx q[0];
rz(1.4627187) q[0];
rz(2.0798551) q[2];
sx q[2];
rz(-0.54414302) q[2];
sx q[2];
rz(-2.6119815) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8569717) q[1];
sx q[1];
rz(-1.0880252) q[1];
sx q[1];
rz(1.8716145) q[1];
rz(1.0542442) q[3];
sx q[3];
rz(-1.787546) q[3];
sx q[3];
rz(2.3473489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68507489) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(-1.7100517) q[2];
rz(0.93959129) q[3];
sx q[3];
rz(-0.51274931) q[3];
sx q[3];
rz(-3.140894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13570304) q[0];
sx q[0];
rz(-0.26517427) q[0];
sx q[0];
rz(1.3294719) q[0];
rz(-1.9895408) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(0.00024814127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86289257) q[0];
sx q[0];
rz(-2.8331146) q[0];
sx q[0];
rz(2.2174979) q[0];
rz(2.1049961) q[2];
sx q[2];
rz(-1.6459668) q[2];
sx q[2];
rz(2.6435564) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50284144) q[1];
sx q[1];
rz(-1.0783245) q[1];
sx q[1];
rz(2.7472463) q[1];
rz(0.080623589) q[3];
sx q[3];
rz(-0.73314694) q[3];
sx q[3];
rz(-2.6997363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0033215) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(0.0058343466) q[2];
rz(0.88998574) q[3];
sx q[3];
rz(-2.4599288) q[3];
sx q[3];
rz(-2.0148923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(1.6538612) q[0];
rz(0.030390175) q[1];
sx q[1];
rz(-2.0149714) q[1];
sx q[1];
rz(0.94246513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93118942) q[0];
sx q[0];
rz(-1.7587587) q[0];
sx q[0];
rz(-0.24954777) q[0];
rz(-pi) q[1];
rz(-0.29877383) q[2];
sx q[2];
rz(-1.818383) q[2];
sx q[2];
rz(-2.5556759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16557517) q[1];
sx q[1];
rz(-1.2846795) q[1];
sx q[1];
rz(-2.0598434) q[1];
rz(-1.0951772) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(-1.7847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.59948644) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(-1.8360809) q[2];
rz(0.54715884) q[3];
sx q[3];
rz(-1.9360417) q[3];
sx q[3];
rz(-2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18043537) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(-2.6590977) q[1];
sx q[1];
rz(-1.1588187) q[1];
sx q[1];
rz(0.85711342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6109638) q[0];
sx q[0];
rz(-0.93994323) q[0];
sx q[0];
rz(-1.1683589) q[0];
x q[1];
rz(-1.6202507) q[2];
sx q[2];
rz(-2.2170545) q[2];
sx q[2];
rz(2.6924999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65322058) q[1];
sx q[1];
rz(-1.3626507) q[1];
sx q[1];
rz(-0.4076165) q[1];
rz(1.9590274) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(1.5807815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7414005) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(-2.8391489) q[2];
rz(-0.97964573) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(-1.2790595) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7804467) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(0.030990344) q[0];
rz(2.567645) q[1];
sx q[1];
rz(-1.4731044) q[1];
sx q[1];
rz(-1.8642289) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325186) q[0];
sx q[0];
rz(-1.2379097) q[0];
sx q[0];
rz(1.9247967) q[0];
rz(2.8988016) q[2];
sx q[2];
rz(-1.0848248) q[2];
sx q[2];
rz(2.9832341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0884235) q[1];
sx q[1];
rz(-2.649247) q[1];
sx q[1];
rz(-1.7264992) q[1];
x q[2];
rz(0.061789024) q[3];
sx q[3];
rz(-1.8553858) q[3];
sx q[3];
rz(-2.9348843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.10451) q[2];
sx q[2];
rz(-0.13585486) q[2];
sx q[2];
rz(0.48745298) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-2.2127377) q[3];
sx q[3];
rz(-3.0776183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071851991) q[0];
sx q[0];
rz(-2.6672279) q[0];
sx q[0];
rz(0.35097861) q[0];
rz(2.9674496) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(1.4016271) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49920666) q[0];
sx q[0];
rz(-2.8642352) q[0];
sx q[0];
rz(-0.46210285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3910037) q[2];
sx q[2];
rz(-1.9074512) q[2];
sx q[2];
rz(2.2845993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3199215) q[1];
sx q[1];
rz(-2.955759) q[1];
sx q[1];
rz(1.7286848) q[1];
x q[2];
rz(0.97934874) q[3];
sx q[3];
rz(-2.2205847) q[3];
sx q[3];
rz(1.5486919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1104687) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-2.5049211) q[2];
rz(-2.0311671) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(0.5184263) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7745895) q[0];
sx q[0];
rz(-0.29357266) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2702613) q[0];
sx q[0];
rz(-3.1360021) q[0];
sx q[0];
rz(-0.63951512) q[0];
rz(-pi) q[1];
rz(-0.27711192) q[2];
sx q[2];
rz(-1.2925576) q[2];
sx q[2];
rz(-2.3269175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0088996) q[1];
sx q[1];
rz(-1.6311967) q[1];
sx q[1];
rz(1.6113043) q[1];
rz(0.070399447) q[3];
sx q[3];
rz(-2.681085) q[3];
sx q[3];
rz(-1.2706336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71558636) q[2];
sx q[2];
rz(-2.6435489) q[2];
sx q[2];
rz(-3.0675724) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(0.35416245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114125) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(-0.62803531) q[2];
sx q[2];
rz(-2.5428094) q[2];
sx q[2];
rz(1.6369292) q[2];
rz(-2.0430123) q[3];
sx q[3];
rz(-0.87971148) q[3];
sx q[3];
rz(1.9426027) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
