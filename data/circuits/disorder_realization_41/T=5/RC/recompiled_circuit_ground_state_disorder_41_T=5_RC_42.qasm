OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5858894) q[0];
sx q[0];
rz(-1.1893505) q[0];
sx q[0];
rz(-1.1638292) q[0];
rz(-3.8883348) q[1];
sx q[1];
rz(0.26489869) q[1];
sx q[1];
rz(12.737552) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88306016) q[0];
sx q[0];
rz(-0.2685606) q[0];
sx q[0];
rz(-2.1826151) q[0];
rz(-pi) q[1];
rz(2.4462893) q[2];
sx q[2];
rz(-2.2334263) q[2];
sx q[2];
rz(-0.60328874) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4894152) q[1];
sx q[1];
rz(-2.846183) q[1];
sx q[1];
rz(1.5728463) q[1];
rz(-pi) q[2];
rz(3.1247507) q[3];
sx q[3];
rz(-1.5762268) q[3];
sx q[3];
rz(0.098244103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.050194) q[2];
sx q[2];
rz(-1.564743) q[2];
sx q[2];
rz(-2.0719299) q[2];
rz(-1.7012677) q[3];
sx q[3];
rz(-1.0245208) q[3];
sx q[3];
rz(1.1721771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9656203) q[0];
sx q[0];
rz(-1.6930641) q[0];
sx q[0];
rz(-2.5536221) q[0];
rz(1.8486456) q[1];
sx q[1];
rz(-0.99423948) q[1];
sx q[1];
rz(0.15486823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706544) q[0];
sx q[0];
rz(-1.0736671) q[0];
sx q[0];
rz(-1.4488104) q[0];
x q[1];
rz(-1.4243519) q[2];
sx q[2];
rz(-1.3322209) q[2];
sx q[2];
rz(-1.5244791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1701291) q[1];
sx q[1];
rz(-1.2358574) q[1];
sx q[1];
rz(-2.0667124) q[1];
rz(-1.7072466) q[3];
sx q[3];
rz(-1.7810571) q[3];
sx q[3];
rz(-2.7790359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64314848) q[2];
sx q[2];
rz(-2.4109106) q[2];
sx q[2];
rz(0.11239642) q[2];
rz(-0.82320881) q[3];
sx q[3];
rz(-1.221849) q[3];
sx q[3];
rz(-0.52197758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68558973) q[0];
sx q[0];
rz(-1.0887479) q[0];
sx q[0];
rz(-2.1734557) q[0];
rz(2.0227506) q[1];
sx q[1];
rz(-1.6066931) q[1];
sx q[1];
rz(1.3333837) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80305014) q[0];
sx q[0];
rz(-2.1678939) q[0];
sx q[0];
rz(0.19495585) q[0];
rz(-pi) q[1];
rz(-2.0393665) q[2];
sx q[2];
rz(-1.0194687) q[2];
sx q[2];
rz(2.9561549) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81840491) q[1];
sx q[1];
rz(-1.0099704) q[1];
sx q[1];
rz(-3.0238999) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9288238) q[3];
sx q[3];
rz(-1.6505401) q[3];
sx q[3];
rz(2.2984576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1916888) q[2];
sx q[2];
rz(-1.0368985) q[2];
sx q[2];
rz(-2.8766768) q[2];
rz(-2.8252937) q[3];
sx q[3];
rz(-0.33360544) q[3];
sx q[3];
rz(1.3914289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3259657) q[0];
sx q[0];
rz(-0.2206603) q[0];
sx q[0];
rz(-1.6478446) q[0];
rz(-1.4119459) q[1];
sx q[1];
rz(-2.1144046) q[1];
sx q[1];
rz(0.6573917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9566204) q[0];
sx q[0];
rz(-1.6746542) q[0];
sx q[0];
rz(1.64479) q[0];
rz(-pi) q[1];
rz(0.78423402) q[2];
sx q[2];
rz(-2.0080697) q[2];
sx q[2];
rz(-1.0904877) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77995357) q[1];
sx q[1];
rz(-1.3758389) q[1];
sx q[1];
rz(-3.1038398) q[1];
x q[2];
rz(-2.3991556) q[3];
sx q[3];
rz(-0.95177878) q[3];
sx q[3];
rz(-0.05832626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9221981) q[2];
sx q[2];
rz(-1.9191091) q[2];
sx q[2];
rz(2.738319) q[2];
rz(-1.203631) q[3];
sx q[3];
rz(-1.5091242) q[3];
sx q[3];
rz(-0.50886124) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4598684) q[0];
sx q[0];
rz(-2.7175856) q[0];
sx q[0];
rz(2.5176609) q[0];
rz(-2.4195747) q[1];
sx q[1];
rz(-1.3474418) q[1];
sx q[1];
rz(-0.10399476) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063833916) q[0];
sx q[0];
rz(-0.88853271) q[0];
sx q[0];
rz(1.7535249) q[0];
x q[1];
rz(-1.2243168) q[2];
sx q[2];
rz(-1.9110381) q[2];
sx q[2];
rz(0.97803309) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7590569) q[1];
sx q[1];
rz(-1.0604572) q[1];
sx q[1];
rz(-2.3063763) q[1];
rz(0.72428836) q[3];
sx q[3];
rz(-0.87832993) q[3];
sx q[3];
rz(-2.9812733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7907052) q[2];
sx q[2];
rz(-2.4217589) q[2];
sx q[2];
rz(2.3525815) q[2];
rz(-1.2645432) q[3];
sx q[3];
rz(-1.0435373) q[3];
sx q[3];
rz(2.1544382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8378545) q[0];
sx q[0];
rz(-0.721295) q[0];
sx q[0];
rz(2.9413132) q[0];
rz(1.6750977) q[1];
sx q[1];
rz(-1.8148345) q[1];
sx q[1];
rz(1.6237578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90676722) q[0];
sx q[0];
rz(-1.2041766) q[0];
sx q[0];
rz(2.5983222) q[0];
rz(-pi) q[1];
rz(-3.0110961) q[2];
sx q[2];
rz(-2.0708212) q[2];
sx q[2];
rz(2.6967626) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1863757) q[1];
sx q[1];
rz(-2.5506335) q[1];
sx q[1];
rz(0.0062828961) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79849859) q[3];
sx q[3];
rz(-1.3991465) q[3];
sx q[3];
rz(2.5018179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5184861) q[2];
sx q[2];
rz(-0.97272626) q[2];
sx q[2];
rz(-2.8002807) q[2];
rz(2.2845279) q[3];
sx q[3];
rz(-1.9144446) q[3];
sx q[3];
rz(-2.6763693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0905867) q[0];
sx q[0];
rz(-2.107928) q[0];
sx q[0];
rz(1.2475913) q[0];
rz(-0.11820758) q[1];
sx q[1];
rz(-1.2756313) q[1];
sx q[1];
rz(-0.035621312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0220862) q[0];
sx q[0];
rz(-1.6883487) q[0];
sx q[0];
rz(-3.1189128) q[0];
rz(-1.1505736) q[2];
sx q[2];
rz(-2.7033119) q[2];
sx q[2];
rz(-0.29017236) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5069711) q[1];
sx q[1];
rz(-1.8160607) q[1];
sx q[1];
rz(0.62758201) q[1];
x q[2];
rz(1.692996) q[3];
sx q[3];
rz(-0.79334337) q[3];
sx q[3];
rz(-2.6650708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15141618) q[2];
sx q[2];
rz(-1.6665062) q[2];
sx q[2];
rz(-2.0659633) q[2];
rz(2.6208124) q[3];
sx q[3];
rz(-2.5213089) q[3];
sx q[3];
rz(-1.5889997) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90466475) q[0];
sx q[0];
rz(-1.0228782) q[0];
sx q[0];
rz(-2.2383595) q[0];
rz(-1.5029933) q[1];
sx q[1];
rz(-1.4475977) q[1];
sx q[1];
rz(1.8797871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98760688) q[0];
sx q[0];
rz(-1.7075451) q[0];
sx q[0];
rz(1.8022949) q[0];
rz(1.1177344) q[2];
sx q[2];
rz(-1.0894962) q[2];
sx q[2];
rz(-2.520176) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8324229) q[1];
sx q[1];
rz(-1.6214402) q[1];
sx q[1];
rz(-3.0502351) q[1];
rz(-2.0715967) q[3];
sx q[3];
rz(-2.1163995) q[3];
sx q[3];
rz(2.9814482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5791851) q[2];
sx q[2];
rz(-0.6520485) q[2];
sx q[2];
rz(-0.11422608) q[2];
rz(-1.3769897) q[3];
sx q[3];
rz(-2.0012794) q[3];
sx q[3];
rz(0.50447869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1219516) q[0];
sx q[0];
rz(-2.1421102) q[0];
sx q[0];
rz(1.9695388) q[0];
rz(-1.7578341) q[1];
sx q[1];
rz(-1.9969321) q[1];
sx q[1];
rz(3.0269472) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31466904) q[0];
sx q[0];
rz(-2.0631587) q[0];
sx q[0];
rz(-2.5040202) q[0];
x q[1];
rz(-1.6646321) q[2];
sx q[2];
rz(-2.0619446) q[2];
sx q[2];
rz(-2.0744086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9108565) q[1];
sx q[1];
rz(-1.641526) q[1];
sx q[1];
rz(-1.8190228) q[1];
rz(-pi) q[2];
rz(-2.3451553) q[3];
sx q[3];
rz(-1.49545) q[3];
sx q[3];
rz(3.0925261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5454448) q[2];
sx q[2];
rz(-2.9344276) q[2];
sx q[2];
rz(0.88085112) q[2];
rz(0.43978575) q[3];
sx q[3];
rz(-1.9763016) q[3];
sx q[3];
rz(-2.0161276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83774829) q[0];
sx q[0];
rz(-1.6551908) q[0];
sx q[0];
rz(3.0434171) q[0];
rz(-0.57442609) q[1];
sx q[1];
rz(-1.6822633) q[1];
sx q[1];
rz(0.59304768) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2132872) q[0];
sx q[0];
rz(-0.79003626) q[0];
sx q[0];
rz(1.1571048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93043987) q[2];
sx q[2];
rz(-2.499492) q[2];
sx q[2];
rz(-0.053089945) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7465311) q[1];
sx q[1];
rz(-1.5840925) q[1];
sx q[1];
rz(-0.011237267) q[1];
rz(-pi) q[2];
rz(1.9022868) q[3];
sx q[3];
rz(-2.176729) q[3];
sx q[3];
rz(-0.87406033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9888931) q[2];
sx q[2];
rz(-1.1836735) q[2];
sx q[2];
rz(2.2031671) q[2];
rz(1.7484131) q[3];
sx q[3];
rz(-1.8311071) q[3];
sx q[3];
rz(-0.2383298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68168454) q[0];
sx q[0];
rz(-0.60612283) q[0];
sx q[0];
rz(-1.7932307) q[0];
rz(1.0472736) q[1];
sx q[1];
rz(-0.92887639) q[1];
sx q[1];
rz(2.1705719) q[1];
rz(3.0964839) q[2];
sx q[2];
rz(-1.5061629) q[2];
sx q[2];
rz(-0.20284222) q[2];
rz(-2.9592414) q[3];
sx q[3];
rz(-1.3538881) q[3];
sx q[3];
rz(2.8543579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
