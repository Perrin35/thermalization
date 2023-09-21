OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2919579) q[0];
sx q[0];
rz(-2.7014974) q[0];
sx q[0];
rz(3.0043998) q[0];
rz(-1.7358915) q[1];
sx q[1];
rz(-1.403221) q[1];
sx q[1];
rz(2.611673) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7512902) q[0];
sx q[0];
rz(-1.4634702) q[0];
sx q[0];
rz(1.1885378) q[0];
rz(0.34420867) q[2];
sx q[2];
rz(-1.1905626) q[2];
sx q[2];
rz(2.3772079) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9143608) q[1];
sx q[1];
rz(-1.1162317) q[1];
sx q[1];
rz(0.55959065) q[1];
rz(-1.3210117) q[3];
sx q[3];
rz(-0.82818177) q[3];
sx q[3];
rz(0.111655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68937504) q[2];
sx q[2];
rz(-1.3000501) q[2];
sx q[2];
rz(-0.33660647) q[2];
rz(-1.6254788) q[3];
sx q[3];
rz(-2.5879526) q[3];
sx q[3];
rz(1.5256933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7933554) q[0];
sx q[0];
rz(-2.0331148) q[0];
sx q[0];
rz(0.021214699) q[0];
rz(1.1938098) q[1];
sx q[1];
rz(-2.1021011) q[1];
sx q[1];
rz(2.3056727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056136925) q[0];
sx q[0];
rz(-1.5888927) q[0];
sx q[0];
rz(1.4319112) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.732227) q[2];
sx q[2];
rz(-1.0592807) q[2];
sx q[2];
rz(2.6057233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0554725) q[1];
sx q[1];
rz(-1.8430084) q[1];
sx q[1];
rz(0.91113669) q[1];
x q[2];
rz(2.4017176) q[3];
sx q[3];
rz(-1.8679108) q[3];
sx q[3];
rz(1.5418996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2772284) q[2];
sx q[2];
rz(-1.9936864) q[2];
sx q[2];
rz(-1.345984) q[2];
rz(0.35955444) q[3];
sx q[3];
rz(-0.94272009) q[3];
sx q[3];
rz(0.49697044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.7132752) q[0];
sx q[0];
rz(-2.064216) q[0];
sx q[0];
rz(-1.0536449) q[0];
rz(-1.9127649) q[1];
sx q[1];
rz(-1.5412953) q[1];
sx q[1];
rz(2.704481) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99840435) q[0];
sx q[0];
rz(-2.8371053) q[0];
sx q[0];
rz(1.8633153) q[0];
rz(-pi) q[1];
x q[1];
rz(1.978546) q[2];
sx q[2];
rz(-1.7951269) q[2];
sx q[2];
rz(-0.11220223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25698173) q[1];
sx q[1];
rz(-1.7092488) q[1];
sx q[1];
rz(-2.715766) q[1];
x q[2];
rz(-1.6821074) q[3];
sx q[3];
rz(-0.49735945) q[3];
sx q[3];
rz(1.5670083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1221216) q[2];
sx q[2];
rz(-0.78142587) q[2];
sx q[2];
rz(1.0220698) q[2];
rz(-1.2381037) q[3];
sx q[3];
rz(-0.3823897) q[3];
sx q[3];
rz(0.4195655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5220752) q[0];
sx q[0];
rz(-1.8958805) q[0];
sx q[0];
rz(-2.1602901) q[0];
rz(-0.13521067) q[1];
sx q[1];
rz(-1.0842666) q[1];
sx q[1];
rz(-2.9503126) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0606196) q[0];
sx q[0];
rz(-1.4220211) q[0];
sx q[0];
rz(3.0540375) q[0];
rz(-pi) q[1];
rz(-0.83696604) q[2];
sx q[2];
rz(-1.3651197) q[2];
sx q[2];
rz(0.93081805) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8297255) q[1];
sx q[1];
rz(-0.79677454) q[1];
sx q[1];
rz(-1.7142678) q[1];
x q[2];
rz(1.220827) q[3];
sx q[3];
rz(-1.063949) q[3];
sx q[3];
rz(-0.99196539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4613351) q[2];
sx q[2];
rz(-2.1562083) q[2];
sx q[2];
rz(2.130924) q[2];
rz(2.3800395) q[3];
sx q[3];
rz(-1.1798309) q[3];
sx q[3];
rz(-0.23553577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33870944) q[0];
sx q[0];
rz(-0.25512472) q[0];
sx q[0];
rz(-0.55661911) q[0];
rz(3.026475) q[1];
sx q[1];
rz(-1.3373673) q[1];
sx q[1];
rz(2.1690878) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92161979) q[0];
sx q[0];
rz(-0.12244206) q[0];
sx q[0];
rz(-2.5570611) q[0];
rz(-pi) q[1];
rz(0.33072492) q[2];
sx q[2];
rz(-0.84825584) q[2];
sx q[2];
rz(1.5028138) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79341187) q[1];
sx q[1];
rz(-0.23985292) q[1];
sx q[1];
rz(2.1972149) q[1];
x q[2];
rz(-3.1049018) q[3];
sx q[3];
rz(-0.9551691) q[3];
sx q[3];
rz(-0.43945593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3187023) q[2];
sx q[2];
rz(-1.0914785) q[2];
sx q[2];
rz(-1.5931607) q[2];
rz(-1.3657773) q[3];
sx q[3];
rz(-2.8184991) q[3];
sx q[3];
rz(-2.1877066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4218629) q[0];
sx q[0];
rz(-1.8780163) q[0];
sx q[0];
rz(1.7156037) q[0];
rz(1.0643719) q[1];
sx q[1];
rz(-2.1247037) q[1];
sx q[1];
rz(-2.7672966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22898856) q[0];
sx q[0];
rz(-1.2545663) q[0];
sx q[0];
rz(-2.7094748) q[0];
rz(-pi) q[1];
rz(2.6475545) q[2];
sx q[2];
rz(-1.3400153) q[2];
sx q[2];
rz(1.5649232) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7172076) q[1];
sx q[1];
rz(-1.3076412) q[1];
sx q[1];
rz(2.5724263) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4738594) q[3];
sx q[3];
rz(-2.1453834) q[3];
sx q[3];
rz(-2.7001065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2465683) q[2];
sx q[2];
rz(-2.4763069) q[2];
sx q[2];
rz(-0.95823112) q[2];
rz(-2.9124177) q[3];
sx q[3];
rz(-1.6848247) q[3];
sx q[3];
rz(2.5206101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36528698) q[0];
sx q[0];
rz(-1.1927274) q[0];
sx q[0];
rz(-2.2348485) q[0];
rz(-1.0892185) q[1];
sx q[1];
rz(-1.4995432) q[1];
sx q[1];
rz(-1.8315171) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.037127) q[0];
sx q[0];
rz(-1.6220399) q[0];
sx q[0];
rz(0.38802223) q[0];
rz(-pi) q[1];
rz(1.9016978) q[2];
sx q[2];
rz(-0.885193) q[2];
sx q[2];
rz(-2.9943525) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9359365) q[1];
sx q[1];
rz(-1.2381136) q[1];
sx q[1];
rz(1.7722305) q[1];
rz(-pi) q[2];
rz(-0.044192627) q[3];
sx q[3];
rz(-2.5857539) q[3];
sx q[3];
rz(-2.1293872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2157796) q[2];
sx q[2];
rz(-0.44162193) q[2];
sx q[2];
rz(-1.4833935) q[2];
rz(-0.27967927) q[3];
sx q[3];
rz(-2.1488583) q[3];
sx q[3];
rz(-3.083995) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9782372) q[0];
sx q[0];
rz(-1.5196479) q[0];
sx q[0];
rz(0.21959198) q[0];
rz(0.50312463) q[1];
sx q[1];
rz(-2.2527835) q[1];
sx q[1];
rz(0.84987744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726639) q[0];
sx q[0];
rz(-1.7964296) q[0];
sx q[0];
rz(-2.9805141) q[0];
x q[1];
rz(0.093896534) q[2];
sx q[2];
rz(-2.2517423) q[2];
sx q[2];
rz(-2.4927793) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97493193) q[1];
sx q[1];
rz(-0.76172511) q[1];
sx q[1];
rz(-1.5207661) q[1];
rz(-pi) q[2];
rz(1.3109342) q[3];
sx q[3];
rz(-1.6901008) q[3];
sx q[3];
rz(-2.3343149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5510817) q[2];
sx q[2];
rz(-1.3122281) q[2];
sx q[2];
rz(-1.760651) q[2];
rz(-2.3855709) q[3];
sx q[3];
rz(-2.9383926) q[3];
sx q[3];
rz(2.7856564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6324156) q[0];
sx q[0];
rz(-2.248705) q[0];
sx q[0];
rz(0.40503043) q[0];
rz(0.45267725) q[1];
sx q[1];
rz(-0.98222268) q[1];
sx q[1];
rz(1.8639494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15690878) q[0];
sx q[0];
rz(-2.589993) q[0];
sx q[0];
rz(0.44996913) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69127609) q[2];
sx q[2];
rz(-1.9553767) q[2];
sx q[2];
rz(-0.28011766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2904418) q[1];
sx q[1];
rz(-2.6944707) q[1];
sx q[1];
rz(-1.6563583) q[1];
rz(0.38655917) q[3];
sx q[3];
rz(-0.87214008) q[3];
sx q[3];
rz(-2.3545635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8032802) q[2];
sx q[2];
rz(-2.7351604) q[2];
sx q[2];
rz(0.35783106) q[2];
rz(1.7221649) q[3];
sx q[3];
rz(-1.8678886) q[3];
sx q[3];
rz(-2.0675802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2232067) q[0];
sx q[0];
rz(-0.077843852) q[0];
sx q[0];
rz(-0.11225587) q[0];
rz(-0.90011251) q[1];
sx q[1];
rz(-2.0745514) q[1];
sx q[1];
rz(-0.21044883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1505517) q[0];
sx q[0];
rz(-1.3897087) q[0];
sx q[0];
rz(0.51245706) q[0];
rz(-1.1591572) q[2];
sx q[2];
rz(-1.7214081) q[2];
sx q[2];
rz(-0.48356907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3887742) q[1];
sx q[1];
rz(-2.7669199) q[1];
sx q[1];
rz(2.545536) q[1];
rz(-2.3841303) q[3];
sx q[3];
rz(-2.8031581) q[3];
sx q[3];
rz(-0.38871845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9615053) q[2];
sx q[2];
rz(-2.4752361) q[2];
sx q[2];
rz(1.5853184) q[2];
rz(-1.8680343) q[3];
sx q[3];
rz(-0.62265101) q[3];
sx q[3];
rz(-2.9343228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5719941) q[0];
sx q[0];
rz(-0.8710237) q[0];
sx q[0];
rz(-1.3652753) q[0];
rz(-0.81644425) q[1];
sx q[1];
rz(-1.2533617) q[1];
sx q[1];
rz(-0.15773699) q[1];
rz(-2.0408761) q[2];
sx q[2];
rz(-1.7190949) q[2];
sx q[2];
rz(0.20863056) q[2];
rz(-2.3002426) q[3];
sx q[3];
rz(-1.311306) q[3];
sx q[3];
rz(-0.73248274) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
