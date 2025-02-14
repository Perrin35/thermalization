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
rz(-2.2633679) q[0];
sx q[0];
rz(-0.72265923) q[0];
sx q[0];
rz(-1.0455796) q[0];
rz(0.0027520952) q[1];
sx q[1];
rz(-1.4581008) q[1];
sx q[1];
rz(2.4207065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2963639) q[0];
sx q[0];
rz(-0.90334821) q[0];
sx q[0];
rz(-0.11752252) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.051250967) q[2];
sx q[2];
rz(-1.379671) q[2];
sx q[2];
rz(-2.4911936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73511368) q[1];
sx q[1];
rz(-1.8476474) q[1];
sx q[1];
rz(-2.7471354) q[1];
x q[2];
rz(-3.0009427) q[3];
sx q[3];
rz(-2.3936134) q[3];
sx q[3];
rz(-1.5281531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8568933) q[2];
sx q[2];
rz(-1.2130986) q[2];
sx q[2];
rz(-2.6846679) q[2];
rz(2.4871155) q[3];
sx q[3];
rz(-0.5623397) q[3];
sx q[3];
rz(-2.7175609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829795) q[0];
sx q[0];
rz(-2.787866) q[0];
sx q[0];
rz(0.60802996) q[0];
rz(-1.9425302) q[1];
sx q[1];
rz(-0.46747318) q[1];
sx q[1];
rz(0.79493585) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1300995) q[0];
sx q[0];
rz(-1.2813632) q[0];
sx q[0];
rz(2.7910146) q[0];
rz(-0.10484597) q[2];
sx q[2];
rz(-2.8730132) q[2];
sx q[2];
rz(-1.9263445) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5175482) q[1];
sx q[1];
rz(-2.4409082) q[1];
sx q[1];
rz(-2.0113883) q[1];
x q[2];
rz(-2.7057418) q[3];
sx q[3];
rz(-2.4727949) q[3];
sx q[3];
rz(-0.28025249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6752211) q[2];
sx q[2];
rz(-1.0543062) q[2];
sx q[2];
rz(-1.6611453) q[2];
rz(0.36888567) q[3];
sx q[3];
rz(-1.2494913) q[3];
sx q[3];
rz(0.7307542) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22496255) q[0];
sx q[0];
rz(-2.1376305) q[0];
sx q[0];
rz(2.4555901) q[0];
rz(-1.2432159) q[1];
sx q[1];
rz(-2.9153283) q[1];
sx q[1];
rz(-0.6023947) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43322726) q[0];
sx q[0];
rz(-2.4658563) q[0];
sx q[0];
rz(0.85322081) q[0];
rz(-pi) q[1];
rz(-1.5302646) q[2];
sx q[2];
rz(-1.6031331) q[2];
sx q[2];
rz(-3.0258816) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2418306) q[1];
sx q[1];
rz(-1.6679523) q[1];
sx q[1];
rz(-1.001557) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32633324) q[3];
sx q[3];
rz(-1.1236533) q[3];
sx q[3];
rz(2.9016419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38597044) q[2];
sx q[2];
rz(-1.8381511) q[2];
sx q[2];
rz(1.8899567) q[2];
rz(2.5899467) q[3];
sx q[3];
rz(-2.9988204) q[3];
sx q[3];
rz(2.2940206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8400693) q[0];
sx q[0];
rz(-1.9020377) q[0];
sx q[0];
rz(-3.0384592) q[0];
rz(-0.81941191) q[1];
sx q[1];
rz(-0.90140072) q[1];
sx q[1];
rz(0.55319667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5846636) q[0];
sx q[0];
rz(-1.9346666) q[0];
sx q[0];
rz(-2.8332024) q[0];
x q[1];
rz(-2.1167088) q[2];
sx q[2];
rz(-1.8736412) q[2];
sx q[2];
rz(-2.5036223) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31732949) q[1];
sx q[1];
rz(-1.1019442) q[1];
sx q[1];
rz(0.15293185) q[1];
rz(-pi) q[2];
rz(3.0402667) q[3];
sx q[3];
rz(-1.1412217) q[3];
sx q[3];
rz(2.5954036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4210522) q[2];
sx q[2];
rz(-2.4002176) q[2];
sx q[2];
rz(2.6067624) q[2];
rz(0.77830166) q[3];
sx q[3];
rz(-2.4141267) q[3];
sx q[3];
rz(2.9261869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4541723) q[0];
sx q[0];
rz(-1.146831) q[0];
sx q[0];
rz(-2.4046894) q[0];
rz(-0.058628254) q[1];
sx q[1];
rz(-0.77719378) q[1];
sx q[1];
rz(-0.58707213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50518513) q[0];
sx q[0];
rz(-0.93707935) q[0];
sx q[0];
rz(-0.039713602) q[0];
x q[1];
rz(1.306702) q[2];
sx q[2];
rz(-2.9042431) q[2];
sx q[2];
rz(-1.3106948) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.67501589) q[1];
sx q[1];
rz(-1.2685923) q[1];
sx q[1];
rz(-0.034053996) q[1];
x q[2];
rz(1.8892509) q[3];
sx q[3];
rz(-1.8618004) q[3];
sx q[3];
rz(-2.165447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8495142) q[2];
sx q[2];
rz(-0.92486113) q[2];
sx q[2];
rz(-0.64797956) q[2];
rz(2.5255711) q[3];
sx q[3];
rz(-1.5024622) q[3];
sx q[3];
rz(-0.047671635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41094008) q[0];
sx q[0];
rz(-2.3766282) q[0];
sx q[0];
rz(0.95482811) q[0];
rz(3.094589) q[1];
sx q[1];
rz(-0.67263293) q[1];
sx q[1];
rz(2.1254553) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3339301) q[0];
sx q[0];
rz(-2.3950845) q[0];
sx q[0];
rz(2.4115415) q[0];
rz(-0.72332763) q[2];
sx q[2];
rz(-1.9942585) q[2];
sx q[2];
rz(-2.9996769) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.816427) q[1];
sx q[1];
rz(-2.0762968) q[1];
sx q[1];
rz(1.5542404) q[1];
rz(1.627907) q[3];
sx q[3];
rz(-2.1942003) q[3];
sx q[3];
rz(0.36950612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5434692) q[2];
sx q[2];
rz(-1.2657961) q[2];
sx q[2];
rz(3.0908965) q[2];
rz(-3.0905881) q[3];
sx q[3];
rz(-1.8404605) q[3];
sx q[3];
rz(0.72371975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26492587) q[0];
sx q[0];
rz(-2.6267138) q[0];
sx q[0];
rz(1.6616954) q[0];
rz(2.0199846) q[1];
sx q[1];
rz(-1.2622967) q[1];
sx q[1];
rz(2.6720537) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8136217) q[0];
sx q[0];
rz(-1.4634973) q[0];
sx q[0];
rz(3.0946428) q[0];
x q[1];
rz(-1.1781349) q[2];
sx q[2];
rz(-1.4265276) q[2];
sx q[2];
rz(-0.96111125) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96315876) q[1];
sx q[1];
rz(-1.3972153) q[1];
sx q[1];
rz(1.0976237) q[1];
x q[2];
rz(-2.363405) q[3];
sx q[3];
rz(-2.0594849) q[3];
sx q[3];
rz(0.60569872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40662128) q[2];
sx q[2];
rz(-1.5903641) q[2];
sx q[2];
rz(1.1173908) q[2];
rz(2.0495074) q[3];
sx q[3];
rz(-1.4039682) q[3];
sx q[3];
rz(2.7058069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0206873) q[0];
sx q[0];
rz(-2.0285323) q[0];
sx q[0];
rz(-0.06509617) q[0];
rz(-0.33196017) q[1];
sx q[1];
rz(-1.8354225) q[1];
sx q[1];
rz(-1.2750221) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.770261) q[0];
sx q[0];
rz(-1.341371) q[0];
sx q[0];
rz(0.0097994402) q[0];
rz(-pi) q[1];
rz(-0.1365928) q[2];
sx q[2];
rz(-1.1316191) q[2];
sx q[2];
rz(0.76132773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.304833) q[1];
sx q[1];
rz(-2.4570434) q[1];
sx q[1];
rz(0.080938653) q[1];
x q[2];
rz(1.3887819) q[3];
sx q[3];
rz(-2.3654147) q[3];
sx q[3];
rz(2.9266199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.28978213) q[2];
sx q[2];
rz(-1.6921356) q[2];
sx q[2];
rz(-0.97274485) q[2];
rz(2.4408477) q[3];
sx q[3];
rz(-2.2795129) q[3];
sx q[3];
rz(0.8249445) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7026611) q[0];
sx q[0];
rz(-2.7362566) q[0];
sx q[0];
rz(-0.40913707) q[0];
rz(1.9955955) q[1];
sx q[1];
rz(-1.0640249) q[1];
sx q[1];
rz(1.8786059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1407023) q[0];
sx q[0];
rz(-1.5789728) q[0];
sx q[0];
rz(-0.00012881669) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97102286) q[2];
sx q[2];
rz(-0.46659887) q[2];
sx q[2];
rz(0.80228892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4597171) q[1];
sx q[1];
rz(-0.72725216) q[1];
sx q[1];
rz(1.2254814) q[1];
rz(-2.5714794) q[3];
sx q[3];
rz(-1.4973204) q[3];
sx q[3];
rz(0.076345293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24141773) q[2];
sx q[2];
rz(-1.2215542) q[2];
sx q[2];
rz(-2.9202666) q[2];
rz(0.48702249) q[3];
sx q[3];
rz(-0.86921391) q[3];
sx q[3];
rz(1.2041436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0601592) q[0];
sx q[0];
rz(-2.8964323) q[0];
sx q[0];
rz(-1.918445) q[0];
rz(3.0606015) q[1];
sx q[1];
rz(-1.6209737) q[1];
sx q[1];
rz(1.5222668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.318383) q[0];
sx q[0];
rz(-1.2391187) q[0];
sx q[0];
rz(-1.7451502) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4143405) q[2];
sx q[2];
rz(-1.1071812) q[2];
sx q[2];
rz(1.9510683) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63357089) q[1];
sx q[1];
rz(-2.7765818) q[1];
sx q[1];
rz(0.99442039) q[1];
x q[2];
rz(-2.6468148) q[3];
sx q[3];
rz(-1.7253381) q[3];
sx q[3];
rz(-3.0734594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0167375) q[2];
sx q[2];
rz(-1.3913466) q[2];
sx q[2];
rz(-0.69046956) q[2];
rz(0.27103439) q[3];
sx q[3];
rz(-2.6764968) q[3];
sx q[3];
rz(-1.5976228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6799714) q[0];
sx q[0];
rz(-2.2053056) q[0];
sx q[0];
rz(2.9619138) q[0];
rz(-0.24196729) q[1];
sx q[1];
rz(-0.89121834) q[1];
sx q[1];
rz(-1.2527087) q[1];
rz(2.7015637) q[2];
sx q[2];
rz(-0.46092214) q[2];
sx q[2];
rz(3.0546247) q[2];
rz(2.8293256) q[3];
sx q[3];
rz(-0.91296416) q[3];
sx q[3];
rz(2.8784646) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
