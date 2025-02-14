OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5724343) q[0];
sx q[0];
rz(7.2814514) q[0];
sx q[0];
rz(9.8132039) q[0];
rz(-1.2216964) q[1];
sx q[1];
rz(-2.1887527) q[1];
sx q[1];
rz(0.0739007) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6361481) q[0];
sx q[0];
rz(-1.0514604) q[0];
sx q[0];
rz(1.8072007) q[0];
rz(-pi) q[1];
rz(2.2227395) q[2];
sx q[2];
rz(-1.5896738) q[2];
sx q[2];
rz(2.9110661) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8732637) q[1];
sx q[1];
rz(-1.2115098) q[1];
sx q[1];
rz(0.84228911) q[1];
x q[2];
rz(-1.8201572) q[3];
sx q[3];
rz(-1.7791181) q[3];
sx q[3];
rz(-1.1492532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4581603) q[2];
sx q[2];
rz(-1.6998123) q[2];
sx q[2];
rz(-1.3275576) q[2];
rz(-2.1761927) q[3];
sx q[3];
rz(-2.5520971) q[3];
sx q[3];
rz(-1.0546257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7205768) q[0];
sx q[0];
rz(-3.133931) q[0];
sx q[0];
rz(-1.3910008) q[0];
rz(1.1945456) q[1];
sx q[1];
rz(-0.60940131) q[1];
sx q[1];
rz(2.4360099) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53853453) q[0];
sx q[0];
rz(-2.0655736) q[0];
sx q[0];
rz(2.323708) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1384355) q[2];
sx q[2];
rz(-2.2213644) q[2];
sx q[2];
rz(-2.7108266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0471474) q[1];
sx q[1];
rz(-1.021968) q[1];
sx q[1];
rz(2.7286163) q[1];
x q[2];
rz(-2.5417305) q[3];
sx q[3];
rz(-1.9981435) q[3];
sx q[3];
rz(-2.7823256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.07448639) q[2];
sx q[2];
rz(-1.4718082) q[2];
sx q[2];
rz(-0.86522317) q[2];
rz(-1.5857961) q[3];
sx q[3];
rz(-0.14388789) q[3];
sx q[3];
rz(1.5998862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14875749) q[0];
sx q[0];
rz(-0.8388297) q[0];
sx q[0];
rz(1.6194153) q[0];
rz(-1.7332227) q[1];
sx q[1];
rz(-2.759628) q[1];
sx q[1];
rz(2.9011889) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0494321) q[0];
sx q[0];
rz(-0.76917523) q[0];
sx q[0];
rz(0.34017684) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9683686) q[2];
sx q[2];
rz(-1.3971739) q[2];
sx q[2];
rz(-2.589203) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8286491) q[1];
sx q[1];
rz(-3.0285081) q[1];
sx q[1];
rz(2.7663016) q[1];
x q[2];
rz(-2.6601276) q[3];
sx q[3];
rz(-1.4507308) q[3];
sx q[3];
rz(-1.039618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1323041) q[2];
sx q[2];
rz(-1.4762286) q[2];
sx q[2];
rz(-2.0557527) q[2];
rz(-3.0911607) q[3];
sx q[3];
rz(-1.4112873) q[3];
sx q[3];
rz(2.1850695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0352935) q[0];
sx q[0];
rz(-1.6989919) q[0];
sx q[0];
rz(0.77044368) q[0];
rz(-0.92535198) q[1];
sx q[1];
rz(-1.4848361) q[1];
sx q[1];
rz(2.3063708) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2352426) q[0];
sx q[0];
rz(-2.6444204) q[0];
sx q[0];
rz(3.1072561) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24150924) q[2];
sx q[2];
rz(-1.4837974) q[2];
sx q[2];
rz(2.3149025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8586352) q[1];
sx q[1];
rz(-1.5834619) q[1];
sx q[1];
rz(2.0974166) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51839197) q[3];
sx q[3];
rz(-0.22991069) q[3];
sx q[3];
rz(-2.0748024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21021065) q[2];
sx q[2];
rz(-1.8385734) q[2];
sx q[2];
rz(0.72969985) q[2];
rz(-2.1424255) q[3];
sx q[3];
rz(-1.3068643) q[3];
sx q[3];
rz(-2.6877747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5900742) q[0];
sx q[0];
rz(-0.22576627) q[0];
sx q[0];
rz(0.99910587) q[0];
rz(-1.0480115) q[1];
sx q[1];
rz(-0.71185714) q[1];
sx q[1];
rz(-2.5968831) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7479046) q[0];
sx q[0];
rz(-1.3917408) q[0];
sx q[0];
rz(0.25435971) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72417132) q[2];
sx q[2];
rz(-2.050638) q[2];
sx q[2];
rz(1.2042696) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5156158) q[1];
sx q[1];
rz(-2.5624609) q[1];
sx q[1];
rz(-2.7939151) q[1];
rz(-pi) q[2];
rz(-1.5673166) q[3];
sx q[3];
rz(-1.9762357) q[3];
sx q[3];
rz(-2.5377948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5576632) q[2];
sx q[2];
rz(-1.1489392) q[2];
sx q[2];
rz(-2.124713) q[2];
rz(-0.66295463) q[3];
sx q[3];
rz(-2.4487977) q[3];
sx q[3];
rz(1.8069161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3589631) q[0];
sx q[0];
rz(-0.49802676) q[0];
sx q[0];
rz(1.1516512) q[0];
rz(2.2478814) q[1];
sx q[1];
rz(-1.7620554) q[1];
sx q[1];
rz(3.025257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3860364) q[0];
sx q[0];
rz(-2.0445637) q[0];
sx q[0];
rz(-0.12065819) q[0];
rz(-pi) q[1];
rz(-1.9151808) q[2];
sx q[2];
rz(-1.55416) q[2];
sx q[2];
rz(-3.0964451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5999122) q[1];
sx q[1];
rz(-2.5697418) q[1];
sx q[1];
rz(1.2979371) q[1];
rz(-3.0553566) q[3];
sx q[3];
rz(-0.75018084) q[3];
sx q[3];
rz(2.3762109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.56670386) q[2];
sx q[2];
rz(-2.7723007) q[2];
sx q[2];
rz(2.134038) q[2];
rz(1.1653853) q[3];
sx q[3];
rz(-2.8663965) q[3];
sx q[3];
rz(0.63024855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4646869) q[0];
sx q[0];
rz(-2.6755264) q[0];
sx q[0];
rz(0.16978547) q[0];
rz(-0.67974293) q[1];
sx q[1];
rz(-0.83234537) q[1];
sx q[1];
rz(-2.2198417) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1257432) q[0];
sx q[0];
rz(-1.4471869) q[0];
sx q[0];
rz(0.25861926) q[0];
rz(-pi) q[1];
rz(-0.1711357) q[2];
sx q[2];
rz(-0.31120069) q[2];
sx q[2];
rz(-2.5301274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55619986) q[1];
sx q[1];
rz(-1.0784282) q[1];
sx q[1];
rz(2.7436168) q[1];
rz(-pi) q[2];
rz(-2.5329221) q[3];
sx q[3];
rz(-2.8117315) q[3];
sx q[3];
rz(0.69486991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0035231) q[2];
sx q[2];
rz(-1.4813083) q[2];
sx q[2];
rz(1.8180234) q[2];
rz(-2.0924163) q[3];
sx q[3];
rz(-1.1314355) q[3];
sx q[3];
rz(-1.3808892) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84999371) q[0];
sx q[0];
rz(-1.1350564) q[0];
sx q[0];
rz(2.3681613) q[0];
rz(0.4666346) q[1];
sx q[1];
rz(-1.0728873) q[1];
sx q[1];
rz(3.1007865) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0590574) q[0];
sx q[0];
rz(-1.5837386) q[0];
sx q[0];
rz(-1.4380161) q[0];
rz(-0.86470072) q[2];
sx q[2];
rz(-1.0664202) q[2];
sx q[2];
rz(2.9838533) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0405214) q[1];
sx q[1];
rz(-0.98993976) q[1];
sx q[1];
rz(-3.0809666) q[1];
rz(-1.0061408) q[3];
sx q[3];
rz(-1.4251544) q[3];
sx q[3];
rz(-0.70051769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9890954) q[2];
sx q[2];
rz(-0.54741198) q[2];
sx q[2];
rz(-3.1067749) q[2];
rz(0.028248938) q[3];
sx q[3];
rz(-1.0476799) q[3];
sx q[3];
rz(-0.69407535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5038274) q[0];
sx q[0];
rz(-0.70264188) q[0];
sx q[0];
rz(1.0850061) q[0];
rz(1.3582235) q[1];
sx q[1];
rz(-2.5085776) q[1];
sx q[1];
rz(2.0620652) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87303783) q[0];
sx q[0];
rz(-1.7760906) q[0];
sx q[0];
rz(-2.1389066) q[0];
x q[1];
rz(-0.44869156) q[2];
sx q[2];
rz(-1.4451705) q[2];
sx q[2];
rz(-0.10754935) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4145045) q[1];
sx q[1];
rz(-1.8782756) q[1];
sx q[1];
rz(0.25924637) q[1];
rz(-pi) q[2];
rz(-0.50864403) q[3];
sx q[3];
rz(-2.2318342) q[3];
sx q[3];
rz(-0.69045541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2435771) q[2];
sx q[2];
rz(-1.54553) q[2];
sx q[2];
rz(0.071694516) q[2];
rz(2.535635) q[3];
sx q[3];
rz(-0.76064435) q[3];
sx q[3];
rz(-3.0683556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83409413) q[0];
sx q[0];
rz(-1.6313169) q[0];
sx q[0];
rz(0.50258762) q[0];
rz(-1.7723627) q[1];
sx q[1];
rz(-1.9160198) q[1];
sx q[1];
rz(0.1756846) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7384199) q[0];
sx q[0];
rz(-1.4851928) q[0];
sx q[0];
rz(-1.9001207) q[0];
rz(-pi) q[1];
rz(2.2403342) q[2];
sx q[2];
rz(-2.2961535) q[2];
sx q[2];
rz(2.8536316) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8524234) q[1];
sx q[1];
rz(-2.2383177) q[1];
sx q[1];
rz(2.4856485) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0629582) q[3];
sx q[3];
rz(-0.70379095) q[3];
sx q[3];
rz(-0.87913556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6548369) q[2];
sx q[2];
rz(-1.882587) q[2];
sx q[2];
rz(2.7969825) q[2];
rz(-1.2279855) q[3];
sx q[3];
rz(-0.69209185) q[3];
sx q[3];
rz(2.1740348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6268613) q[0];
sx q[0];
rz(-1.8296965) q[0];
sx q[0];
rz(-2.1924023) q[0];
rz(1.7687891) q[1];
sx q[1];
rz(-0.95284843) q[1];
sx q[1];
rz(1.3292809) q[1];
rz(-2.8368159) q[2];
sx q[2];
rz(-0.73619107) q[2];
sx q[2];
rz(0.29965055) q[2];
rz(-0.19560858) q[3];
sx q[3];
rz(-0.8630639) q[3];
sx q[3];
rz(0.56474781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
