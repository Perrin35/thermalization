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
rz(-0.93936062) q[0];
sx q[0];
rz(4.2606104) q[0];
sx q[0];
rz(9.4584447) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(1.508498) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6389119) q[0];
sx q[0];
rz(-1.685132) q[0];
sx q[0];
rz(-3.1114716) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.075307651) q[2];
sx q[2];
rz(-1.2774373) q[2];
sx q[2];
rz(1.2874029) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.37225393) q[1];
sx q[1];
rz(-2.4318691) q[1];
sx q[1];
rz(-0.22270568) q[1];
rz(-pi) q[2];
rz(-2.5023828) q[3];
sx q[3];
rz(-0.91851252) q[3];
sx q[3];
rz(-0.88310421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17026751) q[2];
sx q[2];
rz(-2.2496932) q[2];
sx q[2];
rz(-2.3049221) q[2];
rz(-2.3352046) q[3];
sx q[3];
rz(-0.92525768) q[3];
sx q[3];
rz(-2.9562318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065501325) q[0];
sx q[0];
rz(-2.0747023) q[0];
sx q[0];
rz(1.2170894) q[0];
rz(-2.941046) q[1];
sx q[1];
rz(-2.3417818) q[1];
sx q[1];
rz(-2.5977871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4539219) q[0];
sx q[0];
rz(-2.7467496) q[0];
sx q[0];
rz(-2.5107491) q[0];
x q[1];
rz(-0.5949623) q[2];
sx q[2];
rz(-1.5820855) q[2];
sx q[2];
rz(1.6807229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.75419352) q[1];
sx q[1];
rz(-2.6541944) q[1];
sx q[1];
rz(-2.2226366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89558286) q[3];
sx q[3];
rz(-1.5434268) q[3];
sx q[3];
rz(0.41639183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4726938) q[2];
sx q[2];
rz(-0.74287477) q[2];
sx q[2];
rz(1.6061858) q[2];
rz(1.6425447) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(-1.8892052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1010308) q[0];
sx q[0];
rz(-2.5767548) q[0];
sx q[0];
rz(-3.132013) q[0];
rz(-1.0293695) q[1];
sx q[1];
rz(-1.5734438) q[1];
sx q[1];
rz(-1.8052489) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21214616) q[0];
sx q[0];
rz(-0.98185724) q[0];
sx q[0];
rz(-0.94561823) q[0];
rz(-pi) q[1];
rz(0.32149489) q[2];
sx q[2];
rz(-2.0509556) q[2];
sx q[2];
rz(-2.8804422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.0071466999) q[1];
sx q[1];
rz(-1.916051) q[1];
sx q[1];
rz(2.9551639) q[1];
rz(-pi) q[2];
rz(-1.9668749) q[3];
sx q[3];
rz(-1.7273335) q[3];
sx q[3];
rz(-0.96352947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0336527) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(1.3134202) q[2];
rz(-3.1016453) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(-1.1279172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(3.0082598) q[0];
sx q[0];
rz(-1.8513716) q[0];
sx q[0];
rz(0.3966575) q[0];
rz(-2.2465514) q[1];
sx q[1];
rz(-1.4245234) q[1];
sx q[1];
rz(-0.76362124) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39786441) q[0];
sx q[0];
rz(-1.8157122) q[0];
sx q[0];
rz(0.082964925) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12053804) q[2];
sx q[2];
rz(-1.5539031) q[2];
sx q[2];
rz(-2.8617045) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3862003) q[1];
sx q[1];
rz(-1.5441276) q[1];
sx q[1];
rz(-1.3332572) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5337297) q[3];
sx q[3];
rz(-1.8614113) q[3];
sx q[3];
rz(-0.99729482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.414173) q[2];
sx q[2];
rz(-1.6215934) q[2];
sx q[2];
rz(-2.011389) q[2];
rz(1.1388418) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(1.9738919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30911699) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(-2.6600237) q[0];
rz(-2.8140977) q[1];
sx q[1];
rz(-1.7227252) q[1];
sx q[1];
rz(0.96424261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010134546) q[0];
sx q[0];
rz(-1.9572963) q[0];
sx q[0];
rz(-0.23614998) q[0];
x q[1];
rz(0.19179253) q[2];
sx q[2];
rz(-2.1438877) q[2];
sx q[2];
rz(-1.7157784) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3754978) q[1];
sx q[1];
rz(-1.8154486) q[1];
sx q[1];
rz(0.12568139) q[1];
x q[2];
rz(-2.6752969) q[3];
sx q[3];
rz(-2.7096636) q[3];
sx q[3];
rz(-1.4817878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97011906) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(-1.4946651) q[2];
rz(2.1286879) q[3];
sx q[3];
rz(-2.2154112) q[3];
sx q[3];
rz(0.51500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9377015) q[0];
sx q[0];
rz(-2.1370115) q[0];
sx q[0];
rz(1.0585693) q[0];
rz(2.2189498) q[1];
sx q[1];
rz(-1.0170931) q[1];
sx q[1];
rz(-3.0386818) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9227717) q[0];
sx q[0];
rz(-2.0334683) q[0];
sx q[0];
rz(-0.86313049) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40252588) q[2];
sx q[2];
rz(-2.5136097) q[2];
sx q[2];
rz(-1.8875853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6158689) q[1];
sx q[1];
rz(-1.3868887) q[1];
sx q[1];
rz(-2.1083819) q[1];
x q[2];
rz(1.5779805) q[3];
sx q[3];
rz(-2.4251591) q[3];
sx q[3];
rz(-0.35182692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70204488) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(-0.31583819) q[2];
rz(2.1465007) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(0.72527138) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.134326) q[0];
sx q[0];
rz(-0.26679978) q[0];
sx q[0];
rz(2.5079978) q[0];
rz(2.7569356) q[1];
sx q[1];
rz(-1.9622012) q[1];
sx q[1];
rz(0.42919174) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8096247) q[0];
sx q[0];
rz(-1.8242697) q[0];
sx q[0];
rz(1.0868207) q[0];
x q[1];
rz(-1.6056988) q[2];
sx q[2];
rz(-0.88551025) q[2];
sx q[2];
rz(-1.8165782) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1146443) q[1];
sx q[1];
rz(-2.1430795) q[1];
sx q[1];
rz(1.9932822) q[1];
rz(-pi) q[2];
rz(0.46439345) q[3];
sx q[3];
rz(-1.0611609) q[3];
sx q[3];
rz(1.2023752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0523494) q[2];
sx q[2];
rz(-0.77755916) q[2];
sx q[2];
rz(-0.26697978) q[2];
rz(3.0698245) q[3];
sx q[3];
rz(-2.1238582) q[3];
sx q[3];
rz(1.7223541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4778022) q[0];
sx q[0];
rz(-3.0402122) q[0];
sx q[0];
rz(1.7410393) q[0];
rz(1.1446704) q[1];
sx q[1];
rz(-2.0796227) q[1];
sx q[1];
rz(-2.5544419) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42204866) q[0];
sx q[0];
rz(-1.1503771) q[0];
sx q[0];
rz(-2.4499513) q[0];
rz(-1.6875504) q[2];
sx q[2];
rz(-1.1551305) q[2];
sx q[2];
rz(-2.1374201) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.352601) q[1];
sx q[1];
rz(-1.7956211) q[1];
sx q[1];
rz(2.4487752) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8244157) q[3];
sx q[3];
rz(-0.97517386) q[3];
sx q[3];
rz(-0.45444876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8553541) q[2];
sx q[2];
rz(-1.1796395) q[2];
sx q[2];
rz(-0.41909763) q[2];
rz(1.5239483) q[3];
sx q[3];
rz(-0.91846275) q[3];
sx q[3];
rz(-1.6379448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54290259) q[0];
sx q[0];
rz(-1.9200696) q[0];
sx q[0];
rz(2.8322117) q[0];
rz(1.3503374) q[1];
sx q[1];
rz(-0.58809909) q[1];
sx q[1];
rz(-0.5538756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5039936) q[0];
sx q[0];
rz(-1.1553815) q[0];
sx q[0];
rz(1.2211772) q[0];
rz(-pi) q[1];
rz(1.1872841) q[2];
sx q[2];
rz(-1.4330136) q[2];
sx q[2];
rz(-2.1509508) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3024721) q[1];
sx q[1];
rz(-0.97187798) q[1];
sx q[1];
rz(-0.37521404) q[1];
rz(-pi) q[2];
rz(1.8819412) q[3];
sx q[3];
rz(-0.89289819) q[3];
sx q[3];
rz(-2.4571612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.89173335) q[2];
sx q[2];
rz(-1.5313818) q[2];
sx q[2];
rz(1.7943133) q[2];
rz(-2.3740718) q[3];
sx q[3];
rz(-1.9556421) q[3];
sx q[3];
rz(2.2931113) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7193741) q[0];
sx q[0];
rz(-2.8415871) q[0];
sx q[0];
rz(-1.2147709) q[0];
rz(1.0019373) q[1];
sx q[1];
rz(-1.6286214) q[1];
sx q[1];
rz(-0.92438662) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533646) q[0];
sx q[0];
rz(-1.1832675) q[0];
sx q[0];
rz(2.1977682) q[0];
rz(-pi) q[1];
rz(-3.1141823) q[2];
sx q[2];
rz(-2.1152585) q[2];
sx q[2];
rz(1.4278864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50528697) q[1];
sx q[1];
rz(-1.5940548) q[1];
sx q[1];
rz(2.9182698) q[1];
x q[2];
rz(0.89556344) q[3];
sx q[3];
rz(-0.82230869) q[3];
sx q[3];
rz(-3.0117717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7944916) q[2];
sx q[2];
rz(-2.2091986) q[2];
sx q[2];
rz(0.36716983) q[2];
rz(-1.9764887) q[3];
sx q[3];
rz(-2.1737183) q[3];
sx q[3];
rz(0.86618209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1871056) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(2.8754996) q[1];
sx q[1];
rz(-2.2292744) q[1];
sx q[1];
rz(1.6899065) q[1];
rz(-2.8546147) q[2];
sx q[2];
rz(-1.8331883) q[2];
sx q[2];
rz(2.2856648) q[2];
rz(2.9333276) q[3];
sx q[3];
rz(-2.1622445) q[3];
sx q[3];
rz(1.6643659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
