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
rz(-2.1433266) q[0];
sx q[0];
rz(-0.38842595) q[0];
rz(-1.2216964) q[1];
sx q[1];
rz(-2.1887527) q[1];
sx q[1];
rz(-3.067692) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.184351) q[0];
sx q[0];
rz(-1.3660407) q[0];
sx q[0];
rz(0.53157579) q[0];
x q[1];
rz(-0.91885318) q[2];
sx q[2];
rz(-1.5896738) q[2];
sx q[2];
rz(2.9110661) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2683289) q[1];
sx q[1];
rz(-1.2115098) q[1];
sx q[1];
rz(-0.84228911) q[1];
rz(-pi) q[2];
rz(-2.9268215) q[3];
sx q[3];
rz(-1.3269375) q[3];
sx q[3];
rz(-0.36892316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4581603) q[2];
sx q[2];
rz(-1.4417803) q[2];
sx q[2];
rz(-1.8140351) q[2];
rz(-0.96539998) q[3];
sx q[3];
rz(-2.5520971) q[3];
sx q[3];
rz(1.0546257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7205768) q[0];
sx q[0];
rz(-0.0076616658) q[0];
sx q[0];
rz(-1.3910008) q[0];
rz(1.1945456) q[1];
sx q[1];
rz(-0.60940131) q[1];
sx q[1];
rz(2.4360099) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6030581) q[0];
sx q[0];
rz(-2.0655736) q[0];
sx q[0];
rz(2.323708) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4439677) q[2];
sx q[2];
rz(-1.910733) q[2];
sx q[2];
rz(1.4125669) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7930562) q[1];
sx q[1];
rz(-0.67386857) q[1];
sx q[1];
rz(-2.1516031) q[1];
x q[2];
rz(0.59986214) q[3];
sx q[3];
rz(-1.1434492) q[3];
sx q[3];
rz(2.7823256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0671063) q[2];
sx q[2];
rz(-1.4718082) q[2];
sx q[2];
rz(-0.86522317) q[2];
rz(1.5857961) q[3];
sx q[3];
rz(-0.14388789) q[3];
sx q[3];
rz(-1.5998862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14875749) q[0];
sx q[0];
rz(-0.8388297) q[0];
sx q[0];
rz(-1.5221773) q[0];
rz(-1.4083699) q[1];
sx q[1];
rz(-2.759628) q[1];
sx q[1];
rz(-2.9011889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6344389) q[0];
sx q[0];
rz(-0.85570645) q[0];
sx q[0];
rz(-1.2583744) q[0];
rz(0.18797925) q[2];
sx q[2];
rz(-1.1795292) q[2];
sx q[2];
rz(2.1956034) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2061305) q[1];
sx q[1];
rz(-1.67598) q[1];
sx q[1];
rz(-1.5291924) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6601276) q[3];
sx q[3];
rz(-1.6908619) q[3];
sx q[3];
rz(1.039618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1323041) q[2];
sx q[2];
rz(-1.4762286) q[2];
sx q[2];
rz(-1.08584) q[2];
rz(-0.050431937) q[3];
sx q[3];
rz(-1.7303053) q[3];
sx q[3];
rz(2.1850695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10629912) q[0];
sx q[0];
rz(-1.4426008) q[0];
sx q[0];
rz(0.77044368) q[0];
rz(0.92535198) q[1];
sx q[1];
rz(-1.4848361) q[1];
sx q[1];
rz(0.83522183) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94541162) q[0];
sx q[0];
rz(-2.0676488) q[0];
sx q[0];
rz(-1.5894228) q[0];
x q[1];
rz(1.6603819) q[2];
sx q[2];
rz(-1.3302186) q[2];
sx q[2];
rz(2.4188855) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2829574) q[1];
sx q[1];
rz(-1.5581308) q[1];
sx q[1];
rz(1.0441761) q[1];
rz(-pi) q[2];
rz(0.20056574) q[3];
sx q[3];
rz(-1.6839538) q[3];
sx q[3];
rz(0.0029879163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21021065) q[2];
sx q[2];
rz(-1.8385734) q[2];
sx q[2];
rz(-2.4118928) q[2];
rz(0.99916712) q[3];
sx q[3];
rz(-1.3068643) q[3];
sx q[3];
rz(-2.6877747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5900742) q[0];
sx q[0];
rz(-0.22576627) q[0];
sx q[0];
rz(2.1424868) q[0];
rz(2.0935811) q[1];
sx q[1];
rz(-0.71185714) q[1];
sx q[1];
rz(-2.5968831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42349619) q[0];
sx q[0];
rz(-2.831651) q[0];
sx q[0];
rz(-2.518032) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72417132) q[2];
sx q[2];
rz(-2.050638) q[2];
sx q[2];
rz(1.2042696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5156158) q[1];
sx q[1];
rz(-2.5624609) q[1];
sx q[1];
rz(0.34767751) q[1];
x q[2];
rz(0.0081069907) q[3];
sx q[3];
rz(-0.40545344) q[3];
sx q[3];
rz(0.61262006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.58392945) q[2];
sx q[2];
rz(-1.9926535) q[2];
sx q[2];
rz(-2.124713) q[2];
rz(2.478638) q[3];
sx q[3];
rz(-0.69279492) q[3];
sx q[3];
rz(1.3346765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7826295) q[0];
sx q[0];
rz(-0.49802676) q[0];
sx q[0];
rz(1.1516512) q[0];
rz(-2.2478814) q[1];
sx q[1];
rz(-1.7620554) q[1];
sx q[1];
rz(0.11633565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3860364) q[0];
sx q[0];
rz(-1.0970289) q[0];
sx q[0];
rz(3.0209345) q[0];
x q[1];
rz(-0.017673894) q[2];
sx q[2];
rz(-1.9151312) q[2];
sx q[2];
rz(1.6219106) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73977208) q[1];
sx q[1];
rz(-1.7171613) q[1];
sx q[1];
rz(-2.1256281) q[1];
x q[2];
rz(0.74832423) q[3];
sx q[3];
rz(-1.6295506) q[3];
sx q[3];
rz(2.3993381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5748888) q[2];
sx q[2];
rz(-2.7723007) q[2];
sx q[2];
rz(-1.0075547) q[2];
rz(-1.1653853) q[3];
sx q[3];
rz(-0.27519614) q[3];
sx q[3];
rz(0.63024855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6769058) q[0];
sx q[0];
rz(-2.6755264) q[0];
sx q[0];
rz(-2.9718072) q[0];
rz(-2.4618497) q[1];
sx q[1];
rz(-2.3092473) q[1];
sx q[1];
rz(-2.2198417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0089192275) q[0];
sx q[0];
rz(-0.28604315) q[0];
sx q[0];
rz(-2.689365) q[0];
rz(-pi) q[1];
rz(-1.6255195) q[2];
sx q[2];
rz(-1.2642908) q[2];
sx q[2];
rz(2.3505369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55619986) q[1];
sx q[1];
rz(-2.0631644) q[1];
sx q[1];
rz(-0.3979759) q[1];
rz(-pi) q[2];
rz(2.8677651) q[3];
sx q[3];
rz(-1.7570772) q[3];
sx q[3];
rz(-2.8485988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0035231) q[2];
sx q[2];
rz(-1.4813083) q[2];
sx q[2];
rz(1.8180234) q[2];
rz(-1.0491764) q[3];
sx q[3];
rz(-2.0101571) q[3];
sx q[3];
rz(-1.3808892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2915989) q[0];
sx q[0];
rz(-1.1350564) q[0];
sx q[0];
rz(-0.77343136) q[0];
rz(0.4666346) q[1];
sx q[1];
rz(-1.0728873) q[1];
sx q[1];
rz(-0.040806142) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916676) q[0];
sx q[0];
rz(-0.13340575) q[0];
sx q[0];
rz(1.4733423) q[0];
x q[1];
rz(-0.86470072) q[2];
sx q[2];
rz(-2.0751725) q[2];
sx q[2];
rz(0.15773931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0405214) q[1];
sx q[1];
rz(-0.98993976) q[1];
sx q[1];
rz(0.060626027) q[1];
x q[2];
rz(-1.3032622) q[3];
sx q[3];
rz(-2.5604381) q[3];
sx q[3];
rz(-1.095497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1524973) q[2];
sx q[2];
rz(-2.5941807) q[2];
sx q[2];
rz(-3.1067749) q[2];
rz(-0.028248938) q[3];
sx q[3];
rz(-1.0476799) q[3];
sx q[3];
rz(0.69407535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6377653) q[0];
sx q[0];
rz(-2.4389508) q[0];
sx q[0];
rz(1.0850061) q[0];
rz(1.7833692) q[1];
sx q[1];
rz(-2.5085776) q[1];
sx q[1];
rz(1.0795275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5732291) q[0];
sx q[0];
rz(-2.1255593) q[0];
sx q[0];
rz(2.8994096) q[0];
rz(-1.7100542) q[2];
sx q[2];
rz(-1.1258954) q[2];
sx q[2];
rz(-1.5234966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.92380556) q[1];
sx q[1];
rz(-1.3239685) q[1];
sx q[1];
rz(-1.8882165) q[1];
rz(-pi) q[2];
rz(-0.50864403) q[3];
sx q[3];
rz(-2.2318342) q[3];
sx q[3];
rz(-0.69045541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2435771) q[2];
sx q[2];
rz(-1.5960627) q[2];
sx q[2];
rz(0.071694516) q[2];
rz(-0.60595766) q[3];
sx q[3];
rz(-2.3809483) q[3];
sx q[3];
rz(-0.073237091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83409413) q[0];
sx q[0];
rz(-1.6313169) q[0];
sx q[0];
rz(-2.639005) q[0];
rz(-1.7723627) q[1];
sx q[1];
rz(-1.2255729) q[1];
sx q[1];
rz(2.9659081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7384199) q[0];
sx q[0];
rz(-1.4851928) q[0];
sx q[0];
rz(-1.241472) q[0];
rz(-2.530859) q[2];
sx q[2];
rz(-0.94379506) q[2];
sx q[2];
rz(0.58525733) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60424544) q[1];
sx q[1];
rz(-2.2425695) q[1];
sx q[1];
rz(0.91225454) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93264742) q[3];
sx q[3];
rz(-1.89092) q[3];
sx q[3];
rz(0.29039106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.48675576) q[2];
sx q[2];
rz(-1.2590057) q[2];
sx q[2];
rz(0.34461018) q[2];
rz(1.9136072) q[3];
sx q[3];
rz(-2.4495008) q[3];
sx q[3];
rz(0.96755782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5147314) q[0];
sx q[0];
rz(-1.8296965) q[0];
sx q[0];
rz(-2.1924023) q[0];
rz(-1.3728036) q[1];
sx q[1];
rz(-0.95284843) q[1];
sx q[1];
rz(1.3292809) q[1];
rz(-0.30477672) q[2];
sx q[2];
rz(-2.4054016) q[2];
sx q[2];
rz(-2.8419421) q[2];
rz(-2.9459841) q[3];
sx q[3];
rz(-2.2785288) q[3];
sx q[3];
rz(-2.5768448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
