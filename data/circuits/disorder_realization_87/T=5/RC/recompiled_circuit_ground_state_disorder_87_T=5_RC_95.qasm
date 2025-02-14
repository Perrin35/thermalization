OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73622048) q[0];
sx q[0];
rz(-3.1206107) q[0];
sx q[0];
rz(-1.9606645) q[0];
rz(-0.48038545) q[1];
sx q[1];
rz(-1.9863702) q[1];
sx q[1];
rz(-0.46673271) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0281848) q[0];
sx q[0];
rz(-0.81499463) q[0];
sx q[0];
rz(-2.7477376) q[0];
x q[1];
rz(-1.0814904) q[2];
sx q[2];
rz(-1.0439516) q[2];
sx q[2];
rz(2.26109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84604656) q[1];
sx q[1];
rz(-1.8091077) q[1];
sx q[1];
rz(-3.1370276) q[1];
rz(-2.4634697) q[3];
sx q[3];
rz(-1.1145385) q[3];
sx q[3];
rz(0.18517906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0878318) q[2];
sx q[2];
rz(-1.4542645) q[2];
sx q[2];
rz(-1.2254747) q[2];
rz(2.7671704) q[3];
sx q[3];
rz(-2.008308) q[3];
sx q[3];
rz(-0.55330127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9541009) q[0];
sx q[0];
rz(-0.88045374) q[0];
sx q[0];
rz(2.6075897) q[0];
rz(2.7502637) q[1];
sx q[1];
rz(-0.44328872) q[1];
sx q[1];
rz(-0.88723976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75636358) q[0];
sx q[0];
rz(-1.8779715) q[0];
sx q[0];
rz(-0.36184575) q[0];
rz(-pi) q[1];
rz(-0.5440622) q[2];
sx q[2];
rz(-2.3070222) q[2];
sx q[2];
rz(0.53235934) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45184093) q[1];
sx q[1];
rz(-1.4418169) q[1];
sx q[1];
rz(-2.9366628) q[1];
x q[2];
rz(1.0126566) q[3];
sx q[3];
rz(-0.59132708) q[3];
sx q[3];
rz(-1.8048353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54681626) q[2];
sx q[2];
rz(-1.3204601) q[2];
sx q[2];
rz(2.7776264) q[2];
rz(-1.724203) q[3];
sx q[3];
rz(-0.28984362) q[3];
sx q[3];
rz(2.3072306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.904838) q[0];
sx q[0];
rz(-0.098243864) q[0];
sx q[0];
rz(0.33016095) q[0];
rz(-2.5579021) q[1];
sx q[1];
rz(-0.81283641) q[1];
sx q[1];
rz(-0.49461734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8361417) q[0];
sx q[0];
rz(-1.3072805) q[0];
sx q[0];
rz(-1.2900773) q[0];
rz(-pi) q[1];
rz(-2.885899) q[2];
sx q[2];
rz(-1.0937368) q[2];
sx q[2];
rz(3.0050803) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3266625) q[1];
sx q[1];
rz(-1.0644134) q[1];
sx q[1];
rz(0.5270919) q[1];
rz(-pi) q[2];
rz(-1.4026138) q[3];
sx q[3];
rz(-0.59051149) q[3];
sx q[3];
rz(1.8174949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.79564774) q[2];
sx q[2];
rz(-2.4730885) q[2];
sx q[2];
rz(-0.12269679) q[2];
rz(-0.042938558) q[3];
sx q[3];
rz(-2.0624845) q[3];
sx q[3];
rz(-0.19449657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73404679) q[0];
sx q[0];
rz(-0.50739822) q[0];
sx q[0];
rz(-0.88974446) q[0];
rz(-0.45097688) q[1];
sx q[1];
rz(-2.273592) q[1];
sx q[1];
rz(-2.6873592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0182067) q[0];
sx q[0];
rz(-2.7834547) q[0];
sx q[0];
rz(-2.8337155) q[0];
rz(-1.5689911) q[2];
sx q[2];
rz(-1.9981401) q[2];
sx q[2];
rz(1.3722562) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6138623) q[1];
sx q[1];
rz(-2.3017028) q[1];
sx q[1];
rz(1.1682603) q[1];
x q[2];
rz(-2.7866335) q[3];
sx q[3];
rz(-1.6433652) q[3];
sx q[3];
rz(-1.5956958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1526327) q[2];
sx q[2];
rz(-2.5071414) q[2];
sx q[2];
rz(-0.75759849) q[2];
rz(2.9518413) q[3];
sx q[3];
rz(-1.3187783) q[3];
sx q[3];
rz(2.5943713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9309288) q[0];
sx q[0];
rz(-0.78721109) q[0];
sx q[0];
rz(-1.2934562) q[0];
rz(0.58213195) q[1];
sx q[1];
rz(-1.4385834) q[1];
sx q[1];
rz(-0.17939803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.005527786) q[0];
sx q[0];
rz(-1.3910157) q[0];
sx q[0];
rz(1.3639569) q[0];
x q[1];
rz(0.67082884) q[2];
sx q[2];
rz(-1.7636904) q[2];
sx q[2];
rz(-2.2191032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7779667) q[1];
sx q[1];
rz(-2.1898664) q[1];
sx q[1];
rz(2.9113063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2675293) q[3];
sx q[3];
rz(-1.9692326) q[3];
sx q[3];
rz(1.58942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.77357117) q[2];
sx q[2];
rz(-0.37006912) q[2];
sx q[2];
rz(0.5698815) q[2];
rz(-2.701345) q[3];
sx q[3];
rz(-1.2443685) q[3];
sx q[3];
rz(-0.35278916) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2920947) q[0];
sx q[0];
rz(-0.26420132) q[0];
sx q[0];
rz(-2.0000892) q[0];
rz(-1.796272) q[1];
sx q[1];
rz(-0.99015403) q[1];
sx q[1];
rz(-0.17123953) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6692176) q[0];
sx q[0];
rz(-0.49117377) q[0];
sx q[0];
rz(-0.21531658) q[0];
rz(-2.6944567) q[2];
sx q[2];
rz(-2.0512181) q[2];
sx q[2];
rz(2.8775079) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6936612) q[1];
sx q[1];
rz(-1.0180001) q[1];
sx q[1];
rz(-1.3752244) q[1];
rz(3.1290637) q[3];
sx q[3];
rz(-2.1274779) q[3];
sx q[3];
rz(2.0509246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3693927) q[2];
sx q[2];
rz(-0.93595305) q[2];
sx q[2];
rz(-3.0041223) q[2];
rz(1.2482268) q[3];
sx q[3];
rz(-0.56709254) q[3];
sx q[3];
rz(1.221777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0565979) q[0];
sx q[0];
rz(-0.63768142) q[0];
sx q[0];
rz(1.4884663) q[0];
rz(-0.78530637) q[1];
sx q[1];
rz(-1.4005125) q[1];
sx q[1];
rz(-1.8280425) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0037601622) q[0];
sx q[0];
rz(-1.6562471) q[0];
sx q[0];
rz(-1.6944044) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4573757) q[2];
sx q[2];
rz(-1.7608425) q[2];
sx q[2];
rz(-2.6148877) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.25048316) q[1];
sx q[1];
rz(-1.6598071) q[1];
sx q[1];
rz(-1.5989892) q[1];
rz(1.63229) q[3];
sx q[3];
rz(-2.4447421) q[3];
sx q[3];
rz(1.9887672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58275783) q[2];
sx q[2];
rz(-2.1119327) q[2];
sx q[2];
rz(0.47490698) q[2];
rz(-0.20259914) q[3];
sx q[3];
rz(-0.83297268) q[3];
sx q[3];
rz(-1.9995662) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1982034) q[0];
sx q[0];
rz(-0.13309637) q[0];
sx q[0];
rz(2.8357847) q[0];
rz(0.43931475) q[1];
sx q[1];
rz(-2.3190934) q[1];
sx q[1];
rz(1.3154715) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1179349) q[0];
sx q[0];
rz(-0.88738897) q[0];
sx q[0];
rz(1.5562431) q[0];
rz(2.010263) q[2];
sx q[2];
rz(-2.5410598) q[2];
sx q[2];
rz(0.10490049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0076239) q[1];
sx q[1];
rz(-0.79500145) q[1];
sx q[1];
rz(2.9717173) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9813811) q[3];
sx q[3];
rz(-1.2228612) q[3];
sx q[3];
rz(-0.98435005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7266015) q[2];
sx q[2];
rz(-1.5982268) q[2];
sx q[2];
rz(2.7169363) q[2];
rz(-3.1096544) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(2.4922075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53751078) q[0];
sx q[0];
rz(-2.4896955) q[0];
sx q[0];
rz(0.55484581) q[0];
rz(-2.8765053) q[1];
sx q[1];
rz(-1.4684497) q[1];
sx q[1];
rz(1.9404985) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4477312) q[0];
sx q[0];
rz(-2.434512) q[0];
sx q[0];
rz(-0.62753824) q[0];
rz(-pi) q[1];
rz(-2.5910406) q[2];
sx q[2];
rz(-0.32062396) q[2];
sx q[2];
rz(0.52832802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8256067) q[1];
sx q[1];
rz(-2.791643) q[1];
sx q[1];
rz(-3.0341329) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5463806) q[3];
sx q[3];
rz(-2.3204969) q[3];
sx q[3];
rz(1.8745619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.889664) q[2];
sx q[2];
rz(-0.98968518) q[2];
sx q[2];
rz(2.2970693) q[2];
rz(-1.1318413) q[3];
sx q[3];
rz(-2.2364538) q[3];
sx q[3];
rz(1.7822781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1691386) q[0];
sx q[0];
rz(-0.3903946) q[0];
sx q[0];
rz(1.1727232) q[0];
rz(2.4523465) q[1];
sx q[1];
rz(-2.0447562) q[1];
sx q[1];
rz(-0.26395878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.496752) q[0];
sx q[0];
rz(-0.94366108) q[0];
sx q[0];
rz(2.1261313) q[0];
x q[1];
rz(-2.1446635) q[2];
sx q[2];
rz(-2.3341114) q[2];
sx q[2];
rz(-1.6107744) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.86706272) q[1];
sx q[1];
rz(-1.2559051) q[1];
sx q[1];
rz(-0.19537433) q[1];
rz(-pi) q[2];
rz(0.12216025) q[3];
sx q[3];
rz(-2.2111243) q[3];
sx q[3];
rz(-1.8807008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3802203) q[2];
sx q[2];
rz(-1.1470046) q[2];
sx q[2];
rz(-2.0210361) q[2];
rz(-1.5348966) q[3];
sx q[3];
rz(-2.1963019) q[3];
sx q[3];
rz(1.5103316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0705538) q[0];
sx q[0];
rz(-0.65503913) q[0];
sx q[0];
rz(-0.10792637) q[0];
rz(0.72036605) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(-2.3158144) q[2];
sx q[2];
rz(-2.0407531) q[2];
sx q[2];
rz(-2.4498418) q[2];
rz(2.5247305) q[3];
sx q[3];
rz(-1.6322144) q[3];
sx q[3];
rz(-1.4323533) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
