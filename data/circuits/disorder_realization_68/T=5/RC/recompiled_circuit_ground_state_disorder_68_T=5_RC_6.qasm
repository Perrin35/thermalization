OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5691583) q[0];
sx q[0];
rz(-0.99826607) q[0];
sx q[0];
rz(-2.7531667) q[0];
rz(-1.2216964) q[1];
sx q[1];
rz(0.95284) q[1];
sx q[1];
rz(9.3508773) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9572417) q[0];
sx q[0];
rz(-1.3660407) q[0];
sx q[0];
rz(-2.6100169) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5396893) q[2];
sx q[2];
rz(-2.4894161) q[2];
sx q[2];
rz(-1.3155441) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5351341) q[1];
sx q[1];
rz(-0.89790895) q[1];
sx q[1];
rz(2.6752641) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9268215) q[3];
sx q[3];
rz(-1.3269375) q[3];
sx q[3];
rz(-0.36892316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4581603) q[2];
sx q[2];
rz(-1.4417803) q[2];
sx q[2];
rz(1.3275576) q[2];
rz(2.1761927) q[3];
sx q[3];
rz(-0.5894956) q[3];
sx q[3];
rz(-1.0546257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4210159) q[0];
sx q[0];
rz(-3.133931) q[0];
sx q[0];
rz(1.7505919) q[0];
rz(1.947047) q[1];
sx q[1];
rz(-2.5321913) q[1];
sx q[1];
rz(2.4360099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4509038) q[0];
sx q[0];
rz(-0.92508679) q[0];
sx q[0];
rz(0.63668191) q[0];
x q[1];
rz(-2.4439677) q[2];
sx q[2];
rz(-1.910733) q[2];
sx q[2];
rz(-1.4125669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2516219) q[1];
sx q[1];
rz(-1.9202616) q[1];
sx q[1];
rz(0.98213338) q[1];
x q[2];
rz(2.0749749) q[3];
sx q[3];
rz(-1.0312349) q[3];
sx q[3];
rz(-0.93531306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.07448639) q[2];
sx q[2];
rz(-1.4718082) q[2];
sx q[2];
rz(2.2763695) q[2];
rz(-1.5857961) q[3];
sx q[3];
rz(-2.9977048) q[3];
sx q[3];
rz(-1.5998862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0921606) q[0];
sx q[0];
rz(-2.3724174) q[0];
sx q[0];
rz(-2.8014158) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1454606) q[2];
sx q[2];
rz(-0.43197235) q[2];
sx q[2];
rz(-1.4086823) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8286491) q[1];
sx q[1];
rz(-0.11308453) q[1];
sx q[1];
rz(2.7663016) q[1];
rz(1.7060852) q[3];
sx q[3];
rz(-2.0485037) q[3];
sx q[3];
rz(-2.6729134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1323041) q[2];
sx q[2];
rz(-1.4762286) q[2];
sx q[2];
rz(2.0557527) q[2];
rz(-3.0911607) q[3];
sx q[3];
rz(-1.4112873) q[3];
sx q[3];
rz(2.1850695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0352935) q[0];
sx q[0];
rz(-1.6989919) q[0];
sx q[0];
rz(-2.371149) q[0];
rz(2.2162407) q[1];
sx q[1];
rz(-1.4848361) q[1];
sx q[1];
rz(2.3063708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90635008) q[0];
sx q[0];
rz(-2.6444204) q[0];
sx q[0];
rz(-0.034336523) q[0];
x q[1];
rz(2.9000834) q[2];
sx q[2];
rz(-1.4837974) q[2];
sx q[2];
rz(-2.3149025) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8586352) q[1];
sx q[1];
rz(-1.5581308) q[1];
sx q[1];
rz(2.0974166) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20056574) q[3];
sx q[3];
rz(-1.6839538) q[3];
sx q[3];
rz(3.1386047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.931382) q[2];
sx q[2];
rz(-1.8385734) q[2];
sx q[2];
rz(-2.4118928) q[2];
rz(2.1424255) q[3];
sx q[3];
rz(-1.8347284) q[3];
sx q[3];
rz(0.45381799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.5900742) q[0];
sx q[0];
rz(-0.22576627) q[0];
sx q[0];
rz(-0.99910587) q[0];
rz(-1.0480115) q[1];
sx q[1];
rz(-2.4297355) q[1];
sx q[1];
rz(2.5968831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42349619) q[0];
sx q[0];
rz(-2.831651) q[0];
sx q[0];
rz(-0.62356068) q[0];
x q[1];
rz(0.72417132) q[2];
sx q[2];
rz(-1.0909547) q[2];
sx q[2];
rz(-1.9373231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5156158) q[1];
sx q[1];
rz(-0.57913172) q[1];
sx q[1];
rz(0.34767751) q[1];
x q[2];
rz(-1.5673166) q[3];
sx q[3];
rz(-1.165357) q[3];
sx q[3];
rz(-0.60379782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.58392945) q[2];
sx q[2];
rz(-1.1489392) q[2];
sx q[2];
rz(-2.124713) q[2];
rz(0.66295463) q[3];
sx q[3];
rz(-0.69279492) q[3];
sx q[3];
rz(-1.3346765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7826295) q[0];
sx q[0];
rz(-2.6435659) q[0];
sx q[0];
rz(-1.1516512) q[0];
rz(-0.89371124) q[1];
sx q[1];
rz(-1.3795373) q[1];
sx q[1];
rz(-3.025257) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49580816) q[0];
sx q[0];
rz(-2.6538355) q[0];
sx q[0];
rz(-1.3402104) q[0];
rz(1.9151808) q[2];
sx q[2];
rz(-1.5874327) q[2];
sx q[2];
rz(0.045147506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4018206) q[1];
sx q[1];
rz(-1.7171613) q[1];
sx q[1];
rz(-1.0159645) q[1];
x q[2];
rz(1.4907012) q[3];
sx q[3];
rz(-2.3175196) q[3];
sx q[3];
rz(0.88300812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5748888) q[2];
sx q[2];
rz(-2.7723007) q[2];
sx q[2];
rz(-1.0075547) q[2];
rz(1.9762074) q[3];
sx q[3];
rz(-0.27519614) q[3];
sx q[3];
rz(-2.5113441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769058) q[0];
sx q[0];
rz(-0.46606627) q[0];
sx q[0];
rz(0.16978547) q[0];
rz(-0.67974293) q[1];
sx q[1];
rz(-0.83234537) q[1];
sx q[1];
rz(0.92175093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0089192275) q[0];
sx q[0];
rz(-0.28604315) q[0];
sx q[0];
rz(0.45222767) q[0];
x q[1];
rz(1.6255195) q[2];
sx q[2];
rz(-1.2642908) q[2];
sx q[2];
rz(-2.3505369) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8584916) q[1];
sx q[1];
rz(-0.6226317) q[1];
sx q[1];
rz(0.94515349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60867057) q[3];
sx q[3];
rz(-2.8117315) q[3];
sx q[3];
rz(2.4467227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0035231) q[2];
sx q[2];
rz(-1.4813083) q[2];
sx q[2];
rz(1.8180234) q[2];
rz(2.0924163) q[3];
sx q[3];
rz(-2.0101571) q[3];
sx q[3];
rz(1.7607035) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2915989) q[0];
sx q[0];
rz(-2.0065362) q[0];
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
rz(-1.749925) q[0];
sx q[0];
rz(-0.13340575) q[0];
sx q[0];
rz(-1.6682503) q[0];
rz(-pi) q[1];
rz(0.86470072) q[2];
sx q[2];
rz(-1.0664202) q[2];
sx q[2];
rz(-2.9838533) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50302153) q[1];
sx q[1];
rz(-1.5201228) q[1];
sx q[1];
rz(2.1524968) q[1];
rz(-pi) q[2];
rz(2.9696736) q[3];
sx q[3];
rz(-1.0128331) q[3];
sx q[3];
rz(-2.362988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1524973) q[2];
sx q[2];
rz(-0.54741198) q[2];
sx q[2];
rz(3.1067749) q[2];
rz(-3.1133437) q[3];
sx q[3];
rz(-2.0939128) q[3];
sx q[3];
rz(-2.4475173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6377653) q[0];
sx q[0];
rz(-2.4389508) q[0];
sx q[0];
rz(-2.0565865) q[0];
rz(1.7833692) q[1];
sx q[1];
rz(-2.5085776) q[1];
sx q[1];
rz(-2.0620652) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5732291) q[0];
sx q[0];
rz(-1.0160334) q[0];
sx q[0];
rz(-2.8994096) q[0];
rz(-pi) q[1];
rz(1.7100542) q[2];
sx q[2];
rz(-1.1258954) q[2];
sx q[2];
rz(-1.618096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.92380556) q[1];
sx q[1];
rz(-1.3239685) q[1];
sx q[1];
rz(1.8882165) q[1];
x q[2];
rz(-1.0113768) q[3];
sx q[3];
rz(-2.3314948) q[3];
sx q[3];
rz(1.7138077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2435771) q[2];
sx q[2];
rz(-1.5960627) q[2];
sx q[2];
rz(-3.0698981) q[2];
rz(2.535635) q[3];
sx q[3];
rz(-0.76064435) q[3];
sx q[3];
rz(0.073237091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3074985) q[0];
sx q[0];
rz(-1.6313169) q[0];
sx q[0];
rz(0.50258762) q[0];
rz(1.7723627) q[1];
sx q[1];
rz(-1.9160198) q[1];
sx q[1];
rz(-0.1756846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0640819) q[0];
sx q[0];
rz(-0.33987576) q[0];
sx q[0];
rz(-1.3114291) q[0];
rz(-pi) q[1];
rz(-2.2403342) q[2];
sx q[2];
rz(-2.2961535) q[2];
sx q[2];
rz(-2.8536316) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60424544) q[1];
sx q[1];
rz(-0.89902311) q[1];
sx q[1];
rz(0.91225454) q[1];
x q[2];
rz(-0.39145546) q[3];
sx q[3];
rz(-2.1717697) q[3];
sx q[3];
rz(1.5096926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48675576) q[2];
sx q[2];
rz(-1.882587) q[2];
sx q[2];
rz(-0.34461018) q[2];
rz(1.2279855) q[3];
sx q[3];
rz(-2.4495008) q[3];
sx q[3];
rz(2.1740348) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3053038) q[2];
sx q[2];
rz(-2.2659779) q[2];
sx q[2];
rz(-0.1018079) q[2];
rz(0.19560858) q[3];
sx q[3];
rz(-2.2785288) q[3];
sx q[3];
rz(-2.5768448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
