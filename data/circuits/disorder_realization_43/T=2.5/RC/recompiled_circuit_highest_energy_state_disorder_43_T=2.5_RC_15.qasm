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
rz(0.046959538) q[0];
sx q[0];
rz(-1.5488012) q[0];
sx q[0];
rz(-1.7472851) q[0];
rz(2.6810763) q[1];
sx q[1];
rz(-1.2682275) q[1];
sx q[1];
rz(0.4761129) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6980909) q[0];
sx q[0];
rz(-0.68572581) q[0];
sx q[0];
rz(0.94615191) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85501315) q[2];
sx q[2];
rz(-0.35021338) q[2];
sx q[2];
rz(-0.24392715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81836593) q[1];
sx q[1];
rz(-1.4225679) q[1];
sx q[1];
rz(1.5532975) q[1];
x q[2];
rz(0.027276518) q[3];
sx q[3];
rz(-1.6347872) q[3];
sx q[3];
rz(-0.77426739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32690471) q[2];
sx q[2];
rz(-1.9196332) q[2];
sx q[2];
rz(1.2828705) q[2];
rz(-3.0021216) q[3];
sx q[3];
rz(-1.3160416) q[3];
sx q[3];
rz(0.68525806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4271456) q[0];
sx q[0];
rz(-0.35816631) q[0];
sx q[0];
rz(-0.32592475) q[0];
rz(2.4320995) q[1];
sx q[1];
rz(-0.26115099) q[1];
sx q[1];
rz(-2.479898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6174588) q[0];
sx q[0];
rz(-1.5633928) q[0];
sx q[0];
rz(1.5612768) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93962713) q[2];
sx q[2];
rz(-1.1568312) q[2];
sx q[2];
rz(-2.8394073) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1159617) q[1];
sx q[1];
rz(-0.44646663) q[1];
sx q[1];
rz(1.0122416) q[1];
rz(-2.4567408) q[3];
sx q[3];
rz(-2.0461444) q[3];
sx q[3];
rz(-0.97739391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0207396) q[2];
sx q[2];
rz(-1.3384621) q[2];
sx q[2];
rz(0.2600812) q[2];
rz(-0.86380473) q[3];
sx q[3];
rz(-1.443202) q[3];
sx q[3];
rz(-1.1687733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9048555) q[0];
sx q[0];
rz(-2.3598292) q[0];
sx q[0];
rz(-2.8535063) q[0];
rz(-1.2068564) q[1];
sx q[1];
rz(-2.0286045) q[1];
sx q[1];
rz(2.3634214) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93302762) q[0];
sx q[0];
rz(-1.3817245) q[0];
sx q[0];
rz(0.87644942) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3834095) q[2];
sx q[2];
rz(-1.6155525) q[2];
sx q[2];
rz(-2.5999709) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84904387) q[1];
sx q[1];
rz(-2.0660225) q[1];
sx q[1];
rz(-2.0640813) q[1];
rz(-2.7500626) q[3];
sx q[3];
rz(-2.1804296) q[3];
sx q[3];
rz(-0.044332835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.244016) q[2];
sx q[2];
rz(-1.1268758) q[2];
sx q[2];
rz(-0.24813949) q[2];
rz(-1.0330307) q[3];
sx q[3];
rz(-1.4890198) q[3];
sx q[3];
rz(-1.2584794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4397864) q[0];
sx q[0];
rz(-0.93455625) q[0];
sx q[0];
rz(-0.55950657) q[0];
rz(1.7350381) q[1];
sx q[1];
rz(-0.94466698) q[1];
sx q[1];
rz(1.14538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37811144) q[0];
sx q[0];
rz(-1.1147) q[0];
sx q[0];
rz(1.5198288) q[0];
rz(-2.4907095) q[2];
sx q[2];
rz(-2.0325629) q[2];
sx q[2];
rz(-2.8956763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.134702) q[1];
sx q[1];
rz(-2.3242053) q[1];
sx q[1];
rz(-2.7322053) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1405997) q[3];
sx q[3];
rz(-1.583433) q[3];
sx q[3];
rz(2.3848052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8765325) q[2];
sx q[2];
rz(-1.0141076) q[2];
sx q[2];
rz(-1.2476791) q[2];
rz(2.0606591) q[3];
sx q[3];
rz(-1.388988) q[3];
sx q[3];
rz(-1.2601674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0718229) q[0];
sx q[0];
rz(-0.48957303) q[0];
sx q[0];
rz(0.74482942) q[0];
rz(-1.7597594) q[1];
sx q[1];
rz(-1.1188743) q[1];
sx q[1];
rz(-3.0435069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6150655) q[0];
sx q[0];
rz(-1.9221109) q[0];
sx q[0];
rz(-1.243958) q[0];
x q[1];
rz(2.6977243) q[2];
sx q[2];
rz(-1.265268) q[2];
sx q[2];
rz(1.533184) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67691117) q[1];
sx q[1];
rz(-1.5524852) q[1];
sx q[1];
rz(2.4478138) q[1];
x q[2];
rz(0.10048203) q[3];
sx q[3];
rz(-2.049787) q[3];
sx q[3];
rz(2.8707531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7006526) q[2];
sx q[2];
rz(-0.38267371) q[2];
sx q[2];
rz(-1.0901701) q[2];
rz(-0.32554659) q[3];
sx q[3];
rz(-1.6306337) q[3];
sx q[3];
rz(-0.24698273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9149822) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(1.2363303) q[0];
rz(0.63713282) q[1];
sx q[1];
rz(-2.5814711) q[1];
sx q[1];
rz(-2.2083652) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1655344) q[0];
sx q[0];
rz(-1.9325496) q[0];
sx q[0];
rz(0.975859) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1652214) q[2];
sx q[2];
rz(-2.6898807) q[2];
sx q[2];
rz(-3.0690985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5781587) q[1];
sx q[1];
rz(-1.2799885) q[1];
sx q[1];
rz(-1.4842503) q[1];
x q[2];
rz(2.0108443) q[3];
sx q[3];
rz(-0.55956105) q[3];
sx q[3];
rz(2.6236629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8895662) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(1.8760366) q[2];
rz(-1.4243852) q[3];
sx q[3];
rz(-0.5286743) q[3];
sx q[3];
rz(-0.85844794) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76765656) q[0];
sx q[0];
rz(-0.38145426) q[0];
sx q[0];
rz(2.9071627) q[0];
rz(-2.7081721) q[1];
sx q[1];
rz(-1.0973009) q[1];
sx q[1];
rz(-1.1385328) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2018801) q[0];
sx q[0];
rz(-0.98390401) q[0];
sx q[0];
rz(0.71016772) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01485558) q[2];
sx q[2];
rz(-2.6508923) q[2];
sx q[2];
rz(2.4717836) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1416152) q[1];
sx q[1];
rz(-0.50302282) q[1];
sx q[1];
rz(-1.6254237) q[1];
rz(2.906231) q[3];
sx q[3];
rz(-0.70066626) q[3];
sx q[3];
rz(1.2149902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66702691) q[2];
sx q[2];
rz(-0.67956769) q[2];
sx q[2];
rz(-0.95728528) q[2];
rz(-0.75753093) q[3];
sx q[3];
rz(-1.4742955) q[3];
sx q[3];
rz(2.9981414) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677419) q[0];
sx q[0];
rz(-1.1026646) q[0];
sx q[0];
rz(2.962501) q[0];
rz(-0.20009072) q[1];
sx q[1];
rz(-0.6691907) q[1];
sx q[1];
rz(-2.3650513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6497407) q[0];
sx q[0];
rz(-1.0488247) q[0];
sx q[0];
rz(2.9328547) q[0];
x q[1];
rz(-1.1626622) q[2];
sx q[2];
rz(-2.4838349) q[2];
sx q[2];
rz(0.85504261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0697666) q[1];
sx q[1];
rz(-2.7785795) q[1];
sx q[1];
rz(0.20787225) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7257878) q[3];
sx q[3];
rz(-2.8202882) q[3];
sx q[3];
rz(0.50866541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.17579235) q[2];
sx q[2];
rz(-1.0719904) q[2];
sx q[2];
rz(-1.3624066) q[2];
rz(-1.8386748) q[3];
sx q[3];
rz(-1.4785654) q[3];
sx q[3];
rz(-2.1030857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4565444) q[0];
sx q[0];
rz(-1.1948723) q[0];
sx q[0];
rz(3.0273279) q[0];
rz(1.9396797) q[1];
sx q[1];
rz(-2.6004531) q[1];
sx q[1];
rz(-0.094954403) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72391188) q[0];
sx q[0];
rz(-1.0474993) q[0];
sx q[0];
rz(-2.9375034) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6812115) q[2];
sx q[2];
rz(-2.2308439) q[2];
sx q[2];
rz(2.6766863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1324117) q[1];
sx q[1];
rz(-0.55284111) q[1];
sx q[1];
rz(0.10709672) q[1];
rz(-pi) q[2];
rz(-2.0801274) q[3];
sx q[3];
rz(-1.2042787) q[3];
sx q[3];
rz(-1.4118005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3007043) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(0.67187205) q[2];
rz(0.68615174) q[3];
sx q[3];
rz(-2.4785564) q[3];
sx q[3];
rz(-1.2175951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9753863) q[0];
sx q[0];
rz(-2.9500742) q[0];
sx q[0];
rz(1.5331049) q[0];
rz(-2.8168822) q[1];
sx q[1];
rz(-1.277781) q[1];
sx q[1];
rz(0.3130354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1192577) q[0];
sx q[0];
rz(-2.5637163) q[0];
sx q[0];
rz(1.7176601) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3029477) q[2];
sx q[2];
rz(-2.9137879) q[2];
sx q[2];
rz(1.1033664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6021784) q[1];
sx q[1];
rz(-2.115332) q[1];
sx q[1];
rz(-2.3474752) q[1];
rz(-pi) q[2];
rz(0.62449091) q[3];
sx q[3];
rz(-1.9510799) q[3];
sx q[3];
rz(-0.46562574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2031871) q[2];
sx q[2];
rz(-1.1537735) q[2];
sx q[2];
rz(-1.3333092) q[2];
rz(2.5238254) q[3];
sx q[3];
rz(-2.1963162) q[3];
sx q[3];
rz(1.1355404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4278605) q[0];
sx q[0];
rz(-1.8129616) q[0];
sx q[0];
rz(0.34995361) q[0];
rz(-2.2433544) q[1];
sx q[1];
rz(-2.0571092) q[1];
sx q[1];
rz(2.7877997) q[1];
rz(-0.66302115) q[2];
sx q[2];
rz(-2.2793179) q[2];
sx q[2];
rz(0.09003612) q[2];
rz(2.4118123) q[3];
sx q[3];
rz(-1.5530219) q[3];
sx q[3];
rz(-1.1698759) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
