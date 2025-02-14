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
rz(1.5927915) q[0];
sx q[0];
rz(11.172063) q[0];
rz(2.6810763) q[1];
sx q[1];
rz(-1.2682275) q[1];
sx q[1];
rz(-2.6654798) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5232485) q[0];
sx q[0];
rz(-1.1914413) q[0];
sx q[0];
rz(2.1567287) q[0];
rz(-pi) q[1];
rz(-2.9063375) q[2];
sx q[2];
rz(-1.8326743) q[2];
sx q[2];
rz(-0.99083283) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.70041768) q[1];
sx q[1];
rz(-2.9923424) q[1];
sx q[1];
rz(0.11664893) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1143161) q[3];
sx q[3];
rz(-1.5068054) q[3];
sx q[3];
rz(2.3673253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32690471) q[2];
sx q[2];
rz(-1.9196332) q[2];
sx q[2];
rz(-1.8587221) q[2];
rz(0.13947105) q[3];
sx q[3];
rz(-1.8255511) q[3];
sx q[3];
rz(-0.68525806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4271456) q[0];
sx q[0];
rz(-0.35816631) q[0];
sx q[0];
rz(-2.8156679) q[0];
rz(-0.70949316) q[1];
sx q[1];
rz(-2.8804417) q[1];
sx q[1];
rz(2.479898) q[1];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.9847615) q[2];
sx q[2];
rz(2.8394073) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0256309) q[1];
sx q[1];
rz(-2.695126) q[1];
sx q[1];
rz(2.129351) q[1];
rz(-pi) q[2];
rz(0.68485188) q[3];
sx q[3];
rz(-1.0954482) q[3];
sx q[3];
rz(0.97739391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0207396) q[2];
sx q[2];
rz(-1.8031305) q[2];
sx q[2];
rz(-0.2600812) q[2];
rz(-2.2777879) q[3];
sx q[3];
rz(-1.6983906) q[3];
sx q[3];
rz(1.9728194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9048555) q[0];
sx q[0];
rz(-2.3598292) q[0];
sx q[0];
rz(2.8535063) q[0];
rz(1.2068564) q[1];
sx q[1];
rz(-1.1129881) q[1];
sx q[1];
rz(-0.77817121) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818277) q[0];
sx q[0];
rz(-2.4261103) q[0];
sx q[0];
rz(1.2802109) q[0];
rz(-pi) q[1];
rz(0.7581832) q[2];
sx q[2];
rz(-1.5260401) q[2];
sx q[2];
rz(2.5999709) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2925488) q[1];
sx q[1];
rz(-2.0660225) q[1];
sx q[1];
rz(-2.0640813) q[1];
x q[2];
rz(1.0707079) q[3];
sx q[3];
rz(-0.71083655) q[3];
sx q[3];
rz(2.5612166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.244016) q[2];
sx q[2];
rz(-1.1268758) q[2];
sx q[2];
rz(-2.8934532) q[2];
rz(1.0330307) q[3];
sx q[3];
rz(-1.4890198) q[3];
sx q[3];
rz(1.2584794) q[3];
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
rz(-1.7018062) q[0];
sx q[0];
rz(-2.2070364) q[0];
sx q[0];
rz(0.55950657) q[0];
rz(-1.7350381) q[1];
sx q[1];
rz(-2.1969257) q[1];
sx q[1];
rz(1.14538) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7634812) q[0];
sx q[0];
rz(-1.1147) q[0];
sx q[0];
rz(-1.5198288) q[0];
rz(-pi) q[1];
rz(2.4907095) q[2];
sx q[2];
rz(-2.0325629) q[2];
sx q[2];
rz(-0.24591638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.134702) q[1];
sx q[1];
rz(-0.81738735) q[1];
sx q[1];
rz(0.40938739) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5473751) q[3];
sx q[3];
rz(-0.56992793) q[3];
sx q[3];
rz(0.83372926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2650602) q[2];
sx q[2];
rz(-1.0141076) q[2];
sx q[2];
rz(-1.8939135) q[2];
rz(-1.0809336) q[3];
sx q[3];
rz(-1.7526046) q[3];
sx q[3];
rz(1.2601674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0697698) q[0];
sx q[0];
rz(-0.48957303) q[0];
sx q[0];
rz(2.3967632) q[0];
rz(-1.7597594) q[1];
sx q[1];
rz(-1.1188743) q[1];
sx q[1];
rz(-3.0435069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9811889) q[0];
sx q[0];
rz(-1.2645928) q[0];
sx q[0];
rz(-0.36925495) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44386835) q[2];
sx q[2];
rz(-1.265268) q[2];
sx q[2];
rz(-1.6084087) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.225696) q[1];
sx q[1];
rz(-0.69398038) q[1];
sx q[1];
rz(-0.028631239) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3799872) q[3];
sx q[3];
rz(-0.48861438) q[3];
sx q[3];
rz(0.055475108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7006526) q[2];
sx q[2];
rz(-0.38267371) q[2];
sx q[2];
rz(-2.0514226) q[2];
rz(0.32554659) q[3];
sx q[3];
rz(-1.6306337) q[3];
sx q[3];
rz(0.24698273) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9149822) q[0];
sx q[0];
rz(-2.60422) q[0];
sx q[0];
rz(1.2363303) q[0];
rz(2.5044598) q[1];
sx q[1];
rz(-2.5814711) q[1];
sx q[1];
rz(-0.93322745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0649945) q[0];
sx q[0];
rz(-0.68474283) q[0];
sx q[0];
rz(-2.1646654) q[0];
x q[1];
rz(-1.1652214) q[2];
sx q[2];
rz(-0.45171192) q[2];
sx q[2];
rz(-3.0690985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98249032) q[1];
sx q[1];
rz(-1.4878927) q[1];
sx q[1];
rz(-2.8497531) q[1];
rz(-2.8808571) q[3];
sx q[3];
rz(-2.0717295) q[3];
sx q[3];
rz(2.1165533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.25202641) q[2];
sx q[2];
rz(-1.7539975) q[2];
sx q[2];
rz(1.265556) q[2];
rz(-1.4243852) q[3];
sx q[3];
rz(-0.5286743) q[3];
sx q[3];
rz(2.2831447) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76765656) q[0];
sx q[0];
rz(-2.7601384) q[0];
sx q[0];
rz(0.23442991) q[0];
rz(2.7081721) q[1];
sx q[1];
rz(-2.0442918) q[1];
sx q[1];
rz(2.0030599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94110452) q[0];
sx q[0];
rz(-0.88747665) q[0];
sx q[0];
rz(-0.79534689) q[0];
rz(2.6509382) q[2];
sx q[2];
rz(-1.5777967) q[2];
sx q[2];
rz(0.8878844) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5229021) q[1];
sx q[1];
rz(-1.5971208) q[1];
sx q[1];
rz(1.0684038) q[1];
rz(1.7650003) q[3];
sx q[3];
rz(-2.248431) q[3];
sx q[3];
rz(1.6226131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.66702691) q[2];
sx q[2];
rz(-2.462025) q[2];
sx q[2];
rz(-0.95728528) q[2];
rz(0.75753093) q[3];
sx q[3];
rz(-1.6672971) q[3];
sx q[3];
rz(-0.1434513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677419) q[0];
sx q[0];
rz(-1.1026646) q[0];
sx q[0];
rz(0.17909166) q[0];
rz(2.9415019) q[1];
sx q[1];
rz(-0.6691907) q[1];
sx q[1];
rz(-2.3650513) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0901129) q[0];
sx q[0];
rz(-2.5830373) q[0];
sx q[0];
rz(1.9165975) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8440823) q[2];
sx q[2];
rz(-0.97504598) q[2];
sx q[2];
rz(1.7864625) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.69570615) q[1];
sx q[1];
rz(-1.6441455) q[1];
sx q[1];
rz(-0.35584764) q[1];
rz(1.8885047) q[3];
sx q[3];
rz(-1.6195669) q[3];
sx q[3];
rz(1.2093076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17579235) q[2];
sx q[2];
rz(-2.0696023) q[2];
sx q[2];
rz(1.779186) q[2];
rz(-1.3029178) q[3];
sx q[3];
rz(-1.6630273) q[3];
sx q[3];
rz(1.038507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68504828) q[0];
sx q[0];
rz(-1.9467204) q[0];
sx q[0];
rz(-3.0273279) q[0];
rz(-1.9396797) q[1];
sx q[1];
rz(-2.6004531) q[1];
sx q[1];
rz(-3.0466383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72391188) q[0];
sx q[0];
rz(-1.0474993) q[0];
sx q[0];
rz(0.20408922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0906779) q[2];
sx q[2];
rz(-2.3570545) q[2];
sx q[2];
rz(-2.926411) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.11647955) q[1];
sx q[1];
rz(-2.1201061) q[1];
sx q[1];
rz(-1.6366556) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41422959) q[3];
sx q[3];
rz(-1.0981907) q[3];
sx q[3];
rz(2.7850658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8408884) q[2];
sx q[2];
rz(-1.3450832) q[2];
sx q[2];
rz(0.67187205) q[2];
rz(-0.68615174) q[3];
sx q[3];
rz(-2.4785564) q[3];
sx q[3];
rz(-1.9239976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9753863) q[0];
sx q[0];
rz(-2.9500742) q[0];
sx q[0];
rz(1.6084877) q[0];
rz(0.32471049) q[1];
sx q[1];
rz(-1.8638116) q[1];
sx q[1];
rz(2.8285573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4251809) q[0];
sx q[0];
rz(-1.6508174) q[0];
sx q[0];
rz(2.1437313) q[0];
rz(-pi) q[1];
rz(-3.0803142) q[2];
sx q[2];
rz(-1.7903425) q[2];
sx q[2];
rz(-2.3128308) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6021784) q[1];
sx q[1];
rz(-1.0262607) q[1];
sx q[1];
rz(2.3474752) q[1];
rz(-2.5171017) q[3];
sx q[3];
rz(-1.9510799) q[3];
sx q[3];
rz(2.6759669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93840557) q[2];
sx q[2];
rz(-1.9878191) q[2];
sx q[2];
rz(1.3333092) q[2];
rz(2.5238254) q[3];
sx q[3];
rz(-0.94527644) q[3];
sx q[3];
rz(2.0060523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7137322) q[0];
sx q[0];
rz(-1.8129616) q[0];
sx q[0];
rz(0.34995361) q[0];
rz(0.89823828) q[1];
sx q[1];
rz(-2.0571092) q[1];
sx q[1];
rz(2.7877997) q[1];
rz(0.74358616) q[2];
sx q[2];
rz(-1.0844885) q[2];
sx q[2];
rz(1.1906638) q[2];
rz(-2.4118123) q[3];
sx q[3];
rz(-1.5885708) q[3];
sx q[3];
rz(1.9717168) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
