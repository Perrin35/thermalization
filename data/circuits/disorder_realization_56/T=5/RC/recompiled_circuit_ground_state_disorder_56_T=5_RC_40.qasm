OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11256448) q[0];
sx q[0];
rz(-1.636314) q[0];
sx q[0];
rz(0.78328744) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(-2.3496871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0884704) q[0];
sx q[0];
rz(-1.3849918) q[0];
sx q[0];
rz(-2.8775999) q[0];
rz(-pi) q[1];
rz(1.7246805) q[2];
sx q[2];
rz(-1.6072011) q[2];
sx q[2];
rz(-0.77637451) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4152849) q[1];
sx q[1];
rz(-1.0276762) q[1];
sx q[1];
rz(-1.607479) q[1];
rz(-pi) q[2];
rz(0.79444076) q[3];
sx q[3];
rz(-2.0047965) q[3];
sx q[3];
rz(1.7931012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.59114328) q[2];
sx q[2];
rz(-0.14269665) q[2];
sx q[2];
rz(3.0421416) q[2];
rz(-2.8948696) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(2.7439086) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196911) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(0.9285399) q[0];
rz(-1.4506725) q[1];
sx q[1];
rz(-1.2932777) q[1];
sx q[1];
rz(0.64847747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49199864) q[0];
sx q[0];
rz(-1.0186983) q[0];
sx q[0];
rz(-0.92574889) q[0];
rz(0.47467741) q[2];
sx q[2];
rz(-0.62126489) q[2];
sx q[2];
rz(2.9062627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3120756) q[1];
sx q[1];
rz(-1.5551928) q[1];
sx q[1];
rz(-2.9996458) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0233211) q[3];
sx q[3];
rz(-1.857548) q[3];
sx q[3];
rz(-2.3210382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4497946) q[2];
sx q[2];
rz(-0.75642502) q[2];
sx q[2];
rz(-2.2890384) q[2];
rz(-0.84613386) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-2.1107296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(-2.0563828) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.4986821) q[1];
sx q[1];
rz(1.5012213) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417898) q[0];
sx q[0];
rz(-2.9405624) q[0];
sx q[0];
rz(-0.19636671) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9502421) q[2];
sx q[2];
rz(-1.3406959) q[2];
sx q[2];
rz(1.5528581) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2861917) q[1];
sx q[1];
rz(-2.5603189) q[1];
sx q[1];
rz(-1.3894677) q[1];
x q[2];
rz(-0.15764938) q[3];
sx q[3];
rz(-1.356989) q[3];
sx q[3];
rz(0.42663867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40470716) q[2];
sx q[2];
rz(-2.6077304) q[2];
sx q[2];
rz(-2.286818) q[2];
rz(0.9084304) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(-2.0250208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(2.214355) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(-1.9450564) q[0];
rz(1.475097) q[1];
sx q[1];
rz(-2.7207082) q[1];
sx q[1];
rz(1.1845142) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6086585) q[0];
sx q[0];
rz(-2.2581165) q[0];
sx q[0];
rz(-1.88009) q[0];
rz(-2.3886613) q[2];
sx q[2];
rz(-1.4309466) q[2];
sx q[2];
rz(1.953973) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6946751) q[1];
sx q[1];
rz(-1.0694587) q[1];
sx q[1];
rz(2.8828709) q[1];
rz(-pi) q[2];
rz(2.5585737) q[3];
sx q[3];
rz(-2.8172917) q[3];
sx q[3];
rz(-0.7079269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92675942) q[2];
sx q[2];
rz(-0.49761179) q[2];
sx q[2];
rz(-1.8801749) q[2];
rz(0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(-2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.072902) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(3.1121837) q[0];
rz(0.75621653) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(-2.3888033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8653976) q[0];
sx q[0];
rz(-1.5711693) q[0];
sx q[0];
rz(-1.5764109) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2529066) q[2];
sx q[2];
rz(-2.860002) q[2];
sx q[2];
rz(-1.3791305) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1323586) q[1];
sx q[1];
rz(-2.1387324) q[1];
sx q[1];
rz(1.6949937) q[1];
x q[2];
rz(-2.4642508) q[3];
sx q[3];
rz(-0.65653446) q[3];
sx q[3];
rz(-1.3765821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(-2.746554) q[2];
rz(2.6664074) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(2.7648259) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659371) q[0];
sx q[0];
rz(-2.4212615) q[0];
sx q[0];
rz(-2.1790867) q[0];
rz(-2.6267701) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(-0.33822507) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047236) q[0];
sx q[0];
rz(-0.77381182) q[0];
sx q[0];
rz(1.2880727) q[0];
rz(2.9043496) q[2];
sx q[2];
rz(-1.8117935) q[2];
sx q[2];
rz(2.4324696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2185827) q[1];
sx q[1];
rz(-1.5707695) q[1];
sx q[1];
rz(-1.5739417) q[1];
x q[2];
rz(-1.7231483) q[3];
sx q[3];
rz(-2.7151516) q[3];
sx q[3];
rz(1.3063198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.50011355) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(0.55434736) q[2];
rz(2.2902299) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(0.62801492) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3250658) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(0.75041962) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.3651747) q[1];
sx q[1];
rz(2.4837928) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80245012) q[0];
sx q[0];
rz(-0.93451148) q[0];
sx q[0];
rz(-2.0872981) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.051498895) q[2];
sx q[2];
rz(-2.4337075) q[2];
sx q[2];
rz(-2.5619363) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3624576) q[1];
sx q[1];
rz(-1.5541847) q[1];
sx q[1];
rz(-1.8569059) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9522454) q[3];
sx q[3];
rz(-1.7717517) q[3];
sx q[3];
rz(-2.902346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49360069) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(1.6592525) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928612) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(0.29801512) q[0];
rz(2.1030078) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.15377741) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2033635) q[0];
sx q[0];
rz(-1.3171609) q[0];
sx q[0];
rz(-1.4934191) q[0];
x q[1];
rz(-0.91570274) q[2];
sx q[2];
rz(-0.88353744) q[2];
sx q[2];
rz(2.5926431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1250026) q[1];
sx q[1];
rz(-1.224813) q[1];
sx q[1];
rz(1.6660652) q[1];
rz(-1.3078717) q[3];
sx q[3];
rz(-2.701303) q[3];
sx q[3];
rz(2.0844946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.095857233) q[2];
sx q[2];
rz(-2.6726674) q[2];
sx q[2];
rz(1.1886965) q[2];
rz(1.3304322) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(-0.73330283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7981912) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(1.6424302) q[0];
rz(2.5921953) q[1];
sx q[1];
rz(-1.5645942) q[1];
sx q[1];
rz(-0.25064358) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1090013) q[0];
sx q[0];
rz(-1.8792986) q[0];
sx q[0];
rz(0.89621131) q[0];
rz(-pi) q[1];
x q[1];
rz(0.05392404) q[2];
sx q[2];
rz(-0.26130518) q[2];
sx q[2];
rz(0.037029412) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.033015164) q[1];
sx q[1];
rz(-2.3653154) q[1];
sx q[1];
rz(-3.0423711) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9467877) q[3];
sx q[3];
rz(-1.6945632) q[3];
sx q[3];
rz(2.3969802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15829076) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(3.0554092) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30855274) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(-0.41900751) q[0];
rz(0.53681701) q[1];
sx q[1];
rz(-2.5173126) q[1];
sx q[1];
rz(-2.4057665) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8119278) q[0];
sx q[0];
rz(-2.4617534) q[0];
sx q[0];
rz(-2.832704) q[0];
x q[1];
rz(1.0334098) q[2];
sx q[2];
rz(-0.99796406) q[2];
sx q[2];
rz(2.9705641) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.23259737) q[1];
sx q[1];
rz(-2.4687662) q[1];
sx q[1];
rz(-0.35079591) q[1];
rz(-pi) q[2];
rz(-2.5352468) q[3];
sx q[3];
rz(-2.5744573) q[3];
sx q[3];
rz(-2.4727269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1824823) q[2];
sx q[2];
rz(-0.82254326) q[2];
sx q[2];
rz(-0.62270069) q[2];
rz(2.4032118) q[3];
sx q[3];
rz(-1.1850971) q[3];
sx q[3];
rz(-0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.57711346) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(2.4979757) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(2.9109091) q[2];
sx q[2];
rz(-1.9311957) q[2];
sx q[2];
rz(0.25372505) q[2];
rz(1.4293115) q[3];
sx q[3];
rz(-1.8093997) q[3];
sx q[3];
rz(2.3499478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
