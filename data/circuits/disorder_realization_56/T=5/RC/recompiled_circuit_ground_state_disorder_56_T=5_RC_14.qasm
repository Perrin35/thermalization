OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0290282) q[0];
sx q[0];
rz(-1.5052786) q[0];
sx q[0];
rz(2.3583052) q[0];
rz(-2.0593491) q[1];
sx q[1];
rz(-1.6545656) q[1];
sx q[1];
rz(2.3496871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5740252) q[0];
sx q[0];
rz(-1.8301395) q[0];
sx q[0];
rz(1.7631084) q[0];
rz(-pi) q[1];
rz(-1.4169121) q[2];
sx q[2];
rz(-1.5343915) q[2];
sx q[2];
rz(-2.3652181) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3160682) q[1];
sx q[1];
rz(-1.5393942) q[1];
sx q[1];
rz(-2.5981748) q[1];
rz(-pi) q[2];
rz(-0.57609419) q[3];
sx q[3];
rz(-2.2597426) q[3];
sx q[3];
rz(-2.9722633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59114328) q[2];
sx q[2];
rz(-0.14269665) q[2];
sx q[2];
rz(-0.099451065) q[2];
rz(-0.24672306) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1219015) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(1.6909201) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(2.4931152) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4546616) q[0];
sx q[0];
rz(-2.1082011) q[0];
sx q[0];
rz(-2.4848558) q[0];
rz(-pi) q[1];
rz(-0.5669539) q[2];
sx q[2];
rz(-1.3015206) q[2];
sx q[2];
rz(-0.93967162) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3120756) q[1];
sx q[1];
rz(-1.5863998) q[1];
sx q[1];
rz(2.9996458) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4497946) q[2];
sx q[2];
rz(-2.3851676) q[2];
sx q[2];
rz(0.85255426) q[2];
rz(-0.84613386) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-2.1107296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(2.0563828) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(1.6403713) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417898) q[0];
sx q[0];
rz(-2.9405624) q[0];
sx q[0];
rz(-0.19636671) q[0];
rz(-1.8050212) q[2];
sx q[2];
rz(-1.7570436) q[2];
sx q[2];
rz(-3.0795003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(1.752125) q[1];
rz(-2.9839433) q[3];
sx q[3];
rz(-1.7846037) q[3];
sx q[3];
rz(0.42663867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7368855) q[2];
sx q[2];
rz(-2.6077304) q[2];
sx q[2];
rz(-2.286818) q[2];
rz(0.9084304) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(1.1165718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.214355) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(1.9450564) q[0];
rz(1.475097) q[1];
sx q[1];
rz(-2.7207082) q[1];
sx q[1];
rz(1.1845142) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0663721) q[0];
sx q[0];
rz(-2.3983404) q[0];
sx q[0];
rz(-0.35510285) q[0];
rz(-2.3886613) q[2];
sx q[2];
rz(-1.7106461) q[2];
sx q[2];
rz(-1.953973) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8912011) q[1];
sx q[1];
rz(-1.3445092) q[1];
sx q[1];
rz(-2.0865284) q[1];
x q[2];
rz(-2.5585737) q[3];
sx q[3];
rz(-2.8172917) q[3];
sx q[3];
rz(-2.4336658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92675942) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.8801749) q[2];
rz(2.7091806) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.072902) q[0];
sx q[0];
rz(-1.9111159) q[0];
sx q[0];
rz(-0.029408971) q[0];
rz(-0.75621653) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(-0.75278935) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9133224) q[0];
sx q[0];
rz(-3.1359657) q[0];
sx q[0];
rz(1.6371284) q[0];
x q[1];
rz(-0.090172099) q[2];
sx q[2];
rz(-1.3036696) q[2];
sx q[2];
rz(-1.7091441) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.90430561) q[1];
sx q[1];
rz(-0.5798961) q[1];
sx q[1];
rz(-0.19176264) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54075586) q[3];
sx q[3];
rz(-1.1782559) q[3];
sx q[3];
rz(2.768571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14043643) q[2];
sx q[2];
rz(-1.493528) q[2];
sx q[2];
rz(2.746554) q[2];
rz(2.6664074) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(-0.37676677) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756556) q[0];
sx q[0];
rz(-2.4212615) q[0];
sx q[0];
rz(-0.96250594) q[0];
rz(0.51482254) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(0.33822507) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8047236) q[0];
sx q[0];
rz(-0.77381182) q[0];
sx q[0];
rz(1.2880727) q[0];
rz(-pi) q[1];
rz(-1.8184616) q[2];
sx q[2];
rz(-1.3405352) q[2];
sx q[2];
rz(0.80402735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.64778642) q[1];
sx q[1];
rz(-1.5739417) q[1];
sx q[1];
rz(2.6814744e-05) q[1];
rz(0.068840222) q[3];
sx q[3];
rz(-1.9919812) q[3];
sx q[3];
rz(-2.0023458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6414791) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(-0.55434736) q[2];
rz(-0.85136271) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(-0.62801492) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3250658) q[0];
sx q[0];
rz(-0.61282235) q[0];
sx q[0];
rz(-0.75041962) q[0];
rz(2.5665307) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(2.4837928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039666273) q[0];
sx q[0];
rz(-0.79625477) q[0];
sx q[0];
rz(2.5523561) q[0];
rz(2.4343628) q[2];
sx q[2];
rz(-1.5373188) q[2];
sx q[2];
rz(-2.1113124) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
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
rz(0.23924669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.647992) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(1.4823401) q[2];
rz(0.12065398) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7928612) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(-2.8435775) q[0];
rz(2.1030078) q[1];
sx q[1];
rz(-0.54439259) q[1];
sx q[1];
rz(0.15377741) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2033635) q[0];
sx q[0];
rz(-1.8244317) q[0];
sx q[0];
rz(-1.4934191) q[0];
x q[1];
rz(2.3390017) q[2];
sx q[2];
rz(-1.0804515) q[2];
sx q[2];
rz(0.56832321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.01659) q[1];
sx q[1];
rz(-1.224813) q[1];
sx q[1];
rz(1.4755274) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9977536) q[3];
sx q[3];
rz(-1.4597963) q[3];
sx q[3];
rz(-2.3890561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0457354) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.9528961) q[2];
rz(-1.8111604) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(0.73330283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7981912) q[0];
sx q[0];
rz(-2.5372086) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(-2.5921953) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-0.25064358) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9662387) q[0];
sx q[0];
rz(-2.4099775) q[0];
sx q[0];
rz(2.0425955) q[0];
x q[1];
rz(2.8806503) q[2];
sx q[2];
rz(-1.5847209) q[2];
sx q[2];
rz(1.5557289) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.033015164) q[1];
sx q[1];
rz(-2.3653154) q[1];
sx q[1];
rz(-3.0423711) q[1];
x q[2];
rz(0.57119675) q[3];
sx q[3];
rz(-0.2303752) q[3];
sx q[3];
rz(-0.26709589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.15829076) q[2];
sx q[2];
rz(-3.0088708) q[2];
sx q[2];
rz(-2.1251202) q[2];
rz(-3.0554092) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(2.3763954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30855274) q[0];
sx q[0];
rz(-0.85423952) q[0];
sx q[0];
rz(-2.7225851) q[0];
rz(0.53681701) q[1];
sx q[1];
rz(-2.5173126) q[1];
sx q[1];
rz(0.73582617) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3296649) q[0];
sx q[0];
rz(-2.4617534) q[0];
sx q[0];
rz(-0.30888866) q[0];
x q[1];
rz(2.1081829) q[2];
sx q[2];
rz(-2.1436286) q[2];
sx q[2];
rz(-0.17102851) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6169301) q[1];
sx q[1];
rz(-1.3549651) q[1];
sx q[1];
rz(0.642435) q[1];
rz(-pi) q[2];
rz(0.60634585) q[3];
sx q[3];
rz(-0.56713533) q[3];
sx q[3];
rz(2.4727269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1824823) q[2];
sx q[2];
rz(-0.82254326) q[2];
sx q[2];
rz(-2.518892) q[2];
rz(0.73838082) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(2.8625989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644792) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(0.64361698) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(1.9401445) q[2];
sx q[2];
rz(-1.7864173) q[2];
sx q[2];
rz(-1.3997072) q[2];
rz(-0.52538659) q[3];
sx q[3];
rz(-0.27670866) q[3];
sx q[3];
rz(1.8075734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
