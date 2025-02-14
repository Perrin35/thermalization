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
rz(-2.3583052) q[0];
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
rz(0.91814268) q[0];
sx q[0];
rz(-0.32157102) q[0];
sx q[0];
rz(2.5173223) q[0];
rz(-pi) q[1];
rz(1.7246805) q[2];
sx q[2];
rz(-1.6072011) q[2];
sx q[2];
rz(-0.77637451) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72630771) q[1];
sx q[1];
rz(-2.1139164) q[1];
sx q[1];
rz(1.5341137) q[1];
x q[2];
rz(2.1551829) q[3];
sx q[3];
rz(-2.274868) q[3];
sx q[3];
rz(2.5147284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59114328) q[2];
sx q[2];
rz(-0.14269665) q[2];
sx q[2];
rz(3.0421416) q[2];
rz(2.8948696) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-2.7439086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196911) q[0];
sx q[0];
rz(-2.1653403) q[0];
sx q[0];
rz(-0.9285399) q[0];
rz(-1.6909201) q[1];
sx q[1];
rz(-1.2932777) q[1];
sx q[1];
rz(2.4931152) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6715393) q[0];
sx q[0];
rz(-0.82255615) q[0];
sx q[0];
rz(-0.7732735) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47467741) q[2];
sx q[2];
rz(-2.5203278) q[2];
sx q[2];
rz(2.9062627) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14997031) q[1];
sx q[1];
rz(-0.14279616) q[1];
sx q[1];
rz(-3.0317329) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8594555) q[3];
sx q[3];
rz(-1.6842173) q[3];
sx q[3];
rz(0.78383776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4497946) q[2];
sx q[2];
rz(-2.3851676) q[2];
sx q[2];
rz(-2.2890384) q[2];
rz(0.84613386) q[3];
sx q[3];
rz(-1.6328014) q[3];
sx q[3];
rz(-2.1107296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.5071252) q[0];
sx q[0];
rz(-1.9159303) q[0];
sx q[0];
rz(2.0563828) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(1.6403713) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1000993) q[0];
sx q[0];
rz(-1.3736808) q[0];
sx q[0];
rz(-1.6105349) q[0];
rz(-pi) q[1];
rz(2.9502421) q[2];
sx q[2];
rz(-1.3406959) q[2];
sx q[2];
rz(1.5887345) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2861917) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(1.752125) q[1];
rz(-pi) q[2];
rz(-0.94475475) q[3];
sx q[3];
rz(-2.8766655) q[3];
sx q[3];
rz(-1.0696328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7368855) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(0.85477465) q[2];
rz(-0.9084304) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(1.1165718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
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
rz(-1.6664956) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(-1.1845142) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9036983) q[0];
sx q[0];
rz(-1.333295) q[0];
sx q[0];
rz(-2.4303275) q[0];
rz(-2.938561) q[2];
sx q[2];
rz(-0.76329279) q[2];
sx q[2];
rz(-2.610746) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.44691753) q[1];
sx q[1];
rz(-1.0694587) q[1];
sx q[1];
rz(-2.8828709) q[1];
x q[2];
rz(-1.3877901) q[3];
sx q[3];
rz(-1.3015476) q[3];
sx q[3];
rz(-1.8257917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.92675942) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(1.8801749) q[2];
rz(-0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(2.6692218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686907) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(-0.029408971) q[0];
rz(2.3853761) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(-0.75278935) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8653976) q[0];
sx q[0];
rz(-1.5711693) q[0];
sx q[0];
rz(-1.5651817) q[0];
rz(-0.090172099) q[2];
sx q[2];
rz(-1.3036696) q[2];
sx q[2];
rz(1.4324485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.237287) q[1];
sx q[1];
rz(-2.5616966) q[1];
sx q[1];
rz(2.94983) q[1];
rz(2.4642508) q[3];
sx q[3];
rz(-0.65653446) q[3];
sx q[3];
rz(-1.7650106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0011562) q[2];
sx q[2];
rz(-1.493528) q[2];
sx q[2];
rz(-0.39503869) q[2];
rz(-2.6664074) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(0.37676677) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3756556) q[0];
sx q[0];
rz(-0.7203311) q[0];
sx q[0];
rz(-0.96250594) q[0];
rz(2.6267701) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(-0.33822507) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33686906) q[0];
sx q[0];
rz(-2.3677808) q[0];
sx q[0];
rz(1.2880727) q[0];
x q[1];
rz(-2.3338396) q[2];
sx q[2];
rz(-0.33655007) q[2];
sx q[2];
rz(-1.5010264) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64778642) q[1];
sx q[1];
rz(-1.567651) q[1];
sx q[1];
rz(-2.6814744e-05) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.068840222) q[3];
sx q[3];
rz(-1.9919812) q[3];
sx q[3];
rz(2.0023458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.50011355) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(-2.5872453) q[2];
rz(0.85136271) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(0.62801492) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8165269) q[0];
sx q[0];
rz(-0.61282235) q[0];
sx q[0];
rz(-2.391173) q[0];
rz(2.5665307) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(-0.65779984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0938134) q[0];
sx q[0];
rz(-1.1623315) q[0];
sx q[0];
rz(-2.4373217) q[0];
x q[1];
rz(0.051498895) q[2];
sx q[2];
rz(-2.4337075) q[2];
sx q[2];
rz(2.5619363) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7791351) q[1];
sx q[1];
rz(-1.5874079) q[1];
sx q[1];
rz(-1.8569059) q[1];
rz(-pi) q[2];
rz(-2.0714893) q[3];
sx q[3];
rz(-2.7127531) q[3];
sx q[3];
rz(-0.86978149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.647992) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(-1.6592525) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487314) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(0.29801512) q[0];
rz(1.0385849) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-2.9878152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9382291) q[0];
sx q[0];
rz(-1.8244317) q[0];
sx q[0];
rz(-1.4934191) q[0];
rz(-pi) q[1];
rz(0.80259097) q[2];
sx q[2];
rz(-1.0804515) q[2];
sx q[2];
rz(-0.56832321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.01659) q[1];
sx q[1];
rz(-1.224813) q[1];
sx q[1];
rz(-1.6660652) q[1];
rz(-pi) q[2];
rz(-1.3078717) q[3];
sx q[3];
rz(-2.701303) q[3];
sx q[3];
rz(2.0844946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(1.1886965) q[2];
rz(-1.8111604) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(0.73330283) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34340149) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(-0.54939735) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(0.25064358) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3651705) q[0];
sx q[0];
rz(-0.93342268) q[0];
sx q[0];
rz(-0.38743069) q[0];
x q[1];
rz(1.556384) q[2];
sx q[2];
rz(-1.3098798) q[2];
sx q[2];
rz(-0.018785611) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1085775) q[1];
sx q[1];
rz(-2.3653154) q[1];
sx q[1];
rz(3.0423711) q[1];
rz(-pi) q[2];
rz(0.19480494) q[3];
sx q[3];
rz(-1.6945632) q[3];
sx q[3];
rz(-0.74461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-0.13272186) q[2];
sx q[2];
rz(-1.0164725) q[2];
rz(-0.086183444) q[3];
sx q[3];
rz(-2.1250171) q[3];
sx q[3];
rz(2.3763954) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8330399) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(-0.41900751) q[0];
rz(-0.53681701) q[1];
sx q[1];
rz(-2.5173126) q[1];
sx q[1];
rz(-0.73582617) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9402855) q[0];
sx q[0];
rz(-2.2129411) q[0];
sx q[0];
rz(-1.3298195) q[0];
x q[1];
rz(-1.0334098) q[2];
sx q[2];
rz(-2.1436286) q[2];
sx q[2];
rz(-0.17102851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.23259737) q[1];
sx q[1];
rz(-2.4687662) q[1];
sx q[1];
rz(0.35079591) q[1];
x q[2];
rz(2.6594072) q[3];
sx q[3];
rz(-1.8819359) q[3];
sx q[3];
rz(-2.7689214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95911038) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(-2.518892) q[2];
rz(-2.4032118) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(-0.27899376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644792) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(2.4979757) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(1.9401445) q[2];
sx q[2];
rz(-1.7864173) q[2];
sx q[2];
rz(-1.3997072) q[2];
rz(-2.6162061) q[3];
sx q[3];
rz(-2.864884) q[3];
sx q[3];
rz(-1.3340193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
