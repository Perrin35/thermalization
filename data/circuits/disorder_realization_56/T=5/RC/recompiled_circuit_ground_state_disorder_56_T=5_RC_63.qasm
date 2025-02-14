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
rz(-2.0593491) q[1];
sx q[1];
rz(-1.6545656) q[1];
sx q[1];
rz(2.3496871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740252) q[0];
sx q[0];
rz(-1.8301395) q[0];
sx q[0];
rz(-1.7631084) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1047529) q[2];
sx q[2];
rz(-1.7245777) q[2];
sx q[2];
rz(-0.80006727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79719964) q[1];
sx q[1];
rz(-2.5973592) q[1];
sx q[1];
rz(-3.0809156) q[1];
rz(-pi) q[2];
rz(0.57609419) q[3];
sx q[3];
rz(-2.2597426) q[3];
sx q[3];
rz(-0.16932936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5504494) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(3.0421416) q[2];
rz(0.24672306) q[3];
sx q[3];
rz(-1.6171425) q[3];
sx q[3];
rz(-0.39768404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196911) q[0];
sx q[0];
rz(-2.1653403) q[0];
sx q[0];
rz(-2.2130527) q[0];
rz(1.6909201) q[1];
sx q[1];
rz(-1.848315) q[1];
sx q[1];
rz(2.4931152) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47005338) q[0];
sx q[0];
rz(-0.82255615) q[0];
sx q[0];
rz(0.7732735) q[0];
rz(-pi) q[1];
rz(2.6669152) q[2];
sx q[2];
rz(-2.5203278) q[2];
sx q[2];
rz(2.9062627) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26095054) q[1];
sx q[1];
rz(-1.7127258) q[1];
sx q[1];
rz(1.5550343) q[1];
rz(-1.2821372) q[3];
sx q[3];
rz(-1.6842173) q[3];
sx q[3];
rz(2.3577549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4497946) q[2];
sx q[2];
rz(-2.3851676) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6344675) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(-1.0852098) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.4986821) q[1];
sx q[1];
rz(1.5012213) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47848343) q[0];
sx q[0];
rz(-1.609765) q[0];
sx q[0];
rz(-2.9443254) q[0];
x q[1];
rz(-0.19135059) q[2];
sx q[2];
rz(-1.3406959) q[2];
sx q[2];
rz(1.5887345) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13255626) q[1];
sx q[1];
rz(-1.6699797) q[1];
sx q[1];
rz(0.99708359) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15764938) q[3];
sx q[3];
rz(-1.356989) q[3];
sx q[3];
rz(2.714954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.40470716) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(2.2331623) q[3];
sx q[3];
rz(-1.5769438) q[3];
sx q[3];
rz(1.1165718) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.214355) q[0];
sx q[0];
rz(-0.4011918) q[0];
sx q[0];
rz(1.1965363) q[0];
rz(1.475097) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(1.9570785) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0663721) q[0];
sx q[0];
rz(-2.3983404) q[0];
sx q[0];
rz(2.7864898) q[0];
x q[1];
rz(-2.938561) q[2];
sx q[2];
rz(-0.76329279) q[2];
sx q[2];
rz(-2.610746) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2503916) q[1];
sx q[1];
rz(-1.7970835) q[1];
sx q[1];
rz(-2.0865284) q[1];
x q[2];
rz(-2.5585737) q[3];
sx q[3];
rz(-0.32430092) q[3];
sx q[3];
rz(-0.7079269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92675942) q[2];
sx q[2];
rz(-0.49761179) q[2];
sx q[2];
rz(1.8801749) q[2];
rz(0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(0.47237083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686907) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(0.029408971) q[0];
rz(0.75621653) q[1];
sx q[1];
rz(-2.5543946) q[1];
sx q[1];
rz(2.3888033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469893) q[0];
sx q[0];
rz(-1.5764109) q[0];
sx q[0];
rz(0.0003729781) q[0];
x q[1];
rz(-1.3026313) q[2];
sx q[2];
rz(-1.4838306) q[2];
sx q[2];
rz(3.027107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.636105) q[1];
sx q[1];
rz(-1.4661745) q[1];
sx q[1];
rz(2.5701447) q[1];
x q[2];
rz(-2.6008368) q[3];
sx q[3];
rz(-1.9633368) q[3];
sx q[3];
rz(-2.768571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14043643) q[2];
sx q[2];
rz(-1.6480646) q[2];
sx q[2];
rz(2.746554) q[2];
rz(0.47518528) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(0.37676677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659371) q[0];
sx q[0];
rz(-2.4212615) q[0];
sx q[0];
rz(2.1790867) q[0];
rz(-2.6267701) q[1];
sx q[1];
rz(-1.5341026) q[1];
sx q[1];
rz(-2.8033676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7027959) q[0];
sx q[0];
rz(-1.7670175) q[0];
sx q[0];
rz(2.3243734) q[0];
x q[1];
rz(2.3338396) q[2];
sx q[2];
rz(-0.33655007) q[2];
sx q[2];
rz(-1.6405662) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5023313) q[1];
sx q[1];
rz(-0.0031454589) q[1];
sx q[1];
rz(-1.5622713) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0727524) q[3];
sx q[3];
rz(-1.9919812) q[3];
sx q[3];
rz(-2.0023458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.50011355) q[2];
sx q[2];
rz(-1.1621472) q[2];
sx q[2];
rz(0.55434736) q[2];
rz(-0.85136271) q[3];
sx q[3];
rz(-2.7832289) q[3];
sx q[3];
rz(2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3250658) q[0];
sx q[0];
rz(-2.5287703) q[0];
sx q[0];
rz(-0.75041962) q[0];
rz(-2.5665307) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(0.65779984) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039666273) q[0];
sx q[0];
rz(-0.79625477) q[0];
sx q[0];
rz(-2.5523561) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6148241) q[2];
sx q[2];
rz(-2.2775473) q[2];
sx q[2];
rz(2.6296774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7791351) q[1];
sx q[1];
rz(-1.5874079) q[1];
sx q[1];
rz(1.8569059) q[1];
x q[2];
rz(-1.0701033) q[3];
sx q[3];
rz(-0.4288396) q[3];
sx q[3];
rz(-0.86978149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.647992) q[2];
sx q[2];
rz(-1.0085663) q[2];
sx q[2];
rz(-1.4823401) q[2];
rz(-3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(0.83546662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7928612) q[0];
sx q[0];
rz(-0.86935765) q[0];
sx q[0];
rz(2.8435775) q[0];
rz(1.0385849) q[1];
sx q[1];
rz(-0.54439259) q[1];
sx q[1];
rz(2.9878152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9382291) q[0];
sx q[0];
rz(-1.8244317) q[0];
sx q[0];
rz(1.6481736) q[0];
rz(-pi) q[1];
rz(-0.63855448) q[2];
sx q[2];
rz(-2.2306135) q[2];
sx q[2];
rz(1.4294848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3996622) q[1];
sx q[1];
rz(-2.7832418) q[1];
sx q[1];
rz(-2.8835924) q[1];
rz(-pi) q[2];
rz(0.12184398) q[3];
sx q[3];
rz(-1.9949556) q[3];
sx q[3];
rz(0.76790732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(-1.1886965) q[2];
rz(-1.3304322) q[3];
sx q[3];
rz(-1.9537787) q[3];
sx q[3];
rz(2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.7981912) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.6424302) q[0];
rz(0.54939735) q[1];
sx q[1];
rz(-1.5645942) q[1];
sx q[1];
rz(-2.8909491) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77642216) q[0];
sx q[0];
rz(-2.20817) q[0];
sx q[0];
rz(-2.754162) q[0];
x q[1];
rz(1.5852087) q[2];
sx q[2];
rz(-1.3098798) q[2];
sx q[2];
rz(-3.122807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10560606) q[1];
sx q[1];
rz(-0.79933724) q[1];
sx q[1];
rz(1.6677594) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57119675) q[3];
sx q[3];
rz(-0.2303752) q[3];
sx q[3];
rz(-2.8744968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9833019) q[2];
sx q[2];
rz(-0.13272186) q[2];
sx q[2];
rz(2.1251202) q[2];
rz(3.0554092) q[3];
sx q[3];
rz(-2.1250171) q[3];
sx q[3];
rz(-0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30855274) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(-2.7225851) q[0];
rz(-2.6047756) q[1];
sx q[1];
rz(-2.5173126) q[1];
sx q[1];
rz(0.73582617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48437546) q[0];
sx q[0];
rz(-1.3784982) q[0];
sx q[0];
rz(-2.4853287) q[0];
rz(0.64401099) q[2];
sx q[2];
rz(-1.1261054) q[2];
sx q[2];
rz(-1.4294238) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.23259737) q[1];
sx q[1];
rz(-0.67282644) q[1];
sx q[1];
rz(-0.35079591) q[1];
x q[2];
rz(-0.60634585) q[3];
sx q[3];
rz(-2.5744573) q[3];
sx q[3];
rz(2.4727269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95911038) q[2];
sx q[2];
rz(-0.82254326) q[2];
sx q[2];
rz(2.518892) q[2];
rz(-2.4032118) q[3];
sx q[3];
rz(-1.9564956) q[3];
sx q[3];
rz(2.8625989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5644792) q[0];
sx q[0];
rz(-0.27272419) q[0];
sx q[0];
rz(0.96584366) q[0];
rz(-0.64361698) q[1];
sx q[1];
rz(-1.2871965) q[1];
sx q[1];
rz(1.0866477) q[1];
rz(0.23068354) q[2];
sx q[2];
rz(-1.2103969) q[2];
sx q[2];
rz(-2.8878676) q[2];
rz(0.52538659) q[3];
sx q[3];
rz(-2.864884) q[3];
sx q[3];
rz(-1.3340193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
