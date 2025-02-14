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
rz(4.7779067) q[0];
sx q[0];
rz(11.783083) q[0];
rz(1.0822436) q[1];
sx q[1];
rz(-1.487027) q[1];
sx q[1];
rz(0.79190555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.22345) q[0];
sx q[0];
rz(-0.32157102) q[0];
sx q[0];
rz(-0.62427036) q[0];
x q[1];
rz(3.1047529) q[2];
sx q[2];
rz(-1.7245777) q[2];
sx q[2];
rz(0.80006727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3160682) q[1];
sx q[1];
rz(-1.5393942) q[1];
sx q[1];
rz(-0.54341787) q[1];
x q[2];
rz(-0.57609419) q[3];
sx q[3];
rz(-0.8818501) q[3];
sx q[3];
rz(-0.16932936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5504494) q[2];
sx q[2];
rz(-2.998896) q[2];
sx q[2];
rz(-0.099451065) q[2];
rz(0.24672306) q[3];
sx q[3];
rz(-1.5244502) q[3];
sx q[3];
rz(0.39768404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0196911) q[0];
sx q[0];
rz(-0.97625232) q[0];
sx q[0];
rz(-2.2130527) q[0];
rz(-1.6909201) q[1];
sx q[1];
rz(-1.2932777) q[1];
sx q[1];
rz(2.4931152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4546616) q[0];
sx q[0];
rz(-1.0333916) q[0];
sx q[0];
rz(0.65673687) q[0];
rz(1.2546059) q[2];
sx q[2];
rz(-1.0266227) q[2];
sx q[2];
rz(2.3426825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8806421) q[1];
sx q[1];
rz(-1.7127258) q[1];
sx q[1];
rz(-1.5865583) q[1];
x q[2];
rz(-0.11827151) q[3];
sx q[3];
rz(-1.2840446) q[3];
sx q[3];
rz(2.3210382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6917981) q[2];
sx q[2];
rz(-2.3851676) q[2];
sx q[2];
rz(0.85255426) q[2];
rz(0.84613386) q[3];
sx q[3];
rz(-1.5087912) q[3];
sx q[3];
rz(-1.030863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5071252) q[0];
sx q[0];
rz(-1.2256624) q[0];
sx q[0];
rz(-1.0852098) q[0];
rz(2.197544) q[1];
sx q[1];
rz(-1.6429106) q[1];
sx q[1];
rz(-1.5012213) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6631092) q[0];
sx q[0];
rz(-1.609765) q[0];
sx q[0];
rz(2.9443254) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8050212) q[2];
sx q[2];
rz(-1.384549) q[2];
sx q[2];
rz(3.0795003) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.855401) q[1];
sx q[1];
rz(-0.58127379) q[1];
sx q[1];
rz(1.752125) q[1];
rz(-0.15764938) q[3];
sx q[3];
rz(-1.7846037) q[3];
sx q[3];
rz(2.714954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7368855) q[2];
sx q[2];
rz(-0.53386226) q[2];
sx q[2];
rz(-0.85477465) q[2];
rz(-0.9084304) q[3];
sx q[3];
rz(-1.5646489) q[3];
sx q[3];
rz(-1.1165718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92723769) q[0];
sx q[0];
rz(-2.7404009) q[0];
sx q[0];
rz(-1.1965363) q[0];
rz(-1.6664956) q[1];
sx q[1];
rz(-0.42088446) q[1];
sx q[1];
rz(-1.1845142) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0663721) q[0];
sx q[0];
rz(-0.74325222) q[0];
sx q[0];
rz(-0.35510285) q[0];
rz(-pi) q[1];
rz(-1.7613715) q[2];
sx q[2];
rz(-2.3146176) q[2];
sx q[2];
rz(-0.25328749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.44691753) q[1];
sx q[1];
rz(-2.072134) q[1];
sx q[1];
rz(0.25872177) q[1];
rz(-pi) q[2];
rz(0.27359815) q[3];
sx q[3];
rz(-1.7471385) q[3];
sx q[3];
rz(2.8374052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2148332) q[2];
sx q[2];
rz(-2.6439809) q[2];
sx q[2];
rz(-1.2614177) q[2];
rz(-0.43241209) q[3];
sx q[3];
rz(-1.132248) q[3];
sx q[3];
rz(-0.47237083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.072902) q[0];
sx q[0];
rz(-1.2304767) q[0];
sx q[0];
rz(-3.1121837) q[0];
rz(2.3853761) q[1];
sx q[1];
rz(-0.58719802) q[1];
sx q[1];
rz(2.3888033) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2761951) q[0];
sx q[0];
rz(-1.5711693) q[0];
sx q[0];
rz(1.5651817) q[0];
rz(-1.8886861) q[2];
sx q[2];
rz(-2.860002) q[2];
sx q[2];
rz(-1.3791305) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90430561) q[1];
sx q[1];
rz(-2.5616966) q[1];
sx q[1];
rz(-0.19176264) q[1];
x q[2];
rz(-2.4642508) q[3];
sx q[3];
rz(-2.4850582) q[3];
sx q[3];
rz(-1.7650106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.14043643) q[2];
sx q[2];
rz(-1.493528) q[2];
sx q[2];
rz(0.39503869) q[2];
rz(0.47518528) q[3];
sx q[3];
rz(-1.3522215) q[3];
sx q[3];
rz(-2.7648259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.7203311) q[0];
sx q[0];
rz(0.96250594) q[0];
rz(0.51482254) q[1];
sx q[1];
rz(-1.6074901) q[1];
sx q[1];
rz(-0.33822507) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.048934919) q[0];
sx q[0];
rz(-0.83507628) q[0];
sx q[0];
rz(0.2661163) q[0];
rz(-pi) q[1];
rz(0.80775308) q[2];
sx q[2];
rz(-2.8050426) q[2];
sx q[2];
rz(1.5010264) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5023313) q[1];
sx q[1];
rz(-3.1384472) q[1];
sx q[1];
rz(1.5622713) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0727524) q[3];
sx q[3];
rz(-1.1496115) q[3];
sx q[3];
rz(-2.0023458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50011355) q[2];
sx q[2];
rz(-1.9794455) q[2];
sx q[2];
rz(2.5872453) q[2];
rz(2.2902299) q[3];
sx q[3];
rz(-0.35836372) q[3];
sx q[3];
rz(-2.5135777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.8165269) q[0];
sx q[0];
rz(-0.61282235) q[0];
sx q[0];
rz(-0.75041962) q[0];
rz(0.57506192) q[1];
sx q[1];
rz(-1.776418) q[1];
sx q[1];
rz(-2.4837928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039666273) q[0];
sx q[0];
rz(-2.3453379) q[0];
sx q[0];
rz(-0.58923652) q[0];
rz(-0.051498895) q[2];
sx q[2];
rz(-2.4337075) q[2];
sx q[2];
rz(-2.5619363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.15193394) q[1];
sx q[1];
rz(-0.2865782) q[1];
sx q[1];
rz(-1.6295939) q[1];
x q[2];
rz(-2.9255387) q[3];
sx q[3];
rz(-1.1974058) q[3];
sx q[3];
rz(1.4114398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.647992) q[2];
sx q[2];
rz(-2.1330264) q[2];
sx q[2];
rz(1.4823401) q[2];
rz(3.0209387) q[3];
sx q[3];
rz(-1.7574661) q[3];
sx q[3];
rz(2.306126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3487314) q[0];
sx q[0];
rz(-2.272235) q[0];
sx q[0];
rz(-2.8435775) q[0];
rz(2.1030078) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(2.9878152) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9382291) q[0];
sx q[0];
rz(-1.3171609) q[0];
sx q[0];
rz(-1.4934191) q[0];
x q[1];
rz(0.63855448) q[2];
sx q[2];
rz(-0.91097915) q[2];
sx q[2];
rz(1.4294848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3996622) q[1];
sx q[1];
rz(-0.35835086) q[1];
sx q[1];
rz(-2.8835924) q[1];
x q[2];
rz(-1.9977536) q[3];
sx q[3];
rz(-1.4597963) q[3];
sx q[3];
rz(-2.3890561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.095857233) q[2];
sx q[2];
rz(-0.46892527) q[2];
sx q[2];
rz(1.9528961) q[2];
rz(1.3304322) q[3];
sx q[3];
rz(-1.1878139) q[3];
sx q[3];
rz(2.4082898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34340149) q[0];
sx q[0];
rz(-0.60438406) q[0];
sx q[0];
rz(-1.4991624) q[0];
rz(0.54939735) q[1];
sx q[1];
rz(-1.5769985) q[1];
sx q[1];
rz(-0.25064358) q[1];
rz(-pi/2) q[2];
sx q[2];
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
rz(-0.26094238) q[2];
sx q[2];
rz(-1.5568718) q[2];
sx q[2];
rz(-1.5557289) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0359866) q[1];
sx q[1];
rz(-2.3422554) q[1];
sx q[1];
rz(1.6677594) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57119675) q[3];
sx q[3];
rz(-2.9112175) q[3];
sx q[3];
rz(-0.26709589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15829076) q[2];
sx q[2];
rz(-0.13272186) q[2];
sx q[2];
rz(1.0164725) q[2];
rz(0.086183444) q[3];
sx q[3];
rz(-1.0165756) q[3];
sx q[3];
rz(-0.76519722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8330399) q[0];
sx q[0];
rz(-2.2873531) q[0];
sx q[0];
rz(2.7225851) q[0];
rz(-2.6047756) q[1];
sx q[1];
rz(-0.62428004) q[1];
sx q[1];
rz(2.4057665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2013071) q[0];
sx q[0];
rz(-2.2129411) q[0];
sx q[0];
rz(1.8117732) q[0];
rz(-2.4707253) q[2];
sx q[2];
rz(-0.76422526) q[2];
sx q[2];
rz(0.66167458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9365471) q[1];
sx q[1];
rz(-2.1959795) q[1];
sx q[1];
rz(-1.3035098) q[1];
rz(-2.6594072) q[3];
sx q[3];
rz(-1.2596568) q[3];
sx q[3];
rz(0.37267123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95911038) q[2];
sx q[2];
rz(-2.3190494) q[2];
sx q[2];
rz(2.518892) q[2];
rz(-2.4032118) q[3];
sx q[3];
rz(-1.1850971) q[3];
sx q[3];
rz(-2.8625989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57711346) q[0];
sx q[0];
rz(-2.8688685) q[0];
sx q[0];
rz(-2.175749) q[0];
rz(0.64361698) q[1];
sx q[1];
rz(-1.8543961) q[1];
sx q[1];
rz(-2.054945) q[1];
rz(1.2014482) q[2];
sx q[2];
rz(-1.3551753) q[2];
sx q[2];
rz(1.7418855) q[2];
rz(0.24091992) q[3];
sx q[3];
rz(-1.708247) q[3];
sx q[3];
rz(0.74549992) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
