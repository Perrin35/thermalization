OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3747342) q[0];
sx q[0];
rz(-2.2280966) q[0];
sx q[0];
rz(-1.7537533) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(-1.8228276) q[1];
sx q[1];
rz(0.28611046) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494516) q[0];
sx q[0];
rz(-1.1104662) q[0];
sx q[0];
rz(-1.597531) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.181796) q[2];
sx q[2];
rz(-1.5033443) q[2];
sx q[2];
rz(-0.59997257) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1267136) q[1];
sx q[1];
rz(-0.61903799) q[1];
sx q[1];
rz(0.45780413) q[1];
x q[2];
rz(2.0116352) q[3];
sx q[3];
rz(-1.6206656) q[3];
sx q[3];
rz(0.7383371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.046595786) q[2];
sx q[2];
rz(-1.6511714) q[2];
sx q[2];
rz(0.68520927) q[2];
rz(-1.6738711) q[3];
sx q[3];
rz(-2.0592212) q[3];
sx q[3];
rz(-2.8732324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0005223) q[0];
sx q[0];
rz(-1.5474316) q[0];
sx q[0];
rz(-2.1970314) q[0];
rz(0.073307723) q[1];
sx q[1];
rz(-2.3088375) q[1];
sx q[1];
rz(1.2578957) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2908173) q[0];
sx q[0];
rz(-2.3888139) q[0];
sx q[0];
rz(-0.87429177) q[0];
rz(-pi) q[1];
rz(-2.804685) q[2];
sx q[2];
rz(-2.5818129) q[2];
sx q[2];
rz(0.98613662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.98698178) q[1];
sx q[1];
rz(-1.4263337) q[1];
sx q[1];
rz(2.7491436) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1044534) q[3];
sx q[3];
rz(-1.6172774) q[3];
sx q[3];
rz(2.6543736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8153317) q[2];
sx q[2];
rz(-1.7504642) q[2];
sx q[2];
rz(-1.0247256) q[2];
rz(0.68638408) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(-0.91275233) q[3];
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
rz(-0.82655418) q[0];
sx q[0];
rz(-1.8678366) q[0];
sx q[0];
rz(2.3231373) q[0];
rz(-1.9326899) q[1];
sx q[1];
rz(-1.6345638) q[1];
sx q[1];
rz(0.85982972) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6085109) q[0];
sx q[0];
rz(-1.2405335) q[0];
sx q[0];
rz(-0.94655605) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71975885) q[2];
sx q[2];
rz(-1.9406332) q[2];
sx q[2];
rz(-2.134166) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91539327) q[1];
sx q[1];
rz(-2.2186154) q[1];
sx q[1];
rz(-2.1658705) q[1];
rz(-pi) q[2];
rz(0.1158751) q[3];
sx q[3];
rz(-1.795035) q[3];
sx q[3];
rz(1.7622926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7439421) q[2];
sx q[2];
rz(-2.8068779) q[2];
sx q[2];
rz(-3.0261377) q[2];
rz(2.7502821) q[3];
sx q[3];
rz(-2.1152928) q[3];
sx q[3];
rz(1.1439884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7855969) q[0];
sx q[0];
rz(-1.9556671) q[0];
sx q[0];
rz(2.9010229) q[0];
rz(1.3661512) q[1];
sx q[1];
rz(-1.4931449) q[1];
sx q[1];
rz(-1.2496525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6281462) q[0];
sx q[0];
rz(-1.3168939) q[0];
sx q[0];
rz(0.91350072) q[0];
rz(-2.0555105) q[2];
sx q[2];
rz(-2.4668985) q[2];
sx q[2];
rz(-2.7759068) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0084997) q[1];
sx q[1];
rz(-1.805882) q[1];
sx q[1];
rz(-2.4753768) q[1];
rz(-pi) q[2];
rz(1.0650738) q[3];
sx q[3];
rz(-0.95129644) q[3];
sx q[3];
rz(0.92403417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.25632855) q[2];
sx q[2];
rz(-0.54577959) q[2];
sx q[2];
rz(1.3679999) q[2];
rz(-0.3041501) q[3];
sx q[3];
rz(-1.5658028) q[3];
sx q[3];
rz(-2.4450541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2279385) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(-0.53713334) q[0];
rz(0.74983239) q[1];
sx q[1];
rz(-2.0593819) q[1];
sx q[1];
rz(0.73712635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43057399) q[0];
sx q[0];
rz(-1.6665742) q[0];
sx q[0];
rz(-1.6349313) q[0];
rz(-pi) q[1];
rz(3.0546435) q[2];
sx q[2];
rz(-1.3965522) q[2];
sx q[2];
rz(-0.50424313) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3143718) q[1];
sx q[1];
rz(-1.4431568) q[1];
sx q[1];
rz(-2.8807544) q[1];
rz(-pi) q[2];
rz(-0.32325611) q[3];
sx q[3];
rz(-1.4833602) q[3];
sx q[3];
rz(0.78089784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4445112) q[2];
sx q[2];
rz(-0.35327521) q[2];
sx q[2];
rz(2.9949761) q[2];
rz(-1.4491436) q[3];
sx q[3];
rz(-1.2796947) q[3];
sx q[3];
rz(2.6176738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61831063) q[0];
sx q[0];
rz(-2.1250516) q[0];
sx q[0];
rz(1.4215533) q[0];
rz(-1.8824185) q[1];
sx q[1];
rz(-2.4651395) q[1];
sx q[1];
rz(0.4037942) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.161927) q[0];
sx q[0];
rz(-2.296519) q[0];
sx q[0];
rz(0.37958522) q[0];
rz(-pi) q[1];
rz(0.61876671) q[2];
sx q[2];
rz(-2.8293489) q[2];
sx q[2];
rz(1.1340326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67298573) q[1];
sx q[1];
rz(-2.1588209) q[1];
sx q[1];
rz(2.468416) q[1];
x q[2];
rz(-0.29098068) q[3];
sx q[3];
rz(-0.70937362) q[3];
sx q[3];
rz(-3.0090699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.28973618) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(0.16656052) q[2];
rz(-1.5038495) q[3];
sx q[3];
rz(-2.6681191) q[3];
sx q[3];
rz(-3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2378167) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(-0.87752262) q[0];
rz(0.96218836) q[1];
sx q[1];
rz(-1.1233556) q[1];
sx q[1];
rz(2.7872564) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0838881) q[0];
sx q[0];
rz(-1.5885) q[0];
sx q[0];
rz(2.903865) q[0];
rz(-1.4117175) q[2];
sx q[2];
rz(-0.14661486) q[2];
sx q[2];
rz(2.7377812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-7*pi/15) q[1];
sx q[1];
rz(-1.0364658) q[1];
sx q[1];
rz(-2.2855845) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50994344) q[3];
sx q[3];
rz(-1.6367607) q[3];
sx q[3];
rz(-0.13291453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3366036) q[2];
sx q[2];
rz(-0.73770928) q[2];
sx q[2];
rz(2.0965516) q[2];
rz(-2.7392144) q[3];
sx q[3];
rz(-1.9741524) q[3];
sx q[3];
rz(-1.6654525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43110111) q[0];
sx q[0];
rz(-0.11678188) q[0];
sx q[0];
rz(2.2905599) q[0];
rz(2.7878413) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(-0.85465777) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25682043) q[0];
sx q[0];
rz(-1.8516292) q[0];
sx q[0];
rz(-1.7600842) q[0];
x q[1];
rz(0.24979892) q[2];
sx q[2];
rz(-1.7160176) q[2];
sx q[2];
rz(1.2244566) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6254297) q[1];
sx q[1];
rz(-0.99402797) q[1];
sx q[1];
rz(1.7712084) q[1];
x q[2];
rz(0.40290101) q[3];
sx q[3];
rz(-1.4960999) q[3];
sx q[3];
rz(-0.84045974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9510368) q[2];
sx q[2];
rz(-2.4943116) q[2];
sx q[2];
rz(-0.62038842) q[2];
rz(-2.9151211) q[3];
sx q[3];
rz(-1.3969996) q[3];
sx q[3];
rz(-2.5605719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9613551) q[0];
sx q[0];
rz(-1.2018452) q[0];
sx q[0];
rz(-3.1264547) q[0];
rz(-1.0221647) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(1.9415564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7040492) q[0];
sx q[0];
rz(-1.747073) q[0];
sx q[0];
rz(2.2793819) q[0];
x q[1];
rz(2.5483918) q[2];
sx q[2];
rz(-1.5216547) q[2];
sx q[2];
rz(2.1368487) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.697534) q[1];
sx q[1];
rz(-1.6007533) q[1];
sx q[1];
rz(-1.0590963) q[1];
rz(-2.0807812) q[3];
sx q[3];
rz(-2.0037162) q[3];
sx q[3];
rz(2.1149389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6156561) q[2];
sx q[2];
rz(-0.91292149) q[2];
sx q[2];
rz(0.23942648) q[2];
rz(-0.058852363) q[3];
sx q[3];
rz(-0.81164304) q[3];
sx q[3];
rz(-0.27590251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9489768) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(-2.1666727) q[0];
rz(0.67543593) q[1];
sx q[1];
rz(-2.2223739) q[1];
sx q[1];
rz(0.98178896) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47768053) q[0];
sx q[0];
rz(-0.64723158) q[0];
sx q[0];
rz(-2.0114698) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5025595) q[2];
sx q[2];
rz(-1.9815129) q[2];
sx q[2];
rz(1.0194743) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0412933) q[1];
sx q[1];
rz(-2.2791692) q[1];
sx q[1];
rz(-2.7977976) q[1];
rz(-3.1038935) q[3];
sx q[3];
rz(-2.0172098) q[3];
sx q[3];
rz(2.962449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56741095) q[2];
sx q[2];
rz(-0.97102204) q[2];
sx q[2];
rz(1.9166454) q[2];
rz(2.775906) q[3];
sx q[3];
rz(-0.94958011) q[3];
sx q[3];
rz(2.001781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-3.0443307) q[0];
sx q[0];
rz(-1.2428357) q[0];
sx q[0];
rz(-1.6999929) q[0];
rz(0.080009566) q[1];
sx q[1];
rz(-1.3792104) q[1];
sx q[1];
rz(3.1389799) q[1];
rz(1.710113) q[2];
sx q[2];
rz(-1.1291383) q[2];
sx q[2];
rz(0.51869803) q[2];
rz(1.5479709) q[3];
sx q[3];
rz(-1.1183839) q[3];
sx q[3];
rz(-0.090261685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
