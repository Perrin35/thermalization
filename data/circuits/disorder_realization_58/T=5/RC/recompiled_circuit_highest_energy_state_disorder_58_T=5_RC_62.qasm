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
rz(-2.8623924) q[0];
sx q[0];
rz(-1.8912127) q[0];
sx q[0];
rz(-3.1412636) q[0];
rz(-2.012913) q[1];
sx q[1];
rz(-2.0560775) q[1];
sx q[1];
rz(-2.9323199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1827917) q[0];
sx q[0];
rz(-1.3404935) q[0];
sx q[0];
rz(0.15887375) q[0];
rz(-pi) q[1];
rz(-2.5821788) q[2];
sx q[2];
rz(-1.4296247) q[2];
sx q[2];
rz(2.5641297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8289703) q[1];
sx q[1];
rz(-0.43401845) q[1];
sx q[1];
rz(0.75059469) q[1];
x q[2];
rz(-2.3429782) q[3];
sx q[3];
rz(-0.3363401) q[3];
sx q[3];
rz(2.596465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2613232) q[2];
sx q[2];
rz(-1.7797194) q[2];
sx q[2];
rz(-2.9051991) q[2];
rz(-2.7586625) q[3];
sx q[3];
rz(-1.087944) q[3];
sx q[3];
rz(1.6898164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7616538) q[0];
sx q[0];
rz(-0.28872696) q[0];
sx q[0];
rz(-2.5987103) q[0];
rz(0.44034964) q[1];
sx q[1];
rz(-1.1232702) q[1];
sx q[1];
rz(-1.1276668) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.544761) q[0];
sx q[0];
rz(-2.4003599) q[0];
sx q[0];
rz(2.6165589) q[0];
x q[1];
rz(0.93483804) q[2];
sx q[2];
rz(-1.4439772) q[2];
sx q[2];
rz(-0.93709125) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7372514) q[1];
sx q[1];
rz(-1.4228428) q[1];
sx q[1];
rz(0.80929359) q[1];
rz(-2.7064467) q[3];
sx q[3];
rz(-0.75772396) q[3];
sx q[3];
rz(2.0592281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.27309624) q[2];
sx q[2];
rz(-0.79600969) q[2];
sx q[2];
rz(1.2924755) q[2];
rz(-0.83186197) q[3];
sx q[3];
rz(-1.4757194) q[3];
sx q[3];
rz(-2.7003435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.454994) q[0];
sx q[0];
rz(-0.43231493) q[0];
sx q[0];
rz(-1.5469714) q[0];
rz(3.1349685) q[1];
sx q[1];
rz(-2.6057656) q[1];
sx q[1];
rz(2.5033902) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5247045) q[0];
sx q[0];
rz(-1.0959488) q[0];
sx q[0];
rz(3.1010702) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9784967) q[2];
sx q[2];
rz(-0.8434274) q[2];
sx q[2];
rz(-0.076956017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6276826) q[1];
sx q[1];
rz(-1.2177151) q[1];
sx q[1];
rz(-0.38796918) q[1];
rz(-pi) q[2];
rz(-1.9157103) q[3];
sx q[3];
rz(-0.57424358) q[3];
sx q[3];
rz(-1.4019897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39614761) q[2];
sx q[2];
rz(-0.026693176) q[2];
sx q[2];
rz(2.2849582) q[2];
rz(0.73532295) q[3];
sx q[3];
rz(-1.7263128) q[3];
sx q[3];
rz(1.0480688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.624144) q[0];
sx q[0];
rz(-3.090591) q[0];
sx q[0];
rz(-2.2080102) q[0];
rz(2.3816662) q[1];
sx q[1];
rz(-0.81587195) q[1];
sx q[1];
rz(-1.5911721) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.099898) q[0];
sx q[0];
rz(-1.2907488) q[0];
sx q[0];
rz(1.6629908) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7047005) q[2];
sx q[2];
rz(-1.5073188) q[2];
sx q[2];
rz(-1.8885828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62867579) q[1];
sx q[1];
rz(-2.4077291) q[1];
sx q[1];
rz(-0.34725125) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24985866) q[3];
sx q[3];
rz(-1.7267531) q[3];
sx q[3];
rz(0.74710959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7237225) q[2];
sx q[2];
rz(-1.7524065) q[2];
sx q[2];
rz(1.0327445) q[2];
rz(-0.51747733) q[3];
sx q[3];
rz(-1.7101219) q[3];
sx q[3];
rz(-1.8364068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0378157) q[0];
sx q[0];
rz(-0.10194898) q[0];
sx q[0];
rz(-1.284449) q[0];
rz(-0.85353723) q[1];
sx q[1];
rz(-1.6910911) q[1];
sx q[1];
rz(2.6654713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7126851) q[0];
sx q[0];
rz(-1.0700912) q[0];
sx q[0];
rz(2.7963573) q[0];
rz(-2.4867441) q[2];
sx q[2];
rz(-2.691948) q[2];
sx q[2];
rz(-2.8596535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.13083982) q[1];
sx q[1];
rz(-0.71341842) q[1];
sx q[1];
rz(-2.6793007) q[1];
x q[2];
rz(2.3105442) q[3];
sx q[3];
rz(-1.5171192) q[3];
sx q[3];
rz(0.37572467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4173172) q[2];
sx q[2];
rz(-0.83160916) q[2];
sx q[2];
rz(1.6012491) q[2];
rz(-3.0658412) q[3];
sx q[3];
rz(-1.3906393) q[3];
sx q[3];
rz(-2.4116404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041357111) q[0];
sx q[0];
rz(-0.55158177) q[0];
sx q[0];
rz(1.7919737) q[0];
rz(1.8655818) q[1];
sx q[1];
rz(-0.7834692) q[1];
sx q[1];
rz(-2.1003726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43364701) q[0];
sx q[0];
rz(-0.7793847) q[0];
sx q[0];
rz(2.8288366) q[0];
x q[1];
rz(-1.5346933) q[2];
sx q[2];
rz(-2.3182333) q[2];
sx q[2];
rz(2.2359816) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1231844) q[1];
sx q[1];
rz(-2.6861172) q[1];
sx q[1];
rz(-2.5405618) q[1];
rz(2.842343) q[3];
sx q[3];
rz(-0.78093442) q[3];
sx q[3];
rz(2.2921103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.87018806) q[2];
sx q[2];
rz(-1.9726334) q[2];
sx q[2];
rz(-0.99883396) q[2];
rz(-0.53423229) q[3];
sx q[3];
rz(-1.9637008) q[3];
sx q[3];
rz(1.4786725) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6967195) q[0];
sx q[0];
rz(-0.98170009) q[0];
sx q[0];
rz(1.556742) q[0];
rz(2.059767) q[1];
sx q[1];
rz(-1.4993246) q[1];
sx q[1];
rz(2.7781442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7666192) q[0];
sx q[0];
rz(-0.21265499) q[0];
sx q[0];
rz(0.029092832) q[0];
x q[1];
rz(0.29691445) q[2];
sx q[2];
rz(-1.468942) q[2];
sx q[2];
rz(2.7114026) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.195943) q[1];
sx q[1];
rz(-0.87313014) q[1];
sx q[1];
rz(-2.5673668) q[1];
x q[2];
rz(0.15017516) q[3];
sx q[3];
rz(-1.8446823) q[3];
sx q[3];
rz(1.3824579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5536559) q[2];
sx q[2];
rz(-0.16599545) q[2];
sx q[2];
rz(2.0641649) q[2];
rz(-1.59498) q[3];
sx q[3];
rz(-1.8190106) q[3];
sx q[3];
rz(2.4685629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308291) q[0];
sx q[0];
rz(-1.6442693) q[0];
sx q[0];
rz(2.8610435) q[0];
rz(0.55889091) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(1.9299318) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.856177) q[0];
sx q[0];
rz(-0.95024419) q[0];
sx q[0];
rz(0.65914777) q[0];
x q[1];
rz(1.3468387) q[2];
sx q[2];
rz(-1.0887165) q[2];
sx q[2];
rz(0.48165762) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.385541) q[1];
sx q[1];
rz(-1.5935399) q[1];
sx q[1];
rz(2.3616516) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5959963) q[3];
sx q[3];
rz(-0.19807252) q[3];
sx q[3];
rz(-0.36626178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1739379) q[2];
sx q[2];
rz(-1.329198) q[2];
sx q[2];
rz(-0.8052513) q[2];
rz(-1.1236745) q[3];
sx q[3];
rz(-1.1403964) q[3];
sx q[3];
rz(0.32127389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6174018) q[0];
sx q[0];
rz(-2.5179355) q[0];
sx q[0];
rz(0.50530857) q[0];
rz(0.71120039) q[1];
sx q[1];
rz(-0.86831793) q[1];
sx q[1];
rz(-1.1127081) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.959247) q[0];
sx q[0];
rz(-1.2768978) q[0];
sx q[0];
rz(2.6698851) q[0];
rz(-pi) q[1];
rz(-1.4897025) q[2];
sx q[2];
rz(-1.4633274) q[2];
sx q[2];
rz(-1.7484959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1755274) q[1];
sx q[1];
rz(-1.4443782) q[1];
sx q[1];
rz(-1.1299707) q[1];
rz(-pi) q[2];
rz(1.0151583) q[3];
sx q[3];
rz(-1.5446413) q[3];
sx q[3];
rz(-1.2873611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6342371) q[2];
sx q[2];
rz(-1.2398182) q[2];
sx q[2];
rz(-0.235454) q[2];
rz(-2.7638451) q[3];
sx q[3];
rz(-2.7538959) q[3];
sx q[3];
rz(-3.0915251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78558755) q[0];
sx q[0];
rz(-1.8536785) q[0];
sx q[0];
rz(0.054314286) q[0];
rz(2.3771225) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(2.4321709) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13693181) q[0];
sx q[0];
rz(-0.2668902) q[0];
sx q[0];
rz(-1.1679807) q[0];
x q[1];
rz(-0.26086678) q[2];
sx q[2];
rz(-1.4402908) q[2];
sx q[2];
rz(-0.57756391) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5125138) q[1];
sx q[1];
rz(-0.98328749) q[1];
sx q[1];
rz(1.4329989) q[1];
rz(-pi) q[2];
rz(-0.34950236) q[3];
sx q[3];
rz(-2.7785578) q[3];
sx q[3];
rz(0.078840794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8556122) q[2];
sx q[2];
rz(-1.4748272) q[2];
sx q[2];
rz(-0.093185514) q[2];
rz(0.5992254) q[3];
sx q[3];
rz(-2.5803876) q[3];
sx q[3];
rz(0.78527251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4335321) q[0];
sx q[0];
rz(-2.1716433) q[0];
sx q[0];
rz(-0.94950983) q[0];
rz(-0.064619725) q[1];
sx q[1];
rz(-3.0977991) q[1];
sx q[1];
rz(-2.480712) q[1];
rz(2.1223161) q[2];
sx q[2];
rz(-1.5073771) q[2];
sx q[2];
rz(1.3310968) q[2];
rz(-0.6653847) q[3];
sx q[3];
rz(-1.4238085) q[3];
sx q[3];
rz(2.4091638) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
