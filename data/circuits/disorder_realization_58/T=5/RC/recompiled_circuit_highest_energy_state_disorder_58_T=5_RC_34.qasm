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
rz(0.27920029) q[0];
sx q[0];
rz(-1.2503799) q[0];
sx q[0];
rz(-0.00032902349) q[0];
rz(-2.012913) q[1];
sx q[1];
rz(-2.0560775) q[1];
sx q[1];
rz(0.20927277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64855382) q[0];
sx q[0];
rz(-1.7254412) q[0];
sx q[0];
rz(-1.8039319) q[0];
rz(2.5821788) q[2];
sx q[2];
rz(-1.4296247) q[2];
sx q[2];
rz(-2.5641297) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55587253) q[1];
sx q[1];
rz(-1.2798847) q[1];
sx q[1];
rz(2.8147999) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23926034) q[3];
sx q[3];
rz(-1.3321028) q[3];
sx q[3];
rz(0.25588122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2613232) q[2];
sx q[2];
rz(-1.7797194) q[2];
sx q[2];
rz(-0.23639354) q[2];
rz(0.38293019) q[3];
sx q[3];
rz(-1.087944) q[3];
sx q[3];
rz(-1.4517763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7616538) q[0];
sx q[0];
rz(-2.8528657) q[0];
sx q[0];
rz(-2.5987103) q[0];
rz(-0.44034964) q[1];
sx q[1];
rz(-1.1232702) q[1];
sx q[1];
rz(1.1276668) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.544761) q[0];
sx q[0];
rz(-2.4003599) q[0];
sx q[0];
rz(0.5250338) q[0];
rz(-pi) q[1];
rz(-0.93483804) q[2];
sx q[2];
rz(-1.6976154) q[2];
sx q[2];
rz(2.2045014) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1285557) q[1];
sx q[1];
rz(-0.77289509) q[1];
sx q[1];
rz(1.3580639) q[1];
rz(-0.43514593) q[3];
sx q[3];
rz(-0.75772396) q[3];
sx q[3];
rz(1.0823646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27309624) q[2];
sx q[2];
rz(-0.79600969) q[2];
sx q[2];
rz(1.8491171) q[2];
rz(-2.3097307) q[3];
sx q[3];
rz(-1.6658733) q[3];
sx q[3];
rz(-2.7003435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.454994) q[0];
sx q[0];
rz(-2.7092777) q[0];
sx q[0];
rz(1.5946213) q[0];
rz(3.1349685) q[1];
sx q[1];
rz(-0.5358271) q[1];
sx q[1];
rz(-2.5033902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4362559) q[0];
sx q[0];
rz(-0.47644189) q[0];
sx q[0];
rz(1.492155) q[0];
rz(-0.41902991) q[2];
sx q[2];
rz(-2.3263676) q[2];
sx q[2];
rz(-2.4885675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51391006) q[1];
sx q[1];
rz(-1.9238775) q[1];
sx q[1];
rz(-0.38796918) q[1];
rz(-pi) q[2];
rz(-1.2258824) q[3];
sx q[3];
rz(-0.57424358) q[3];
sx q[3];
rz(1.4019897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.745445) q[2];
sx q[2];
rz(-0.026693176) q[2];
sx q[2];
rz(2.2849582) q[2];
rz(2.4062697) q[3];
sx q[3];
rz(-1.7263128) q[3];
sx q[3];
rz(2.0935238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.624144) q[0];
sx q[0];
rz(-0.051001661) q[0];
sx q[0];
rz(0.93358246) q[0];
rz(-0.75992641) q[1];
sx q[1];
rz(-2.3257207) q[1];
sx q[1];
rz(1.5911721) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5869414) q[0];
sx q[0];
rz(-1.6593895) q[0];
sx q[0];
rz(2.8604126) q[0];
rz(-pi) q[1];
rz(2.9924852) q[2];
sx q[2];
rz(-2.7004074) q[2];
sx q[2];
rz(-0.18277482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67956808) q[1];
sx q[1];
rz(-1.3408519) q[1];
sx q[1];
rz(-2.438222) q[1];
rz(-pi) q[2];
rz(2.5752034) q[3];
sx q[3];
rz(-0.29367125) q[3];
sx q[3];
rz(1.7711414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4178702) q[2];
sx q[2];
rz(-1.3891862) q[2];
sx q[2];
rz(1.0327445) q[2];
rz(2.6241153) q[3];
sx q[3];
rz(-1.7101219) q[3];
sx q[3];
rz(1.3051858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1037769) q[0];
sx q[0];
rz(-3.0396437) q[0];
sx q[0];
rz(-1.8571437) q[0];
rz(-2.2880554) q[1];
sx q[1];
rz(-1.6910911) q[1];
sx q[1];
rz(0.4761214) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697107) q[0];
sx q[0];
rz(-0.5997385) q[0];
sx q[0];
rz(2.1246596) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7760162) q[2];
sx q[2];
rz(-1.8387059) q[2];
sx q[2];
rz(0.6838201) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13083982) q[1];
sx q[1];
rz(-2.4281742) q[1];
sx q[1];
rz(-2.6793007) q[1];
rz(-pi) q[2];
rz(-3.0689803) q[3];
sx q[3];
rz(-2.3092306) q[3];
sx q[3];
rz(1.8975964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1002355) q[0];
sx q[0];
rz(-0.55158177) q[0];
sx q[0];
rz(1.7919737) q[0];
rz(1.8655818) q[1];
sx q[1];
rz(-2.3581235) q[1];
sx q[1];
rz(-1.0412201) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2305377) q[0];
sx q[0];
rz(-1.7887702) q[0];
sx q[0];
rz(-0.75454069) q[0];
rz(-pi) q[1];
rz(3.1026672) q[2];
sx q[2];
rz(-0.74813977) q[2];
sx q[2];
rz(0.85252658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7752616) q[1];
sx q[1];
rz(-1.9420673) q[1];
sx q[1];
rz(-1.300578) q[1];
rz(0.29924969) q[3];
sx q[3];
rz(-2.3606582) q[3];
sx q[3];
rz(-0.84948236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87018806) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(2.1427587) q[2];
rz(2.6073604) q[3];
sx q[3];
rz(-1.1778919) q[3];
sx q[3];
rz(1.6629201) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6967195) q[0];
sx q[0];
rz(-0.98170009) q[0];
sx q[0];
rz(1.5848507) q[0];
rz(-1.0818256) q[1];
sx q[1];
rz(-1.4993246) q[1];
sx q[1];
rz(2.7781442) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666192) q[0];
sx q[0];
rz(-2.9289377) q[0];
sx q[0];
rz(0.029092832) q[0];
x q[1];
rz(-0.29691445) q[2];
sx q[2];
rz(-1.468942) q[2];
sx q[2];
rz(-2.7114026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.195943) q[1];
sx q[1];
rz(-2.2684625) q[1];
sx q[1];
rz(-2.5673668) q[1];
x q[2];
rz(2.9914175) q[3];
sx q[3];
rz(-1.2969103) q[3];
sx q[3];
rz(-1.7591348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5536559) q[2];
sx q[2];
rz(-0.16599545) q[2];
sx q[2];
rz(1.0774277) q[2];
rz(-1.59498) q[3];
sx q[3];
rz(-1.3225821) q[3];
sx q[3];
rz(-2.4685629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107635) q[0];
sx q[0];
rz(-1.6442693) q[0];
sx q[0];
rz(-2.8610435) q[0];
rz(0.55889091) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(-1.2116609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794811) q[0];
sx q[0];
rz(-2.0923847) q[0];
sx q[0];
rz(0.83569367) q[0];
rz(-pi) q[1];
rz(1.3468387) q[2];
sx q[2];
rz(-2.0528762) q[2];
sx q[2];
rz(-0.48165762) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9338464) q[1];
sx q[1];
rz(-2.3504815) q[1];
sx q[1];
rz(1.6027811) q[1];
x q[2];
rz(0.0050571613) q[3];
sx q[3];
rz(-1.3727875) q[3];
sx q[3];
rz(-0.34055948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9676548) q[2];
sx q[2];
rz(-1.329198) q[2];
sx q[2];
rz(-2.3363414) q[2];
rz(2.0179181) q[3];
sx q[3];
rz(-1.1403964) q[3];
sx q[3];
rz(-2.8203188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6174018) q[0];
sx q[0];
rz(-2.5179355) q[0];
sx q[0];
rz(-0.50530857) q[0];
rz(2.4303923) q[1];
sx q[1];
rz(-0.86831793) q[1];
sx q[1];
rz(1.1127081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.959247) q[0];
sx q[0];
rz(-1.8646949) q[0];
sx q[0];
rz(-2.6698851) q[0];
rz(-pi) q[1];
rz(-1.6518902) q[2];
sx q[2];
rz(-1.6782653) q[2];
sx q[2];
rz(-1.7484959) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6869097) q[1];
sx q[1];
rz(-1.1337323) q[1];
sx q[1];
rz(-3.0019762) q[1];
x q[2];
rz(-1.0151583) q[3];
sx q[3];
rz(-1.5446413) q[3];
sx q[3];
rz(1.2873611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6342371) q[2];
sx q[2];
rz(-1.2398182) q[2];
sx q[2];
rz(2.9061387) q[2];
rz(-0.37774751) q[3];
sx q[3];
rz(-0.38769671) q[3];
sx q[3];
rz(-3.0915251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78558755) q[0];
sx q[0];
rz(-1.2879141) q[0];
sx q[0];
rz(-0.054314286) q[0];
rz(-0.76447019) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(-0.70942172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13693181) q[0];
sx q[0];
rz(-2.8747025) q[0];
sx q[0];
rz(-1.1679807) q[0];
rz(-pi) q[1];
rz(1.7058175) q[2];
sx q[2];
rz(-1.3121989) q[2];
sx q[2];
rz(-2.1136346) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1349984) q[1];
sx q[1];
rz(-1.6853764) q[1];
sx q[1];
rz(-2.5496818) q[1];
x q[2];
rz(2.798779) q[3];
sx q[3];
rz(-1.6926994) q[3];
sx q[3];
rz(1.3212622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8556122) q[2];
sx q[2];
rz(-1.6667655) q[2];
sx q[2];
rz(-0.093185514) q[2];
rz(0.5992254) q[3];
sx q[3];
rz(-0.56120509) q[3];
sx q[3];
rz(-0.78527251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7080606) q[0];
sx q[0];
rz(-2.1716433) q[0];
sx q[0];
rz(-0.94950983) q[0];
rz(3.0769729) q[1];
sx q[1];
rz(-3.0977991) q[1];
sx q[1];
rz(-2.480712) q[1];
rz(-1.6914037) q[2];
sx q[2];
rz(-2.5868138) q[2];
sx q[2];
rz(-0.13704337) q[2];
rz(-2.906218) q[3];
sx q[3];
rz(-2.4625851) q[3];
sx q[3];
rz(-2.4877297) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
