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
rz(1.9361629) q[0];
sx q[0];
rz(-0.054447629) q[0];
sx q[0];
rz(0.86583889) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(3.08334) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6893495) q[0];
sx q[0];
rz(-1.5238683) q[0];
sx q[0];
rz(1.7353135) q[0];
rz(2.0664332) q[2];
sx q[2];
rz(-1.6358536) q[2];
sx q[2];
rz(-3.0954454) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4076947) q[1];
sx q[1];
rz(-1.0419894) q[1];
sx q[1];
rz(2.4771792) q[1];
rz(2.3184003) q[3];
sx q[3];
rz(-0.45837742) q[3];
sx q[3];
rz(-2.2708054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8286579) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(-2.177296) q[2];
rz(-0.079484552) q[3];
sx q[3];
rz(-1.8389523) q[3];
sx q[3];
rz(0.77118072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013414772) q[0];
sx q[0];
rz(-0.34341136) q[0];
sx q[0];
rz(1.6012023) q[0];
rz(-2.3004498) q[1];
sx q[1];
rz(-1.4079739) q[1];
sx q[1];
rz(-2.2841563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090912772) q[0];
sx q[0];
rz(-1.5653531) q[0];
sx q[0];
rz(2.0093657) q[0];
x q[1];
rz(0.93480627) q[2];
sx q[2];
rz(-1.1804712) q[2];
sx q[2];
rz(1.9644511) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5799478) q[1];
sx q[1];
rz(-1.0542634) q[1];
sx q[1];
rz(-2.0793316) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27910618) q[3];
sx q[3];
rz(-1.6630904) q[3];
sx q[3];
rz(1.5971605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88605827) q[2];
sx q[2];
rz(-2.8958791) q[2];
sx q[2];
rz(1.9319755) q[2];
rz(-1.7602734) q[3];
sx q[3];
rz(-1.8807024) q[3];
sx q[3];
rz(-2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7968314) q[0];
sx q[0];
rz(-1.3153356) q[0];
sx q[0];
rz(1.3321846) q[0];
rz(-1.6650797) q[1];
sx q[1];
rz(-1.8107199) q[1];
sx q[1];
rz(0.61980334) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9136494) q[0];
sx q[0];
rz(-1.648128) q[0];
sx q[0];
rz(-1.5864775) q[0];
x q[1];
rz(0.99360433) q[2];
sx q[2];
rz(-2.8393778) q[2];
sx q[2];
rz(2.4156918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1683783) q[1];
sx q[1];
rz(-1.4992979) q[1];
sx q[1];
rz(-0.31798239) q[1];
x q[2];
rz(-1.920082) q[3];
sx q[3];
rz(-0.6699282) q[3];
sx q[3];
rz(2.0921897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0009813112) q[2];
sx q[2];
rz(-0.32537127) q[2];
sx q[2];
rz(0.78835431) q[2];
rz(-1.2761448) q[3];
sx q[3];
rz(-1.2272464) q[3];
sx q[3];
rz(0.86281323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514423) q[0];
sx q[0];
rz(-1.7553512) q[0];
sx q[0];
rz(-2.0748806) q[0];
rz(-1.3534631) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(-1.1563168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24329127) q[0];
sx q[0];
rz(-1.4292846) q[0];
sx q[0];
rz(3.1226052) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5300445) q[2];
sx q[2];
rz(-1.5446071) q[2];
sx q[2];
rz(-1.9350236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.621832) q[1];
sx q[1];
rz(-1.1982188) q[1];
sx q[1];
rz(2.8308059) q[1];
x q[2];
rz(1.5045936) q[3];
sx q[3];
rz(-0.90747031) q[3];
sx q[3];
rz(-2.5732793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68916965) q[2];
sx q[2];
rz(-0.66340476) q[2];
sx q[2];
rz(3.0827674) q[2];
rz(-2.5022653) q[3];
sx q[3];
rz(-1.9202193) q[3];
sx q[3];
rz(-0.85734573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(0.095379742) q[0];
sx q[0];
rz(-1.4154499) q[0];
sx q[0];
rz(0.010490622) q[0];
rz(1.563021) q[1];
sx q[1];
rz(-2.450921) q[1];
sx q[1];
rz(0.48935997) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78016475) q[0];
sx q[0];
rz(-2.1892244) q[0];
sx q[0];
rz(-1.7649107) q[0];
x q[1];
rz(2.0496558) q[2];
sx q[2];
rz(-0.74981028) q[2];
sx q[2];
rz(0.15287003) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53392473) q[1];
sx q[1];
rz(-1.9235538) q[1];
sx q[1];
rz(-2.4466848) q[1];
x q[2];
rz(1.5813835) q[3];
sx q[3];
rz(-1.540031) q[3];
sx q[3];
rz(-0.72584541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8420777) q[2];
sx q[2];
rz(-2.044951) q[2];
sx q[2];
rz(0.00046029885) q[2];
rz(2.4680303) q[3];
sx q[3];
rz(-2.2431777) q[3];
sx q[3];
rz(0.7091929) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29086581) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(1.3036183) q[0];
rz(0.29779008) q[1];
sx q[1];
rz(-0.89556634) q[1];
sx q[1];
rz(-0.29388014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93730629) q[0];
sx q[0];
rz(-2.1171452) q[0];
sx q[0];
rz(2.2541304) q[0];
rz(0.43115012) q[2];
sx q[2];
rz(-2.2047055) q[2];
sx q[2];
rz(2.0946787) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3658947) q[1];
sx q[1];
rz(-1.130321) q[1];
sx q[1];
rz(1.0431402) q[1];
x q[2];
rz(-0.63346432) q[3];
sx q[3];
rz(-2.7567299) q[3];
sx q[3];
rz(0.72609392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46875724) q[2];
sx q[2];
rz(-1.4361278) q[2];
sx q[2];
rz(2.920816) q[2];
rz(-1.8686434) q[3];
sx q[3];
rz(-1.7468529) q[3];
sx q[3];
rz(1.9934191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64977589) q[0];
sx q[0];
rz(-0.62115541) q[0];
sx q[0];
rz(2.5671) q[0];
rz(-0.54221398) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(2.7508459) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0076535) q[0];
sx q[0];
rz(-1.7152733) q[0];
sx q[0];
rz(-0.59916117) q[0];
rz(-1.0263881) q[2];
sx q[2];
rz(-1.5261391) q[2];
sx q[2];
rz(0.80151973) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71902117) q[1];
sx q[1];
rz(-0.41677058) q[1];
sx q[1];
rz(-3.1279492) q[1];
rz(0.096116738) q[3];
sx q[3];
rz(-0.36727723) q[3];
sx q[3];
rz(-2.8678081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6812402) q[2];
sx q[2];
rz(-1.9443024) q[2];
sx q[2];
rz(0.61275068) q[2];
rz(-2.1593306) q[3];
sx q[3];
rz(-1.8270315) q[3];
sx q[3];
rz(0.46510348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.9923582) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(-0.055334844) q[0];
rz(-1.5845567) q[1];
sx q[1];
rz(-2.0535856) q[1];
sx q[1];
rz(0.43992821) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377559) q[0];
sx q[0];
rz(-1.3246264) q[0];
sx q[0];
rz(0.47292821) q[0];
rz(-pi) q[1];
rz(-0.49100809) q[2];
sx q[2];
rz(-1.7517118) q[2];
sx q[2];
rz(2.343585) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4127369) q[1];
sx q[1];
rz(-1.2555033) q[1];
sx q[1];
rz(2.5607177) q[1];
rz(-2.1371577) q[3];
sx q[3];
rz(-1.3310588) q[3];
sx q[3];
rz(2.9842003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7679193) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(-2.8295753) q[2];
rz(-0.9196552) q[3];
sx q[3];
rz(-1.7731526) q[3];
sx q[3];
rz(-2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298518) q[0];
sx q[0];
rz(-0.98167247) q[0];
sx q[0];
rz(0.028129015) q[0];
rz(1.6955388) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(0.45949724) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7336373) q[0];
sx q[0];
rz(-1.6064294) q[0];
sx q[0];
rz(-3.1130457) q[0];
rz(-1.9295646) q[2];
sx q[2];
rz(-0.78132403) q[2];
sx q[2];
rz(1.5272905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64134899) q[1];
sx q[1];
rz(-0.8952221) q[1];
sx q[1];
rz(-0.18570034) q[1];
rz(-pi) q[2];
rz(1.2351843) q[3];
sx q[3];
rz(-0.40235717) q[3];
sx q[3];
rz(-1.5188007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10099899) q[2];
sx q[2];
rz(-1.3446292) q[2];
sx q[2];
rz(-0.25516587) q[2];
rz(-0.98662871) q[3];
sx q[3];
rz(-0.41472236) q[3];
sx q[3];
rz(-1.7184947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7607018) q[0];
sx q[0];
rz(-2.07169) q[0];
sx q[0];
rz(1.3421407) q[0];
rz(-0.0019207151) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(-2.0916746) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59219071) q[0];
sx q[0];
rz(-1.4225385) q[0];
sx q[0];
rz(-1.7853944) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7467612) q[2];
sx q[2];
rz(-1.5696042) q[2];
sx q[2];
rz(-0.81874412) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.109595) q[1];
sx q[1];
rz(-1.2032615) q[1];
sx q[1];
rz(-1.8046741) q[1];
x q[2];
rz(-1.330659) q[3];
sx q[3];
rz(-0.74455816) q[3];
sx q[3];
rz(1.9743686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2125825) q[2];
sx q[2];
rz(-1.9276103) q[2];
sx q[2];
rz(-0.49599656) q[2];
rz(1.4604733) q[3];
sx q[3];
rz(-1.7998327) q[3];
sx q[3];
rz(-1.1869441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7620508) q[0];
sx q[0];
rz(-1.1045781) q[0];
sx q[0];
rz(-1.7212575) q[0];
rz(-0.63915359) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(2.6668702) q[2];
sx q[2];
rz(-1.3480777) q[2];
sx q[2];
rz(2.011275) q[2];
rz(2.3231376) q[3];
sx q[3];
rz(-2.0764805) q[3];
sx q[3];
rz(0.71972328) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
