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
rz(-0.24212317) q[0];
sx q[0];
rz(5.557856) q[0];
sx q[0];
rz(10.331487) q[0];
rz(-0.46696219) q[1];
sx q[1];
rz(-2.2388206) q[1];
sx q[1];
rz(-0.24388193) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8236602) q[0];
sx q[0];
rz(-2.3711304) q[0];
sx q[0];
rz(2.5567358) q[0];
rz(-pi) q[1];
rz(0.73286956) q[2];
sx q[2];
rz(-0.83031619) q[2];
sx q[2];
rz(-0.24093369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4441057) q[1];
sx q[1];
rz(-0.19679697) q[1];
sx q[1];
rz(2.7641731) q[1];
rz(-2.9296845) q[3];
sx q[3];
rz(-2.3125907) q[3];
sx q[3];
rz(0.96860368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4563518) q[2];
sx q[2];
rz(-0.46619236) q[2];
sx q[2];
rz(2.1166128) q[2];
rz(0.13925615) q[3];
sx q[3];
rz(-1.1956297) q[3];
sx q[3];
rz(0.16417424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76992947) q[0];
sx q[0];
rz(-1.0924082) q[0];
sx q[0];
rz(1.9364233) q[0];
rz(1.6734164) q[1];
sx q[1];
rz(-1.2638777) q[1];
sx q[1];
rz(1.42043) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.842671) q[0];
sx q[0];
rz(-0.76220977) q[0];
sx q[0];
rz(-1.803714) q[0];
rz(3.0864079) q[2];
sx q[2];
rz(-1.8351438) q[2];
sx q[2];
rz(0.73433876) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0542595) q[1];
sx q[1];
rz(-1.9078322) q[1];
sx q[1];
rz(-0.23566206) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4010664) q[3];
sx q[3];
rz(-1.2051925) q[3];
sx q[3];
rz(1.4359563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18985192) q[2];
sx q[2];
rz(-0.18988374) q[2];
sx q[2];
rz(-0.80113775) q[2];
rz(2.7028911) q[3];
sx q[3];
rz(-1.2480241) q[3];
sx q[3];
rz(-1.1130822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7640215) q[0];
sx q[0];
rz(-0.29733297) q[0];
sx q[0];
rz(-2.0024894) q[0];
rz(-2.8746919) q[1];
sx q[1];
rz(-1.2511988) q[1];
sx q[1];
rz(-2.7762754) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26317715) q[0];
sx q[0];
rz(-0.52367822) q[0];
sx q[0];
rz(-0.19864638) q[0];
x q[1];
rz(0.92411218) q[2];
sx q[2];
rz(-1.1050997) q[2];
sx q[2];
rz(2.6135824) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2899785) q[1];
sx q[1];
rz(-1.9383926) q[1];
sx q[1];
rz(-0.062003597) q[1];
rz(-0.64008681) q[3];
sx q[3];
rz(-2.1800506) q[3];
sx q[3];
rz(-2.9971788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.9109362) q[2];
sx q[2];
rz(-0.74030423) q[2];
sx q[2];
rz(-0.12348565) q[2];
rz(0.89928818) q[3];
sx q[3];
rz(-1.4495918) q[3];
sx q[3];
rz(-1.9062769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3673636) q[0];
sx q[0];
rz(-1.6772567) q[0];
sx q[0];
rz(2.2520219) q[0];
rz(-0.74596897) q[1];
sx q[1];
rz(-1.6791226) q[1];
sx q[1];
rz(-1.7568582) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33147062) q[0];
sx q[0];
rz(-1.5359914) q[0];
sx q[0];
rz(-1.6286544) q[0];
rz(-0.64176527) q[2];
sx q[2];
rz(-1.7671529) q[2];
sx q[2];
rz(-1.3349443) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1390723) q[1];
sx q[1];
rz(-2.5031282) q[1];
sx q[1];
rz(-2.8248766) q[1];
x q[2];
rz(2.607658) q[3];
sx q[3];
rz(-1.5980532) q[3];
sx q[3];
rz(-1.7202924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7437462) q[2];
sx q[2];
rz(-0.83796871) q[2];
sx q[2];
rz(2.1261334) q[2];
rz(3.1189392) q[3];
sx q[3];
rz(-1.7706784) q[3];
sx q[3];
rz(-3.0034351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7669693) q[0];
sx q[0];
rz(-1.8648819) q[0];
sx q[0];
rz(-2.936777) q[0];
rz(2.9264033) q[1];
sx q[1];
rz(-0.22061017) q[1];
sx q[1];
rz(0.87517103) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9143863) q[0];
sx q[0];
rz(-1.7191186) q[0];
sx q[0];
rz(-1.877584) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3059101) q[2];
sx q[2];
rz(-1.4561404) q[2];
sx q[2];
rz(-2.3029652) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4103247) q[1];
sx q[1];
rz(-1.7909697) q[1];
sx q[1];
rz(2.8360044) q[1];
rz(-pi) q[2];
rz(-1.3761843) q[3];
sx q[3];
rz(-2.5014957) q[3];
sx q[3];
rz(2.3605337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.63336664) q[2];
sx q[2];
rz(-0.88819155) q[2];
sx q[2];
rz(1.7249031) q[2];
rz(-2.21777) q[3];
sx q[3];
rz(-0.97771907) q[3];
sx q[3];
rz(2.0098604) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1205587) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(2.293204) q[0];
rz(1.5757164) q[1];
sx q[1];
rz(-1.3278241) q[1];
sx q[1];
rz(-2.1956992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3065328) q[0];
sx q[0];
rz(-0.9047821) q[0];
sx q[0];
rz(1.205797) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81056548) q[2];
sx q[2];
rz(-1.7195576) q[2];
sx q[2];
rz(-2.8523977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1523244) q[1];
sx q[1];
rz(-2.3890308) q[1];
sx q[1];
rz(-2.4965733) q[1];
x q[2];
rz(0.49094851) q[3];
sx q[3];
rz(-1.3958389) q[3];
sx q[3];
rz(-3.0096731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1762323) q[2];
sx q[2];
rz(-1.2364028) q[2];
sx q[2];
rz(-1.7410834) q[2];
rz(-1.4817574) q[3];
sx q[3];
rz(-1.7912495) q[3];
sx q[3];
rz(-0.19937936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7401212) q[0];
sx q[0];
rz(-0.49982247) q[0];
sx q[0];
rz(-1.9299782) q[0];
rz(-2.6648193) q[1];
sx q[1];
rz(-1.8649273) q[1];
sx q[1];
rz(1.5062987) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7619373) q[0];
sx q[0];
rz(-0.64122824) q[0];
sx q[0];
rz(2.0564709) q[0];
rz(-1.2477307) q[2];
sx q[2];
rz(-1.808721) q[2];
sx q[2];
rz(2.0428773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16852779) q[1];
sx q[1];
rz(-2.2650411) q[1];
sx q[1];
rz(1.7315052) q[1];
rz(-pi) q[2];
rz(-2.659306) q[3];
sx q[3];
rz(-1.3083754) q[3];
sx q[3];
rz(-0.69604128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8022884) q[2];
sx q[2];
rz(-1.6243287) q[2];
sx q[2];
rz(-0.18292546) q[2];
rz(-0.66169345) q[3];
sx q[3];
rz(-1.2122093) q[3];
sx q[3];
rz(0.70484149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6905007) q[0];
sx q[0];
rz(-1.4649614) q[0];
sx q[0];
rz(0.14725421) q[0];
rz(2.9946949) q[1];
sx q[1];
rz(-2.2203827) q[1];
sx q[1];
rz(-1.1796835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3203293) q[0];
sx q[0];
rz(-2.593145) q[0];
sx q[0];
rz(0.68446496) q[0];
rz(0.4381035) q[2];
sx q[2];
rz(-2.0102276) q[2];
sx q[2];
rz(0.319744) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1547566) q[1];
sx q[1];
rz(-0.97797003) q[1];
sx q[1];
rz(1.4378217) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3875119) q[3];
sx q[3];
rz(-1.6508266) q[3];
sx q[3];
rz(2.4802993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94072056) q[2];
sx q[2];
rz(-2.2042037) q[2];
sx q[2];
rz(-2.3737523) q[2];
rz(2.614894) q[3];
sx q[3];
rz(-1.0517164) q[3];
sx q[3];
rz(-2.5832978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(1.4621157) q[0];
sx q[0];
rz(-0.18590346) q[0];
sx q[0];
rz(2.0515077) q[0];
rz(1.7915626) q[1];
sx q[1];
rz(-1.7410024) q[1];
sx q[1];
rz(0.42492351) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23608828) q[0];
sx q[0];
rz(-2.6668735) q[0];
sx q[0];
rz(1.6951872) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0504136) q[2];
sx q[2];
rz(-2.489739) q[2];
sx q[2];
rz(-0.21287795) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7817345) q[1];
sx q[1];
rz(-1.1175851) q[1];
sx q[1];
rz(-1.9537326) q[1];
rz(-2.0941689) q[3];
sx q[3];
rz(-2.2096388) q[3];
sx q[3];
rz(2.1372319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9099884) q[2];
sx q[2];
rz(-2.6308172) q[2];
sx q[2];
rz(0.52499047) q[2];
rz(1.3548343) q[3];
sx q[3];
rz(-0.89559186) q[3];
sx q[3];
rz(1.3625328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43596426) q[0];
sx q[0];
rz(-2.4810677) q[0];
sx q[0];
rz(0.10667644) q[0];
rz(1.0542997) q[1];
sx q[1];
rz(-1.7091457) q[1];
sx q[1];
rz(1.7342825) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6330235) q[0];
sx q[0];
rz(-2.3481524) q[0];
sx q[0];
rz(0.88280789) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58511795) q[2];
sx q[2];
rz(-0.74293929) q[2];
sx q[2];
rz(-1.9415426) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.046316) q[1];
sx q[1];
rz(-1.3477948) q[1];
sx q[1];
rz(-0.15671244) q[1];
rz(-pi) q[2];
rz(2.3913283) q[3];
sx q[3];
rz(-2.0553053) q[3];
sx q[3];
rz(-1.110525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4878896) q[2];
sx q[2];
rz(-0.50450486) q[2];
sx q[2];
rz(2.8686236) q[2];
rz(-0.87131396) q[3];
sx q[3];
rz(-1.4600236) q[3];
sx q[3];
rz(-3.0070846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7731666) q[0];
sx q[0];
rz(-1.1933403) q[0];
sx q[0];
rz(-2.2884952) q[0];
rz(2.7586965) q[1];
sx q[1];
rz(-1.8538414) q[1];
sx q[1];
rz(-0.014898653) q[1];
rz(2.0128925) q[2];
sx q[2];
rz(-1.0539891) q[2];
sx q[2];
rz(-2.5602387) q[2];
rz(-1.6712345) q[3];
sx q[3];
rz(-0.77340579) q[3];
sx q[3];
rz(0.24271942) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
