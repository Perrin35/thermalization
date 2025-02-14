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
rz(2.5166729) q[0];
sx q[0];
rz(-1.807037) q[0];
sx q[0];
rz(-3.0883489) q[0];
rz(-2.6553395) q[1];
sx q[1];
rz(-3.0142205) q[1];
sx q[1];
rz(1.4349586) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.931728) q[0];
sx q[0];
rz(-1.7057944) q[0];
sx q[0];
rz(-0.28688669) q[0];
rz(1.6900914) q[2];
sx q[2];
rz(-1.8362459) q[2];
sx q[2];
rz(-2.1448898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26027825) q[1];
sx q[1];
rz(-1.5173755) q[1];
sx q[1];
rz(-2.747606) q[1];
rz(-pi) q[2];
rz(-2.1987183) q[3];
sx q[3];
rz(-0.74340313) q[3];
sx q[3];
rz(1.8441895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3367553) q[2];
sx q[2];
rz(-1.7982322) q[2];
sx q[2];
rz(0.27377823) q[2];
rz(1.9233507) q[3];
sx q[3];
rz(-2.4680586) q[3];
sx q[3];
rz(-0.28606733) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022920595) q[0];
sx q[0];
rz(-0.76347041) q[0];
sx q[0];
rz(2.3619695) q[0];
rz(0.66462213) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(1.2453311) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095663505) q[0];
sx q[0];
rz(-1.6324658) q[0];
sx q[0];
rz(1.685919) q[0];
x q[1];
rz(1.4636366) q[2];
sx q[2];
rz(-2.4273211) q[2];
sx q[2];
rz(2.799054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9960963) q[1];
sx q[1];
rz(-0.22095535) q[1];
sx q[1];
rz(1.851382) q[1];
x q[2];
rz(0.20540463) q[3];
sx q[3];
rz(-1.0213791) q[3];
sx q[3];
rz(-2.7618264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43210426) q[2];
sx q[2];
rz(-0.69807845) q[2];
sx q[2];
rz(-3.091605) q[2];
rz(1.4738119) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(0.93562359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899984) q[0];
sx q[0];
rz(-2.279156) q[0];
sx q[0];
rz(1.9631901) q[0];
rz(-1.1475457) q[1];
sx q[1];
rz(-2.9541364) q[1];
sx q[1];
rz(0.60290927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1499304) q[0];
sx q[0];
rz(-1.4105182) q[0];
sx q[0];
rz(-0.10326881) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1772867) q[2];
sx q[2];
rz(-0.78769257) q[2];
sx q[2];
rz(2.380013) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29290043) q[1];
sx q[1];
rz(-2.22977) q[1];
sx q[1];
rz(1.0408569) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1985377) q[3];
sx q[3];
rz(-0.78305093) q[3];
sx q[3];
rz(0.32367009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9511562) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(-2.7583165) q[2];
rz(1.8303309) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(-3.0200628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7168147) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(2.2130261) q[0];
rz(-0.83944744) q[1];
sx q[1];
rz(-1.3176368) q[1];
sx q[1];
rz(-2.3323257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1195212) q[0];
sx q[0];
rz(-0.68953994) q[0];
sx q[0];
rz(2.1697696) q[0];
rz(-pi) q[1];
rz(1.4014612) q[2];
sx q[2];
rz(-1.1462829) q[2];
sx q[2];
rz(-0.2178387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88621186) q[1];
sx q[1];
rz(-1.6396697) q[1];
sx q[1];
rz(-1.8167472) q[1];
rz(-pi) q[2];
rz(0.92692848) q[3];
sx q[3];
rz(-1.2371105) q[3];
sx q[3];
rz(-2.9179171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3889918) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(1.4468225) q[2];
rz(-2.0145448) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(0.62895044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9160354) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(1.7115364) q[0];
rz(0.23722181) q[1];
sx q[1];
rz(-2.1784541) q[1];
sx q[1];
rz(2.846834) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5200318) q[0];
sx q[0];
rz(-2.6803117) q[0];
sx q[0];
rz(-0.3256021) q[0];
x q[1];
rz(-2.0343658) q[2];
sx q[2];
rz(-2.4245302) q[2];
sx q[2];
rz(-2.0328558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.39043104) q[1];
sx q[1];
rz(-1.1220349) q[1];
sx q[1];
rz(-0.46350355) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5676457) q[3];
sx q[3];
rz(-1.6960295) q[3];
sx q[3];
rz(2.4698911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99990591) q[2];
sx q[2];
rz(-2.9408231) q[2];
sx q[2];
rz(3.0086009) q[2];
rz(0.76987949) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(2.8225115) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8517476) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.7293365) q[0];
rz(0.051636592) q[1];
sx q[1];
rz(-0.62283555) q[1];
sx q[1];
rz(0.19270611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4225005) q[0];
sx q[0];
rz(-2.0845045) q[0];
sx q[0];
rz(-1.5490319) q[0];
rz(2.4102845) q[2];
sx q[2];
rz(-1.0885065) q[2];
sx q[2];
rz(2.363229) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6289857) q[1];
sx q[1];
rz(-2.4977105) q[1];
sx q[1];
rz(-2.3888247) q[1];
x q[2];
rz(1.9514378) q[3];
sx q[3];
rz(-2.1771113) q[3];
sx q[3];
rz(-2.9553138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0685588) q[2];
sx q[2];
rz(-1.8517588) q[2];
sx q[2];
rz(1.7603091) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-2.2541855) q[3];
sx q[3];
rz(-1.380434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16159049) q[0];
sx q[0];
rz(-1.1450293) q[0];
sx q[0];
rz(-0.34647754) q[0];
rz(-1.0193635) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(1.4541218) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2812735) q[0];
sx q[0];
rz(-2.3695282) q[0];
sx q[0];
rz(0.66769974) q[0];
x q[1];
rz(2.2242111) q[2];
sx q[2];
rz(-1.1251984) q[2];
sx q[2];
rz(-0.080597045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2994308) q[1];
sx q[1];
rz(-0.60363942) q[1];
sx q[1];
rz(1.0697212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4164657) q[3];
sx q[3];
rz(-1.4714421) q[3];
sx q[3];
rz(1.9973444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5399897) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(0.23923624) q[2];
rz(-0.38419497) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(-0.95170963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20188986) q[0];
sx q[0];
rz(-1.4677784) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(-1.2205623) q[1];
sx q[1];
rz(-2.2790597) q[1];
sx q[1];
rz(-2.9753704) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3275422) q[0];
sx q[0];
rz(-2.6084427) q[0];
sx q[0];
rz(-2.0943805) q[0];
rz(-1.2532089) q[2];
sx q[2];
rz(-0.86195213) q[2];
sx q[2];
rz(-0.52121693) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4835755) q[1];
sx q[1];
rz(-1.3596351) q[1];
sx q[1];
rz(1.5394475) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9907624) q[3];
sx q[3];
rz(-1.8059429) q[3];
sx q[3];
rz(1.3104035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.347747) q[2];
sx q[2];
rz(-2.2913427) q[2];
sx q[2];
rz(-1.0469077) q[2];
rz(1.2547803) q[3];
sx q[3];
rz(-1.0898277) q[3];
sx q[3];
rz(1.2966398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60085249) q[0];
sx q[0];
rz(-2.7678601) q[0];
sx q[0];
rz(1.1131713) q[0];
rz(2.3241849) q[1];
sx q[1];
rz(-1.6956885) q[1];
sx q[1];
rz(2.7730952) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21462003) q[0];
sx q[0];
rz(-1.9153739) q[0];
sx q[0];
rz(-2.8520682) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13780528) q[2];
sx q[2];
rz(-0.74713444) q[2];
sx q[2];
rz(-1.8835619) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52942136) q[1];
sx q[1];
rz(-2.2670806) q[1];
sx q[1];
rz(-1.9208292) q[1];
rz(2.1300674) q[3];
sx q[3];
rz(-0.91010053) q[3];
sx q[3];
rz(1.504577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9937399) q[2];
sx q[2];
rz(-2.5233614) q[2];
sx q[2];
rz(-1.7842133) q[2];
rz(-2.3944858) q[3];
sx q[3];
rz(-0.96304572) q[3];
sx q[3];
rz(0.67556206) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52637446) q[0];
sx q[0];
rz(-2.6850061) q[0];
sx q[0];
rz(1.0427465) q[0];
rz(0.77049795) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(2.3811293) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9077936) q[0];
sx q[0];
rz(-1.2408419) q[0];
sx q[0];
rz(-0.25262649) q[0];
rz(-1.3622401) q[2];
sx q[2];
rz(-2.957666) q[2];
sx q[2];
rz(-0.49801258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14899602) q[1];
sx q[1];
rz(-1.9933874) q[1];
sx q[1];
rz(1.484297) q[1];
x q[2];
rz(-2.4861927) q[3];
sx q[3];
rz(-1.2714579) q[3];
sx q[3];
rz(-1.4531524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.20327917) q[2];
sx q[2];
rz(-2.5310897) q[2];
sx q[2];
rz(2.7321613) q[2];
rz(-1.5583386) q[3];
sx q[3];
rz(-2.6810472) q[3];
sx q[3];
rz(-0.62622768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64774491) q[0];
sx q[0];
rz(-1.0418325) q[0];
sx q[0];
rz(-2.7000725) q[0];
rz(-1.5589177) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(0.51519338) q[2];
sx q[2];
rz(-0.11820785) q[2];
sx q[2];
rz(-1.6109656) q[2];
rz(-1.4382451) q[3];
sx q[3];
rz(-1.3302531) q[3];
sx q[3];
rz(2.1147685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
