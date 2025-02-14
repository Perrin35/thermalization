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
rz(2.1287542) q[0];
sx q[0];
rz(-1.1829809) q[0];
sx q[0];
rz(0.24721375) q[0];
rz(1.5565058) q[1];
sx q[1];
rz(-1.0143919) q[1];
sx q[1];
rz(-0.95451626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.010095) q[0];
sx q[0];
rz(-2.7662377) q[0];
sx q[0];
rz(1.6502871) q[0];
x q[1];
rz(-2.688301) q[2];
sx q[2];
rz(-1.5201836) q[2];
sx q[2];
rz(2.5091189) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2672429) q[1];
sx q[1];
rz(-1.8863719) q[1];
sx q[1];
rz(-0.79488956) q[1];
rz(-2.3684816) q[3];
sx q[3];
rz(-1.3476117) q[3];
sx q[3];
rz(0.091835819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1859833) q[2];
sx q[2];
rz(-2.535203) q[2];
sx q[2];
rz(-0.28140226) q[2];
rz(1.227281) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(-1.9277771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1948552) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(2.5545004) q[0];
rz(0.52014822) q[1];
sx q[1];
rz(-1.9303493) q[1];
sx q[1];
rz(-2.928226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30946975) q[0];
sx q[0];
rz(-1.7209527) q[0];
sx q[0];
rz(0.062376745) q[0];
x q[1];
rz(0.93307067) q[2];
sx q[2];
rz(-0.75680774) q[2];
sx q[2];
rz(-2.840586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50700106) q[1];
sx q[1];
rz(-1.6640955) q[1];
sx q[1];
rz(-3.0970645) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1150241) q[3];
sx q[3];
rz(-2.5682862) q[3];
sx q[3];
rz(0.70070964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0352036) q[2];
sx q[2];
rz(-0.28855244) q[2];
sx q[2];
rz(-0.2717379) q[2];
rz(-0.64618293) q[3];
sx q[3];
rz(-1.313442) q[3];
sx q[3];
rz(-1.2077695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81994098) q[0];
sx q[0];
rz(-0.8701179) q[0];
sx q[0];
rz(-0.2680378) q[0];
rz(2.9053814) q[1];
sx q[1];
rz(-1.3583438) q[1];
sx q[1];
rz(1.4150298) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9504323) q[0];
sx q[0];
rz(-1.6163905) q[0];
sx q[0];
rz(-0.10204236) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7790952) q[2];
sx q[2];
rz(-1.4115184) q[2];
sx q[2];
rz(0.25937072) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1620924) q[1];
sx q[1];
rz(-2.0071173) q[1];
sx q[1];
rz(-1.8772582) q[1];
rz(-pi) q[2];
rz(2.4684942) q[3];
sx q[3];
rz(-1.3521225) q[3];
sx q[3];
rz(-0.84450561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1872824) q[2];
sx q[2];
rz(-0.3287181) q[2];
sx q[2];
rz(1.7598565) q[2];
rz(-0.36661822) q[3];
sx q[3];
rz(-1.4724052) q[3];
sx q[3];
rz(3.0439175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5496552) q[0];
sx q[0];
rz(-0.20579919) q[0];
sx q[0];
rz(0.016121443) q[0];
rz(1.6751809) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(2.256567) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8649053) q[0];
sx q[0];
rz(-1.7430796) q[0];
sx q[0];
rz(0.40934632) q[0];
rz(2.9134404) q[2];
sx q[2];
rz(-1.9288262) q[2];
sx q[2];
rz(-1.6659074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.64513834) q[1];
sx q[1];
rz(-2.0296603) q[1];
sx q[1];
rz(-1.9978844) q[1];
rz(-pi) q[2];
rz(-1.0543941) q[3];
sx q[3];
rz(-1.9266911) q[3];
sx q[3];
rz(-0.10117029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5625988) q[2];
sx q[2];
rz(-2.6851974) q[2];
sx q[2];
rz(1.5835416) q[2];
rz(1.7700178) q[3];
sx q[3];
rz(-1.8669502) q[3];
sx q[3];
rz(0.19545999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1564932) q[0];
sx q[0];
rz(-1.766196) q[0];
sx q[0];
rz(-2.9920355) q[0];
rz(-2.0984446) q[1];
sx q[1];
rz(-2.1688192) q[1];
sx q[1];
rz(-2.3057888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59398932) q[0];
sx q[0];
rz(-2.7271252) q[0];
sx q[0];
rz(2.2311121) q[0];
x q[1];
rz(2.2611375) q[2];
sx q[2];
rz(-1.5780851) q[2];
sx q[2];
rz(-1.7178423) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0975102) q[1];
sx q[1];
rz(-0.20803504) q[1];
sx q[1];
rz(1.5745622) q[1];
rz(-pi) q[2];
rz(-2.1981524) q[3];
sx q[3];
rz(-0.30590484) q[3];
sx q[3];
rz(2.4750575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.122637) q[2];
sx q[2];
rz(-1.6654207) q[2];
sx q[2];
rz(2.7254851) q[2];
rz(0.082402669) q[3];
sx q[3];
rz(-2.1722138) q[3];
sx q[3];
rz(-0.20127067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06269726) q[0];
sx q[0];
rz(-2.1128928) q[0];
sx q[0];
rz(2.0507574) q[0];
rz(1.1385607) q[1];
sx q[1];
rz(-1.1999612) q[1];
sx q[1];
rz(0.60637766) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2045468) q[0];
sx q[0];
rz(-1.8228476) q[0];
sx q[0];
rz(0.65752794) q[0];
rz(-pi) q[1];
rz(-0.80198432) q[2];
sx q[2];
rz(-1.238566) q[2];
sx q[2];
rz(1.3114245) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5910037) q[1];
sx q[1];
rz(-2.2865851) q[1];
sx q[1];
rz(0.71859931) q[1];
x q[2];
rz(0.13952556) q[3];
sx q[3];
rz(-2.4184113) q[3];
sx q[3];
rz(0.65825247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.828317) q[2];
sx q[2];
rz(-2.5930391) q[2];
sx q[2];
rz(0.11087785) q[2];
rz(1.5634792) q[3];
sx q[3];
rz(-0.21136798) q[3];
sx q[3];
rz(1.8160688) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9368847) q[0];
sx q[0];
rz(-2.3218343) q[0];
sx q[0];
rz(-2.4844266) q[0];
rz(-1.8183297) q[1];
sx q[1];
rz(-2.3902049) q[1];
sx q[1];
rz(2.0164067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60356451) q[0];
sx q[0];
rz(-1.5549994) q[0];
sx q[0];
rz(-1.5558916) q[0];
rz(-0.31382896) q[2];
sx q[2];
rz(-1.6633361) q[2];
sx q[2];
rz(1.0864663) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4309498) q[1];
sx q[1];
rz(-2.498327) q[1];
sx q[1];
rz(-0.67420723) q[1];
x q[2];
rz(-2.4666305) q[3];
sx q[3];
rz(-1.3276023) q[3];
sx q[3];
rz(-1.4362818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.40921673) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(0.96119514) q[2];
rz(0.0088648908) q[3];
sx q[3];
rz(-0.88518393) q[3];
sx q[3];
rz(-1.9302906) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0601785) q[0];
sx q[0];
rz(-3.0280805) q[0];
sx q[0];
rz(2.7139582) q[0];
rz(-1.112452) q[1];
sx q[1];
rz(-1.3149657) q[1];
sx q[1];
rz(0.99666673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8450981) q[0];
sx q[0];
rz(-1.0834435) q[0];
sx q[0];
rz(-2.4250406) q[0];
rz(-0.76297595) q[2];
sx q[2];
rz(-1.7691649) q[2];
sx q[2];
rz(-1.6301483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3370875) q[1];
sx q[1];
rz(-1.8356016) q[1];
sx q[1];
rz(2.1534285) q[1];
rz(-pi) q[2];
rz(-0.31672041) q[3];
sx q[3];
rz(-2.5473928) q[3];
sx q[3];
rz(-2.7946928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6796278) q[2];
sx q[2];
rz(-0.5373911) q[2];
sx q[2];
rz(-2.4930387) q[2];
rz(2.8280761) q[3];
sx q[3];
rz(-0.88839141) q[3];
sx q[3];
rz(-3.0891109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.736883) q[0];
sx q[0];
rz(-2.9209904) q[0];
sx q[0];
rz(0.96694651) q[0];
rz(2.1838358) q[1];
sx q[1];
rz(-1.0916595) q[1];
sx q[1];
rz(-2.0584094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5388018) q[0];
sx q[0];
rz(-1.1447009) q[0];
sx q[0];
rz(-1.8417131) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4114001) q[2];
sx q[2];
rz(-1.0456523) q[2];
sx q[2];
rz(0.60827916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5450648) q[1];
sx q[1];
rz(-1.178243) q[1];
sx q[1];
rz(-0.82253455) q[1];
rz(0.47934909) q[3];
sx q[3];
rz(-1.5707301) q[3];
sx q[3];
rz(-2.655055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.12072418) q[2];
sx q[2];
rz(-2.3126297) q[2];
sx q[2];
rz(1.4618358) q[2];
rz(2.7802137) q[3];
sx q[3];
rz(-1.3249818) q[3];
sx q[3];
rz(2.1593275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(0.35230377) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(-0.58255449) q[0];
rz(1.7763058) q[1];
sx q[1];
rz(-1.0368232) q[1];
sx q[1];
rz(2.9127311) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6270646) q[0];
sx q[0];
rz(-0.43958966) q[0];
sx q[0];
rz(-2.702781) q[0];
rz(-pi) q[1];
rz(-0.17915704) q[2];
sx q[2];
rz(-1.3115885) q[2];
sx q[2];
rz(1.955206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46581546) q[1];
sx q[1];
rz(-2.9699616) q[1];
sx q[1];
rz(-1.1824498) q[1];
rz(-pi) q[2];
rz(-1.1353605) q[3];
sx q[3];
rz(-1.9276016) q[3];
sx q[3];
rz(2.8716068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87054306) q[2];
sx q[2];
rz(-0.55176631) q[2];
sx q[2];
rz(-1.40847) q[2];
rz(-0.72077858) q[3];
sx q[3];
rz(-1.1495138) q[3];
sx q[3];
rz(1.518998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0558753) q[0];
sx q[0];
rz(-1.7039104) q[0];
sx q[0];
rz(-1.3088551) q[0];
rz(1.5720639) q[1];
sx q[1];
rz(-0.61620284) q[1];
sx q[1];
rz(-2.9175704) q[1];
rz(-2.1355676) q[2];
sx q[2];
rz(-1.7823162) q[2];
sx q[2];
rz(1.6020365) q[2];
rz(-1.4103945) q[3];
sx q[3];
rz(-1.7216202) q[3];
sx q[3];
rz(-1.3826821) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
