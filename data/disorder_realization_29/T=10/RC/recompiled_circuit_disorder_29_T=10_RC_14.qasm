OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0260789) q[0];
sx q[0];
rz(-1.6576515) q[0];
sx q[0];
rz(-2.8154362) q[0];
rz(-1.1905319) q[1];
sx q[1];
rz(-1.3500554) q[1];
sx q[1];
rz(-1.5989369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9261949) q[0];
sx q[0];
rz(-0.30456671) q[0];
sx q[0];
rz(0.79415168) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60249451) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(2.9162625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5988679) q[1];
sx q[1];
rz(-1.3398783) q[1];
sx q[1];
rz(0.24778194) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1278586) q[3];
sx q[3];
rz(-0.78679689) q[3];
sx q[3];
rz(2.9247583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(-0.88511434) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(-1.0990748) q[0];
rz(0.66501578) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(-0.87759334) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646334) q[0];
sx q[0];
rz(-1.7517125) q[0];
sx q[0];
rz(0.69676708) q[0];
rz(-pi) q[1];
rz(-0.36466937) q[2];
sx q[2];
rz(-1.4450577) q[2];
sx q[2];
rz(-1.0276577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.5223915) q[1];
sx q[1];
rz(-0.11208216) q[1];
sx q[1];
rz(-1.2680608) q[1];
rz(-pi) q[2];
rz(0.86032805) q[3];
sx q[3];
rz(-0.48947696) q[3];
sx q[3];
rz(-1.5599172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.39891222) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(-2.0111283) q[2];
rz(-1.8418664) q[3];
sx q[3];
rz(-1.2239417) q[3];
sx q[3];
rz(-1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14625064) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(0.28999844) q[0];
rz(2.4747804) q[1];
sx q[1];
rz(-1.0338444) q[1];
sx q[1];
rz(3.0677632) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4742972) q[0];
sx q[0];
rz(-1.6601666) q[0];
sx q[0];
rz(1.5090293) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29088144) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(1.3128624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5362629) q[1];
sx q[1];
rz(-1.7483835) q[1];
sx q[1];
rz(-0.0568584) q[1];
rz(-pi) q[2];
rz(1.4521493) q[3];
sx q[3];
rz(-1.2216179) q[3];
sx q[3];
rz(2.196764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4905711) q[2];
sx q[2];
rz(-1.1940424) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(-2.0329287) q[3];
sx q[3];
rz(-2.3587148) q[3];
sx q[3];
rz(-1.1289319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17764238) q[0];
sx q[0];
rz(-0.17340604) q[0];
sx q[0];
rz(-1.1886764) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-0.97266346) q[1];
sx q[1];
rz(1.7046938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80285145) q[0];
sx q[0];
rz(-2.2995298) q[0];
sx q[0];
rz(1.4472423) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2595348) q[2];
sx q[2];
rz(-1.2584104) q[2];
sx q[2];
rz(2.600008) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.23197933) q[1];
sx q[1];
rz(-2.5464499) q[1];
sx q[1];
rz(-0.91722782) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8021637) q[3];
sx q[3];
rz(-0.18076104) q[3];
sx q[3];
rz(-2.4651431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58549762) q[2];
sx q[2];
rz(-1.4336339) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(0.99075738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595903) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(2.7686152) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.9517027) q[1];
sx q[1];
rz(-1.8251098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81985695) q[0];
sx q[0];
rz(-2.5279191) q[0];
sx q[0];
rz(0.86921285) q[0];
rz(-pi) q[1];
rz(2.3773642) q[2];
sx q[2];
rz(-2.7856305) q[2];
sx q[2];
rz(-2.228235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97536918) q[1];
sx q[1];
rz(-1.6227286) q[1];
sx q[1];
rz(3.1293392) q[1];
rz(2.1042215) q[3];
sx q[3];
rz(-2.0539527) q[3];
sx q[3];
rz(2.0869568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10107723) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(2.3357847) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(-2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7508115) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(3.0506296) q[0];
rz(0.85995752) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(1.3202753) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4754007) q[0];
sx q[0];
rz(-2.0087419) q[0];
sx q[0];
rz(1.1613208) q[0];
rz(-2.5713021) q[2];
sx q[2];
rz(-1.4207134) q[2];
sx q[2];
rz(-2.1652086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7652313) q[1];
sx q[1];
rz(-2.3665641) q[1];
sx q[1];
rz(-2.2750957) q[1];
rz(1.9404066) q[3];
sx q[3];
rz(-2.4513621) q[3];
sx q[3];
rz(1.6604916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0992574) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(2.6565334) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(1.4453567) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27959529) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(-0.55554187) q[0];
rz(0.034596054) q[1];
sx q[1];
rz(-2.3831773) q[1];
sx q[1];
rz(-1.7506036) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7000274) q[0];
sx q[0];
rz(-1.0374984) q[0];
sx q[0];
rz(-1.1839068) q[0];
x q[1];
rz(-2.7752635) q[2];
sx q[2];
rz(-1.3372256) q[2];
sx q[2];
rz(-0.69586588) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.890306) q[1];
sx q[1];
rz(-1.6975926) q[1];
sx q[1];
rz(2.2433953) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4963385) q[3];
sx q[3];
rz(-2.3450608) q[3];
sx q[3];
rz(1.9394685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81327072) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(-0.39548809) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46273461) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(1.6495552) q[0];
rz(-0.95343268) q[1];
sx q[1];
rz(-1.1089193) q[1];
sx q[1];
rz(-1.7038201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9490818) q[0];
sx q[0];
rz(-1.7532945) q[0];
sx q[0];
rz(-1.048868) q[0];
rz(-pi) q[1];
rz(2.4687924) q[2];
sx q[2];
rz(-0.2163419) q[2];
sx q[2];
rz(-1.0018444) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9287712) q[1];
sx q[1];
rz(-2.050424) q[1];
sx q[1];
rz(0.70478435) q[1];
x q[2];
rz(1.5290456) q[3];
sx q[3];
rz(-2.6937727) q[3];
sx q[3];
rz(-2.8988422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4961204) q[2];
sx q[2];
rz(-2.7068553) q[2];
sx q[2];
rz(0.17871857) q[2];
rz(2.2802165) q[3];
sx q[3];
rz(-1.2025611) q[3];
sx q[3];
rz(-0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9361967) q[0];
sx q[0];
rz(-1.9367171) q[0];
sx q[0];
rz(1.0937011) q[0];
rz(-0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(1.0587943) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5878384) q[0];
sx q[0];
rz(-2.5231579) q[0];
sx q[0];
rz(-1.0879602) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86973127) q[2];
sx q[2];
rz(-2.5517002) q[2];
sx q[2];
rz(1.9260977) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0000856) q[1];
sx q[1];
rz(-0.39751378) q[1];
sx q[1];
rz(1.8915528) q[1];
rz(-pi) q[2];
rz(1.9731673) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(-2.6228867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8686691) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(0.41040928) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(-0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-2.8503382) q[0];
rz(-2.5323396) q[1];
sx q[1];
rz(-1.4657425) q[1];
sx q[1];
rz(1.4321009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0390022) q[0];
sx q[0];
rz(-2.0615091) q[0];
sx q[0];
rz(-1.2105932) q[0];
x q[1];
rz(2.6851995) q[2];
sx q[2];
rz(-0.49955873) q[2];
sx q[2];
rz(1.6299786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4869625) q[1];
sx q[1];
rz(-1.527206) q[1];
sx q[1];
rz(-0.55200465) q[1];
rz(1.7581975) q[3];
sx q[3];
rz(-2.7029576) q[3];
sx q[3];
rz(2.8965829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46135205) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(1.1414026) q[2];
rz(-1.5348148) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(0.46943584) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5466945) q[0];
sx q[0];
rz(-1.5190769) q[0];
sx q[0];
rz(1.4357823) q[0];
rz(2.7728511) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(1.776406) q[2];
sx q[2];
rz(-1.3484577) q[2];
sx q[2];
rz(1.2586013) q[2];
rz(-0.42320078) q[3];
sx q[3];
rz(-1.3907708) q[3];
sx q[3];
rz(1.6650865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
