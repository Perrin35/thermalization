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
rz(4.6255339) q[0];
sx q[0];
rz(12.892527) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(-1.7915373) q[1];
sx q[1];
rz(1.5989369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5393028) q[0];
sx q[0];
rz(-1.7825583) q[0];
sx q[0];
rz(1.35023) q[0];
rz(0.60249451) q[2];
sx q[2];
rz(-1.3598816) q[2];
sx q[2];
rz(-2.9162625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1118288) q[1];
sx q[1];
rz(-1.3297237) q[1];
sx q[1];
rz(-1.3328711) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78674973) q[3];
sx q[3];
rz(-1.5610715) q[3];
sx q[3];
rz(1.3442638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2797543) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(0.88511434) q[2];
rz(-0.72201133) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(-0.0074145934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-0.95887029) q[0];
sx q[0];
rz(2.0425178) q[0];
rz(-2.4765769) q[1];
sx q[1];
rz(-1.4140833) q[1];
sx q[1];
rz(-0.87759334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1359522) q[0];
sx q[0];
rz(-0.71605039) q[0];
sx q[0];
rz(2.8639249) q[0];
rz(-2.7769233) q[2];
sx q[2];
rz(-1.4450577) q[2];
sx q[2];
rz(1.0276577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6192012) q[1];
sx q[1];
rz(-0.11208216) q[1];
sx q[1];
rz(1.2680608) q[1];
rz(2.8072076) q[3];
sx q[3];
rz(-1.2063724) q[3];
sx q[3];
rz(0.78727608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39891222) q[2];
sx q[2];
rz(-1.30554) q[2];
sx q[2];
rz(1.1304643) q[2];
rz(1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(-1.4484423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.995342) q[0];
sx q[0];
rz(-1.2820219) q[0];
sx q[0];
rz(0.28999844) q[0];
rz(-2.4747804) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(-0.07382948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2732669) q[0];
sx q[0];
rz(-0.10859117) q[0];
sx q[0];
rz(-0.6032087) q[0];
rz(-2.8507112) q[2];
sx q[2];
rz(-1.8996432) q[2];
sx q[2];
rz(-1.8287303) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2245582) q[1];
sx q[1];
rz(-2.9552166) q[1];
sx q[1];
rz(-1.8774377) q[1];
rz(-0.35145268) q[3];
sx q[3];
rz(-1.4593399) q[3];
sx q[3];
rz(-0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(2.4948965) q[2];
rz(2.0329287) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(2.0126608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(-2.1229318) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(-1.4368988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68543363) q[0];
sx q[0];
rz(-1.478727) q[0];
sx q[0];
rz(2.409056) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3827299) q[2];
sx q[2];
rz(-0.43735158) q[2];
sx q[2];
rz(-0.26668374) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.77364591) q[1];
sx q[1];
rz(-1.918643) q[1];
sx q[1];
rz(1.0775954) q[1];
x q[2];
rz(-1.4401682) q[3];
sx q[3];
rz(-1.4454953) q[3];
sx q[3];
rz(-3.0076722) q[3];
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
rz(1.1222703) q[2];
rz(1.026011) q[3];
sx q[3];
rz(-2.3882073) q[3];
sx q[3];
rz(-0.99075738) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68200237) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(-2.7686152) q[0];
rz(2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(-1.3164828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7861159) q[0];
sx q[0];
rz(-1.9516203) q[0];
sx q[0];
rz(1.0771846) q[0];
x q[1];
rz(-0.76422845) q[2];
sx q[2];
rz(-0.35596213) q[2];
sx q[2];
rz(2.228235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.97536918) q[1];
sx q[1];
rz(-1.6227286) q[1];
sx q[1];
rz(0.012253461) q[1];
rz(-pi) q[2];
rz(0.54721197) q[3];
sx q[3];
rz(-1.1037165) q[3];
sx q[3];
rz(0.78391778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-2.310029) q[2];
sx q[2];
rz(2.3357847) q[2];
rz(0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(1.0114975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7508115) q[0];
sx q[0];
rz(-1.3178786) q[0];
sx q[0];
rz(-3.0506296) q[0];
rz(-0.85995752) q[1];
sx q[1];
rz(-2.0188315) q[1];
sx q[1];
rz(1.8213173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0108171) q[0];
sx q[0];
rz(-2.5512619) q[0];
sx q[0];
rz(-2.4369795) q[0];
rz(-2.5713021) q[2];
sx q[2];
rz(-1.7208793) q[2];
sx q[2];
rz(-0.97638408) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.37636138) q[1];
sx q[1];
rz(-2.3665641) q[1];
sx q[1];
rz(-2.2750957) q[1];
x q[2];
rz(-1.201186) q[3];
sx q[3];
rz(-0.69023057) q[3];
sx q[3];
rz(-1.6604916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0992574) q[2];
sx q[2];
rz(-2.1741185) q[2];
sx q[2];
rz(-0.60097224) q[2];
rz(0.48505923) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-2.8619974) q[0];
sx q[0];
rz(-1.9626564) q[0];
sx q[0];
rz(-2.5860508) q[0];
rz(-3.1069966) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.7506036) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7000274) q[0];
sx q[0];
rz(-1.0374984) q[0];
sx q[0];
rz(-1.1839068) q[0];
rz(-pi) q[1];
rz(1.820307) q[2];
sx q[2];
rz(-1.2148641) q[2];
sx q[2];
rz(-0.78636679) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3607189) q[1];
sx q[1];
rz(-0.90457537) q[1];
sx q[1];
rz(-0.16155508) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77565149) q[3];
sx q[3];
rz(-1.5175879) q[3];
sx q[3];
rz(-2.720811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.81327072) q[2];
sx q[2];
rz(-0.44181028) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.8528806) q[3];
sx q[3];
rz(-1.5356531) q[3];
sx q[3];
rz(-2.4718463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.678858) q[0];
sx q[0];
rz(-2.8023219) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(-0.95343268) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(-1.4377726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4576365) q[0];
sx q[0];
rz(-2.5914765) q[0];
sx q[0];
rz(1.9253299) q[0];
rz(-pi) q[1];
rz(-2.9713694) q[2];
sx q[2];
rz(-1.7049689) q[2];
sx q[2];
rz(1.2302878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21282141) q[1];
sx q[1];
rz(-1.0911687) q[1];
sx q[1];
rz(0.70478435) q[1];
x q[2];
rz(1.5290456) q[3];
sx q[3];
rz(-2.6937727) q[3];
sx q[3];
rz(0.24275045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64547223) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(2.9628741) q[2];
rz(2.2802165) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(0.3716968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.20539595) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(2.0478915) q[0];
rz(0.73668346) q[1];
sx q[1];
rz(-1.271558) q[1];
sx q[1];
rz(-1.0587943) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5537542) q[0];
sx q[0];
rz(-2.5231579) q[0];
sx q[0];
rz(2.0536325) q[0];
rz(-0.40760298) q[2];
sx q[2];
rz(-2.0098445) q[2];
sx q[2];
rz(-1.1328732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72653786) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(1.9499669) q[1];
rz(-1.9731673) q[3];
sx q[3];
rz(-1.1306922) q[3];
sx q[3];
rz(0.518706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.27292353) q[2];
sx q[2];
rz(-2.4464641) q[2];
sx q[2];
rz(1.4808562) q[2];
rz(0.41040928) q[3];
sx q[3];
rz(-1.4199665) q[3];
sx q[3];
rz(-0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8695628) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-0.29125443) q[0];
rz(0.60925305) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.7094918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8490484) q[0];
sx q[0];
rz(-1.2546854) q[0];
sx q[0];
rz(-2.6228117) q[0];
x q[1];
rz(0.45551331) q[2];
sx q[2];
rz(-1.3580772) q[2];
sx q[2];
rz(-0.46609512) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0309279) q[1];
sx q[1];
rz(-1.0193766) q[1];
sx q[1];
rz(1.621978) q[1];
rz(-pi) q[2];
rz(0.087177353) q[3];
sx q[3];
rz(-2.0012337) q[3];
sx q[3];
rz(-3.103053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6802406) q[2];
sx q[2];
rz(-2.0246918) q[2];
sx q[2];
rz(2.0001901) q[2];
rz(1.5348148) q[3];
sx q[3];
rz(-1.9635868) q[3];
sx q[3];
rz(-2.6721568) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948982) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(-0.36874157) q[1];
sx q[1];
rz(-1.2422961) q[1];
sx q[1];
rz(3.0098343) q[1];
rz(-2.4070807) q[2];
sx q[2];
rz(-2.8399158) q[2];
sx q[2];
rz(0.5010571) q[2];
rz(2.7183919) q[3];
sx q[3];
rz(-1.3907708) q[3];
sx q[3];
rz(1.6650865) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];