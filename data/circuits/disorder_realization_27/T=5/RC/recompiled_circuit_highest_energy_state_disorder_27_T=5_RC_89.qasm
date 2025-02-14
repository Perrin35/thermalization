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
rz(1.2394387) q[0];
sx q[0];
rz(-1.8128938) q[0];
sx q[0];
rz(-0.21188307) q[0];
rz(-1.4172685) q[1];
sx q[1];
rz(-2.6091726) q[1];
sx q[1];
rz(0.37791696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.657562) q[0];
sx q[0];
rz(-0.0095417984) q[0];
sx q[0];
rz(-1.5242759) q[0];
rz(-pi) q[1];
rz(-2.7013999) q[2];
sx q[2];
rz(-2.4131916) q[2];
sx q[2];
rz(0.0022526646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18641414) q[1];
sx q[1];
rz(-2.4354773) q[1];
sx q[1];
rz(2.7104055) q[1];
x q[2];
rz(-2.3842281) q[3];
sx q[3];
rz(-0.93776449) q[3];
sx q[3];
rz(2.4773134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2892896) q[2];
sx q[2];
rz(-2.0630344) q[2];
sx q[2];
rz(-0.079785384) q[2];
rz(2.1763109) q[3];
sx q[3];
rz(-1.4502757) q[3];
sx q[3];
rz(-2.4108346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47736436) q[0];
sx q[0];
rz(-1.0034765) q[0];
sx q[0];
rz(3.0526414) q[0];
rz(-1.3720007) q[1];
sx q[1];
rz(-1.4274495) q[1];
sx q[1];
rz(1.1044097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1386344) q[0];
sx q[0];
rz(-1.3413652) q[0];
sx q[0];
rz(-2.5087439) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3588793) q[2];
sx q[2];
rz(-1.9316976) q[2];
sx q[2];
rz(-0.95907839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1762878) q[1];
sx q[1];
rz(-1.4215648) q[1];
sx q[1];
rz(2.01841) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9761451) q[3];
sx q[3];
rz(-2.5962127) q[3];
sx q[3];
rz(-0.65217962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.264512) q[2];
sx q[2];
rz(-2.0161714) q[2];
sx q[2];
rz(-1.6748927) q[2];
rz(0.029021164) q[3];
sx q[3];
rz(-2.1252188) q[3];
sx q[3];
rz(0.69028729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90917176) q[0];
sx q[0];
rz(-0.066815289) q[0];
sx q[0];
rz(2.8636041) q[0];
rz(1.7104644) q[1];
sx q[1];
rz(-0.93228308) q[1];
sx q[1];
rz(0.40036449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6872266) q[0];
sx q[0];
rz(-0.80601826) q[0];
sx q[0];
rz(0.63712742) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8258926) q[2];
sx q[2];
rz(-0.9005024) q[2];
sx q[2];
rz(-0.75302659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4440205) q[1];
sx q[1];
rz(-0.50067798) q[1];
sx q[1];
rz(0.33659192) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8002162) q[3];
sx q[3];
rz(-1.2814643) q[3];
sx q[3];
rz(-0.067726243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4572738) q[2];
sx q[2];
rz(-1.532734) q[2];
sx q[2];
rz(-2.686783) q[2];
rz(1.2416035) q[3];
sx q[3];
rz(-2.1865215) q[3];
sx q[3];
rz(-0.72030592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.95551816) q[0];
sx q[0];
rz(-1.8121413) q[0];
sx q[0];
rz(0.00051001471) q[0];
rz(-0.60091248) q[1];
sx q[1];
rz(-0.83952236) q[1];
sx q[1];
rz(0.14437637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2287284) q[0];
sx q[0];
rz(-1.0236386) q[0];
sx q[0];
rz(-0.24258258) q[0];
x q[1];
rz(-1.4207441) q[2];
sx q[2];
rz(-1.320082) q[2];
sx q[2];
rz(-1.5423519) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7971273) q[1];
sx q[1];
rz(-2.4542311) q[1];
sx q[1];
rz(2.7784154) q[1];
rz(-pi) q[2];
rz(-0.5471121) q[3];
sx q[3];
rz(-1.4792076) q[3];
sx q[3];
rz(1.509059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1986177) q[2];
sx q[2];
rz(-1.4451507) q[2];
sx q[2];
rz(-2.1790738) q[2];
rz(1.6019542) q[3];
sx q[3];
rz(-1.736015) q[3];
sx q[3];
rz(-0.34077728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2365504) q[0];
sx q[0];
rz(-1.4555229) q[0];
sx q[0];
rz(-1.9566253) q[0];
rz(-2.9153337) q[1];
sx q[1];
rz(-0.87892756) q[1];
sx q[1];
rz(-2.102898) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5796367) q[0];
sx q[0];
rz(-0.59369864) q[0];
sx q[0];
rz(0.14621347) q[0];
rz(-2.4657927) q[2];
sx q[2];
rz(-0.86188176) q[2];
sx q[2];
rz(-1.0647286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9925527) q[1];
sx q[1];
rz(-0.91616183) q[1];
sx q[1];
rz(-2.4813985) q[1];
x q[2];
rz(-1.709278) q[3];
sx q[3];
rz(-2.6977728) q[3];
sx q[3];
rz(-2.6276772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.431939) q[2];
sx q[2];
rz(-2.4434872) q[2];
sx q[2];
rz(2.8254438) q[2];
rz(-1.6759253) q[3];
sx q[3];
rz(-2.0114055) q[3];
sx q[3];
rz(-1.3795615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50452152) q[0];
sx q[0];
rz(-0.9032473) q[0];
sx q[0];
rz(1.1035408) q[0];
rz(-0.48031131) q[1];
sx q[1];
rz(-0.66245285) q[1];
sx q[1];
rz(2.5441817) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7561853) q[0];
sx q[0];
rz(-0.28235897) q[0];
sx q[0];
rz(2.307111) q[0];
x q[1];
rz(1.3296606) q[2];
sx q[2];
rz(-0.84000194) q[2];
sx q[2];
rz(0.75474778) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0874071) q[1];
sx q[1];
rz(-0.30255908) q[1];
sx q[1];
rz(2.5564479) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0496063) q[3];
sx q[3];
rz(-1.2013271) q[3];
sx q[3];
rz(-1.2928158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9395113) q[2];
sx q[2];
rz(-1.5152405) q[2];
sx q[2];
rz(-1.0763947) q[2];
rz(-1.1427897) q[3];
sx q[3];
rz(-2.3484774) q[3];
sx q[3];
rz(2.6514261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48653212) q[0];
sx q[0];
rz(-1.158411) q[0];
sx q[0];
rz(3.1031188) q[0];
rz(3.0768652) q[1];
sx q[1];
rz(-1.7579301) q[1];
sx q[1];
rz(-0.23385349) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9486987) q[0];
sx q[0];
rz(-1.3854376) q[0];
sx q[0];
rz(-2.9167487) q[0];
x q[1];
rz(0.87217561) q[2];
sx q[2];
rz(-0.46469122) q[2];
sx q[2];
rz(-1.2810436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.95215248) q[1];
sx q[1];
rz(-2.2169211) q[1];
sx q[1];
rz(1.6829856) q[1];
rz(-3.0117118) q[3];
sx q[3];
rz(-1.6280773) q[3];
sx q[3];
rz(-0.86806017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.48876277) q[2];
sx q[2];
rz(-1.5583928) q[2];
sx q[2];
rz(-0.59824198) q[2];
rz(0.11387842) q[3];
sx q[3];
rz(-1.3860393) q[3];
sx q[3];
rz(2.2940476) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7365731) q[0];
sx q[0];
rz(-0.84252715) q[0];
sx q[0];
rz(-0.2463499) q[0];
rz(1.8440638) q[1];
sx q[1];
rz(-2.0228701) q[1];
sx q[1];
rz(-0.49682239) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56820541) q[0];
sx q[0];
rz(-2.3620689) q[0];
sx q[0];
rz(-2.5318092) q[0];
rz(2.970817) q[2];
sx q[2];
rz(-2.429212) q[2];
sx q[2];
rz(0.26584372) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0441974) q[1];
sx q[1];
rz(-0.5438416) q[1];
sx q[1];
rz(1.103765) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15487352) q[3];
sx q[3];
rz(-1.0836156) q[3];
sx q[3];
rz(-0.39184141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42314998) q[2];
sx q[2];
rz(-2.6079874) q[2];
sx q[2];
rz(1.9971087) q[2];
rz(-1.9874969) q[3];
sx q[3];
rz(-1.4242947) q[3];
sx q[3];
rz(-2.9437039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4474739) q[0];
sx q[0];
rz(-0.23451528) q[0];
sx q[0];
rz(3.0754454) q[0];
rz(1.9937531) q[1];
sx q[1];
rz(-1.7639953) q[1];
sx q[1];
rz(-0.60417169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8291289) q[0];
sx q[0];
rz(-1.365571) q[0];
sx q[0];
rz(-0.28136307) q[0];
rz(-pi) q[1];
rz(2.6148817) q[2];
sx q[2];
rz(-1.4350495) q[2];
sx q[2];
rz(1.900224) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79163247) q[1];
sx q[1];
rz(-2.5436328) q[1];
sx q[1];
rz(2.7954742) q[1];
rz(-pi) q[2];
rz(-2.484453) q[3];
sx q[3];
rz(-2.1517188) q[3];
sx q[3];
rz(2.7681153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26312795) q[2];
sx q[2];
rz(-0.85004127) q[2];
sx q[2];
rz(-0.43295941) q[2];
rz(1.2290907) q[3];
sx q[3];
rz(-1.9357598) q[3];
sx q[3];
rz(-1.8142726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94806725) q[0];
sx q[0];
rz(-1.075241) q[0];
sx q[0];
rz(1.2731592) q[0];
rz(-0.46514568) q[1];
sx q[1];
rz(-1.3659313) q[1];
sx q[1];
rz(2.926631) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33782712) q[0];
sx q[0];
rz(-1.318232) q[0];
sx q[0];
rz(-2.4830677) q[0];
rz(2.0029055) q[2];
sx q[2];
rz(-1.0800763) q[2];
sx q[2];
rz(-3.1068206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4196288) q[1];
sx q[1];
rz(-0.89744324) q[1];
sx q[1];
rz(-2.2339905) q[1];
rz(-pi) q[2];
rz(1.8278024) q[3];
sx q[3];
rz(-2.1177835) q[3];
sx q[3];
rz(-1.5316602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2994069) q[2];
sx q[2];
rz(-0.81373787) q[2];
sx q[2];
rz(0.95883933) q[2];
rz(-1.1622608) q[3];
sx q[3];
rz(-1.4872888) q[3];
sx q[3];
rz(1.6963262) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5398298) q[0];
sx q[0];
rz(-1.5400664) q[0];
sx q[0];
rz(-1.6590317) q[0];
rz(-2.4330347) q[1];
sx q[1];
rz(-2.9500912) q[1];
sx q[1];
rz(-0.74831829) q[1];
rz(2.6513349) q[2];
sx q[2];
rz(-0.62050642) q[2];
sx q[2];
rz(-1.4668589) q[2];
rz(-0.21227588) q[3];
sx q[3];
rz(-2.4083781) q[3];
sx q[3];
rz(2.0177672) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
