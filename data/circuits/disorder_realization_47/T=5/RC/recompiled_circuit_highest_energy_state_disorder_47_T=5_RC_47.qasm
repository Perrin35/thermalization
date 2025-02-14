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
rz(1.6506305) q[0];
sx q[0];
rz(-1.4181674) q[0];
sx q[0];
rz(11.081063) q[0];
rz(1.0506884) q[1];
sx q[1];
rz(-1.7424072) q[1];
sx q[1];
rz(-0.73837003) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0169058) q[0];
sx q[0];
rz(-1.7217741) q[0];
sx q[0];
rz(0.24589234) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9638871) q[2];
sx q[2];
rz(-1.7121268) q[2];
sx q[2];
rz(-1.7766952) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.112632) q[1];
sx q[1];
rz(-2.6765209) q[1];
sx q[1];
rz(-2.8264224) q[1];
rz(-pi) q[2];
rz(-0.61951903) q[3];
sx q[3];
rz(-1.3196919) q[3];
sx q[3];
rz(0.34003231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-0.28306511) q[2];
sx q[2];
rz(2.7281249) q[2];
rz(2.5773898) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(1.9570501) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7269932) q[0];
sx q[0];
rz(-2.0818384) q[0];
sx q[0];
rz(3.0376036) q[0];
rz(-3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(-0.41735059) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433421) q[0];
sx q[0];
rz(-0.62055991) q[0];
sx q[0];
rz(2.1055566) q[0];
x q[1];
rz(-0.012243791) q[2];
sx q[2];
rz(-1.2070388) q[2];
sx q[2];
rz(3.0573483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0318533) q[1];
sx q[1];
rz(-1.268034) q[1];
sx q[1];
rz(-3.0979373) q[1];
rz(2.1363229) q[3];
sx q[3];
rz(-1.1124753) q[3];
sx q[3];
rz(-1.895973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61791164) q[2];
sx q[2];
rz(-2.1126426) q[2];
sx q[2];
rz(-0.20393142) q[2];
rz(1.6265053) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(-1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5390891) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(-2.7432826) q[0];
rz(-1.0771982) q[1];
sx q[1];
rz(-0.64882433) q[1];
sx q[1];
rz(2.5881252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22455118) q[0];
sx q[0];
rz(-1.5812435) q[0];
sx q[0];
rz(0.73442572) q[0];
rz(-2.4301694) q[2];
sx q[2];
rz(-1.4069948) q[2];
sx q[2];
rz(1.244551) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0981507) q[1];
sx q[1];
rz(-0.69313184) q[1];
sx q[1];
rz(-2.5767703) q[1];
rz(-pi) q[2];
rz(2.3970277) q[3];
sx q[3];
rz(-2.8855763) q[3];
sx q[3];
rz(-2.721173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.074097721) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(-1.925776) q[2];
rz(1.0007693) q[3];
sx q[3];
rz(-1.3832904) q[3];
sx q[3];
rz(-1.8892939) q[3];
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
rz(-pi/2) q[3];
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
rz(-2.2860586) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(1.4372987) q[0];
rz(1.8057436) q[1];
sx q[1];
rz(-0.62594405) q[1];
sx q[1];
rz(0.45509532) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0741827) q[0];
sx q[0];
rz(-0.95522308) q[0];
sx q[0];
rz(1.4727946) q[0];
rz(-pi) q[1];
rz(1.0735157) q[2];
sx q[2];
rz(-0.4182564) q[2];
sx q[2];
rz(-2.4474395) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1170821) q[1];
sx q[1];
rz(-2.7221537) q[1];
sx q[1];
rz(-0.6371577) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4761127) q[3];
sx q[3];
rz(-1.3631184) q[3];
sx q[3];
rz(0.39765795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2745634) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(-2.3095798) q[2];
rz(-2.4882107) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298544) q[0];
sx q[0];
rz(-1.8269704) q[0];
sx q[0];
rz(-0.6889371) q[0];
rz(0.2116994) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(-0.45417085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1566166) q[0];
sx q[0];
rz(-2.1097221) q[0];
sx q[0];
rz(2.6205553) q[0];
rz(-pi) q[1];
rz(-1.0866223) q[2];
sx q[2];
rz(-2.4459029) q[2];
sx q[2];
rz(-1.5351968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8555357) q[1];
sx q[1];
rz(-2.2814676) q[1];
sx q[1];
rz(-2.2362806) q[1];
x q[2];
rz(2.1319785) q[3];
sx q[3];
rz(-2.8343763) q[3];
sx q[3];
rz(3.0209783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6058558) q[2];
sx q[2];
rz(-1.3593295) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(-1.2153252) q[3];
sx q[3];
rz(-1.5769703) q[3];
sx q[3];
rz(2.5069265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11584347) q[0];
sx q[0];
rz(-2.8307493) q[0];
sx q[0];
rz(-2.6281443) q[0];
rz(0.91649857) q[1];
sx q[1];
rz(-1.1767574) q[1];
sx q[1];
rz(-1.7880012) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3037864) q[0];
sx q[0];
rz(-0.45188658) q[0];
sx q[0];
rz(2.9511098) q[0];
rz(-pi) q[1];
rz(1.1613834) q[2];
sx q[2];
rz(-1.231603) q[2];
sx q[2];
rz(1.3320751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1648851) q[1];
sx q[1];
rz(-2.0583377) q[1];
sx q[1];
rz(-1.5406002) q[1];
rz(-0.065033536) q[3];
sx q[3];
rz(-1.659593) q[3];
sx q[3];
rz(2.9119583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4765656) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(-3.108007) q[2];
rz(0.50470662) q[3];
sx q[3];
rz(-1.4387771) q[3];
sx q[3];
rz(-2.3869042) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228834) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(-1.1705742) q[0];
rz(0.19078828) q[1];
sx q[1];
rz(-2.6759594) q[1];
sx q[1];
rz(-0.64291397) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87827728) q[0];
sx q[0];
rz(-1.8366792) q[0];
sx q[0];
rz(-3.0356867) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91251164) q[2];
sx q[2];
rz(-0.40900074) q[2];
sx q[2];
rz(2.8308979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1509779) q[1];
sx q[1];
rz(-2.7508368) q[1];
sx q[1];
rz(2.4713466) q[1];
rz(-pi) q[2];
rz(2.7369381) q[3];
sx q[3];
rz(-1.4611562) q[3];
sx q[3];
rz(3.0325923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.99870318) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(1.9996803) q[2];
rz(-3.111908) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(-2.3619385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2340045) q[0];
sx q[0];
rz(-1.9588082) q[0];
sx q[0];
rz(1.8335861) q[0];
rz(-1.1169149) q[1];
sx q[1];
rz(-1.4968548) q[1];
sx q[1];
rz(2.9294779) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9418535) q[0];
sx q[0];
rz(-2.3638569) q[0];
sx q[0];
rz(2.0286125) q[0];
rz(-pi) q[1];
rz(1.3364974) q[2];
sx q[2];
rz(-1.8673225) q[2];
sx q[2];
rz(1.4594452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6740283) q[1];
sx q[1];
rz(-2.0851567) q[1];
sx q[1];
rz(-0.5089203) q[1];
x q[2];
rz(-0.66069938) q[3];
sx q[3];
rz(-0.96008077) q[3];
sx q[3];
rz(0.92044059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9461296) q[2];
sx q[2];
rz(-1.3529494) q[2];
sx q[2];
rz(1.8606261) q[2];
rz(1.7383176) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(-2.4173071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4487149) q[0];
sx q[0];
rz(-0.71664482) q[0];
sx q[0];
rz(2.2591059) q[0];
rz(2.8485883) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(-1.1262456) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4868797) q[0];
sx q[0];
rz(-2.3695787) q[0];
sx q[0];
rz(2.604542) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5416667) q[2];
sx q[2];
rz(-1.5675822) q[2];
sx q[2];
rz(-0.41157882) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3346904) q[1];
sx q[1];
rz(-0.77065361) q[1];
sx q[1];
rz(1.4937964) q[1];
rz(-pi) q[2];
rz(0.4057986) q[3];
sx q[3];
rz(-1.8631336) q[3];
sx q[3];
rz(-2.1617011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61665159) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(2.9743527) q[2];
rz(-0.031115726) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(2.4955366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328218) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(-2.3311145) q[0];
rz(-1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(0.097361758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7435849) q[0];
sx q[0];
rz(-0.79774374) q[0];
sx q[0];
rz(-1.8037075) q[0];
x q[1];
rz(1.951959) q[2];
sx q[2];
rz(-1.544017) q[2];
sx q[2];
rz(-1.8130515) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78079911) q[1];
sx q[1];
rz(-2.153984) q[1];
sx q[1];
rz(1.5128283) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11730365) q[3];
sx q[3];
rz(-1.3051093) q[3];
sx q[3];
rz(-0.22646389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.29791609) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(-1.8664912) q[2];
rz(-0.44025931) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(0.27537235) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-0.62200017) q[1];
sx q[1];
rz(-1.2983464) q[1];
sx q[1];
rz(1.3141528) q[1];
rz(-1.6820108) q[2];
sx q[2];
rz(-1.6485452) q[2];
sx q[2];
rz(-2.9948276) q[2];
rz(2.9138184) q[3];
sx q[3];
rz(-2.4854599) q[3];
sx q[3];
rz(1.5067185) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
