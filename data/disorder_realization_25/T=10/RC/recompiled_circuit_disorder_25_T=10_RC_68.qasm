OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(-2.3526469) q[1];
sx q[1];
rz(0.087892428) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1577507) q[0];
sx q[0];
rz(-1.112126) q[0];
sx q[0];
rz(1.0755324) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42638875) q[2];
sx q[2];
rz(-1.083056) q[2];
sx q[2];
rz(-2.3418155) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.76074782) q[1];
sx q[1];
rz(-1.5686085) q[1];
sx q[1];
rz(-0.50427498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5929234) q[3];
sx q[3];
rz(-2.8828354) q[3];
sx q[3];
rz(2.4165137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(-1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.3585842) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.233261) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9995995) q[0];
sx q[0];
rz(-0.90786952) q[0];
sx q[0];
rz(-0.87447383) q[0];
rz(-2.9821175) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(-2.8603539) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43014363) q[1];
sx q[1];
rz(-1.5549545) q[1];
sx q[1];
rz(1.364691) q[1];
x q[2];
rz(-0.29970701) q[3];
sx q[3];
rz(-2.4518659) q[3];
sx q[3];
rz(-1.0563861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1374986) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(0.0022350524) q[2];
rz(2.3114752) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925791) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(-0.85154831) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(-3.070389) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.203513) q[0];
sx q[0];
rz(-0.31103125) q[0];
sx q[0];
rz(-0.87840338) q[0];
rz(-pi) q[1];
rz(-2.0257646) q[2];
sx q[2];
rz(-1.3340545) q[2];
sx q[2];
rz(-0.67982212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60656602) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(1.7391298) q[1];
x q[2];
rz(0.76020469) q[3];
sx q[3];
rz(-1.01902) q[3];
sx q[3];
rz(-2.9387568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3123793) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(-1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(1.4473787) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81340504) q[0];
sx q[0];
rz(-1.7338599) q[0];
sx q[0];
rz(-1.3379315) q[0];
rz(1.5201735) q[2];
sx q[2];
rz(-0.17855893) q[2];
sx q[2];
rz(-0.52390097) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4923258) q[1];
sx q[1];
rz(-2.4267303) q[1];
sx q[1];
rz(1.8148755) q[1];
rz(-pi) q[2];
rz(-1.8131687) q[3];
sx q[3];
rz(-1.216785) q[3];
sx q[3];
rz(-2.6995475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(-0.14492598) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-0.01468006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-2.4575535) q[0];
rz(1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(1.9285944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1319259) q[0];
sx q[0];
rz(-1.3103232) q[0];
sx q[0];
rz(-0.38513222) q[0];
rz(-2.331779) q[2];
sx q[2];
rz(-2.3286616) q[2];
sx q[2];
rz(0.44250689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1274384) q[1];
sx q[1];
rz(-1.6886097) q[1];
sx q[1];
rz(0.25728667) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9641987) q[3];
sx q[3];
rz(-1.3131485) q[3];
sx q[3];
rz(-0.28211668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(-2.2563289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(2.5568331) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8557381) q[0];
sx q[0];
rz(-1.8830839) q[0];
sx q[0];
rz(2.5166442) q[0];
rz(-pi) q[1];
rz(-3.1044664) q[2];
sx q[2];
rz(-0.6147487) q[2];
sx q[2];
rz(0.45072134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0271921) q[1];
sx q[1];
rz(-0.90585867) q[1];
sx q[1];
rz(-2.7726638) q[1];
x q[2];
rz(1.5201969) q[3];
sx q[3];
rz(-1.4046602) q[3];
sx q[3];
rz(-1.43464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16453234) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(-1.6284846) q[2];
rz(2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-0.62455463) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(-3.0896297) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(-0.60639492) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.279778) q[0];
sx q[0];
rz(-1.3347515) q[0];
sx q[0];
rz(2.1104913) q[0];
rz(-pi) q[1];
rz(-1.1029926) q[2];
sx q[2];
rz(-0.71873795) q[2];
sx q[2];
rz(-1.4008092) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7402621) q[1];
sx q[1];
rz(-1.6832502) q[1];
sx q[1];
rz(-0.80660352) q[1];
rz(-pi) q[2];
rz(2.2783979) q[3];
sx q[3];
rz(-2.6226225) q[3];
sx q[3];
rz(-0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(2.0243747) q[3];
sx q[3];
rz(-0.17337392) q[3];
sx q[3];
rz(1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.8023119) q[0];
sx q[0];
rz(-0.34582129) q[0];
sx q[0];
rz(-2.3876277) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65093016) q[0];
sx q[0];
rz(-2.7224446) q[0];
sx q[0];
rz(1.7463669) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67543244) q[2];
sx q[2];
rz(-0.54985148) q[2];
sx q[2];
rz(-2.7477802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7714872) q[1];
sx q[1];
rz(-0.58413726) q[1];
sx q[1];
rz(2.7258956) q[1];
rz(2.5826107) q[3];
sx q[3];
rz(-1.6056656) q[3];
sx q[3];
rz(-0.99136409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7446328) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(2.0869758) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(1.9931591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52590695) q[0];
sx q[0];
rz(-1.6427759) q[0];
sx q[0];
rz(-3.1176438) q[0];
x q[1];
rz(2.7788413) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(2.9438058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.941274) q[1];
sx q[1];
rz(-0.66220821) q[1];
sx q[1];
rz(2.3108285) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0461379) q[3];
sx q[3];
rz(-0.3534375) q[3];
sx q[3];
rz(-0.27775882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0917197) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-0.77825528) q[2];
rz(0.19566472) q[3];
sx q[3];
rz(-2.1019432) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8675999) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(0.25319779) q[0];
rz(2.4018438) q[1];
sx q[1];
rz(-0.7363798) q[1];
sx q[1];
rz(2.443312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4822599) q[0];
sx q[0];
rz(-1.2685304) q[0];
sx q[0];
rz(-1.6460101) q[0];
rz(-pi) q[1];
rz(3.1168078) q[2];
sx q[2];
rz(-1.4966655) q[2];
sx q[2];
rz(2.5075846) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89430289) q[1];
sx q[1];
rz(-1.786788) q[1];
sx q[1];
rz(-2.0386019) q[1];
rz(-pi) q[2];
rz(-1.9395589) q[3];
sx q[3];
rz(-0.59951111) q[3];
sx q[3];
rz(2.5332019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(-2.3013766) q[2];
rz(0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(-3.1124658) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-1.2226979) q[2];
sx q[2];
rz(-2.1445027) q[2];
sx q[2];
rz(2.4377844) q[2];
rz(2.2610353) q[3];
sx q[3];
rz(-1.4546483) q[3];
sx q[3];
rz(1.6402257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
