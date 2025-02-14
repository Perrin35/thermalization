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
rz(-2.0909042) q[1];
sx q[1];
rz(4.8839999) q[1];
sx q[1];
rz(13.304741) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90621072) q[0];
sx q[0];
rz(-0.28774187) q[0];
sx q[0];
rz(-2.5830027) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9264439) q[2];
sx q[2];
rz(-2.7251149) q[2];
sx q[2];
rz(-0.53336018) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.112632) q[1];
sx q[1];
rz(-0.46507177) q[1];
sx q[1];
rz(-2.8264224) q[1];
rz(0.61951903) q[3];
sx q[3];
rz(-1.8219007) q[3];
sx q[3];
rz(-2.8015603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38306132) q[2];
sx q[2];
rz(-0.28306511) q[2];
sx q[2];
rz(2.7281249) q[2];
rz(2.5773898) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(-1.1845425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41459945) q[0];
sx q[0];
rz(-1.0597543) q[0];
sx q[0];
rz(-3.0376036) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(0.41735059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433421) q[0];
sx q[0];
rz(-0.62055991) q[0];
sx q[0];
rz(-1.0360361) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9345788) q[2];
sx q[2];
rz(-1.5593537) q[2];
sx q[2];
rz(-1.6506843) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2552143) q[1];
sx q[1];
rz(-2.8357949) q[1];
sx q[1];
rz(1.4319819) q[1];
rz(-pi) q[2];
rz(0.82666918) q[3];
sx q[3];
rz(-2.429768) q[3];
sx q[3];
rz(-2.8579808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.61791164) q[2];
sx q[2];
rz(-1.0289501) q[2];
sx q[2];
rz(0.20393142) q[2];
rz(-1.5150874) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(-1.6319857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5390891) q[0];
sx q[0];
rz(-1.4668377) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(-1.0771982) q[1];
sx q[1];
rz(-2.4927683) q[1];
sx q[1];
rz(0.55346742) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9170415) q[0];
sx q[0];
rz(-1.5603491) q[0];
sx q[0];
rz(-2.4071669) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24793779) q[2];
sx q[2];
rz(-0.72681475) q[2];
sx q[2];
rz(3.0023129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.043442) q[1];
sx q[1];
rz(-0.69313184) q[1];
sx q[1];
rz(0.56482238) q[1];
x q[2];
rz(2.9514246) q[3];
sx q[3];
rz(-1.7432508) q[3];
sx q[3];
rz(1.2631388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0674949) q[2];
sx q[2];
rz(-0.40253887) q[2];
sx q[2];
rz(-1.2158166) q[2];
rz(1.0007693) q[3];
sx q[3];
rz(-1.3832904) q[3];
sx q[3];
rz(1.2522987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(0.85553402) q[0];
sx q[0];
rz(-1.0581886) q[0];
sx q[0];
rz(1.704294) q[0];
rz(-1.8057436) q[1];
sx q[1];
rz(-0.62594405) q[1];
sx q[1];
rz(2.6864973) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9055331) q[0];
sx q[0];
rz(-0.62232557) q[0];
sx q[0];
rz(0.13747352) q[0];
rz(1.1983775) q[2];
sx q[2];
rz(-1.7657868) q[2];
sx q[2];
rz(0.41620987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1170821) q[1];
sx q[1];
rz(-2.7221537) q[1];
sx q[1];
rz(2.504435) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7198445) q[3];
sx q[3];
rz(-2.9136326) q[3];
sx q[3];
rz(2.3123119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2745634) q[2];
sx q[2];
rz(-0.44963351) q[2];
sx q[2];
rz(-2.3095798) q[2];
rz(2.4882107) q[3];
sx q[3];
rz(-0.91975206) q[3];
sx q[3];
rz(0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5117383) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(2.4526556) q[0];
rz(-0.2116994) q[1];
sx q[1];
rz(-1.4392821) q[1];
sx q[1];
rz(0.45417085) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4409597) q[0];
sx q[0];
rz(-2.0122177) q[0];
sx q[0];
rz(0.9671797) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93438678) q[2];
sx q[2];
rz(-1.8737405) q[2];
sx q[2];
rz(-2.7222939) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28605697) q[1];
sx q[1];
rz(-2.2814676) q[1];
sx q[1];
rz(-2.2362806) q[1];
x q[2];
rz(1.8332043) q[3];
sx q[3];
rz(-1.7324362) q[3];
sx q[3];
rz(1.9899881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5357369) q[2];
sx q[2];
rz(-1.7822632) q[2];
sx q[2];
rz(-0.5109171) q[2];
rz(1.9262675) q[3];
sx q[3];
rz(-1.5769703) q[3];
sx q[3];
rz(-0.63466614) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257492) q[0];
sx q[0];
rz(-0.31084335) q[0];
sx q[0];
rz(-2.6281443) q[0];
rz(0.91649857) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(1.7880012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3037864) q[0];
sx q[0];
rz(-0.45188658) q[0];
sx q[0];
rz(-2.9511098) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36717461) q[2];
sx q[2];
rz(-1.1859787) q[2];
sx q[2];
rz(-0.095331017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97670758) q[1];
sx q[1];
rz(-2.0583377) q[1];
sx q[1];
rz(-1.5406002) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6597802) q[3];
sx q[3];
rz(-1.5060194) q[3];
sx q[3];
rz(1.806206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6650271) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(2.636886) q[3];
sx q[3];
rz(-1.4387771) q[3];
sx q[3];
rz(-0.75468841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1228834) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(-1.1705742) q[0];
rz(-0.19078828) q[1];
sx q[1];
rz(-2.6759594) q[1];
sx q[1];
rz(-2.4986787) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.64775) q[0];
sx q[0];
rz(-0.28573418) q[0];
sx q[0];
rz(1.2005376) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9011074) q[2];
sx q[2];
rz(-1.3250371) q[2];
sx q[2];
rz(1.877223) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9289405) q[1];
sx q[1];
rz(-1.3319322) q[1];
sx q[1];
rz(-0.3122621) q[1];
rz(-pi) q[2];
rz(-0.40465458) q[3];
sx q[3];
rz(-1.4611562) q[3];
sx q[3];
rz(-0.10900036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-0.10991749) q[2];
sx q[2];
rz(1.9996803) q[2];
rz(-3.111908) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(0.77965411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90758816) q[0];
sx q[0];
rz(-1.1827844) q[0];
sx q[0];
rz(-1.8335861) q[0];
rz(-1.1169149) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(-2.9294779) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19973913) q[0];
sx q[0];
rz(-0.77773577) q[0];
sx q[0];
rz(1.1129802) q[0];
rz(-pi) q[1];
rz(0.30435698) q[2];
sx q[2];
rz(-1.3469014) q[2];
sx q[2];
rz(0.18098132) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6740283) q[1];
sx q[1];
rz(-2.0851567) q[1];
sx q[1];
rz(-0.5089203) q[1];
rz(-pi) q[2];
rz(-0.85101012) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(-1.855576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1954631) q[2];
sx q[2];
rz(-1.7886432) q[2];
sx q[2];
rz(1.8606261) q[2];
rz(1.403275) q[3];
sx q[3];
rz(-0.91126982) q[3];
sx q[3];
rz(-0.72428552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69287777) q[0];
sx q[0];
rz(-2.4249478) q[0];
sx q[0];
rz(0.88248673) q[0];
rz(0.29300434) q[1];
sx q[1];
rz(-0.92331433) q[1];
sx q[1];
rz(-2.015347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.654713) q[0];
sx q[0];
rz(-2.3695787) q[0];
sx q[0];
rz(2.604542) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5669021) q[2];
sx q[2];
rz(-0.97087395) q[2];
sx q[2];
rz(1.1614161) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80690224) q[1];
sx q[1];
rz(-2.370939) q[1];
sx q[1];
rz(-1.6477963) q[1];
x q[2];
rz(-0.65139095) q[3];
sx q[3];
rz(-2.6462501) q[3];
sx q[3];
rz(-1.9598531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5249411) q[2];
sx q[2];
rz(-2.5280819) q[2];
sx q[2];
rz(-2.9743527) q[2];
rz(0.031115726) q[3];
sx q[3];
rz(-1.9626706) q[3];
sx q[3];
rz(-2.4955366) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50877082) q[0];
sx q[0];
rz(-1.8080067) q[0];
sx q[0];
rz(2.3311145) q[0];
rz(-1.3088016) q[1];
sx q[1];
rz(-0.93593132) q[1];
sx q[1];
rz(-3.0442309) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33695147) q[0];
sx q[0];
rz(-1.7367678) q[0];
sx q[0];
rz(-2.354855) q[0];
rz(-pi) q[1];
rz(0.028848666) q[2];
sx q[2];
rz(-1.9518153) q[2];
sx q[2];
rz(2.888607) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.88579455) q[1];
sx q[1];
rz(-0.58572873) q[1];
sx q[1];
rz(-3.0540007) q[1];
rz(-pi) q[2];
rz(1.8382357) q[3];
sx q[3];
rz(-1.6839661) q[3];
sx q[3];
rz(-1.8281931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8436766) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(-1.2751014) q[2];
rz(-0.44025931) q[3];
sx q[3];
rz(-2.2831254) q[3];
sx q[3];
rz(-0.27537235) q[3];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020582747) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(0.62200017) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(-0.95876454) q[2];
sx q[2];
rz(-0.13560451) q[2];
sx q[2];
rz(2.3252631) q[2];
rz(1.7429327) q[3];
sx q[3];
rz(-2.2071916) q[3];
sx q[3];
rz(-1.9194736) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
