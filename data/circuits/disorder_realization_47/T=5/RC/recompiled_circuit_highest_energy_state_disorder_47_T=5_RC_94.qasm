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
rz(-1.4909622) q[0];
sx q[0];
rz(-1.7234252) q[0];
sx q[0];
rz(-1.6562847) q[0];
rz(1.0506884) q[1];
sx q[1];
rz(-1.7424072) q[1];
sx q[1];
rz(-0.73837003) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90621072) q[0];
sx q[0];
rz(-2.8538508) q[0];
sx q[0];
rz(-0.55858992) q[0];
rz(-pi) q[1];
rz(-1.2151488) q[2];
sx q[2];
rz(-2.7251149) q[2];
sx q[2];
rz(2.6082325) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.112632) q[1];
sx q[1];
rz(-2.6765209) q[1];
sx q[1];
rz(-0.3151703) q[1];
rz(-pi) q[2];
rz(1.876023) q[3];
sx q[3];
rz(-2.1681227) q[3];
sx q[3];
rz(1.0553774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7585313) q[2];
sx q[2];
rz(-2.8585275) q[2];
sx q[2];
rz(-0.4134678) q[2];
rz(0.56420285) q[3];
sx q[3];
rz(-1.6744303) q[3];
sx q[3];
rz(1.1845425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7269932) q[0];
sx q[0];
rz(-2.0818384) q[0];
sx q[0];
rz(-0.10398908) q[0];
rz(3.0005786) q[1];
sx q[1];
rz(-2.4555989) q[1];
sx q[1];
rz(0.41735059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.418103) q[0];
sx q[0];
rz(-1.8716629) q[0];
sx q[0];
rz(2.1221493) q[0];
rz(1.6029458) q[2];
sx q[2];
rz(-2.7776383) q[2];
sx q[2];
rz(0.049843069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1097393) q[1];
sx q[1];
rz(-1.268034) q[1];
sx q[1];
rz(-0.043655386) q[1];
x q[2];
rz(-1.0052698) q[3];
sx q[3];
rz(-2.0291174) q[3];
sx q[3];
rz(-1.2456196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.61791164) q[2];
sx q[2];
rz(-1.0289501) q[2];
sx q[2];
rz(-0.20393142) q[2];
rz(-1.6265053) q[3];
sx q[3];
rz(-0.47658673) q[3];
sx q[3];
rz(-1.509607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6025036) q[0];
sx q[0];
rz(-1.6747549) q[0];
sx q[0];
rz(0.3983101) q[0];
rz(-2.0643945) q[1];
sx q[1];
rz(-2.4927683) q[1];
sx q[1];
rz(2.5881252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.80478) q[0];
sx q[0];
rz(-2.3051728) q[0];
sx q[0];
rz(-1.5848716) q[0];
rz(-pi) q[1];
rz(-0.24793779) q[2];
sx q[2];
rz(-2.4147779) q[2];
sx q[2];
rz(3.0023129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.043442) q[1];
sx q[1];
rz(-2.4484608) q[1];
sx q[1];
rz(2.5767703) q[1];
rz(1.7463527) q[3];
sx q[3];
rz(-1.3834828) q[3];
sx q[3];
rz(-2.8009149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.074097721) q[2];
sx q[2];
rz(-2.7390538) q[2];
sx q[2];
rz(-1.925776) q[2];
rz(-1.0007693) q[3];
sx q[3];
rz(-1.3832904) q[3];
sx q[3];
rz(1.8892939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85553402) q[0];
sx q[0];
rz(-2.0834041) q[0];
sx q[0];
rz(1.4372987) q[0];
rz(-1.3358491) q[1];
sx q[1];
rz(-0.62594405) q[1];
sx q[1];
rz(0.45509532) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44667654) q[0];
sx q[0];
rz(-1.4908264) q[0];
sx q[0];
rz(0.6178426) q[0];
rz(2.0680769) q[2];
sx q[2];
rz(-2.7233363) q[2];
sx q[2];
rz(-2.4474395) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1170821) q[1];
sx q[1];
rz(-0.41943892) q[1];
sx q[1];
rz(0.6371577) q[1];
rz(-pi) q[2];
rz(-0.20858553) q[3];
sx q[3];
rz(-1.6634395) q[3];
sx q[3];
rz(1.988033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2745634) q[2];
sx q[2];
rz(-2.6919591) q[2];
sx q[2];
rz(-2.3095798) q[2];
rz(0.653382) q[3];
sx q[3];
rz(-2.2218406) q[3];
sx q[3];
rz(0.67460361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298544) q[0];
sx q[0];
rz(-1.3146223) q[0];
sx q[0];
rz(-2.4526556) q[0];
rz(-0.2116994) q[1];
sx q[1];
rz(-1.7023106) q[1];
sx q[1];
rz(2.6874218) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70063299) q[0];
sx q[0];
rz(-2.0122177) q[0];
sx q[0];
rz(-2.174413) q[0];
rz(-pi) q[1];
rz(-0.93438678) q[2];
sx q[2];
rz(-1.8737405) q[2];
sx q[2];
rz(-2.7222939) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28605697) q[1];
sx q[1];
rz(-2.2814676) q[1];
sx q[1];
rz(0.90531207) q[1];
x q[2];
rz(-2.1319785) q[3];
sx q[3];
rz(-2.8343763) q[3];
sx q[3];
rz(-3.0209783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6058558) q[2];
sx q[2];
rz(-1.7822632) q[2];
sx q[2];
rz(-2.6306756) q[2];
rz(1.2153252) q[3];
sx q[3];
rz(-1.5769703) q[3];
sx q[3];
rz(0.63466614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11584347) q[0];
sx q[0];
rz(-0.31084335) q[0];
sx q[0];
rz(-2.6281443) q[0];
rz(-0.91649857) q[1];
sx q[1];
rz(-1.9648353) q[1];
sx q[1];
rz(-1.7880012) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62666639) q[0];
sx q[0];
rz(-2.0139222) q[0];
sx q[0];
rz(-1.6624381) q[0];
x q[1];
rz(-0.36717461) q[2];
sx q[2];
rz(-1.9556139) q[2];
sx q[2];
rz(3.0462616) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1648851) q[1];
sx q[1];
rz(-2.0583377) q[1];
sx q[1];
rz(1.5406002) q[1];
rz(-pi) q[2];
rz(1.6597802) q[3];
sx q[3];
rz(-1.6355733) q[3];
sx q[3];
rz(1.3353867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6650271) q[2];
sx q[2];
rz(-2.5883784) q[2];
sx q[2];
rz(0.033585699) q[2];
rz(0.50470662) q[3];
sx q[3];
rz(-1.7028156) q[3];
sx q[3];
rz(-0.75468841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.018709239) q[0];
sx q[0];
rz(-1.2241192) q[0];
sx q[0];
rz(-1.9710185) q[0];
rz(0.19078828) q[1];
sx q[1];
rz(-0.46563322) q[1];
sx q[1];
rz(-2.4986787) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49384269) q[0];
sx q[0];
rz(-0.28573418) q[0];
sx q[0];
rz(-1.2005376) q[0];
rz(-pi) q[1];
rz(1.9011074) q[2];
sx q[2];
rz(-1.3250371) q[2];
sx q[2];
rz(1.2643697) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1509779) q[1];
sx q[1];
rz(-2.7508368) q[1];
sx q[1];
rz(-2.4713466) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6899818) q[3];
sx q[3];
rz(-1.9728807) q[3];
sx q[3];
rz(1.5086255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1428895) q[2];
sx q[2];
rz(-3.0316752) q[2];
sx q[2];
rz(1.1419123) q[2];
rz(-3.111908) q[3];
sx q[3];
rz(-1.5342865) q[3];
sx q[3];
rz(-2.3619385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(2.0246778) q[1];
sx q[1];
rz(-1.6447379) q[1];
sx q[1];
rz(-2.9294779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19973913) q[0];
sx q[0];
rz(-2.3638569) q[0];
sx q[0];
rz(-2.0286125) q[0];
x q[1];
rz(1.8050952) q[2];
sx q[2];
rz(-1.2742701) q[2];
sx q[2];
rz(1.4594452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4675644) q[1];
sx q[1];
rz(-2.0851567) q[1];
sx q[1];
rz(2.6326724) q[1];
x q[2];
rz(0.85101012) q[3];
sx q[3];
rz(-2.27423) q[3];
sx q[3];
rz(1.855576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9461296) q[2];
sx q[2];
rz(-1.3529494) q[2];
sx q[2];
rz(-1.2809666) q[2];
rz(-1.403275) q[3];
sx q[3];
rz(-2.2303228) q[3];
sx q[3];
rz(-0.72428552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
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
rz(-2.2182783) q[1];
sx q[1];
rz(2.015347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6544271) q[0];
sx q[0];
rz(-1.2058655) q[0];
sx q[0];
rz(-0.69661822) q[0];
rz(1.5669021) q[2];
sx q[2];
rz(-0.97087395) q[2];
sx q[2];
rz(1.9801766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3224015) q[1];
sx q[1];
rz(-1.5171851) q[1];
sx q[1];
rz(0.80162572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4902017) q[3];
sx q[3];
rz(-2.6462501) q[3];
sx q[3];
rz(-1.1817396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61665159) q[2];
sx q[2];
rz(-0.61351073) q[2];
sx q[2];
rz(-0.16723995) q[2];
rz(-0.031115726) q[3];
sx q[3];
rz(-1.1789221) q[3];
sx q[3];
rz(-2.4955366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50877082) q[0];
sx q[0];
rz(-1.333586) q[0];
sx q[0];
rz(0.81047812) q[0];
rz(-1.8327911) q[1];
sx q[1];
rz(-2.2056613) q[1];
sx q[1];
rz(0.097361758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705436) q[0];
sx q[0];
rz(-0.80035058) q[0];
sx q[0];
rz(0.23231028) q[0];
x q[1];
rz(-0.028848666) q[2];
sx q[2];
rz(-1.1897773) q[2];
sx q[2];
rz(2.888607) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75805) q[1];
sx q[1];
rz(-1.6191747) q[1];
sx q[1];
rz(0.58396062) q[1];
rz(3.024289) q[3];
sx q[3];
rz(-1.3051093) q[3];
sx q[3];
rz(-2.9151288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29791609) q[2];
sx q[2];
rz(-1.0138136) q[2];
sx q[2];
rz(-1.8664912) q[2];
rz(-2.7013333) q[3];
sx q[3];
rz(-0.85846725) q[3];
sx q[3];
rz(2.8662203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1210099) q[0];
sx q[0];
rz(-2.5204211) q[0];
sx q[0];
rz(-0.68516635) q[0];
rz(-2.5195925) q[1];
sx q[1];
rz(-1.8432462) q[1];
sx q[1];
rz(-1.8274399) q[1];
rz(3.0633625) q[2];
sx q[2];
rz(-1.4599192) q[2];
sx q[2];
rz(-1.4327049) q[2];
rz(-2.9138184) q[3];
sx q[3];
rz(-0.65613272) q[3];
sx q[3];
rz(-1.6348742) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
