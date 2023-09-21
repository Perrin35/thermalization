OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(1.9265494) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44170609) q[0];
sx q[0];
rz(-0.21472782) q[0];
sx q[0];
rz(1.0617274) q[0];
x q[1];
rz(-2.7772929) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(0.42687624) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.58303761) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(2.1117044) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9925572) q[3];
sx q[3];
rz(-1.8334853) q[3];
sx q[3];
rz(-1.5955773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3540196) q[2];
sx q[2];
rz(-0.95280567) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(-0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(-2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(2.8027957) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3829271) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(-1.7090319) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73866663) q[2];
sx q[2];
rz(-1.7324442) q[2];
sx q[2];
rz(-0.48036286) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6290366) q[1];
sx q[1];
rz(-1.3965544) q[1];
sx q[1];
rz(2.0577355) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48861309) q[3];
sx q[3];
rz(-2.6693137) q[3];
sx q[3];
rz(2.3924535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(2.4831333) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(0.69491274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(0.43310305) q[0];
rz(1.9494879) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(2.5862397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-2.0534678) q[0];
sx q[0];
rz(2.32248) q[0];
x q[1];
rz(-2.8279282) q[2];
sx q[2];
rz(-2.1137538) q[2];
sx q[2];
rz(-2.0344337) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66545031) q[1];
sx q[1];
rz(-0.94201554) q[1];
sx q[1];
rz(-2.9412342) q[1];
rz(-pi) q[2];
rz(-1.025612) q[3];
sx q[3];
rz(-1.3979997) q[3];
sx q[3];
rz(1.2196513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4042523) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(-2.897443) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(-0.81480169) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6738621) q[0];
sx q[0];
rz(-0.47463372) q[0];
sx q[0];
rz(-0.69068308) q[0];
rz(2.7093676) q[2];
sx q[2];
rz(-2.4556293) q[2];
sx q[2];
rz(1.9908817) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5863122) q[1];
sx q[1];
rz(-1.781207) q[1];
sx q[1];
rz(-1.8754688) q[1];
x q[2];
rz(-1.049794) q[3];
sx q[3];
rz(-1.5265326) q[3];
sx q[3];
rz(-1.1417768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-2.9684084) q[2];
rz(-2.611768) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.859905) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(1.2777404) q[1];
sx q[1];
rz(-2.3294096) q[1];
sx q[1];
rz(-3.0854991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1144975) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(-2.5278805) q[0];
rz(-pi) q[1];
rz(-2.7943139) q[2];
sx q[2];
rz(-1.1355073) q[2];
sx q[2];
rz(-1.6821282) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5406815) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(-1.5555698) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.089245307) q[3];
sx q[3];
rz(-2.1292994) q[3];
sx q[3];
rz(2.9104779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-0.43219217) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(-1.4661219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0284001) q[0];
sx q[0];
rz(-0.88554651) q[0];
sx q[0];
rz(-2.4940441) q[0];
rz(1.2619069) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19676767) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(1.0374271) q[0];
rz(-0.94021057) q[2];
sx q[2];
rz(-1.5988837) q[2];
sx q[2];
rz(2.0036151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8050025) q[1];
sx q[1];
rz(-1.0953566) q[1];
sx q[1];
rz(-2.5617983) q[1];
rz(-2.9783863) q[3];
sx q[3];
rz(-1.2833793) q[3];
sx q[3];
rz(-2.4678469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59297562) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-2.9086869) q[0];
sx q[0];
rz(0.74321157) q[0];
rz(-1.6339533) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(2.5315703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7689432) q[0];
sx q[0];
rz(-2.5805051) q[0];
sx q[0];
rz(0.48604301) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4201944) q[2];
sx q[2];
rz(-2.4820231) q[2];
sx q[2];
rz(1.9285551) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.89027379) q[1];
sx q[1];
rz(-0.90369019) q[1];
sx q[1];
rz(0.031884738) q[1];
x q[2];
rz(1.1915019) q[3];
sx q[3];
rz(-1.3790352) q[3];
sx q[3];
rz(-2.3281043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.6992133) q[2];
sx q[2];
rz(-0.91840333) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-0.38890719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.780705) q[0];
sx q[0];
rz(-0.66910678) q[0];
sx q[0];
rz(-1.5135182) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37305957) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(-1.2241227) q[0];
rz(-pi) q[1];
rz(-1.893115) q[2];
sx q[2];
rz(-1.467448) q[2];
sx q[2];
rz(-1.5600187) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.089162) q[1];
sx q[1];
rz(-2.6658635) q[1];
sx q[1];
rz(0.22389852) q[1];
rz(-pi) q[2];
rz(1.1235808) q[3];
sx q[3];
rz(-1.8641324) q[3];
sx q[3];
rz(-0.48188996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(0.68391189) q[2];
rz(1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(-1.7470523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972315) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(-0.99705237) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-0.7448147) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574402) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-1.1293344) q[0];
rz(-pi) q[1];
rz(-0.3955598) q[2];
sx q[2];
rz(-2.308508) q[2];
sx q[2];
rz(-1.8849444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77764952) q[1];
sx q[1];
rz(-1.3476106) q[1];
sx q[1];
rz(1.9294444) q[1];
rz(-0.78299384) q[3];
sx q[3];
rz(-1.8564312) q[3];
sx q[3];
rz(1.0234969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1054489) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.9469117) q[2];
rz(2.1448994) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(-0.99635807) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-2.7375896) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.4844091) q[1];
sx q[1];
rz(1.1709447) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02689657) q[0];
sx q[0];
rz(-0.42744246) q[0];
sx q[0];
rz(-0.052274152) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77565907) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(1.0722216) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3552637) q[1];
sx q[1];
rz(-1.5621645) q[1];
sx q[1];
rz(1.5934056) q[1];
x q[2];
rz(1.2788494) q[3];
sx q[3];
rz(-1.1598831) q[3];
sx q[3];
rz(1.3565855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6955473) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(0.5919624) q[2];
rz(2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.6177572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(0.099427632) q[1];
sx q[1];
rz(-1.8933404) q[1];
sx q[1];
rz(1.0642687) q[1];
rz(-0.912491) q[2];
sx q[2];
rz(-2.1952663) q[2];
sx q[2];
rz(-1.2637539) q[2];
rz(1.2373274) q[3];
sx q[3];
rz(-1.5860535) q[3];
sx q[3];
rz(0.64762583) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
