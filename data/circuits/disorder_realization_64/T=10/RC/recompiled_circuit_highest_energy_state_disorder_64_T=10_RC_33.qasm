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
rz(-1.0567226) q[0];
sx q[0];
rz(-1.280008) q[0];
sx q[0];
rz(-0.7260538) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(-1.0322303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3772349) q[0];
sx q[0];
rz(-2.4774158) q[0];
sx q[0];
rz(-2.9773877) q[0];
x q[1];
rz(-0.039041877) q[2];
sx q[2];
rz(-2.0882975) q[2];
sx q[2];
rz(1.3618748) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0192248) q[1];
sx q[1];
rz(-1.4488954) q[1];
sx q[1];
rz(2.8278964) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.80268245) q[3];
sx q[3];
rz(-1.9415874) q[3];
sx q[3];
rz(0.30457531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7008179) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(0.087892858) q[2];
rz(1.4676189) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(2.1427593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(1.4922967) q[0];
rz(2.3536033) q[1];
sx q[1];
rz(-1.9580656) q[1];
sx q[1];
rz(2.1764887) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1247281) q[0];
sx q[0];
rz(-1.6627321) q[0];
sx q[0];
rz(-0.73081907) q[0];
rz(-0.1568753) q[2];
sx q[2];
rz(-2.3986536) q[2];
sx q[2];
rz(-0.76240957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.31307855) q[1];
sx q[1];
rz(-1.113849) q[1];
sx q[1];
rz(0.20707737) q[1];
rz(-pi) q[2];
rz(1.6613668) q[3];
sx q[3];
rz(-1.1027059) q[3];
sx q[3];
rz(0.15453574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.013177055) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.817912) q[2];
rz(-0.88661083) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(2.9429341) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24432467) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(-2.4460728) q[0];
rz(2.3059402) q[1];
sx q[1];
rz(-1.7213768) q[1];
sx q[1];
rz(2.4998891) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7135237) q[0];
sx q[0];
rz(-1.7155353) q[0];
sx q[0];
rz(2.3109396) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19410816) q[2];
sx q[2];
rz(-0.72065565) q[2];
sx q[2];
rz(-1.0247165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5776331) q[1];
sx q[1];
rz(-2.0490987) q[1];
sx q[1];
rz(0.35525124) q[1];
rz(-pi) q[2];
rz(-2.4536127) q[3];
sx q[3];
rz(-1.1065346) q[3];
sx q[3];
rz(1.9339428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9782763) q[2];
sx q[2];
rz(-1.1288319) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(2.002772) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9982933) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(3.0553014) q[1];
sx q[1];
rz(-0.50420612) q[1];
sx q[1];
rz(0.48826826) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241656) q[0];
sx q[0];
rz(-2.0323951) q[0];
sx q[0];
rz(1.6469) q[0];
rz(-pi) q[1];
rz(2.5163842) q[2];
sx q[2];
rz(-0.80179384) q[2];
sx q[2];
rz(1.2259353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.77535532) q[1];
sx q[1];
rz(-1.7960494) q[1];
sx q[1];
rz(2.3083592) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9776439) q[3];
sx q[3];
rz(-1.9449807) q[3];
sx q[3];
rz(-1.0519303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6733072) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(-2.9735273) q[2];
rz(-0.42260653) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(2.9882123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76822686) q[0];
sx q[0];
rz(-2.234937) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(0.52781421) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6884228) q[0];
sx q[0];
rz(-0.97081447) q[0];
sx q[0];
rz(2.1196516) q[0];
rz(-pi) q[1];
rz(1.509769) q[2];
sx q[2];
rz(-1.5102855) q[2];
sx q[2];
rz(2.4946545) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76705383) q[1];
sx q[1];
rz(-0.77195814) q[1];
sx q[1];
rz(2.1943387) q[1];
rz(-pi) q[2];
rz(0.47021535) q[3];
sx q[3];
rz(-0.99567669) q[3];
sx q[3];
rz(1.0660389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8268383) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(2.4763988) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(0.78072602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219532) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(-0.39805472) q[0];
rz(1.4705426) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(-0.7801396) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3624293) q[0];
sx q[0];
rz(-2.0423959) q[0];
sx q[0];
rz(2.7410024) q[0];
rz(-pi) q[1];
rz(1.3684811) q[2];
sx q[2];
rz(-1.8087808) q[2];
sx q[2];
rz(-1.1954824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4698339) q[1];
sx q[1];
rz(-0.90171725) q[1];
sx q[1];
rz(-1.2059187) q[1];
rz(1.9530746) q[3];
sx q[3];
rz(-2.4197289) q[3];
sx q[3];
rz(-0.75306276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2842497) q[2];
sx q[2];
rz(-1.1025068) q[2];
sx q[2];
rz(-0.16172376) q[2];
rz(0.11180793) q[3];
sx q[3];
rz(-2.4155152) q[3];
sx q[3];
rz(0.73981729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1424471) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(0.69851843) q[0];
rz(2.0493719) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(0.05365595) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87854496) q[0];
sx q[0];
rz(-1.8302655) q[0];
sx q[0];
rz(-0.83570133) q[0];
x q[1];
rz(2.0094107) q[2];
sx q[2];
rz(-1.918186) q[2];
sx q[2];
rz(-1.6288822) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0567017) q[1];
sx q[1];
rz(-2.6680601) q[1];
sx q[1];
rz(-3.0242306) q[1];
rz(-pi) q[2];
rz(1.3732304) q[3];
sx q[3];
rz(-0.65614349) q[3];
sx q[3];
rz(0.11696091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1412389) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(-2.9290283) q[2];
rz(2.8138748) q[3];
sx q[3];
rz(-0.98723427) q[3];
sx q[3];
rz(1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(0.32972202) q[0];
rz(0.7715191) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(-0.82569295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59300789) q[0];
sx q[0];
rz(-2.4720044) q[0];
sx q[0];
rz(0.14051814) q[0];
rz(-3.0631376) q[2];
sx q[2];
rz(-2.1657073) q[2];
sx q[2];
rz(1.1537781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4976501) q[1];
sx q[1];
rz(-1.6794599) q[1];
sx q[1];
rz(2.3791299) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4903622) q[3];
sx q[3];
rz(-2.8222601) q[3];
sx q[3];
rz(-2.8791219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9416435) q[2];
sx q[2];
rz(-1.4120833) q[2];
sx q[2];
rz(-1.6597718) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.900443) q[0];
sx q[0];
rz(-0.41541442) q[0];
sx q[0];
rz(0.39733091) q[0];
rz(0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(1.0362524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5118714) q[0];
sx q[0];
rz(-0.40490926) q[0];
sx q[0];
rz(-0.029442336) q[0];
rz(2.7888227) q[2];
sx q[2];
rz(-1.4401541) q[2];
sx q[2];
rz(-0.58707159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81779382) q[1];
sx q[1];
rz(-1.7897871) q[1];
sx q[1];
rz(-1.4794255) q[1];
rz(-pi) q[2];
rz(-2.0898706) q[3];
sx q[3];
rz(-1.241893) q[3];
sx q[3];
rz(-2.2954706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1194666) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(-2.8890166) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(2.778497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0040358) q[0];
sx q[0];
rz(-0.79020774) q[0];
sx q[0];
rz(-2.9055415) q[0];
rz(-3.0421742) q[1];
sx q[1];
rz(-2.3853018) q[1];
sx q[1];
rz(-1.258446) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1374219) q[0];
sx q[0];
rz(-1.6485456) q[0];
sx q[0];
rz(0.56044062) q[0];
rz(-pi) q[1];
rz(-2.4781371) q[2];
sx q[2];
rz(-1.9773615) q[2];
sx q[2];
rz(-2.6321509) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26791456) q[1];
sx q[1];
rz(-2.2482052) q[1];
sx q[1];
rz(1.4235086) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7332337) q[3];
sx q[3];
rz(-1.1757188) q[3];
sx q[3];
rz(1.4407106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9749757) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(1.0413337) q[2];
rz(-0.58639041) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-2.7753684) q[0];
sx q[0];
rz(-2.0851705) q[0];
sx q[0];
rz(2.0023517) q[0];
rz(-2.1495023) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(-0.31789657) q[2];
sx q[2];
rz(-1.7902725) q[2];
sx q[2];
rz(-0.27223311) q[2];
rz(-0.59320368) q[3];
sx q[3];
rz(-1.9402258) q[3];
sx q[3];
rz(2.6711834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
