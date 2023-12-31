OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(1.3442858) q[1];
sx q[1];
rz(-2.9872515) q[1];
sx q[1];
rz(0.57758346) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87969765) q[0];
sx q[0];
rz(-1.4059773) q[0];
sx q[0];
rz(-0.66189712) q[0];
x q[1];
rz(2.0787813) q[2];
sx q[2];
rz(-1.0928109) q[2];
sx q[2];
rz(-2.9302772) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56596245) q[1];
sx q[1];
rz(-1.168025) q[1];
sx q[1];
rz(2.8494542) q[1];
rz(-pi) q[2];
x q[2];
rz(2.203381) q[3];
sx q[3];
rz(-2.1198366) q[3];
sx q[3];
rz(-3.0905746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.43705964) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(-2.4543767) q[2];
rz(-1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(-0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17094831) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(-1.2600391) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0802404) q[0];
sx q[0];
rz(-0.70525673) q[0];
sx q[0];
rz(2.4387226) q[0];
rz(-1.4405865) q[2];
sx q[2];
rz(-1.4381593) q[2];
sx q[2];
rz(-0.33078937) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0963124) q[1];
sx q[1];
rz(-1.7573866) q[1];
sx q[1];
rz(0.50218302) q[1];
rz(2.1917079) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(-1.3148395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.1876748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1903494) q[0];
sx q[0];
rz(-2.9028063) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(-0.13126016) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(-2.0551596) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7032996) q[0];
sx q[0];
rz(-0.46143954) q[0];
sx q[0];
rz(-1.567054) q[0];
rz(1.7227206) q[2];
sx q[2];
rz(-0.65290367) q[2];
sx q[2];
rz(2.8002847) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74124386) q[1];
sx q[1];
rz(-1.2819918) q[1];
sx q[1];
rz(-1.5762394) q[1];
rz(-pi) q[2];
rz(1.1658737) q[3];
sx q[3];
rz(-1.7305264) q[3];
sx q[3];
rz(-2.3111642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0125668) q[2];
sx q[2];
rz(-1.4386703) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.8846735) q[3];
sx q[3];
rz(2.3390521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6894158) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(-3.0932328) q[0];
rz(-0.16391779) q[1];
sx q[1];
rz(-0.36968958) q[1];
sx q[1];
rz(-1.4455459) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7955129) q[0];
sx q[0];
rz(-2.4582986) q[0];
sx q[0];
rz(0.28738316) q[0];
rz(0.21638685) q[2];
sx q[2];
rz(-1.6587703) q[2];
sx q[2];
rz(0.83425922) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85019894) q[1];
sx q[1];
rz(-1.7899917) q[1];
sx q[1];
rz(-1.5047969) q[1];
rz(-2.8533832) q[3];
sx q[3];
rz(-2.0734348) q[3];
sx q[3];
rz(0.99311815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0754898) q[2];
sx q[2];
rz(-1.7843856) q[2];
sx q[2];
rz(2.1172822) q[2];
rz(1.5284437) q[3];
sx q[3];
rz(-1.6201092) q[3];
sx q[3];
rz(-2.9083692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7016474) q[0];
sx q[0];
rz(-0.82413903) q[0];
sx q[0];
rz(1.8540927) q[0];
rz(2.8225186) q[1];
sx q[1];
rz(-1.5998452) q[1];
sx q[1];
rz(0.85420001) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6720649) q[0];
sx q[0];
rz(-2.2846691) q[0];
sx q[0];
rz(-3.1307427) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2703083) q[2];
sx q[2];
rz(-1.9565906) q[2];
sx q[2];
rz(-2.3988349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2367868) q[1];
sx q[1];
rz(-2.2854837) q[1];
sx q[1];
rz(2.3805815) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1339995) q[3];
sx q[3];
rz(-1.0922722) q[3];
sx q[3];
rz(-2.0625045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-2.5482735) q[2];
sx q[2];
rz(-0.577315) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(-2.234941) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0034870738) q[0];
sx q[0];
rz(-1.071799) q[0];
sx q[0];
rz(-3.0694718) q[0];
rz(1.1068608) q[1];
sx q[1];
rz(-0.51263428) q[1];
sx q[1];
rz(-0.12621005) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5835411) q[0];
sx q[0];
rz(-1.3582894) q[0];
sx q[0];
rz(-1.538518) q[0];
rz(-pi) q[1];
rz(-2.3665479) q[2];
sx q[2];
rz(-1.8473052) q[2];
sx q[2];
rz(1.0724049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0969442) q[1];
sx q[1];
rz(-1.7587887) q[1];
sx q[1];
rz(-0.85) q[1];
rz(-pi) q[2];
x q[2];
rz(1.24228) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(-1.3086705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(2.3136247) q[3];
sx q[3];
rz(-2.8505846) q[3];
sx q[3];
rz(1.1221788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5489952) q[0];
sx q[0];
rz(-2.7153375) q[0];
sx q[0];
rz(0.098408498) q[0];
rz(1.1920284) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0347621) q[0];
sx q[0];
rz(-2.1666399) q[0];
sx q[0];
rz(1.4164657) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9497629) q[2];
sx q[2];
rz(-0.26608135) q[2];
sx q[2];
rz(1.9907469) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.050484867) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(2.1662103) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8758994) q[3];
sx q[3];
rz(-2.4943647) q[3];
sx q[3];
rz(-0.94097394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6825535) q[2];
sx q[2];
rz(-1.295853) q[2];
sx q[2];
rz(-2.7977978) q[2];
rz(-2.5750459) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(2.6678273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7664117) q[0];
sx q[0];
rz(-1.8090929) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(-0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(0.87160814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9974737) q[0];
sx q[0];
rz(-1.561957) q[0];
sx q[0];
rz(-0.075450443) q[0];
x q[1];
rz(0.34198728) q[2];
sx q[2];
rz(-2.7139398) q[2];
sx q[2];
rz(-0.74117408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.034681319) q[1];
sx q[1];
rz(-0.92202631) q[1];
sx q[1];
rz(-1.1120863) q[1];
rz(0.96846795) q[3];
sx q[3];
rz(-1.6864711) q[3];
sx q[3];
rz(-1.0761716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72620755) q[2];
sx q[2];
rz(-1.0691079) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.6428044) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.9198445) q[0];
rz(-2.9755759) q[1];
sx q[1];
rz(-1.32042) q[1];
sx q[1];
rz(1.5244012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0696408) q[0];
sx q[0];
rz(-1.4884293) q[0];
sx q[0];
rz(0.029043341) q[0];
rz(-1.379307) q[2];
sx q[2];
rz(-2.0404437) q[2];
sx q[2];
rz(-2.7188403) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.081454885) q[1];
sx q[1];
rz(-2.7739035) q[1];
sx q[1];
rz(2.0245488) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41096656) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.021585492) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-2.7900556) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(-0.74469152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979424) q[0];
sx q[0];
rz(-0.90181667) q[0];
sx q[0];
rz(1.836401) q[0];
rz(-2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(2.8881853) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2828335) q[0];
sx q[0];
rz(-1.4677605) q[0];
sx q[0];
rz(1.1742924) q[0];
rz(-1.9099405) q[2];
sx q[2];
rz(-1.412743) q[2];
sx q[2];
rz(-0.16976419) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0903783) q[1];
sx q[1];
rz(-1.9131294) q[1];
sx q[1];
rz(1.8263837) q[1];
rz(-pi) q[2];
rz(-1.1425584) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(2.464307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5499251) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(-0.5029451) q[2];
rz(-2.2425966) q[3];
sx q[3];
rz(-1.2939724) q[3];
sx q[3];
rz(-1.9780654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1702561) q[0];
sx q[0];
rz(-1.5383056) q[0];
sx q[0];
rz(-2.8785895) q[0];
rz(2.4304216) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(-2.860643) q[2];
sx q[2];
rz(-1.8200257) q[2];
sx q[2];
rz(-0.65489468) q[2];
rz(-0.9234225) q[3];
sx q[3];
rz(-2.2150726) q[3];
sx q[3];
rz(1.1005145) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
