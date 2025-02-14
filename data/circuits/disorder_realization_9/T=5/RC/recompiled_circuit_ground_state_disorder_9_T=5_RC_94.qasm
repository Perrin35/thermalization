OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5300712) q[0];
sx q[0];
rz(-2.0687215) q[0];
sx q[0];
rz(1.4273509) q[0];
rz(0.99675769) q[1];
sx q[1];
rz(3.7550959) q[1];
sx q[1];
rz(9.4902314) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6320541) q[0];
sx q[0];
rz(-2.6231714) q[0];
sx q[0];
rz(-2.753756) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10721389) q[2];
sx q[2];
rz(-2.7709024) q[2];
sx q[2];
rz(3.0840741) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0330645) q[1];
sx q[1];
rz(-2.4208596) q[1];
sx q[1];
rz(2.0374551) q[1];
rz(-pi) q[2];
rz(-0.61566004) q[3];
sx q[3];
rz(-1.1634852) q[3];
sx q[3];
rz(3.1074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2685711) q[2];
sx q[2];
rz(-0.34653386) q[2];
sx q[2];
rz(-1.2968501) q[2];
rz(-1.1067363) q[3];
sx q[3];
rz(-0.96702558) q[3];
sx q[3];
rz(-1.4510252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035037128) q[0];
sx q[0];
rz(-0.7283926) q[0];
sx q[0];
rz(1.534071) q[0];
rz(-0.84397498) q[1];
sx q[1];
rz(-1.5183828) q[1];
sx q[1];
rz(0.27110505) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26626884) q[0];
sx q[0];
rz(-1.3014587) q[0];
sx q[0];
rz(-1.5311509) q[0];
x q[1];
rz(-0.098625318) q[2];
sx q[2];
rz(-1.319241) q[2];
sx q[2];
rz(-1.7401623) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7162937) q[1];
sx q[1];
rz(-1.4746086) q[1];
sx q[1];
rz(0.41751087) q[1];
rz(-pi) q[2];
rz(-0.45174349) q[3];
sx q[3];
rz(-2.9752033) q[3];
sx q[3];
rz(1.1903669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8035182) q[2];
sx q[2];
rz(-0.98793554) q[2];
sx q[2];
rz(2.3482077) q[2];
rz(-3.0986943) q[3];
sx q[3];
rz(-1.9148613) q[3];
sx q[3];
rz(2.4734917) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0778097) q[0];
sx q[0];
rz(-2.6580647) q[0];
sx q[0];
rz(0.10312816) q[0];
rz(0.25281301) q[1];
sx q[1];
rz(-1.4272855) q[1];
sx q[1];
rz(-0.3140744) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2677339) q[0];
sx q[0];
rz(-1.1815869) q[0];
sx q[0];
rz(-0.53797526) q[0];
rz(-pi) q[1];
rz(-3.0843749) q[2];
sx q[2];
rz(-2.0524745) q[2];
sx q[2];
rz(-1.6149947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.16437) q[1];
sx q[1];
rz(-0.80406351) q[1];
sx q[1];
rz(0.82694808) q[1];
x q[2];
rz(1.7400916) q[3];
sx q[3];
rz(-1.3162656) q[3];
sx q[3];
rz(-1.0485378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.126943) q[2];
sx q[2];
rz(-0.64121556) q[2];
sx q[2];
rz(2.3250697) q[2];
rz(-2.4294295) q[3];
sx q[3];
rz(-1.687259) q[3];
sx q[3];
rz(-0.85675353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6461058) q[0];
sx q[0];
rz(-0.093558885) q[0];
sx q[0];
rz(0.58864546) q[0];
rz(2.7694287) q[1];
sx q[1];
rz(-1.8446422) q[1];
sx q[1];
rz(0.94474307) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78658453) q[0];
sx q[0];
rz(-1.2517126) q[0];
sx q[0];
rz(-0.068885013) q[0];
rz(0.2619602) q[2];
sx q[2];
rz(-1.0925781) q[2];
sx q[2];
rz(-0.062104696) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0081629) q[1];
sx q[1];
rz(-1.1436698) q[1];
sx q[1];
rz(-2.5168946) q[1];
rz(-2.0382763) q[3];
sx q[3];
rz(-0.40492461) q[3];
sx q[3];
rz(-1.5011464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2302236) q[2];
sx q[2];
rz(-1.5350124) q[2];
sx q[2];
rz(-1.6564507) q[2];
rz(2.7459775) q[3];
sx q[3];
rz(-1.6499358) q[3];
sx q[3];
rz(0.0609456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7570801) q[0];
sx q[0];
rz(-0.39329305) q[0];
sx q[0];
rz(-1.242189) q[0];
rz(3.0643265) q[1];
sx q[1];
rz(-1.9237513) q[1];
sx q[1];
rz(-0.088931106) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1031694) q[0];
sx q[0];
rz(-2.9154202) q[0];
sx q[0];
rz(0.75407501) q[0];
rz(-pi) q[1];
rz(-2.9776666) q[2];
sx q[2];
rz(-2.5897615) q[2];
sx q[2];
rz(1.1140119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4684852) q[1];
sx q[1];
rz(-1.846828) q[1];
sx q[1];
rz(-1.0795781) q[1];
rz(1.4375646) q[3];
sx q[3];
rz(-1.4733088) q[3];
sx q[3];
rz(0.19791244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.58101216) q[2];
sx q[2];
rz(-1.3835013) q[2];
sx q[2];
rz(-1.9182473) q[2];
rz(1.2286202) q[3];
sx q[3];
rz(-2.2556428) q[3];
sx q[3];
rz(0.74738735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05518797) q[0];
sx q[0];
rz(-2.3265657) q[0];
sx q[0];
rz(-0.95274693) q[0];
rz(1.3505666) q[1];
sx q[1];
rz(-1.2341276) q[1];
sx q[1];
rz(1.1635253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19630963) q[0];
sx q[0];
rz(-2.3988814) q[0];
sx q[0];
rz(-0.16603827) q[0];
rz(-pi) q[1];
rz(0.38062747) q[2];
sx q[2];
rz(-2.5109595) q[2];
sx q[2];
rz(0.45898123) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7147594) q[1];
sx q[1];
rz(-1.3275255) q[1];
sx q[1];
rz(-3.0882278) q[1];
x q[2];
rz(2.5361193) q[3];
sx q[3];
rz(-1.9289199) q[3];
sx q[3];
rz(0.78152657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20840883) q[2];
sx q[2];
rz(-1.9219425) q[2];
sx q[2];
rz(1.5677933) q[2];
rz(0.21582223) q[3];
sx q[3];
rz(-0.62041557) q[3];
sx q[3];
rz(-2.9847434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1167574) q[0];
sx q[0];
rz(-0.76501608) q[0];
sx q[0];
rz(1.5417954) q[0];
rz(-2.722591) q[1];
sx q[1];
rz(-2.0332789) q[1];
sx q[1];
rz(-1.7760743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5657264) q[0];
sx q[0];
rz(-2.49263) q[0];
sx q[0];
rz(1.8945694) q[0];
rz(0.78132479) q[2];
sx q[2];
rz(-1.1622687) q[2];
sx q[2];
rz(0.21121001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1394386) q[1];
sx q[1];
rz(-1.4055058) q[1];
sx q[1];
rz(-2.2846245) q[1];
x q[2];
rz(1.0700052) q[3];
sx q[3];
rz(-0.65088256) q[3];
sx q[3];
rz(-0.18903519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8968203) q[2];
sx q[2];
rz(-0.12910566) q[2];
sx q[2];
rz(-1.6884035) q[2];
rz(-1.3990654) q[3];
sx q[3];
rz(-1.1996256) q[3];
sx q[3];
rz(-2.3340268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52706611) q[0];
sx q[0];
rz(-2.2699321) q[0];
sx q[0];
rz(0.51323071) q[0];
rz(-0.97663438) q[1];
sx q[1];
rz(-2.4726424) q[1];
sx q[1];
rz(-1.4248779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34951613) q[0];
sx q[0];
rz(-2.2648025) q[0];
sx q[0];
rz(-2.9727139) q[0];
rz(0.11995671) q[2];
sx q[2];
rz(-0.55177125) q[2];
sx q[2];
rz(0.063869501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2234474) q[1];
sx q[1];
rz(-0.53358101) q[1];
sx q[1];
rz(2.1672638) q[1];
x q[2];
rz(1.7376457) q[3];
sx q[3];
rz(-1.4544974) q[3];
sx q[3];
rz(-0.97613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49154115) q[2];
sx q[2];
rz(-1.0074793) q[2];
sx q[2];
rz(2.7752191) q[2];
rz(2.5036687) q[3];
sx q[3];
rz(-1.3831235) q[3];
sx q[3];
rz(1.5095507) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76846182) q[0];
sx q[0];
rz(-2.8167384) q[0];
sx q[0];
rz(-1.1573867) q[0];
rz(0.53818446) q[1];
sx q[1];
rz(-1.8073852) q[1];
sx q[1];
rz(0.023712637) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0932281) q[0];
sx q[0];
rz(-2.1914838) q[0];
sx q[0];
rz(1.4976504) q[0];
x q[1];
rz(-0.41466434) q[2];
sx q[2];
rz(-0.86277308) q[2];
sx q[2];
rz(-1.9239349) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7598724) q[1];
sx q[1];
rz(-2.4982735) q[1];
sx q[1];
rz(-2.5030959) q[1];
x q[2];
rz(1.0162418) q[3];
sx q[3];
rz(-2.9062677) q[3];
sx q[3];
rz(0.91431844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9862765) q[2];
sx q[2];
rz(-1.3331058) q[2];
sx q[2];
rz(-0.39815608) q[2];
rz(0.46477535) q[3];
sx q[3];
rz(-2.0538752) q[3];
sx q[3];
rz(-1.8006511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5584797) q[0];
sx q[0];
rz(-2.8093503) q[0];
sx q[0];
rz(2.1930021) q[0];
rz(-2.1743383) q[1];
sx q[1];
rz(-1.9890669) q[1];
sx q[1];
rz(-2.1327877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44701496) q[0];
sx q[0];
rz(-2.3409746) q[0];
sx q[0];
rz(0.009616212) q[0];
rz(-2.1494184) q[2];
sx q[2];
rz(-2.1632901) q[2];
sx q[2];
rz(-0.078607056) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12647835) q[1];
sx q[1];
rz(-2.6228944) q[1];
sx q[1];
rz(0.37094231) q[1];
x q[2];
rz(-2.4835577) q[3];
sx q[3];
rz(-1.74538) q[3];
sx q[3];
rz(3.0327219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41445109) q[2];
sx q[2];
rz(-2.264617) q[2];
sx q[2];
rz(-3.1205175) q[2];
rz(2.6575139) q[3];
sx q[3];
rz(-1.1295854) q[3];
sx q[3];
rz(2.4626203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2912343) q[0];
sx q[0];
rz(-1.6073011) q[0];
sx q[0];
rz(-1.4144443) q[0];
rz(0.5961295) q[1];
sx q[1];
rz(-2.3565751) q[1];
sx q[1];
rz(-0.48859488) q[1];
rz(-0.1487371) q[2];
sx q[2];
rz(-1.3631459) q[2];
sx q[2];
rz(2.7447328) q[2];
rz(1.3652741) q[3];
sx q[3];
rz(-1.4649656) q[3];
sx q[3];
rz(-3.0595915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
