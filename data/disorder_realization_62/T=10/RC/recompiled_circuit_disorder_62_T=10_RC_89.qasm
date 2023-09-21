OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8695495) q[0];
sx q[0];
rz(-0.46177319) q[0];
sx q[0];
rz(0.053134993) q[0];
rz(-2.7825836) q[2];
sx q[2];
rz(-2.3308672) q[2];
sx q[2];
rz(0.63149482) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.692688) q[1];
sx q[1];
rz(-1.2346039) q[1];
sx q[1];
rz(2.8640139) q[1];
rz(-pi) q[2];
rz(1.9840368) q[3];
sx q[3];
rz(-2.8765656) q[3];
sx q[3];
rz(2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(-2.5845394) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(-0.59659514) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81235028) q[0];
sx q[0];
rz(-1.3717522) q[0];
sx q[0];
rz(-0.45254405) q[0];
rz(2.3833582) q[2];
sx q[2];
rz(-0.97823921) q[2];
sx q[2];
rz(0.50298467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.021124161) q[1];
sx q[1];
rz(-1.7906584) q[1];
sx q[1];
rz(-0.0004417858) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5561043) q[3];
sx q[3];
rz(-0.5726632) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(-0.80667574) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.249041) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-2.1121315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91988504) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(2.1973781) q[0];
rz(-pi) q[1];
rz(0.37000334) q[2];
sx q[2];
rz(-2.4577933) q[2];
sx q[2];
rz(1.5918658) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6368235) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(-1.7117281) q[1];
rz(-pi) q[2];
rz(-0.13266487) q[3];
sx q[3];
rz(-0.74436114) q[3];
sx q[3];
rz(-0.45385195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(-1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(-1.625659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3665873) q[0];
sx q[0];
rz(-2.2457652) q[0];
sx q[0];
rz(1.3875899) q[0];
rz(-1.0117202) q[2];
sx q[2];
rz(-1.6435197) q[2];
sx q[2];
rz(0.59414547) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1277395) q[1];
sx q[1];
rz(-1.8426367) q[1];
sx q[1];
rz(-0.66687648) q[1];
rz(0.58873119) q[3];
sx q[3];
rz(-1.7541459) q[3];
sx q[3];
rz(0.4962894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76379124) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.8466922) q[2];
rz(0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(-2.5095818) q[1];
sx q[1];
rz(-0.72223392) q[1];
sx q[1];
rz(-2.246726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7475815) q[0];
sx q[0];
rz(-2.3590439) q[0];
sx q[0];
rz(-0.6015908) q[0];
rz(-pi) q[1];
rz(0.61770265) q[2];
sx q[2];
rz(-2.3139944) q[2];
sx q[2];
rz(0.4862116) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27069651) q[1];
sx q[1];
rz(-1.1333229) q[1];
sx q[1];
rz(-2.2938964) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5785061) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(-2.8849998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.9990702) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
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
rz(-2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-2.6766052) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87100055) q[0];
sx q[0];
rz(-0.56068476) q[0];
sx q[0];
rz(-0.90765783) q[0];
rz(-pi) q[1];
rz(-2.8000185) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(0.14788936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8919233) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(0.93530099) q[1];
rz(-0.24352169) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(-0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(0.4450376) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-0.26708189) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(0.15596095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6316846) q[0];
sx q[0];
rz(-1.8242867) q[0];
sx q[0];
rz(-1.7476728) q[0];
rz(-1.9975125) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(2.0063426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8487726) q[1];
sx q[1];
rz(-2.9584868) q[1];
sx q[1];
rz(1.8715026) q[1];
rz(-pi) q[2];
rz(-2.6050623) q[3];
sx q[3];
rz(-1.1324258) q[3];
sx q[3];
rz(2.8560672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(-2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9496562) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(-1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(-0.24857323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6876858) q[0];
sx q[0];
rz(-1.0605863) q[0];
sx q[0];
rz(1.2841671) q[0];
x q[1];
rz(0.82069223) q[2];
sx q[2];
rz(-0.65260115) q[2];
sx q[2];
rz(-1.5936268) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8057115) q[1];
sx q[1];
rz(-0.98166775) q[1];
sx q[1];
rz(-2.5475403) q[1];
rz(-pi) q[2];
rz(-1.1484654) q[3];
sx q[3];
rz(-2.9557807) q[3];
sx q[3];
rz(1.8728135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865737) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(1.3185906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628172) q[0];
sx q[0];
rz(-1.6512198) q[0];
sx q[0];
rz(2.0128987) q[0];
rz(-0.73762383) q[2];
sx q[2];
rz(-2.3361623) q[2];
sx q[2];
rz(-1.9464303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9642155) q[1];
sx q[1];
rz(-0.69987684) q[1];
sx q[1];
rz(-1.7978976) q[1];
x q[2];
rz(0.41140822) q[3];
sx q[3];
rz(-2.3477926) q[3];
sx q[3];
rz(1.0994764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(2.7424157) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(-1.5426853) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(1.261196) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1051837) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(-1.6856134) q[0];
x q[1];
rz(-1.5485498) q[2];
sx q[2];
rz(-0.97264475) q[2];
sx q[2];
rz(2.1200137) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2640472) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(1.643814) q[1];
rz(-pi) q[2];
rz(0.68393771) q[3];
sx q[3];
rz(-1.6947019) q[3];
sx q[3];
rz(2.370196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(0.56979257) q[2];
rz(1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-0.026253168) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(-0.59109296) q[3];
sx q[3];
rz(-1.4553087) q[3];
sx q[3];
rz(1.6607264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
