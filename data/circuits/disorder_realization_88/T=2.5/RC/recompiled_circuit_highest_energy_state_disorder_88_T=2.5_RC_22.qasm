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
rz(-3.1336194) q[0];
sx q[0];
rz(-2.70533) q[0];
sx q[0];
rz(1.0101969) q[0];
rz(-2.9695192) q[1];
sx q[1];
rz(-0.80614027) q[1];
sx q[1];
rz(-2.1176718) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4181241) q[0];
sx q[0];
rz(-2.0586657) q[0];
sx q[0];
rz(0.0077631891) q[0];
rz(-1.802581) q[2];
sx q[2];
rz(-1.9531774) q[2];
sx q[2];
rz(-0.2222375) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0144452) q[1];
sx q[1];
rz(-1.5513485) q[1];
sx q[1];
rz(-1.4102742) q[1];
rz(0.33396592) q[3];
sx q[3];
rz(-2.8779038) q[3];
sx q[3];
rz(2.7602344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1930834) q[2];
sx q[2];
rz(-2.812959) q[2];
sx q[2];
rz(-0.23552093) q[2];
rz(-1.0282907) q[3];
sx q[3];
rz(-1.2582658) q[3];
sx q[3];
rz(1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4368206) q[0];
sx q[0];
rz(-1.3587767) q[0];
sx q[0];
rz(0.92192465) q[0];
rz(0.1164662) q[1];
sx q[1];
rz(-1.4215697) q[1];
sx q[1];
rz(-2.2540653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015892643) q[0];
sx q[0];
rz(-1.3845446) q[0];
sx q[0];
rz(-0.13293477) q[0];
x q[1];
rz(0.15917425) q[2];
sx q[2];
rz(-0.35458699) q[2];
sx q[2];
rz(-1.872242) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4749467) q[1];
sx q[1];
rz(-1.4578867) q[1];
sx q[1];
rz(-0.571588) q[1];
rz(-pi) q[2];
rz(-0.87882249) q[3];
sx q[3];
rz(-1.678595) q[3];
sx q[3];
rz(-2.381058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.55887115) q[2];
sx q[2];
rz(-2.116394) q[2];
sx q[2];
rz(-2.9026046) q[2];
rz(-0.72238266) q[3];
sx q[3];
rz(-2.3014849) q[3];
sx q[3];
rz(1.1649789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.8257985) q[0];
sx q[0];
rz(-2.1738985) q[0];
sx q[0];
rz(-2.353299) q[0];
rz(-1.7514508) q[1];
sx q[1];
rz(-0.43594053) q[1];
sx q[1];
rz(1.699126) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9787131) q[0];
sx q[0];
rz(-1.7121234) q[0];
sx q[0];
rz(-3.0706329) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2295775) q[2];
sx q[2];
rz(-2.680353) q[2];
sx q[2];
rz(-1.7597424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6953814) q[1];
sx q[1];
rz(-1.3331116) q[1];
sx q[1];
rz(2.2878245) q[1];
rz(-pi) q[2];
rz(1.0679108) q[3];
sx q[3];
rz(-1.4137795) q[3];
sx q[3];
rz(0.63105052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.73551377) q[2];
sx q[2];
rz(-1.9903851) q[2];
sx q[2];
rz(0.25685143) q[2];
rz(-0.86366051) q[3];
sx q[3];
rz(-2.05859) q[3];
sx q[3];
rz(0.5736205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094512016) q[0];
sx q[0];
rz(-0.032945078) q[0];
sx q[0];
rz(-1.5091913) q[0];
rz(-2.8250561) q[1];
sx q[1];
rz(-1.5455952) q[1];
sx q[1];
rz(-2.1260156) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2779396) q[0];
sx q[0];
rz(-2.4653593) q[0];
sx q[0];
rz(2.3689224) q[0];
x q[1];
rz(-1.2945097) q[2];
sx q[2];
rz(-2.3450548) q[2];
sx q[2];
rz(2.403402) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8311743) q[1];
sx q[1];
rz(-1.948658) q[1];
sx q[1];
rz(-2.8997594) q[1];
x q[2];
rz(-1.7523132) q[3];
sx q[3];
rz(-2.3861775) q[3];
sx q[3];
rz(2.8528573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29869002) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(-1.451937) q[2];
rz(-1.9076094) q[3];
sx q[3];
rz(-1.0943639) q[3];
sx q[3];
rz(0.88172495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7538309) q[0];
sx q[0];
rz(-1.8330638) q[0];
sx q[0];
rz(2.0552788) q[0];
rz(1.7408675) q[1];
sx q[1];
rz(-0.81175214) q[1];
sx q[1];
rz(3.1382255) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28890507) q[0];
sx q[0];
rz(-0.59453911) q[0];
sx q[0];
rz(-3.1311228) q[0];
rz(-pi) q[1];
rz(-2.0515767) q[2];
sx q[2];
rz(-0.64489105) q[2];
sx q[2];
rz(-0.99586785) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0778179) q[1];
sx q[1];
rz(-1.6269725) q[1];
sx q[1];
rz(0.50775524) q[1];
rz(-0.60008757) q[3];
sx q[3];
rz(-1.5494124) q[3];
sx q[3];
rz(0.96259538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1668261) q[2];
sx q[2];
rz(-2.4415015) q[2];
sx q[2];
rz(0.33494803) q[2];
rz(2.8454928) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(-0.12224841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5122546) q[0];
sx q[0];
rz(-2.6423995) q[0];
sx q[0];
rz(-0.41803023) q[0];
rz(2.4755075) q[1];
sx q[1];
rz(-1.8981551) q[1];
sx q[1];
rz(-2.8935249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9739415) q[0];
sx q[0];
rz(-2.535916) q[0];
sx q[0];
rz(2.1781237) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0204629) q[2];
sx q[2];
rz(-1.5213335) q[2];
sx q[2];
rz(-0.63724697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0937178) q[1];
sx q[1];
rz(-1.1823726) q[1];
sx q[1];
rz(-1.018143) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94132386) q[3];
sx q[3];
rz(-1.3622614) q[3];
sx q[3];
rz(-2.1311614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.69998133) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(1.9786394) q[2];
rz(-0.59398389) q[3];
sx q[3];
rz(-1.5409527) q[3];
sx q[3];
rz(-0.99015132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.537259) q[0];
sx q[0];
rz(-2.8626677) q[0];
sx q[0];
rz(3.0766686) q[0];
rz(-0.36554947) q[1];
sx q[1];
rz(-1.7100916) q[1];
sx q[1];
rz(2.4117267) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.192963) q[0];
sx q[0];
rz(-2.2092487) q[0];
sx q[0];
rz(-3.0045511) q[0];
rz(1.7634298) q[2];
sx q[2];
rz(-0.81708497) q[2];
sx q[2];
rz(-1.1258923) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87903467) q[1];
sx q[1];
rz(-1.015519) q[1];
sx q[1];
rz(-3.0295275) q[1];
x q[2];
rz(-1.1564674) q[3];
sx q[3];
rz(-0.92770139) q[3];
sx q[3];
rz(-0.82655686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2713752) q[2];
sx q[2];
rz(-1.7729746) q[2];
sx q[2];
rz(-1.082487) q[2];
rz(-1.6535951) q[3];
sx q[3];
rz(-2.7797785) q[3];
sx q[3];
rz(2.7910119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7545886) q[0];
sx q[0];
rz(-2.258774) q[0];
sx q[0];
rz(2.8104695) q[0];
rz(-1.6638727) q[1];
sx q[1];
rz(-2.176087) q[1];
sx q[1];
rz(-2.8211735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77003682) q[0];
sx q[0];
rz(-1.5965684) q[0];
sx q[0];
rz(-3.0427264) q[0];
rz(-pi) q[1];
rz(0.1542145) q[2];
sx q[2];
rz(-1.7153828) q[2];
sx q[2];
rz(0.54942451) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3255405) q[1];
sx q[1];
rz(-2.7150702) q[1];
sx q[1];
rz(-0.037530516) q[1];
rz(-pi) q[2];
rz(3.0842264) q[3];
sx q[3];
rz(-0.88729492) q[3];
sx q[3];
rz(2.7025616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7122345) q[2];
sx q[2];
rz(-1.4171436) q[2];
sx q[2];
rz(2.7799535) q[2];
rz(2.2278191) q[3];
sx q[3];
rz(-2.8445966) q[3];
sx q[3];
rz(-2.698212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906032) q[0];
sx q[0];
rz(-2.3539703) q[0];
sx q[0];
rz(-3.0255764) q[0];
rz(-1.8138255) q[1];
sx q[1];
rz(-1.9783744) q[1];
sx q[1];
rz(1.2001002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3108771) q[0];
sx q[0];
rz(-0.19234622) q[0];
sx q[0];
rz(-0.098523454) q[0];
rz(-0.16559439) q[2];
sx q[2];
rz(-2.2496918) q[2];
sx q[2];
rz(2.0904581) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0018275) q[1];
sx q[1];
rz(-1.5854302) q[1];
sx q[1];
rz(-0.56165265) q[1];
x q[2];
rz(1.3312134) q[3];
sx q[3];
rz(-1.4093883) q[3];
sx q[3];
rz(0.39612285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3001331) q[2];
sx q[2];
rz(-1.7068784) q[2];
sx q[2];
rz(-0.33162281) q[2];
rz(-2.8578109) q[3];
sx q[3];
rz(-2.0399358) q[3];
sx q[3];
rz(-2.7200429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3109741) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(2.9750138) q[0];
rz(-0.84959787) q[1];
sx q[1];
rz(-0.46509898) q[1];
sx q[1];
rz(0.073624484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66921418) q[0];
sx q[0];
rz(-1.8249236) q[0];
sx q[0];
rz(2.4082802) q[0];
rz(-pi) q[1];
rz(0.4771941) q[2];
sx q[2];
rz(-1.0721803) q[2];
sx q[2];
rz(2.3248364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5391748) q[1];
sx q[1];
rz(-2.5218369) q[1];
sx q[1];
rz(2.5386993) q[1];
rz(-pi) q[2];
rz(2.0517379) q[3];
sx q[3];
rz(-2.1409799) q[3];
sx q[3];
rz(-0.4057623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2354043) q[2];
sx q[2];
rz(-2.3222458) q[2];
sx q[2];
rz(0.6582312) q[2];
rz(1.8002347) q[3];
sx q[3];
rz(-2.1737289) q[3];
sx q[3];
rz(2.2906176) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8258719) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(2.6201164) q[1];
sx q[1];
rz(-1.5670525) q[1];
sx q[1];
rz(-1.570931) q[1];
rz(2.1076219) q[2];
sx q[2];
rz(-0.85672094) q[2];
sx q[2];
rz(-2.232205) q[2];
rz(-1.7375199) q[3];
sx q[3];
rz(-2.1862595) q[3];
sx q[3];
rz(-1.4487472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
