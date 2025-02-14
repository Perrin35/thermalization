OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.705536) q[0];
sx q[0];
rz(-2.2826865) q[0];
sx q[0];
rz(-0.98281759) q[0];
rz(4.3217826) q[1];
sx q[1];
rz(5.5569841) q[1];
sx q[1];
rz(10.32294) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8456062) q[0];
sx q[0];
rz(-1.4665415) q[0];
sx q[0];
rz(-0.3180252) q[0];
rz(-2.6001789) q[2];
sx q[2];
rz(-1.2595508) q[2];
sx q[2];
rz(-1.0203066) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31831989) q[1];
sx q[1];
rz(-1.1538236) q[1];
sx q[1];
rz(0.61064536) q[1];
rz(-2.996436) q[3];
sx q[3];
rz(-1.4719799) q[3];
sx q[3];
rz(1.0001118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0183705) q[2];
sx q[2];
rz(-1.8223338) q[2];
sx q[2];
rz(-2.3940864) q[2];
rz(-1.2627164) q[3];
sx q[3];
rz(-1.5371753) q[3];
sx q[3];
rz(0.023716299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8338407) q[0];
sx q[0];
rz(-0.045578651) q[0];
sx q[0];
rz(-2.0763092) q[0];
rz(2.4899958) q[1];
sx q[1];
rz(-0.35531303) q[1];
sx q[1];
rz(1.5078872) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9621443) q[0];
sx q[0];
rz(-1.0217669) q[0];
sx q[0];
rz(0.35513504) q[0];
rz(2.6670611) q[2];
sx q[2];
rz(-1.2515278) q[2];
sx q[2];
rz(-0.87488562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2245141) q[1];
sx q[1];
rz(-1.1076756) q[1];
sx q[1];
rz(-0.45189894) q[1];
rz(-pi) q[2];
rz(1.5298858) q[3];
sx q[3];
rz(-0.61167756) q[3];
sx q[3];
rz(3.1168695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4773341) q[2];
sx q[2];
rz(-1.4724255) q[2];
sx q[2];
rz(-0.092546917) q[2];
rz(-2.7009098) q[3];
sx q[3];
rz(-0.40658545) q[3];
sx q[3];
rz(-1.1455166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6630845) q[0];
sx q[0];
rz(-2.9496851) q[0];
sx q[0];
rz(-1.2364016) q[0];
rz(2.9036486) q[1];
sx q[1];
rz(-1.6704208) q[1];
sx q[1];
rz(-2.1740289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1205115) q[0];
sx q[0];
rz(-0.97645611) q[0];
sx q[0];
rz(-2.3011776) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6009623) q[2];
sx q[2];
rz(-1.8701359) q[2];
sx q[2];
rz(-0.23330748) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6131276) q[1];
sx q[1];
rz(-1.2674164) q[1];
sx q[1];
rz(0.19618285) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7015412) q[3];
sx q[3];
rz(-2.7249911) q[3];
sx q[3];
rz(-0.9804014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6588916) q[2];
sx q[2];
rz(-2.7693558) q[2];
sx q[2];
rz(2.2230395) q[2];
rz(2.3679768) q[3];
sx q[3];
rz(-2.2548803) q[3];
sx q[3];
rz(0.95019597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2660265) q[0];
sx q[0];
rz(-0.73930621) q[0];
sx q[0];
rz(2.8493122) q[0];
rz(2.267011) q[1];
sx q[1];
rz(-0.62868172) q[1];
sx q[1];
rz(-0.94327092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2349512) q[0];
sx q[0];
rz(-2.4138193) q[0];
sx q[0];
rz(2.0208738) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9246632) q[2];
sx q[2];
rz(-2.0299021) q[2];
sx q[2];
rz(0.87068671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5041793) q[1];
sx q[1];
rz(-2.4054745) q[1];
sx q[1];
rz(2.0448684) q[1];
rz(-pi) q[2];
rz(0.99340474) q[3];
sx q[3];
rz(-1.3311989) q[3];
sx q[3];
rz(1.8519608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5248519) q[2];
sx q[2];
rz(-1.1139694) q[2];
sx q[2];
rz(0.88391602) q[2];
rz(2.0058477) q[3];
sx q[3];
rz(-0.92848778) q[3];
sx q[3];
rz(-0.69042027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8166872) q[0];
sx q[0];
rz(-0.08520928) q[0];
sx q[0];
rz(-2.9499522) q[0];
rz(1.2958255) q[1];
sx q[1];
rz(-1.192966) q[1];
sx q[1];
rz(-2.6909018) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5553868) q[0];
sx q[0];
rz(-1.4795917) q[0];
sx q[0];
rz(2.2074998) q[0];
rz(-0.24709267) q[2];
sx q[2];
rz(-0.28045248) q[2];
sx q[2];
rz(-2.7541898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6932842) q[1];
sx q[1];
rz(-2.5599883) q[1];
sx q[1];
rz(-1.4395836) q[1];
rz(2.8487318) q[3];
sx q[3];
rz(-1.5080875) q[3];
sx q[3];
rz(-0.69578275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4638764) q[2];
sx q[2];
rz(-1.8666942) q[2];
sx q[2];
rz(-2.1837168) q[2];
rz(0.052834474) q[3];
sx q[3];
rz(-2.5962679) q[3];
sx q[3];
rz(2.7132645) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65474725) q[0];
sx q[0];
rz(-1.0423648) q[0];
sx q[0];
rz(-0.47252193) q[0];
rz(2.3782702) q[1];
sx q[1];
rz(-2.4210052) q[1];
sx q[1];
rz(-0.29019132) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1520574) q[0];
sx q[0];
rz(-0.88071918) q[0];
sx q[0];
rz(-1.8581273) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0792208) q[2];
sx q[2];
rz(-2.5334458) q[2];
sx q[2];
rz(-1.7660994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7136894) q[1];
sx q[1];
rz(-2.2285502) q[1];
sx q[1];
rz(-1.7005928) q[1];
rz(-pi) q[2];
rz(0.84110705) q[3];
sx q[3];
rz(-0.26454138) q[3];
sx q[3];
rz(-3.0244248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7758238) q[2];
sx q[2];
rz(-0.86161986) q[2];
sx q[2];
rz(0.59930581) q[2];
rz(1.4186836) q[3];
sx q[3];
rz(-1.8214046) q[3];
sx q[3];
rz(2.1883709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391984) q[0];
sx q[0];
rz(-1.9559487) q[0];
sx q[0];
rz(1.8307357) q[0];
rz(-1.5015548) q[1];
sx q[1];
rz(-0.40032598) q[1];
sx q[1];
rz(-2.5379873) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5561605) q[0];
sx q[0];
rz(-1.7045665) q[0];
sx q[0];
rz(1.1473535) q[0];
rz(-pi) q[1];
rz(-0.55159388) q[2];
sx q[2];
rz(-1.3836129) q[2];
sx q[2];
rz(2.1647705) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7052066) q[1];
sx q[1];
rz(-1.5028186) q[1];
sx q[1];
rz(-1.8118565) q[1];
rz(-pi) q[2];
rz(1.7619753) q[3];
sx q[3];
rz(-2.303225) q[3];
sx q[3];
rz(-2.5922056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32460585) q[2];
sx q[2];
rz(-0.88674712) q[2];
sx q[2];
rz(0.19503197) q[2];
rz(-2.4845691) q[3];
sx q[3];
rz(-1.421509) q[3];
sx q[3];
rz(-0.11369625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3312155) q[0];
sx q[0];
rz(-1.4137784) q[0];
sx q[0];
rz(-3.069416) q[0];
rz(0.78308925) q[1];
sx q[1];
rz(-0.63242811) q[1];
sx q[1];
rz(-2.2236688) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24195053) q[0];
sx q[0];
rz(-0.52142087) q[0];
sx q[0];
rz(-0.42556406) q[0];
rz(-0.62064865) q[2];
sx q[2];
rz(-1.7908515) q[2];
sx q[2];
rz(1.4043851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5707122) q[1];
sx q[1];
rz(-1.2868795) q[1];
sx q[1];
rz(1.8181807) q[1];
rz(-pi) q[2];
rz(2.0948671) q[3];
sx q[3];
rz(-2.6597866) q[3];
sx q[3];
rz(0.039299358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5607295) q[2];
sx q[2];
rz(-1.9274638) q[2];
sx q[2];
rz(-1.5011935) q[2];
rz(0.8864657) q[3];
sx q[3];
rz(-1.3366046) q[3];
sx q[3];
rz(3.0361573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60370541) q[0];
sx q[0];
rz(-2.7925346) q[0];
sx q[0];
rz(2.6339997) q[0];
rz(2.4137068) q[1];
sx q[1];
rz(-1.4898841) q[1];
sx q[1];
rz(-2.0354039) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9637601) q[0];
sx q[0];
rz(-0.42362693) q[0];
sx q[0];
rz(-0.92064853) q[0];
rz(-pi) q[1];
x q[1];
rz(1.395325) q[2];
sx q[2];
rz(-1.9242491) q[2];
sx q[2];
rz(-2.8025016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.111208) q[1];
sx q[1];
rz(-1.2413352) q[1];
sx q[1];
rz(-1.9558286) q[1];
rz(-pi) q[2];
rz(2.15739) q[3];
sx q[3];
rz(-1.1569258) q[3];
sx q[3];
rz(-2.4936287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.20405208) q[2];
sx q[2];
rz(-2.0378518) q[2];
sx q[2];
rz(2.6943915) q[2];
rz(-1.4372545) q[3];
sx q[3];
rz(-0.81366003) q[3];
sx q[3];
rz(-1.6566431) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60780418) q[0];
sx q[0];
rz(-2.0543126) q[0];
sx q[0];
rz(0.5150038) q[0];
rz(1.9683413) q[1];
sx q[1];
rz(-1.8517905) q[1];
sx q[1];
rz(2.692093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7106626) q[0];
sx q[0];
rz(-1.8966513) q[0];
sx q[0];
rz(-0.023825721) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1859288) q[2];
sx q[2];
rz(-1.1329831) q[2];
sx q[2];
rz(-2.4085338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.044552) q[1];
sx q[1];
rz(-2.1637193) q[1];
sx q[1];
rz(2.135599) q[1];
rz(-pi) q[2];
rz(1.2925384) q[3];
sx q[3];
rz(-0.738916) q[3];
sx q[3];
rz(0.41058879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0476734) q[2];
sx q[2];
rz(-2.8670222) q[2];
sx q[2];
rz(-2.3982415) q[2];
rz(2.7798233) q[3];
sx q[3];
rz(-1.0310562) q[3];
sx q[3];
rz(-0.37775347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.34306985) q[0];
sx q[0];
rz(-1.1953851) q[0];
sx q[0];
rz(0.60085798) q[0];
rz(0.15380225) q[1];
sx q[1];
rz(-1.8228795) q[1];
sx q[1];
rz(0.22620329) q[1];
rz(1.9696196) q[2];
sx q[2];
rz(-0.33179596) q[2];
sx q[2];
rz(-1.1992762) q[2];
rz(2.5220925) q[3];
sx q[3];
rz(-1.8826857) q[3];
sx q[3];
rz(2.7553325) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
