OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(3.9225188) q[0];
sx q[0];
rz(9.6315686) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(1.0607399) q[1];
sx q[1];
rz(12.386204) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313523) q[0];
sx q[0];
rz(-0.58824476) q[0];
sx q[0];
rz(2.504185) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51282672) q[2];
sx q[2];
rz(-1.4564118) q[2];
sx q[2];
rz(-2.0778542) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.441592) q[1];
sx q[1];
rz(-2.4543426) q[1];
sx q[1];
rz(-1.4541461) q[1];
rz(-2.9163625) q[3];
sx q[3];
rz(-0.69437829) q[3];
sx q[3];
rz(0.089689342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(1.2940548) q[2];
rz(0.40575746) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(-0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(3.0157715) q[0];
rz(2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(-1.7696101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67324154) q[0];
sx q[0];
rz(-1.5090319) q[0];
sx q[0];
rz(-0.67586918) q[0];
rz(1.4470909) q[2];
sx q[2];
rz(-1.0281841) q[2];
sx q[2];
rz(-1.6306842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11916493) q[1];
sx q[1];
rz(-0.25401527) q[1];
sx q[1];
rz(2.959842) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6729229) q[3];
sx q[3];
rz(-1.3523003) q[3];
sx q[3];
rz(3.1359429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8877318) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(-0.99622336) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(-1.0590142) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3740736) q[0];
sx q[0];
rz(-0.78157434) q[0];
sx q[0];
rz(2.563345) q[0];
x q[1];
rz(1.1946482) q[2];
sx q[2];
rz(-2.7333439) q[2];
sx q[2];
rz(0.86366913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3036084) q[1];
sx q[1];
rz(-1.260699) q[1];
sx q[1];
rz(-0.012197818) q[1];
rz(-pi) q[2];
rz(-0.16626658) q[3];
sx q[3];
rz(-0.85321745) q[3];
sx q[3];
rz(2.1239514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(-0.26322571) q[2];
rz(2.0227382) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(-2.0746453) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.432503) q[0];
sx q[0];
rz(-1.6095918) q[0];
sx q[0];
rz(-1.0409271) q[0];
x q[1];
rz(-2.8550451) q[2];
sx q[2];
rz(-2.8612125) q[2];
sx q[2];
rz(2.7718411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.74944118) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(0.72026003) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9846109) q[3];
sx q[3];
rz(-0.51895751) q[3];
sx q[3];
rz(-1.8116236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.84714326) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(-0.28953826) q[2];
rz(-2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(-2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510058) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(2.3068413) q[0];
rz(-1.897215) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(-1.7117737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1236876) q[0];
sx q[0];
rz(-0.062739685) q[0];
sx q[0];
rz(0.040266589) q[0];
rz(-pi) q[1];
rz(-2.447117) q[2];
sx q[2];
rz(-2.777213) q[2];
sx q[2];
rz(0.37730544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5606219) q[1];
sx q[1];
rz(-0.79172687) q[1];
sx q[1];
rz(-0.03246275) q[1];
rz(-pi) q[2];
rz(-0.91246446) q[3];
sx q[3];
rz(-1.4294799) q[3];
sx q[3];
rz(-0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4804068) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.8583813) q[2];
rz(-3.0094106) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(-2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-2.5710411) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30734277) q[0];
sx q[0];
rz(-1.9295921) q[0];
sx q[0];
rz(-0.3105727) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.19679) q[2];
sx q[2];
rz(-1.9960969) q[2];
sx q[2];
rz(2.5380295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.693336) q[1];
sx q[1];
rz(-2.2168471) q[1];
sx q[1];
rz(2.4962884) q[1];
rz(1.5930575) q[3];
sx q[3];
rz(-2.2189757) q[3];
sx q[3];
rz(0.19433403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657848) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(2.7835223) q[0];
rz(-0.2886731) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(1.3105062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34328) q[0];
sx q[0];
rz(-0.57050059) q[0];
sx q[0];
rz(2.7891854) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4438379) q[2];
sx q[2];
rz(-2.0332094) q[2];
sx q[2];
rz(-2.4042839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9999034) q[1];
sx q[1];
rz(-2.018376) q[1];
sx q[1];
rz(-2.5738641) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2566124) q[3];
sx q[3];
rz(-1.0776057) q[3];
sx q[3];
rz(-2.5964338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(-0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.7060446) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26577935) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(1.2365201) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(1.1901201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8461788) q[0];
sx q[0];
rz(-0.23970397) q[0];
sx q[0];
rz(-1.833605) q[0];
rz(-pi) q[1];
rz(-0.868452) q[2];
sx q[2];
rz(-0.10570279) q[2];
sx q[2];
rz(2.5603574) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5684143) q[1];
sx q[1];
rz(-1.8548994) q[1];
sx q[1];
rz(-1.2969639) q[1];
rz(-pi) q[2];
x q[2];
rz(0.45300071) q[3];
sx q[3];
rz(-1.060033) q[3];
sx q[3];
rz(2.9651027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3999346) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(-1.3828145) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(-0.32489052) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(2.9898306) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230726) q[0];
sx q[0];
rz(-1.1306445) q[0];
sx q[0];
rz(-1.7258304) q[0];
rz(-0.62990909) q[2];
sx q[2];
rz(-1.3917149) q[2];
sx q[2];
rz(1.5392898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4703428) q[1];
sx q[1];
rz(-1.8103231) q[1];
sx q[1];
rz(-2.4688979) q[1];
rz(-pi) q[2];
rz(2.1612694) q[3];
sx q[3];
rz(-1.1579517) q[3];
sx q[3];
rz(2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(2.7091743) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(2.7684257) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-0.7235136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31691566) q[0];
sx q[0];
rz(-2.4002889) q[0];
sx q[0];
rz(2.9497428) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11156545) q[2];
sx q[2];
rz(-2.4590883) q[2];
sx q[2];
rz(2.0276558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53612667) q[1];
sx q[1];
rz(-1.3007174) q[1];
sx q[1];
rz(-2.3675766) q[1];
rz(2.8711653) q[3];
sx q[3];
rz(-1.575483) q[3];
sx q[3];
rz(-0.0029431012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2202806) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-2.5881361) q[2];
rz(2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217459) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(0.28221054) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(0.7278022) q[2];
sx q[2];
rz(-2.5143378) q[2];
sx q[2];
rz(0.85744748) q[2];
rz(-1.2890733) q[3];
sx q[3];
rz(-1.1829794) q[3];
sx q[3];
rz(1.1869528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];