OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73206168) q[0];
sx q[0];
rz(-1.7763897) q[0];
sx q[0];
rz(2.1172297) q[0];
rz(0.60511869) q[1];
sx q[1];
rz(-0.53202283) q[1];
sx q[1];
rz(-1.9722809) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3496075) q[0];
sx q[0];
rz(-1.2149997) q[0];
sx q[0];
rz(-2.4277359) q[0];
rz(-pi) q[1];
rz(0.77180441) q[2];
sx q[2];
rz(-2.3183841) q[2];
sx q[2];
rz(-0.31847218) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3633903) q[1];
sx q[1];
rz(-1.9323903) q[1];
sx q[1];
rz(1.1671288) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23006769) q[3];
sx q[3];
rz(-2.821273) q[3];
sx q[3];
rz(0.39428082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8756276) q[2];
sx q[2];
rz(-0.83845323) q[2];
sx q[2];
rz(-1.3226091) q[2];
rz(-0.30098513) q[3];
sx q[3];
rz(-0.61166489) q[3];
sx q[3];
rz(1.7606364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.009636119) q[0];
sx q[0];
rz(-0.29254237) q[0];
sx q[0];
rz(-2.6665376) q[0];
rz(-1.7430199) q[1];
sx q[1];
rz(-2.1865632) q[1];
sx q[1];
rz(-1.0377201) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25181928) q[0];
sx q[0];
rz(-0.60476859) q[0];
sx q[0];
rz(2.7483727) q[0];
rz(0.77387626) q[2];
sx q[2];
rz(-0.69303382) q[2];
sx q[2];
rz(-0.80640031) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.84320074) q[1];
sx q[1];
rz(-1.1644662) q[1];
sx q[1];
rz(0.028224736) q[1];
rz(-2.5492937) q[3];
sx q[3];
rz(-0.90511887) q[3];
sx q[3];
rz(0.75331068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66118801) q[2];
sx q[2];
rz(-1.3385237) q[2];
sx q[2];
rz(-0.084687106) q[2];
rz(-2.7627913) q[3];
sx q[3];
rz(-2.8642604) q[3];
sx q[3];
rz(-1.144073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.5304853) q[0];
sx q[0];
rz(-2.0407016) q[0];
sx q[0];
rz(-2.1858922) q[0];
rz(2.7509007) q[1];
sx q[1];
rz(-0.57084584) q[1];
sx q[1];
rz(-0.57317615) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5432376) q[0];
sx q[0];
rz(-2.7063473) q[0];
sx q[0];
rz(-1.4198562) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3689752) q[2];
sx q[2];
rz(-1.3396016) q[2];
sx q[2];
rz(-0.85341838) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1307615) q[1];
sx q[1];
rz(-1.8567994) q[1];
sx q[1];
rz(-2.7421013) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8543386) q[3];
sx q[3];
rz(-0.58218282) q[3];
sx q[3];
rz(-0.69609387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0009784) q[2];
sx q[2];
rz(-0.30423519) q[2];
sx q[2];
rz(-0.19392459) q[2];
rz(-3.0443232) q[3];
sx q[3];
rz(-1.8563742) q[3];
sx q[3];
rz(-2.9320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(0.28213421) q[0];
sx q[0];
rz(-0.54444805) q[0];
sx q[0];
rz(2.5909246) q[0];
rz(2.0129054) q[1];
sx q[1];
rz(-1.0602602) q[1];
sx q[1];
rz(-2.7788924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.761128) q[0];
sx q[0];
rz(-1.3306381) q[0];
sx q[0];
rz(-2.2855177) q[0];
rz(-pi) q[1];
rz(2.0998459) q[2];
sx q[2];
rz(-3.1051271) q[2];
sx q[2];
rz(1.0066102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94707045) q[1];
sx q[1];
rz(-2.2010989) q[1];
sx q[1];
rz(-1.6987726) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6710715) q[3];
sx q[3];
rz(-0.62697151) q[3];
sx q[3];
rz(-2.7659622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68391934) q[2];
sx q[2];
rz(-1.0689015) q[2];
sx q[2];
rz(0.0030227946) q[2];
rz(2.4827042) q[3];
sx q[3];
rz(-0.342841) q[3];
sx q[3];
rz(-0.80250424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9534849) q[0];
sx q[0];
rz(-2.4664724) q[0];
sx q[0];
rz(-3.127393) q[0];
rz(3.1242127) q[1];
sx q[1];
rz(-2.1936369) q[1];
sx q[1];
rz(-1.4594706) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5495816) q[0];
sx q[0];
rz(-1.3322543) q[0];
sx q[0];
rz(2.9907988) q[0];
rz(2.2541788) q[2];
sx q[2];
rz(-1.9878584) q[2];
sx q[2];
rz(2.8487157) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7534605) q[1];
sx q[1];
rz(-1.7963444) q[1];
sx q[1];
rz(-1.7727586) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0586117) q[3];
sx q[3];
rz(-1.2851464) q[3];
sx q[3];
rz(0.80174996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3061299) q[2];
sx q[2];
rz(-0.43734044) q[2];
sx q[2];
rz(-2.2996976) q[2];
rz(-1.016559) q[3];
sx q[3];
rz(-2.026365) q[3];
sx q[3];
rz(-1.564933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.013997812) q[0];
sx q[0];
rz(-2.4348149) q[0];
sx q[0];
rz(2.5573964) q[0];
rz(-1.2305413) q[1];
sx q[1];
rz(-2.0293593) q[1];
sx q[1];
rz(3.0029283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16620557) q[0];
sx q[0];
rz(-1.5674601) q[0];
sx q[0];
rz(-0.020394527) q[0];
rz(-pi) q[1];
rz(2.9406592) q[2];
sx q[2];
rz(-1.4954508) q[2];
sx q[2];
rz(-0.88694015) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12339679) q[1];
sx q[1];
rz(-2.1235848) q[1];
sx q[1];
rz(0.68560302) q[1];
x q[2];
rz(-1.9964553) q[3];
sx q[3];
rz(-2.2666551) q[3];
sx q[3];
rz(-1.3130207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77506322) q[2];
sx q[2];
rz(-1.083192) q[2];
sx q[2];
rz(1.3343875) q[2];
rz(1.9813609) q[3];
sx q[3];
rz(-1.3637873) q[3];
sx q[3];
rz(3.0464723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60550624) q[0];
sx q[0];
rz(-1.1331929) q[0];
sx q[0];
rz(-2.6050674) q[0];
rz(-0.58553186) q[1];
sx q[1];
rz(-0.1383055) q[1];
sx q[1];
rz(2.5172863) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7621967) q[0];
sx q[0];
rz(-1.3610098) q[0];
sx q[0];
rz(2.4136153) q[0];
rz(-pi) q[1];
rz(2.1820071) q[2];
sx q[2];
rz(-0.72003905) q[2];
sx q[2];
rz(2.6401273) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.149951) q[1];
sx q[1];
rz(-2.2551564) q[1];
sx q[1];
rz(-2.409694) q[1];
x q[2];
rz(0.47847139) q[3];
sx q[3];
rz(-2.1395184) q[3];
sx q[3];
rz(2.3346321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6162993) q[2];
sx q[2];
rz(-2.0234183) q[2];
sx q[2];
rz(-0.38254151) q[2];
rz(-3.110102) q[3];
sx q[3];
rz(-0.68920207) q[3];
sx q[3];
rz(0.1077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-0.87576762) q[0];
sx q[0];
rz(-0.28755292) q[0];
sx q[0];
rz(-0.13993046) q[0];
rz(-1.4639927) q[1];
sx q[1];
rz(-0.96698562) q[1];
sx q[1];
rz(0.12891842) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0807954) q[0];
sx q[0];
rz(-1.1757438) q[0];
sx q[0];
rz(-0.11443826) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.244197) q[2];
sx q[2];
rz(-2.348263) q[2];
sx q[2];
rz(-0.71133864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8184549) q[1];
sx q[1];
rz(-1.8433246) q[1];
sx q[1];
rz(0.44169359) q[1];
x q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-2.401712) q[3];
sx q[3];
rz(3.0144514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.504618) q[2];
sx q[2];
rz(-2.0998173) q[2];
sx q[2];
rz(-2.5174985) q[2];
rz(-2.9028153) q[3];
sx q[3];
rz(-1.5437361) q[3];
sx q[3];
rz(1.0673267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36214608) q[0];
sx q[0];
rz(-1.0405259) q[0];
sx q[0];
rz(1.2497485) q[0];
rz(-0.016013913) q[1];
sx q[1];
rz(-2.3857954) q[1];
sx q[1];
rz(-1.790766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46247813) q[0];
sx q[0];
rz(-0.20972855) q[0];
sx q[0];
rz(2.1049343) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58198858) q[2];
sx q[2];
rz(-2.4798923) q[2];
sx q[2];
rz(2.2616507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1180601) q[1];
sx q[1];
rz(-2.2242821) q[1];
sx q[1];
rz(-0.94783028) q[1];
rz(-pi) q[2];
rz(-1.8072855) q[3];
sx q[3];
rz(-2.3097976) q[3];
sx q[3];
rz(-0.47154271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0525557) q[2];
sx q[2];
rz(-2.3904843) q[2];
sx q[2];
rz(3.0272711) q[2];
rz(1.8814686) q[3];
sx q[3];
rz(-2.114664) q[3];
sx q[3];
rz(1.982622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2262912) q[0];
sx q[0];
rz(-1.4878595) q[0];
sx q[0];
rz(0.21324883) q[0];
rz(-0.419871) q[1];
sx q[1];
rz(-2.1397736) q[1];
sx q[1];
rz(0.54668033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28323805) q[0];
sx q[0];
rz(-0.75728098) q[0];
sx q[0];
rz(-1.8961294) q[0];
rz(-pi) q[1];
rz(0.47537739) q[2];
sx q[2];
rz(-1.7065085) q[2];
sx q[2];
rz(-3.140608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.039779546) q[1];
sx q[1];
rz(-1.8598078) q[1];
sx q[1];
rz(1.6092369) q[1];
x q[2];
rz(-2.91594) q[3];
sx q[3];
rz(-1.435558) q[3];
sx q[3];
rz(1.7547363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13835779) q[2];
sx q[2];
rz(-1.582575) q[2];
sx q[2];
rz(-3.1372916) q[2];
rz(0.99758482) q[3];
sx q[3];
rz(-2.6514566) q[3];
sx q[3];
rz(-0.51013851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1417086) q[0];
sx q[0];
rz(-2.0337491) q[0];
sx q[0];
rz(0.98325892) q[0];
rz(0.44395631) q[1];
sx q[1];
rz(-2.8580491) q[1];
sx q[1];
rz(-1.8681189) q[1];
rz(-3.1267358) q[2];
sx q[2];
rz(-1.5236293) q[2];
sx q[2];
rz(-2.2894494) q[2];
rz(-1.0989582) q[3];
sx q[3];
rz(-1.8660587) q[3];
sx q[3];
rz(2.3793424) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
