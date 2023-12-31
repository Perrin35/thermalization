OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0922001) q[0];
sx q[0];
rz(-0.28486262) q[0];
sx q[0];
rz(2.0869135) q[0];
rz(-1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(1.7681233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76525926) q[1];
sx q[1];
rz(-1.9783124) q[1];
sx q[1];
rz(2.5024662) q[1];
x q[2];
rz(-0.7198556) q[3];
sx q[3];
rz(-2.8350283) q[3];
sx q[3];
rz(-2.7048064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0063643) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-2.485086) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(-2.7242993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1332557) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(-0.061925109) q[1];
sx q[1];
rz(-1.2281111) q[1];
sx q[1];
rz(-2.6541236) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66747626) q[0];
sx q[0];
rz(-1.6352909) q[0];
sx q[0];
rz(-1.6518031) q[0];
rz(-pi) q[1];
rz(0.42178085) q[2];
sx q[2];
rz(-1.1676719) q[2];
sx q[2];
rz(-0.031907206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5408052) q[1];
sx q[1];
rz(-1.4122737) q[1];
sx q[1];
rz(-1.7608789) q[1];
rz(1.3120193) q[3];
sx q[3];
rz(-0.32734713) q[3];
sx q[3];
rz(-0.60206383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9618824) q[2];
sx q[2];
rz(-1.2625182) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(-3.1029207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9449126) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-2.9512067) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-2.7540837) q[1];
sx q[1];
rz(0.22274676) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4653141) q[0];
sx q[0];
rz(-1.197581) q[0];
sx q[0];
rz(2.7200384) q[0];
rz(-0.22521714) q[2];
sx q[2];
rz(-1.9334963) q[2];
sx q[2];
rz(2.1622216) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1258771) q[1];
sx q[1];
rz(-1.5029969) q[1];
sx q[1];
rz(-2.5281639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.907903) q[3];
sx q[3];
rz(-1.8030231) q[3];
sx q[3];
rz(0.043134886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7524183) q[0];
sx q[0];
rz(-2.8635633) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(2.2213675) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5266787) q[0];
sx q[0];
rz(-2.3290312) q[0];
sx q[0];
rz(2.148669) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4843066) q[2];
sx q[2];
rz(-1.930634) q[2];
sx q[2];
rz(2.6038225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72318075) q[1];
sx q[1];
rz(-2.3310117) q[1];
sx q[1];
rz(0.65631887) q[1];
x q[2];
rz(-2.6357189) q[3];
sx q[3];
rz(-0.40588356) q[3];
sx q[3];
rz(-2.1959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.197864) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(0.81930339) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6827877) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(2.6884902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83349193) q[0];
sx q[0];
rz(-1.6080603) q[0];
sx q[0];
rz(1.6444071) q[0];
rz(-3.1066936) q[2];
sx q[2];
rz(-1.4279832) q[2];
sx q[2];
rz(-2.3226483) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52757712) q[1];
sx q[1];
rz(-1.3506538) q[1];
sx q[1];
rz(3.0575391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.971332) q[3];
sx q[3];
rz(-1.1820275) q[3];
sx q[3];
rz(-2.142981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12864628) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(0.72367469) q[3];
sx q[3];
rz(-1.8975763) q[3];
sx q[3];
rz(2.2912912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0742764) q[0];
sx q[0];
rz(-2.5551676) q[0];
sx q[0];
rz(0.61074722) q[0];
rz(1.4432625) q[2];
sx q[2];
rz(-2.6886534) q[2];
sx q[2];
rz(2.4085101) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5634798) q[1];
sx q[1];
rz(-2.0260749) q[1];
sx q[1];
rz(-0.77225765) q[1];
rz(-pi) q[2];
rz(-0.38576491) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-0.5849134) q[2];
sx q[2];
rz(-1.9980105) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(2.7440199) q[0];
rz(1.827084) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.1869173) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2448954) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(-2.0969735) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.600012) q[2];
sx q[2];
rz(-1.6299106) q[2];
sx q[2];
rz(-0.37994775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4759051) q[1];
sx q[1];
rz(-2.396282) q[1];
sx q[1];
rz(1.3728632) q[1];
rz(-pi) q[2];
rz(0.62992923) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.8354841) q[2];
sx q[2];
rz(-1.7857893) q[2];
rz(1.6342182) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-2.142895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81925201) q[0];
sx q[0];
rz(-1.1871908) q[0];
sx q[0];
rz(2.2858455) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(2.1112679) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5949769) q[0];
sx q[0];
rz(-0.82406509) q[0];
sx q[0];
rz(-1.0975518) q[0];
x q[1];
rz(-1.6542997) q[2];
sx q[2];
rz(-2.6132772) q[2];
sx q[2];
rz(-0.63937843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47146637) q[1];
sx q[1];
rz(-1.867749) q[1];
sx q[1];
rz(-0.15836497) q[1];
rz(0.86352591) q[3];
sx q[3];
rz(-2.415495) q[3];
sx q[3];
rz(1.2442148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(-3.0714152) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.4708054) q[3];
sx q[3];
rz(-0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.4940015) q[0];
sx q[0];
rz(-1.7140472) q[0];
sx q[0];
rz(-2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(1.1351599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5495957) q[0];
sx q[0];
rz(-1.9599008) q[0];
sx q[0];
rz(-1.7745716) q[0];
x q[1];
rz(1.4831545) q[2];
sx q[2];
rz(-2.8200375) q[2];
sx q[2];
rz(-1.9784387) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70442048) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(-1.0874332) q[1];
rz(3.1177403) q[3];
sx q[3];
rz(-0.94787041) q[3];
sx q[3];
rz(0.21104392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51745522) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5448145) q[2];
rz(-2.4638713) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.90821663) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(-2.7897575) q[0];
rz(-0.31967638) q[1];
sx q[1];
rz(-2.7668178) q[1];
sx q[1];
rz(0.19616729) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10712121) q[0];
sx q[0];
rz(-0.84566085) q[0];
sx q[0];
rz(0.79191533) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7061383) q[2];
sx q[2];
rz(-1.6207098) q[2];
sx q[2];
rz(0.70336715) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-2.4073283) q[1];
sx q[1];
rz(1.5478565) q[1];
rz(-pi) q[2];
rz(-2.4862643) q[3];
sx q[3];
rz(-0.95643759) q[3];
sx q[3];
rz(2.4650246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7508042) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022973013) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(1.8267869) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(1.5224456) q[2];
sx q[2];
rz(-1.804525) q[2];
sx q[2];
rz(-1.3934025) q[2];
rz(3.0951981) q[3];
sx q[3];
rz(-1.3413324) q[3];
sx q[3];
rz(-3.0975292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
