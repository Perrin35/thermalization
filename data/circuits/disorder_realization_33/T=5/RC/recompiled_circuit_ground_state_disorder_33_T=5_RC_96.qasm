OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44242087) q[0];
sx q[0];
rz(-2.3306263) q[0];
sx q[0];
rz(2.6851658) q[0];
rz(2.2189848) q[1];
sx q[1];
rz(-0.95093095) q[1];
sx q[1];
rz(2.9589597) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62902495) q[0];
sx q[0];
rz(-2.1898666) q[0];
sx q[0];
rz(-2.7946212) q[0];
rz(0.47477291) q[2];
sx q[2];
rz(-0.92657303) q[2];
sx q[2];
rz(-0.34851532) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.071242407) q[1];
sx q[1];
rz(-1.6082941) q[1];
sx q[1];
rz(2.9686023) q[1];
x q[2];
rz(1.1785281) q[3];
sx q[3];
rz(-1.0844106) q[3];
sx q[3];
rz(1.148759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7543588) q[2];
sx q[2];
rz(-2.0626455) q[2];
sx q[2];
rz(-0.31164247) q[2];
rz(-1.8841057) q[3];
sx q[3];
rz(-0.25185549) q[3];
sx q[3];
rz(0.18865147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1455014) q[0];
sx q[0];
rz(-1.4459193) q[0];
sx q[0];
rz(1.2392932) q[0];
rz(2.7481825) q[1];
sx q[1];
rz(-2.0937803) q[1];
sx q[1];
rz(1.0308824) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6604583) q[0];
sx q[0];
rz(-1.0054614) q[0];
sx q[0];
rz(-2.7271366) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8003045) q[2];
sx q[2];
rz(-1.82845) q[2];
sx q[2];
rz(2.9837556) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0075025) q[1];
sx q[1];
rz(-2.1975027) q[1];
sx q[1];
rz(0.43372633) q[1];
x q[2];
rz(0.88406422) q[3];
sx q[3];
rz(-2.2746448) q[3];
sx q[3];
rz(-0.91020179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82765141) q[2];
sx q[2];
rz(-0.88259077) q[2];
sx q[2];
rz(0.35448709) q[2];
rz(0.83141023) q[3];
sx q[3];
rz(-2.3737213) q[3];
sx q[3];
rz(-1.6262511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671234) q[0];
sx q[0];
rz(-0.18018436) q[0];
sx q[0];
rz(2.6572976) q[0];
rz(2.2593185) q[1];
sx q[1];
rz(-2.2703998) q[1];
sx q[1];
rz(-2.5994515) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4216293) q[0];
sx q[0];
rz(-2.6609592) q[0];
sx q[0];
rz(2.292657) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9464821) q[2];
sx q[2];
rz(-0.32719041) q[2];
sx q[2];
rz(-1.3766152) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0391991) q[1];
sx q[1];
rz(-2.3619283) q[1];
sx q[1];
rz(-1.3179661) q[1];
x q[2];
rz(-3.077742) q[3];
sx q[3];
rz(-1.333916) q[3];
sx q[3];
rz(-2.3781535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3279646) q[2];
sx q[2];
rz(-0.70983228) q[2];
sx q[2];
rz(-2.1072809) q[2];
rz(-1.3616925) q[3];
sx q[3];
rz(-1.1797649) q[3];
sx q[3];
rz(-2.5907607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4807602) q[0];
sx q[0];
rz(-1.3897422) q[0];
sx q[0];
rz(-1.047629) q[0];
rz(1.3230336) q[1];
sx q[1];
rz(-2.5593457) q[1];
sx q[1];
rz(-0.38527647) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.229363) q[0];
sx q[0];
rz(-1.1426539) q[0];
sx q[0];
rz(-2.1332801) q[0];
rz(-3.0550256) q[2];
sx q[2];
rz(-0.7184775) q[2];
sx q[2];
rz(-0.95667808) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.99285728) q[1];
sx q[1];
rz(-0.22051935) q[1];
sx q[1];
rz(0.33111568) q[1];
rz(3.0317467) q[3];
sx q[3];
rz(-0.89021909) q[3];
sx q[3];
rz(-2.3504013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.72383991) q[2];
sx q[2];
rz(-0.94261348) q[2];
sx q[2];
rz(0.27457944) q[2];
rz(2.5802021) q[3];
sx q[3];
rz(-2.0761108) q[3];
sx q[3];
rz(1.6339462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1356337) q[0];
sx q[0];
rz(-2.5649286) q[0];
sx q[0];
rz(0.86719257) q[0];
rz(2.7658956) q[1];
sx q[1];
rz(-2.5521894) q[1];
sx q[1];
rz(-1.9120749) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50841416) q[0];
sx q[0];
rz(-1.408218) q[0];
sx q[0];
rz(0.44277097) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0528238) q[2];
sx q[2];
rz(-1.8395632) q[2];
sx q[2];
rz(2.4160699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.873949) q[1];
sx q[1];
rz(-2.3726844) q[1];
sx q[1];
rz(-2.257009) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8939871) q[3];
sx q[3];
rz(-0.5115307) q[3];
sx q[3];
rz(-2.3681896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1162954) q[2];
sx q[2];
rz(-1.2325492) q[2];
sx q[2];
rz(-0.06289014) q[2];
rz(-0.40999117) q[3];
sx q[3];
rz(-2.4153109) q[3];
sx q[3];
rz(0.15967742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9922239) q[0];
sx q[0];
rz(-0.069644444) q[0];
sx q[0];
rz(2.3058291) q[0];
rz(1.1831076) q[1];
sx q[1];
rz(-1.2703905) q[1];
sx q[1];
rz(2.3979208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1069458) q[0];
sx q[0];
rz(-1.8726882) q[0];
sx q[0];
rz(0.56580122) q[0];
rz(1.9299149) q[2];
sx q[2];
rz(-0.15378498) q[2];
sx q[2];
rz(-1.502047) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2433155) q[1];
sx q[1];
rz(-0.87496266) q[1];
sx q[1];
rz(-2.9030373) q[1];
rz(-0.25571574) q[3];
sx q[3];
rz(-1.6220495) q[3];
sx q[3];
rz(1.071614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8379197) q[2];
sx q[2];
rz(-0.49771365) q[2];
sx q[2];
rz(0.79353235) q[2];
rz(0.41905904) q[3];
sx q[3];
rz(-1.4895118) q[3];
sx q[3];
rz(-0.52217531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1917052) q[0];
sx q[0];
rz(-0.97923034) q[0];
sx q[0];
rz(1.9784084) q[0];
rz(-0.78041068) q[1];
sx q[1];
rz(-0.33640948) q[1];
sx q[1];
rz(-1.6132678) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53188092) q[0];
sx q[0];
rz(-1.5059294) q[0];
sx q[0];
rz(2.2676629) q[0];
rz(-0.18272551) q[2];
sx q[2];
rz(-1.2264681) q[2];
sx q[2];
rz(1.8553714) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0039499) q[1];
sx q[1];
rz(-2.9423135) q[1];
sx q[1];
rz(-1.4635221) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27604528) q[3];
sx q[3];
rz(-1.9791043) q[3];
sx q[3];
rz(0.21430548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10716001) q[2];
sx q[2];
rz(-1.1792504) q[2];
sx q[2];
rz(3.072928) q[2];
rz(-2.5717403) q[3];
sx q[3];
rz(-2.6611501) q[3];
sx q[3];
rz(0.3705875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625075) q[0];
sx q[0];
rz(-0.67254368) q[0];
sx q[0];
rz(1.8261209) q[0];
rz(-1.2830118) q[1];
sx q[1];
rz(-0.43671572) q[1];
sx q[1];
rz(0.049093094) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9271761) q[0];
sx q[0];
rz(-0.11888725) q[0];
sx q[0];
rz(2.2772339) q[0];
rz(-pi) q[1];
rz(-2.2678864) q[2];
sx q[2];
rz(-2.3915136) q[2];
sx q[2];
rz(1.9421792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9994574) q[1];
sx q[1];
rz(-1.189636) q[1];
sx q[1];
rz(2.3232404) q[1];
x q[2];
rz(-3.0144948) q[3];
sx q[3];
rz(-0.93080257) q[3];
sx q[3];
rz(2.618263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38356885) q[2];
sx q[2];
rz(-2.263676) q[2];
sx q[2];
rz(-0.54086584) q[2];
rz(-2.0567549) q[3];
sx q[3];
rz(-0.70298755) q[3];
sx q[3];
rz(1.7613523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2329907) q[0];
sx q[0];
rz(-2.1451696) q[0];
sx q[0];
rz(-3.0365699) q[0];
rz(-2.5999293) q[1];
sx q[1];
rz(-2.2553406) q[1];
sx q[1];
rz(-0.35596102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5830688) q[0];
sx q[0];
rz(-1.6083058) q[0];
sx q[0];
rz(1.4315228) q[0];
x q[1];
rz(2.2952406) q[2];
sx q[2];
rz(-1.0484107) q[2];
sx q[2];
rz(1.236793) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7193422) q[1];
sx q[1];
rz(-2.1394891) q[1];
sx q[1];
rz(0.45180068) q[1];
rz(-1.2139236) q[3];
sx q[3];
rz(-0.3398474) q[3];
sx q[3];
rz(1.7879037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9291222) q[2];
sx q[2];
rz(-1.4996303) q[2];
sx q[2];
rz(-1.8222202) q[2];
rz(1.3245373) q[3];
sx q[3];
rz(-1.5732485) q[3];
sx q[3];
rz(2.1738079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20973715) q[0];
sx q[0];
rz(-3.1071438) q[0];
sx q[0];
rz(1.4631648) q[0];
rz(-1.4596918) q[1];
sx q[1];
rz(-1.1594783) q[1];
sx q[1];
rz(0.34585888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1347156) q[0];
sx q[0];
rz(-1.6505989) q[0];
sx q[0];
rz(0.047264506) q[0];
rz(-pi) q[1];
rz(1.5058636) q[2];
sx q[2];
rz(-1.254515) q[2];
sx q[2];
rz(2.8313178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74259669) q[1];
sx q[1];
rz(-2.1861939) q[1];
sx q[1];
rz(2.1324498) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8000051) q[3];
sx q[3];
rz(-0.95529592) q[3];
sx q[3];
rz(3.0581491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86089245) q[2];
sx q[2];
rz(-1.2714551) q[2];
sx q[2];
rz(0.26091179) q[2];
rz(-1.4704618) q[3];
sx q[3];
rz(-1.4025531) q[3];
sx q[3];
rz(-1.7117333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3093001) q[0];
sx q[0];
rz(-2.0335048) q[0];
sx q[0];
rz(-0.19620398) q[0];
rz(0.55083864) q[1];
sx q[1];
rz(-1.6549587) q[1];
sx q[1];
rz(0.5400198) q[1];
rz(-0.2395128) q[2];
sx q[2];
rz(-0.85776599) q[2];
sx q[2];
rz(0.89606482) q[2];
rz(-2.8367219) q[3];
sx q[3];
rz(-0.75406995) q[3];
sx q[3];
rz(-2.0537805) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
