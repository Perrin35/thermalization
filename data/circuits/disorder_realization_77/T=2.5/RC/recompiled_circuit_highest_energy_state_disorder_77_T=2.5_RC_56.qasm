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
rz(-0.4975118) q[0];
sx q[0];
rz(-1.8026135) q[0];
sx q[0];
rz(-2.8259377) q[0];
rz(-1.9714126) q[1];
sx q[1];
rz(-2.7538731) q[1];
sx q[1];
rz(2.5375836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8020129) q[0];
sx q[0];
rz(-2.8348742) q[0];
sx q[0];
rz(1.7392327) q[0];
x q[1];
rz(-0.96678218) q[2];
sx q[2];
rz(-0.76672115) q[2];
sx q[2];
rz(-1.3441835) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7249704) q[1];
sx q[1];
rz(-1.23471) q[1];
sx q[1];
rz(-1.027847) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2420869) q[3];
sx q[3];
rz(-2.8162873) q[3];
sx q[3];
rz(1.5055552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9903367) q[2];
sx q[2];
rz(-1.1948816) q[2];
sx q[2];
rz(-1.7830431) q[2];
rz(-1.5938866) q[3];
sx q[3];
rz(-0.93021506) q[3];
sx q[3];
rz(-1.6319298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5556521) q[0];
sx q[0];
rz(-2.5037615) q[0];
sx q[0];
rz(1.1974539) q[0];
rz(-2.6400631) q[1];
sx q[1];
rz(-1.5634368) q[1];
sx q[1];
rz(1.7053568) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2510808) q[0];
sx q[0];
rz(-0.62411004) q[0];
sx q[0];
rz(-1.1680383) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0244963) q[2];
sx q[2];
rz(-2.4352322) q[2];
sx q[2];
rz(-1.7603086) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3489248) q[1];
sx q[1];
rz(-1.7670005) q[1];
sx q[1];
rz(-1.6755107) q[1];
x q[2];
rz(-1.3615492) q[3];
sx q[3];
rz(-1.1236785) q[3];
sx q[3];
rz(-0.9613049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2531835) q[2];
sx q[2];
rz(-1.9056355) q[2];
sx q[2];
rz(-3.1186228) q[2];
rz(0.90211558) q[3];
sx q[3];
rz(-1.918856) q[3];
sx q[3];
rz(-0.42660108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70812923) q[0];
sx q[0];
rz(-1.7925649) q[0];
sx q[0];
rz(-0.9915114) q[0];
rz(-1.0333215) q[1];
sx q[1];
rz(-0.28305498) q[1];
sx q[1];
rz(1.2352157) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2184705) q[0];
sx q[0];
rz(-1.8184796) q[0];
sx q[0];
rz(-2.5367303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0472544) q[2];
sx q[2];
rz(-1.0152349) q[2];
sx q[2];
rz(-0.68651344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0924401) q[1];
sx q[1];
rz(-1.3235281) q[1];
sx q[1];
rz(-2.33806) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3091812) q[3];
sx q[3];
rz(-0.73736546) q[3];
sx q[3];
rz(-0.6399782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11423763) q[2];
sx q[2];
rz(-0.75216746) q[2];
sx q[2];
rz(0.99119622) q[2];
rz(1.2973971) q[3];
sx q[3];
rz(-1.8810561) q[3];
sx q[3];
rz(-1.4170925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39795136) q[0];
sx q[0];
rz(-2.209111) q[0];
sx q[0];
rz(0.81746307) q[0];
rz(-2.815333) q[1];
sx q[1];
rz(-1.1810818) q[1];
sx q[1];
rz(-2.356333) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1251251) q[0];
sx q[0];
rz(-2.7332691) q[0];
sx q[0];
rz(0.84618469) q[0];
rz(2.8815124) q[2];
sx q[2];
rz(-2.4818015) q[2];
sx q[2];
rz(1.545797) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98475115) q[1];
sx q[1];
rz(-1.9923216) q[1];
sx q[1];
rz(-2.0578029) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4912075) q[3];
sx q[3];
rz(-1.2168988) q[3];
sx q[3];
rz(-2.6864664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0205445) q[2];
sx q[2];
rz(-2.497017) q[2];
sx q[2];
rz(-2.5784967) q[2];
rz(2.2480887) q[3];
sx q[3];
rz(-1.1734791) q[3];
sx q[3];
rz(-1.9068498) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9017482) q[0];
sx q[0];
rz(-2.8185066) q[0];
sx q[0];
rz(-1.595994) q[0];
rz(-1.0218989) q[1];
sx q[1];
rz(-1.1232168) q[1];
sx q[1];
rz(0.18128577) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2342628) q[0];
sx q[0];
rz(-1.0815623) q[0];
sx q[0];
rz(-0.2184043) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0490169) q[2];
sx q[2];
rz(-0.47738722) q[2];
sx q[2];
rz(0.30711781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4605356) q[1];
sx q[1];
rz(-1.3486282) q[1];
sx q[1];
rz(-0.028789) q[1];
rz(-0.44293176) q[3];
sx q[3];
rz(-1.4287523) q[3];
sx q[3];
rz(1.0012817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.01866092) q[2];
sx q[2];
rz(-0.67492008) q[2];
sx q[2];
rz(-1.2459416) q[2];
rz(-2.5490226) q[3];
sx q[3];
rz(-1.6962681) q[3];
sx q[3];
rz(-0.84111253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5331921) q[0];
sx q[0];
rz(-2.1493981) q[0];
sx q[0];
rz(0.30884185) q[0];
rz(-2.8187075) q[1];
sx q[1];
rz(-0.37728089) q[1];
sx q[1];
rz(-0.13993941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9062146) q[0];
sx q[0];
rz(-1.8188634) q[0];
sx q[0];
rz(0.95309044) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0467911) q[2];
sx q[2];
rz(-1.1166414) q[2];
sx q[2];
rz(2.4081782) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6778826) q[1];
sx q[1];
rz(-1.8720798) q[1];
sx q[1];
rz(-2.3150285) q[1];
rz(-pi) q[2];
rz(-2.841137) q[3];
sx q[3];
rz(-2.2078035) q[3];
sx q[3];
rz(2.9954122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4796925) q[2];
sx q[2];
rz(-2.4883344) q[2];
sx q[2];
rz(0.24197401) q[2];
rz(0.74609977) q[3];
sx q[3];
rz(-1.7753764) q[3];
sx q[3];
rz(-1.411875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11874966) q[0];
sx q[0];
rz(-2.0976837) q[0];
sx q[0];
rz(-1.8705077) q[0];
rz(-0.56308833) q[1];
sx q[1];
rz(-0.92910281) q[1];
sx q[1];
rz(-1.44106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898687) q[0];
sx q[0];
rz(-1.2996724) q[0];
sx q[0];
rz(0.95377484) q[0];
rz(-pi) q[1];
rz(-1.4323498) q[2];
sx q[2];
rz(-1.9488153) q[2];
sx q[2];
rz(-1.0919309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5477834) q[1];
sx q[1];
rz(-1.154222) q[1];
sx q[1];
rz(1.3273456) q[1];
rz(-pi) q[2];
rz(2.4164124) q[3];
sx q[3];
rz(-1.4075052) q[3];
sx q[3];
rz(2.0189328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.21961221) q[2];
sx q[2];
rz(-1.2954804) q[2];
sx q[2];
rz(-1.8249576) q[2];
rz(0.5611788) q[3];
sx q[3];
rz(-2.1771274) q[3];
sx q[3];
rz(-1.5285899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.104326) q[0];
sx q[0];
rz(-2.9236561) q[0];
sx q[0];
rz(-0.88687801) q[0];
rz(0.2001702) q[1];
sx q[1];
rz(-1.472241) q[1];
sx q[1];
rz(1.0221457) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6157142) q[0];
sx q[0];
rz(-1.7019042) q[0];
sx q[0];
rz(0.025630533) q[0];
rz(-pi) q[1];
rz(-2.3365302) q[2];
sx q[2];
rz(-2.2230122) q[2];
sx q[2];
rz(2.7307939) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8274668) q[1];
sx q[1];
rz(-2.3955401) q[1];
sx q[1];
rz(-1.57342) q[1];
rz(1.2194013) q[3];
sx q[3];
rz(-2.1179652) q[3];
sx q[3];
rz(-1.6148293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48978051) q[2];
sx q[2];
rz(-2.688789) q[2];
sx q[2];
rz(0.81857267) q[2];
rz(-0.98322785) q[3];
sx q[3];
rz(-2.1849617) q[3];
sx q[3];
rz(2.6035068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1061358) q[0];
sx q[0];
rz(-1.9945194) q[0];
sx q[0];
rz(1.5863093) q[0];
rz(-2.4447794) q[1];
sx q[1];
rz(-2.9882444) q[1];
sx q[1];
rz(-2.3568025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12778388) q[0];
sx q[0];
rz(-2.0079566) q[0];
sx q[0];
rz(-0.60370914) q[0];
x q[1];
rz(-1.3680787) q[2];
sx q[2];
rz(-1.7200816) q[2];
sx q[2];
rz(1.9307856) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0052629) q[1];
sx q[1];
rz(-1.6044093) q[1];
sx q[1];
rz(2.7064347) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56702964) q[3];
sx q[3];
rz(-2.4483557) q[3];
sx q[3];
rz(0.45613134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67184225) q[2];
sx q[2];
rz(-1.8639114) q[2];
sx q[2];
rz(-0.81076852) q[2];
rz(-1.9226711) q[3];
sx q[3];
rz(-2.6007077) q[3];
sx q[3];
rz(2.0764652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588381) q[0];
sx q[0];
rz(-2.8420119) q[0];
sx q[0];
rz(1.4240356) q[0];
rz(0.77761039) q[1];
sx q[1];
rz(-0.97336665) q[1];
sx q[1];
rz(2.3861859) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7292405) q[0];
sx q[0];
rz(-1.5753813) q[0];
sx q[0];
rz(1.3608749) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.099775) q[2];
sx q[2];
rz(-0.39678573) q[2];
sx q[2];
rz(2.7355742) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3302119) q[1];
sx q[1];
rz(-0.63618681) q[1];
sx q[1];
rz(2.5274171) q[1];
x q[2];
rz(2.741119) q[3];
sx q[3];
rz(-2.4254834) q[3];
sx q[3];
rz(-3.0587089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5648254) q[2];
sx q[2];
rz(-0.29594031) q[2];
sx q[2];
rz(0.15288615) q[2];
rz(1.2711924) q[3];
sx q[3];
rz(-1.8891687) q[3];
sx q[3];
rz(2.474474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85970238) q[0];
sx q[0];
rz(-2.6489881) q[0];
sx q[0];
rz(-0.76982605) q[0];
rz(1.9181171) q[1];
sx q[1];
rz(-1.2212831) q[1];
sx q[1];
rz(2.4881359) q[1];
rz(2.5443947) q[2];
sx q[2];
rz(-1.2033249) q[2];
sx q[2];
rz(-2.0133599) q[2];
rz(1.9742692) q[3];
sx q[3];
rz(-1.8971761) q[3];
sx q[3];
rz(-1.182233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
