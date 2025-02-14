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
rz(0.97419089) q[0];
sx q[0];
rz(-2.2926957) q[0];
sx q[0];
rz(1.726024) q[0];
rz(-0.8610791) q[1];
sx q[1];
rz(1.5134892) q[1];
sx q[1];
rz(8.6585915) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7216033) q[0];
sx q[0];
rz(-3.0675305) q[0];
sx q[0];
rz(2.6377399) q[0];
rz(-3.1076961) q[2];
sx q[2];
rz(-1.4501415) q[2];
sx q[2];
rz(-1.1519037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76152245) q[1];
sx q[1];
rz(-0.20343154) q[1];
sx q[1];
rz(0.037271715) q[1];
rz(-pi) q[2];
rz(1.3283861) q[3];
sx q[3];
rz(-2.7178353) q[3];
sx q[3];
rz(-1.3668982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2506931) q[2];
sx q[2];
rz(-0.74960274) q[2];
sx q[2];
rz(0.062798373) q[2];
rz(1.7167669) q[3];
sx q[3];
rz(-1.7664322) q[3];
sx q[3];
rz(0.73960251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0733136) q[0];
sx q[0];
rz(-0.24389076) q[0];
sx q[0];
rz(2.0252315) q[0];
rz(-1.5721389) q[1];
sx q[1];
rz(-2.8469323) q[1];
sx q[1];
rz(2.0223298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9452764) q[0];
sx q[0];
rz(-0.8633365) q[0];
sx q[0];
rz(1.5697126) q[0];
rz(-pi) q[1];
rz(-2.767973) q[2];
sx q[2];
rz(-1.116215) q[2];
sx q[2];
rz(-1.0475243) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5183603) q[1];
sx q[1];
rz(-2.6122852) q[1];
sx q[1];
rz(3.0119986) q[1];
rz(-pi) q[2];
rz(2.4517566) q[3];
sx q[3];
rz(-0.82131006) q[3];
sx q[3];
rz(-0.46508712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.08673) q[2];
sx q[2];
rz(-0.77069288) q[2];
sx q[2];
rz(-1.9697624) q[2];
rz(2.150676) q[3];
sx q[3];
rz(-2.1358229) q[3];
sx q[3];
rz(1.2539554) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1427926) q[0];
sx q[0];
rz(-0.87738377) q[0];
sx q[0];
rz(0.073609322) q[0];
rz(0.40120801) q[1];
sx q[1];
rz(-2.7707272) q[1];
sx q[1];
rz(2.4295095) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3129666) q[0];
sx q[0];
rz(-0.12477984) q[0];
sx q[0];
rz(-1.1484042) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0019348423) q[2];
sx q[2];
rz(-0.65604416) q[2];
sx q[2];
rz(1.7334138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.69171762) q[1];
sx q[1];
rz(-0.32539168) q[1];
sx q[1];
rz(-2.8001333) q[1];
rz(-pi) q[2];
rz(-0.058319326) q[3];
sx q[3];
rz(-2.1033896) q[3];
sx q[3];
rz(0.67067702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.74505836) q[2];
sx q[2];
rz(-1.859954) q[2];
sx q[2];
rz(0.38953951) q[2];
rz(0.39858308) q[3];
sx q[3];
rz(-1.4094718) q[3];
sx q[3];
rz(-1.7194974) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1280549) q[0];
sx q[0];
rz(-0.062352926) q[0];
sx q[0];
rz(-0.41559497) q[0];
rz(1.8889113) q[1];
sx q[1];
rz(-1.1377734) q[1];
sx q[1];
rz(0.91928732) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.141826) q[0];
sx q[0];
rz(-1.6117038) q[0];
sx q[0];
rz(1.4072988) q[0];
rz(-pi) q[1];
rz(-0.083468632) q[2];
sx q[2];
rz(-0.26390227) q[2];
sx q[2];
rz(2.6696799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34621908) q[1];
sx q[1];
rz(-0.88118964) q[1];
sx q[1];
rz(2.3350689) q[1];
rz(-2.6773964) q[3];
sx q[3];
rz(-2.3168543) q[3];
sx q[3];
rz(-0.27955258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84819841) q[2];
sx q[2];
rz(-1.6323117) q[2];
sx q[2];
rz(-2.1634114) q[2];
rz(0.50655043) q[3];
sx q[3];
rz(-1.4975486) q[3];
sx q[3];
rz(-2.1250561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070234805) q[0];
sx q[0];
rz(-0.6986343) q[0];
sx q[0];
rz(-0.51682669) q[0];
rz(-2.0961854) q[1];
sx q[1];
rz(-1.005859) q[1];
sx q[1];
rz(0.4019092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43278431) q[0];
sx q[0];
rz(-1.1371181) q[0];
sx q[0];
rz(-2.3545803) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0974196) q[2];
sx q[2];
rz(-2.8136471) q[2];
sx q[2];
rz(-1.8484268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.13788183) q[1];
sx q[1];
rz(-1.5750843) q[1];
sx q[1];
rz(1.3715308) q[1];
rz(0.95540445) q[3];
sx q[3];
rz(-1.6699353) q[3];
sx q[3];
rz(-2.4038199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7261293) q[2];
sx q[2];
rz(-1.9620506) q[2];
sx q[2];
rz(-0.22280517) q[2];
rz(2.2807617) q[3];
sx q[3];
rz(-1.5102883) q[3];
sx q[3];
rz(-1.2739325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.083542682) q[0];
sx q[0];
rz(-1.2419751) q[0];
sx q[0];
rz(-0.85839957) q[0];
rz(0.79943132) q[1];
sx q[1];
rz(-1.0170499) q[1];
sx q[1];
rz(1.3851091) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0388865) q[0];
sx q[0];
rz(-1.2974076) q[0];
sx q[0];
rz(0.1929184) q[0];
rz(-pi) q[1];
rz(1.7313174) q[2];
sx q[2];
rz(-0.33000144) q[2];
sx q[2];
rz(2.3397719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3333545) q[1];
sx q[1];
rz(-0.97443354) q[1];
sx q[1];
rz(1.4247231) q[1];
x q[2];
rz(-1.9989789) q[3];
sx q[3];
rz(-0.52985668) q[3];
sx q[3];
rz(2.4827905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4110306) q[2];
sx q[2];
rz(-2.4081814) q[2];
sx q[2];
rz(0.15156558) q[2];
rz(-1.9568806) q[3];
sx q[3];
rz(-1.1368753) q[3];
sx q[3];
rz(0.90004164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93722349) q[0];
sx q[0];
rz(-1.7974412) q[0];
sx q[0];
rz(-2.3580661) q[0];
rz(0.20768684) q[1];
sx q[1];
rz(-0.91778225) q[1];
sx q[1];
rz(-1.2088561) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1087813) q[0];
sx q[0];
rz(-2.0447123) q[0];
sx q[0];
rz(-1.117871) q[0];
x q[1];
rz(0.98548205) q[2];
sx q[2];
rz(-0.83067465) q[2];
sx q[2];
rz(-2.5565846) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5294339) q[1];
sx q[1];
rz(-2.4116257) q[1];
sx q[1];
rz(1.5117701) q[1];
rz(-pi) q[2];
rz(1.1999628) q[3];
sx q[3];
rz(-1.4238024) q[3];
sx q[3];
rz(-1.6021836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2850767) q[2];
sx q[2];
rz(-2.3975942) q[2];
sx q[2];
rz(0.12459717) q[2];
rz(-1.3397269) q[3];
sx q[3];
rz(-0.41346574) q[3];
sx q[3];
rz(-1.193803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0143921) q[0];
sx q[0];
rz(-1.4407644) q[0];
sx q[0];
rz(0.097323962) q[0];
rz(-2.4864181) q[1];
sx q[1];
rz(-0.56071463) q[1];
sx q[1];
rz(2.4526144) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5118588) q[0];
sx q[0];
rz(-1.0847158) q[0];
sx q[0];
rz(3.0575471) q[0];
x q[1];
rz(-1.6332878) q[2];
sx q[2];
rz(-2.7890076) q[2];
sx q[2];
rz(-3.0931713) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1199159) q[1];
sx q[1];
rz(-0.6890474) q[1];
sx q[1];
rz(-2.9400548) q[1];
x q[2];
rz(0.19676713) q[3];
sx q[3];
rz(-0.69878687) q[3];
sx q[3];
rz(-3.1270032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24136209) q[2];
sx q[2];
rz(-0.9598383) q[2];
sx q[2];
rz(1.0660561) q[2];
rz(-2.7212932) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(-1.6828116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578388) q[0];
sx q[0];
rz(-2.3729615) q[0];
sx q[0];
rz(2.2572857) q[0];
rz(-0.090944313) q[1];
sx q[1];
rz(-0.93946409) q[1];
sx q[1];
rz(-1.9627862) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6323157) q[0];
sx q[0];
rz(-0.17296436) q[0];
sx q[0];
rz(2.7559682) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5680204) q[2];
sx q[2];
rz(-0.66589543) q[2];
sx q[2];
rz(-1.3022873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.072047) q[1];
sx q[1];
rz(-0.58696771) q[1];
sx q[1];
rz(-2.8295641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9409845) q[3];
sx q[3];
rz(-1.4535913) q[3];
sx q[3];
rz(2.8497118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74943298) q[2];
sx q[2];
rz(-1.9955111) q[2];
sx q[2];
rz(-0.47449365) q[2];
rz(3.0992881) q[3];
sx q[3];
rz(-1.5174805) q[3];
sx q[3];
rz(-1.8208985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6831191) q[0];
sx q[0];
rz(-1.1011769) q[0];
sx q[0];
rz(-1.3018357) q[0];
rz(0.63109541) q[1];
sx q[1];
rz(-1.0429691) q[1];
sx q[1];
rz(-1.016681) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4563718) q[0];
sx q[0];
rz(-2.4349762) q[0];
sx q[0];
rz(-2.0232692) q[0];
x q[1];
rz(2.0687556) q[2];
sx q[2];
rz(-3.0291621) q[2];
sx q[2];
rz(0.87804483) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70003033) q[1];
sx q[1];
rz(-0.71876457) q[1];
sx q[1];
rz(1.584254) q[1];
x q[2];
rz(2.437463) q[3];
sx q[3];
rz(-2.1065154) q[3];
sx q[3];
rz(2.2975722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3116698) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(2.1666918) q[2];
rz(-1.7389065) q[3];
sx q[3];
rz(-0.97167531) q[3];
sx q[3];
rz(-1.5909125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0065895157) q[0];
sx q[0];
rz(-0.74151825) q[0];
sx q[0];
rz(-2.1666574) q[0];
rz(-1.7368175) q[1];
sx q[1];
rz(-1.4092696) q[1];
sx q[1];
rz(-0.41913941) q[1];
rz(-2.1640833) q[2];
sx q[2];
rz(-2.0089663) q[2];
sx q[2];
rz(0.51539863) q[2];
rz(1.6831123) q[3];
sx q[3];
rz(-1.1426123) q[3];
sx q[3];
rz(-1.860426) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
