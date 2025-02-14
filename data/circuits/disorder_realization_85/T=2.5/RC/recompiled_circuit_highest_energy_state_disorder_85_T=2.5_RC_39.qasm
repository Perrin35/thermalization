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
rz(-2.1674018) q[0];
sx q[0];
rz(-0.84889698) q[0];
sx q[0];
rz(1.4155686) q[0];
rz(2.2805136) q[1];
sx q[1];
rz(-1.5134892) q[1];
sx q[1];
rz(2.3754062) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9250037) q[0];
sx q[0];
rz(-1.6356409) q[0];
sx q[0];
rz(-1.606604) q[0];
rz(-pi) q[1];
rz(-1.6915198) q[2];
sx q[2];
rz(-1.6044464) q[2];
sx q[2];
rz(-0.41481123) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3688226) q[1];
sx q[1];
rz(-1.5632679) q[1];
sx q[1];
rz(-0.20329411) q[1];
x q[2];
rz(0.10786003) q[3];
sx q[3];
rz(-1.1601891) q[3];
sx q[3];
rz(1.6317898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8908995) q[2];
sx q[2];
rz(-2.3919899) q[2];
sx q[2];
rz(-3.0787943) q[2];
rz(-1.7167669) q[3];
sx q[3];
rz(-1.3751605) q[3];
sx q[3];
rz(-2.4019901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0733136) q[0];
sx q[0];
rz(-2.8977019) q[0];
sx q[0];
rz(2.0252315) q[0];
rz(1.5694537) q[1];
sx q[1];
rz(-2.8469323) q[1];
sx q[1];
rz(2.0223298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.946944) q[0];
sx q[0];
rz(-0.70746052) q[0];
sx q[0];
rz(0.0012673541) q[0];
x q[1];
rz(-1.0874027) q[2];
sx q[2];
rz(-1.2366939) q[2];
sx q[2];
rz(2.7887864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77308622) q[1];
sx q[1];
rz(-2.0952053) q[1];
sx q[1];
rz(-1.4953411) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4495115) q[3];
sx q[3];
rz(-1.0861703) q[3];
sx q[3];
rz(-1.6177819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0548627) q[2];
sx q[2];
rz(-2.3708998) q[2];
sx q[2];
rz(1.1718303) q[2];
rz(-2.150676) q[3];
sx q[3];
rz(-1.0057697) q[3];
sx q[3];
rz(1.2539554) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9988001) q[0];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9799412) q[0];
sx q[0];
rz(-1.5197541) q[0];
sx q[0];
rz(-1.4568842) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5722857) q[2];
sx q[2];
rz(-2.226839) q[2];
sx q[2];
rz(1.4106205) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5873196) q[1];
sx q[1];
rz(-1.6780507) q[1];
sx q[1];
rz(0.3078021) q[1];
rz(1.6693657) q[3];
sx q[3];
rz(-0.53547066) q[3];
sx q[3];
rz(2.5853973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.74505836) q[2];
sx q[2];
rz(-1.2816387) q[2];
sx q[2];
rz(2.7520531) q[2];
rz(2.7430096) q[3];
sx q[3];
rz(-1.7321209) q[3];
sx q[3];
rz(-1.7194974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.013537708) q[0];
sx q[0];
rz(-0.062352926) q[0];
sx q[0];
rz(0.41559497) q[0];
rz(1.2526814) q[1];
sx q[1];
rz(-2.0038192) q[1];
sx q[1];
rz(-2.2223053) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18596953) q[0];
sx q[0];
rz(-0.16849314) q[0];
sx q[0];
rz(1.3244434) q[0];
x q[1];
rz(2.8785673) q[2];
sx q[2];
rz(-1.5925455) q[2];
sx q[2];
rz(-1.1794752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7953736) q[1];
sx q[1];
rz(-2.260403) q[1];
sx q[1];
rz(-2.3350689) q[1];
x q[2];
rz(2.3727338) q[3];
sx q[3];
rz(-1.9058085) q[3];
sx q[3];
rz(1.5227536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2933942) q[2];
sx q[2];
rz(-1.5092809) q[2];
sx q[2];
rz(-2.1634114) q[2];
rz(2.6350422) q[3];
sx q[3];
rz(-1.6440441) q[3];
sx q[3];
rz(1.0165366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070234805) q[0];
sx q[0];
rz(-0.6986343) q[0];
sx q[0];
rz(0.51682669) q[0];
rz(1.0454073) q[1];
sx q[1];
rz(-1.005859) q[1];
sx q[1];
rz(-2.7396835) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73905009) q[0];
sx q[0];
rz(-2.2687904) q[0];
sx q[0];
rz(-2.151346) q[0];
rz(-1.8568618) q[2];
sx q[2];
rz(-1.4081883) q[2];
sx q[2];
rz(0.78078496) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13788183) q[1];
sx q[1];
rz(-1.5665083) q[1];
sx q[1];
rz(-1.3715308) q[1];
rz(-pi) q[2];
rz(-1.4001717) q[3];
sx q[3];
rz(-2.5192891) q[3];
sx q[3];
rz(0.69392747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7261293) q[2];
sx q[2];
rz(-1.1795421) q[2];
sx q[2];
rz(0.22280517) q[2];
rz(0.86083096) q[3];
sx q[3];
rz(-1.5102883) q[3];
sx q[3];
rz(-1.8676602) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.05805) q[0];
sx q[0];
rz(-1.2419751) q[0];
sx q[0];
rz(-0.85839957) q[0];
rz(2.3421613) q[1];
sx q[1];
rz(-2.1245427) q[1];
sx q[1];
rz(1.3851091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0388865) q[0];
sx q[0];
rz(-1.2974076) q[0];
sx q[0];
rz(-2.9486743) q[0];
x q[1];
rz(1.8968514) q[2];
sx q[2];
rz(-1.5189803) q[2];
sx q[2];
rz(-2.5246132) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68011682) q[1];
sx q[1];
rz(-1.6915186) q[1];
sx q[1];
rz(0.6013479) q[1];
rz(2.9030209) q[3];
sx q[3];
rz(-2.0485463) q[3];
sx q[3];
rz(2.9693317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4110306) q[2];
sx q[2];
rz(-2.4081814) q[2];
sx q[2];
rz(-0.15156558) q[2];
rz(1.1847121) q[3];
sx q[3];
rz(-1.1368753) q[3];
sx q[3];
rz(0.90004164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93722349) q[0];
sx q[0];
rz(-1.3441514) q[0];
sx q[0];
rz(0.7835266) q[0];
rz(-0.20768684) q[1];
sx q[1];
rz(-0.91778225) q[1];
sx q[1];
rz(1.2088561) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75654001) q[0];
sx q[0];
rz(-1.970743) q[0];
sx q[0];
rz(0.51838309) q[0];
rz(-pi) q[1];
rz(0.54401347) q[2];
sx q[2];
rz(-2.233783) q[2];
sx q[2];
rz(-1.3617409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4502984) q[1];
sx q[1];
rz(-2.2992059) q[1];
sx q[1];
rz(-0.052740514) q[1];
rz(-1.9586669) q[3];
sx q[3];
rz(-0.3976477) q[3];
sx q[3];
rz(-0.39163799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2850767) q[2];
sx q[2];
rz(-2.3975942) q[2];
sx q[2];
rz(3.0169955) q[2];
rz(1.8018657) q[3];
sx q[3];
rz(-2.7281269) q[3];
sx q[3];
rz(-1.9477897) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1272005) q[0];
sx q[0];
rz(-1.4407644) q[0];
sx q[0];
rz(3.0442687) q[0];
rz(-0.65517455) q[1];
sx q[1];
rz(-2.580878) q[1];
sx q[1];
rz(-0.68897828) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6902721) q[0];
sx q[0];
rz(-0.49271944) q[0];
sx q[0];
rz(1.413233) q[0];
rz(1.9227486) q[2];
sx q[2];
rz(-1.5492288) q[2];
sx q[2];
rz(-1.5810313) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.6072907) q[1];
sx q[1];
rz(-1.698415) q[1];
sx q[1];
rz(-2.4625596) q[1];
x q[2];
rz(1.407988) q[3];
sx q[3];
rz(-2.2534784) q[3];
sx q[3];
rz(-2.9014661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.24136209) q[2];
sx q[2];
rz(-2.1817544) q[2];
sx q[2];
rz(-2.0755365) q[2];
rz(2.7212932) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(-1.458781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0578388) q[0];
sx q[0];
rz(-2.3729615) q[0];
sx q[0];
rz(0.88430697) q[0];
rz(-0.090944313) q[1];
sx q[1];
rz(-0.93946409) q[1];
sx q[1];
rz(-1.9627862) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509277) q[0];
sx q[0];
rz(-0.17296436) q[0];
sx q[0];
rz(0.38562445) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2366899) q[2];
sx q[2];
rz(-1.5690815) q[2];
sx q[2];
rz(0.26632613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2388666) q[1];
sx q[1];
rz(-1.3999434) q[1];
sx q[1];
rz(0.56436964) q[1];
rz(-1.4512153) q[3];
sx q[3];
rz(-1.371583) q[3];
sx q[3];
rz(-1.8389033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3921597) q[2];
sx q[2];
rz(-1.9955111) q[2];
sx q[2];
rz(0.47449365) q[2];
rz(3.0992881) q[3];
sx q[3];
rz(-1.5174805) q[3];
sx q[3];
rz(-1.8208985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6831191) q[0];
sx q[0];
rz(-2.0404158) q[0];
sx q[0];
rz(-1.3018357) q[0];
rz(-2.5104972) q[1];
sx q[1];
rz(-2.0986235) q[1];
sx q[1];
rz(-2.1249117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6852208) q[0];
sx q[0];
rz(-0.70661649) q[0];
sx q[0];
rz(-2.0232692) q[0];
rz(-3.0877168) q[2];
sx q[2];
rz(-1.6695256) q[2];
sx q[2];
rz(-1.3786664) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70003033) q[1];
sx q[1];
rz(-2.4228281) q[1];
sx q[1];
rz(1.584254) q[1];
rz(0.74211699) q[3];
sx q[3];
rz(-2.285503) q[3];
sx q[3];
rz(1.8736739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8299228) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(0.97490087) q[2];
rz(-1.4026862) q[3];
sx q[3];
rz(-2.1699173) q[3];
sx q[3];
rz(-1.5909125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1350031) q[0];
sx q[0];
rz(-2.4000744) q[0];
sx q[0];
rz(0.97493521) q[0];
rz(1.4047752) q[1];
sx q[1];
rz(-1.4092696) q[1];
sx q[1];
rz(-0.41913941) q[1];
rz(-2.2683138) q[2];
sx q[2];
rz(-2.4200405) q[2];
sx q[2];
rz(2.6478052) q[2];
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
