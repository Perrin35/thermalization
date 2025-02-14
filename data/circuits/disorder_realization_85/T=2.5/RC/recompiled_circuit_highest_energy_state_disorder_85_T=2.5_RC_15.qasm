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
rz(-1.6281035) q[1];
sx q[1];
rz(-2.3754062) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.216589) q[0];
sx q[0];
rz(-1.6356409) q[0];
sx q[0];
rz(1.5349887) q[0];
rz(-pi) q[1];
rz(-1.2982326) q[2];
sx q[2];
rz(-3.0162891) q[2];
sx q[2];
rz(-1.426515) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3420145) q[1];
sx q[1];
rz(-1.7740846) q[1];
sx q[1];
rz(1.5631097) q[1];
rz(-pi) q[2];
rz(3.0337326) q[3];
sx q[3];
rz(-1.9814035) q[3];
sx q[3];
rz(1.6317898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2506931) q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06827908) q[0];
sx q[0];
rz(-2.8977019) q[0];
sx q[0];
rz(-2.0252315) q[0];
rz(-1.5694537) q[1];
sx q[1];
rz(-2.8469323) q[1];
sx q[1];
rz(1.1192628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1946487) q[0];
sx q[0];
rz(-0.70746052) q[0];
sx q[0];
rz(0.0012673541) q[0];
rz(0.92932911) q[2];
sx q[2];
rz(-2.5615938) q[2];
sx q[2];
rz(1.7763688) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.77308622) q[1];
sx q[1];
rz(-1.0463873) q[1];
sx q[1];
rz(1.4953411) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4517566) q[3];
sx q[3];
rz(-0.82131006) q[3];
sx q[3];
rz(0.46508712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.08673) q[2];
sx q[2];
rz(-2.3708998) q[2];
sx q[2];
rz(-1.9697624) q[2];
rz(0.99091667) q[3];
sx q[3];
rz(-1.0057697) q[3];
sx q[3];
rz(-1.8876373) q[3];
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
rz(-1.1427926) q[0];
sx q[0];
rz(-0.87738377) q[0];
sx q[0];
rz(-0.073609322) q[0];
rz(-2.7403846) q[1];
sx q[1];
rz(-0.3708655) q[1];
sx q[1];
rz(-2.4295095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1616514) q[0];
sx q[0];
rz(-1.5197541) q[0];
sx q[0];
rz(1.6847085) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4855494) q[2];
sx q[2];
rz(-1.5696161) q[2];
sx q[2];
rz(0.16108433) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.69171762) q[1];
sx q[1];
rz(-0.32539168) q[1];
sx q[1];
rz(2.8001333) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0832733) q[3];
sx q[3];
rz(-1.0382031) q[3];
sx q[3];
rz(-2.4709156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3965343) q[2];
sx q[2];
rz(-1.2816387) q[2];
sx q[2];
rz(0.38953951) q[2];
rz(-2.7430096) q[3];
sx q[3];
rz(-1.4094718) q[3];
sx q[3];
rz(1.4220953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1280549) q[0];
sx q[0];
rz(-3.0792397) q[0];
sx q[0];
rz(0.41559497) q[0];
rz(1.2526814) q[1];
sx q[1];
rz(-2.0038192) q[1];
sx q[1];
rz(0.91928732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556231) q[0];
sx q[0];
rz(-0.16849314) q[0];
sx q[0];
rz(1.3244434) q[0];
rz(3.058124) q[2];
sx q[2];
rz(-2.8776904) q[2];
sx q[2];
rz(0.47191274) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.34621908) q[1];
sx q[1];
rz(-2.260403) q[1];
sx q[1];
rz(2.3350689) q[1];
rz(-pi) q[2];
rz(0.76885884) q[3];
sx q[3];
rz(-1.2357841) q[3];
sx q[3];
rz(1.5227536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84819841) q[2];
sx q[2];
rz(-1.6323117) q[2];
sx q[2];
rz(0.97818127) q[2];
rz(2.6350422) q[3];
sx q[3];
rz(-1.6440441) q[3];
sx q[3];
rz(1.0165366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-0.070234805) q[0];
sx q[0];
rz(-0.6986343) q[0];
sx q[0];
rz(2.624766) q[0];
rz(-1.0454073) q[1];
sx q[1];
rz(-1.005859) q[1];
sx q[1];
rz(-0.4019092) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5346679) q[0];
sx q[0];
rz(-0.87550301) q[0];
sx q[0];
rz(2.5625227) q[0];
rz(-0.16936765) q[2];
sx q[2];
rz(-1.2886088) q[2];
sx q[2];
rz(-2.3991632) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.13788183) q[1];
sx q[1];
rz(-1.5665083) q[1];
sx q[1];
rz(1.7700618) q[1];
rz(0.95540445) q[3];
sx q[3];
rz(-1.4716574) q[3];
sx q[3];
rz(-0.73777273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4154633) q[2];
sx q[2];
rz(-1.1795421) q[2];
sx q[2];
rz(2.9187875) q[2];
rz(-2.2807617) q[3];
sx q[3];
rz(-1.6313044) q[3];
sx q[3];
rz(-1.2739325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.05805) q[0];
sx q[0];
rz(-1.8996176) q[0];
sx q[0];
rz(0.85839957) q[0];
rz(-2.3421613) q[1];
sx q[1];
rz(-2.1245427) q[1];
sx q[1];
rz(1.7564836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0388865) q[0];
sx q[0];
rz(-1.2974076) q[0];
sx q[0];
rz(0.1929184) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4102752) q[2];
sx q[2];
rz(-0.33000144) q[2];
sx q[2];
rz(0.80182072) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3333545) q[1];
sx q[1];
rz(-0.97443354) q[1];
sx q[1];
rz(-1.4247231) q[1];
x q[2];
rz(-1.9989789) q[3];
sx q[3];
rz(-2.611736) q[3];
sx q[3];
rz(0.65880218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73056209) q[2];
sx q[2];
rz(-2.4081814) q[2];
sx q[2];
rz(-2.9900271) q[2];
rz(-1.1847121) q[3];
sx q[3];
rz(-1.1368753) q[3];
sx q[3];
rz(-0.90004164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93722349) q[0];
sx q[0];
rz(-1.3441514) q[0];
sx q[0];
rz(2.3580661) q[0];
rz(-0.20768684) q[1];
sx q[1];
rz(-0.91778225) q[1];
sx q[1];
rz(-1.9327365) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3850526) q[0];
sx q[0];
rz(-1.970743) q[0];
sx q[0];
rz(2.6232096) q[0];
x q[1];
rz(2.3105588) q[2];
sx q[2];
rz(-1.990982) q[2];
sx q[2];
rz(-0.56545602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6121587) q[1];
sx q[1];
rz(-2.4116257) q[1];
sx q[1];
rz(1.6298226) q[1];
x q[2];
rz(2.9840488) q[3];
sx q[3];
rz(-1.9374401) q[3];
sx q[3];
rz(-3.1160924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2850767) q[2];
sx q[2];
rz(-0.74399844) q[2];
sx q[2];
rz(-3.0169955) q[2];
rz(-1.3397269) q[3];
sx q[3];
rz(-0.41346574) q[3];
sx q[3];
rz(1.9477897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1272005) q[0];
sx q[0];
rz(-1.7008282) q[0];
sx q[0];
rz(-0.097323962) q[0];
rz(-0.65517455) q[1];
sx q[1];
rz(-2.580878) q[1];
sx q[1];
rz(2.4526144) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6902721) q[0];
sx q[0];
rz(-0.49271944) q[0];
sx q[0];
rz(-1.7283597) q[0];
x q[1];
rz(3.1186173) q[2];
sx q[2];
rz(-1.2189294) q[2];
sx q[2];
rz(-0.018154649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1199159) q[1];
sx q[1];
rz(-2.4525453) q[1];
sx q[1];
rz(-0.20153789) q[1];
rz(-pi) q[2];
x q[2];
rz(1.407988) q[3];
sx q[3];
rz(-0.88811425) q[3];
sx q[3];
rz(2.9014661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9002306) q[2];
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
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083753839) q[0];
sx q[0];
rz(-0.76863113) q[0];
sx q[0];
rz(-0.88430697) q[0];
rz(3.0506483) q[1];
sx q[1];
rz(-2.2021286) q[1];
sx q[1];
rz(1.9627862) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0232078) q[0];
sx q[0];
rz(-1.410648) q[0];
sx q[0];
rz(1.5051756) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0021807533) q[2];
sx q[2];
rz(-2.2366887) q[2];
sx q[2];
rz(-1.8357753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7024421) q[1];
sx q[1];
rz(-2.1259754) q[1];
sx q[1];
rz(-1.3693643) q[1];
x q[2];
rz(0.53369681) q[3];
sx q[3];
rz(-2.9096536) q[3];
sx q[3];
rz(0.75702778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74943298) q[2];
sx q[2];
rz(-1.9955111) q[2];
sx q[2];
rz(-2.667099) q[2];
rz(3.0992881) q[3];
sx q[3];
rz(-1.5174805) q[3];
sx q[3];
rz(1.3206941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6831191) q[0];
sx q[0];
rz(-1.1011769) q[0];
sx q[0];
rz(1.839757) q[0];
rz(0.63109541) q[1];
sx q[1];
rz(-1.0429691) q[1];
sx q[1];
rz(2.1249117) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0250769) q[0];
sx q[0];
rz(-2.1943551) q[0];
sx q[0];
rz(0.35720346) q[0];
rz(-pi) q[1];
rz(2.0687556) q[2];
sx q[2];
rz(-3.0291621) q[2];
sx q[2];
rz(0.87804483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.423682) q[1];
sx q[1];
rz(-0.85211098) q[1];
sx q[1];
rz(0.011773034) q[1];
x q[2];
rz(0.70412962) q[3];
sx q[3];
rz(-2.1065154) q[3];
sx q[3];
rz(0.84402049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3116698) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(-2.1666918) q[2];
rz(-1.4026862) q[3];
sx q[3];
rz(-2.1699173) q[3];
sx q[3];
rz(1.5506802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1350031) q[0];
sx q[0];
rz(-2.4000744) q[0];
sx q[0];
rz(0.97493521) q[0];
rz(1.7368175) q[1];
sx q[1];
rz(-1.732323) q[1];
sx q[1];
rz(2.7224532) q[1];
rz(-0.97750935) q[2];
sx q[2];
rz(-1.1326264) q[2];
sx q[2];
rz(-2.626194) q[2];
rz(-1.4584803) q[3];
sx q[3];
rz(-1.1426123) q[3];
sx q[3];
rz(-1.860426) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
