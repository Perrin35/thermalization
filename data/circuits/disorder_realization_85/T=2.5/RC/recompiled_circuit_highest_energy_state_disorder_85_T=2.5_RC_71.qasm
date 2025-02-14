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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4199894) q[0];
sx q[0];
rz(-0.074062183) q[0];
sx q[0];
rz(-0.5038528) q[0];
x q[1];
rz(1.6915198) q[2];
sx q[2];
rz(-1.5371463) q[2];
sx q[2];
rz(-0.41481123) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3800702) q[1];
sx q[1];
rz(-2.9381611) q[1];
sx q[1];
rz(-0.037271715) q[1];
rz(1.1580519) q[3];
sx q[3];
rz(-1.6696602) q[3];
sx q[3];
rz(3.1237941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8908995) q[2];
sx q[2];
rz(-2.3919899) q[2];
sx q[2];
rz(-3.0787943) q[2];
rz(1.7167669) q[3];
sx q[3];
rz(-1.3751605) q[3];
sx q[3];
rz(2.4019901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06827908) q[0];
sx q[0];
rz(-2.8977019) q[0];
sx q[0];
rz(1.1163611) q[0];
rz(1.5694537) q[1];
sx q[1];
rz(-2.8469323) q[1];
sx q[1];
rz(-1.1192628) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7664082) q[0];
sx q[0];
rz(-1.5699727) q[0];
sx q[0];
rz(2.4341325) q[0];
rz(-pi) q[1];
rz(2.2122635) q[2];
sx q[2];
rz(-0.57999883) q[2];
sx q[2];
rz(-1.3652238) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62323233) q[1];
sx q[1];
rz(-2.6122852) q[1];
sx q[1];
rz(0.12959403) q[1];
x q[2];
rz(0.6898361) q[3];
sx q[3];
rz(-2.3202826) q[3];
sx q[3];
rz(-0.46508712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0548627) q[2];
sx q[2];
rz(-2.3708998) q[2];
sx q[2];
rz(-1.1718303) q[2];
rz(-0.99091667) q[3];
sx q[3];
rz(-2.1358229) q[3];
sx q[3];
rz(1.2539554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9988001) q[0];
sx q[0];
rz(-2.2642089) q[0];
sx q[0];
rz(0.073609322) q[0];
rz(2.7403846) q[1];
sx q[1];
rz(-2.7707272) q[1];
sx q[1];
rz(0.71208316) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7382848) q[0];
sx q[0];
rz(-1.6845595) q[0];
sx q[0];
rz(3.0902181) q[0];
rz(-pi) q[1];
rz(-2.4855494) q[2];
sx q[2];
rz(-1.5719766) q[2];
sx q[2];
rz(-0.16108433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.449875) q[1];
sx q[1];
rz(-2.816201) q[1];
sx q[1];
rz(-2.8001333) q[1];
rz(-0.058319326) q[3];
sx q[3];
rz(-2.1033896) q[3];
sx q[3];
rz(-2.4709156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74505836) q[2];
sx q[2];
rz(-1.2816387) q[2];
sx q[2];
rz(0.38953951) q[2];
rz(0.39858308) q[3];
sx q[3];
rz(-1.7321209) q[3];
sx q[3];
rz(1.7194974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013537708) q[0];
sx q[0];
rz(-0.062352926) q[0];
sx q[0];
rz(-2.7259977) q[0];
rz(-1.2526814) q[1];
sx q[1];
rz(-2.0038192) q[1];
sx q[1];
rz(-0.91928732) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556231) q[0];
sx q[0];
rz(-2.9730995) q[0];
sx q[0];
rz(-1.8171492) q[0];
rz(2.8785673) q[2];
sx q[2];
rz(-1.5925455) q[2];
sx q[2];
rz(-1.1794752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.369345) q[1];
sx q[1];
rz(-1.0076081) q[1];
sx q[1];
rz(2.289829) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46419628) q[3];
sx q[3];
rz(-2.3168543) q[3];
sx q[3];
rz(2.8620401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84819841) q[2];
sx q[2];
rz(-1.5092809) q[2];
sx q[2];
rz(2.1634114) q[2];
rz(0.50655043) q[3];
sx q[3];
rz(-1.4975486) q[3];
sx q[3];
rz(-2.1250561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070234805) q[0];
sx q[0];
rz(-2.4429584) q[0];
sx q[0];
rz(2.624766) q[0];
rz(-2.0961854) q[1];
sx q[1];
rz(-1.005859) q[1];
sx q[1];
rz(-2.7396835) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6069247) q[0];
sx q[0];
rz(-2.2660896) q[0];
sx q[0];
rz(-0.57906998) q[0];
rz(-1.8568618) q[2];
sx q[2];
rz(-1.7334043) q[2];
sx q[2];
rz(2.3608077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4116844) q[1];
sx q[1];
rz(-2.9422816) q[1];
sx q[1];
rz(1.5924551) q[1];
rz(-pi) q[2];
rz(-1.7414209) q[3];
sx q[3];
rz(-0.62230357) q[3];
sx q[3];
rz(0.69392747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7261293) q[2];
sx q[2];
rz(-1.9620506) q[2];
sx q[2];
rz(0.22280517) q[2];
rz(0.86083096) q[3];
sx q[3];
rz(-1.6313044) q[3];
sx q[3];
rz(-1.2739325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.083542682) q[0];
sx q[0];
rz(-1.8996176) q[0];
sx q[0];
rz(-0.85839957) q[0];
rz(-0.79943132) q[1];
sx q[1];
rz(-1.0170499) q[1];
sx q[1];
rz(1.7564836) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41539587) q[0];
sx q[0];
rz(-1.3851278) q[0];
sx q[0];
rz(-1.2925005) q[0];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-2.0771602) q[1];
sx q[1];
rz(-0.61187498) q[1];
sx q[1];
rz(0.21122698) q[1];
rz(-pi) q[2];
rz(-1.0812182) q[3];
sx q[3];
rz(-1.7822232) q[3];
sx q[3];
rz(-1.8544153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4110306) q[2];
sx q[2];
rz(-2.4081814) q[2];
sx q[2];
rz(-2.9900271) q[2];
rz(-1.9568806) q[3];
sx q[3];
rz(-2.0047174) q[3];
sx q[3];
rz(2.241551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93722349) q[0];
sx q[0];
rz(-1.7974412) q[0];
sx q[0];
rz(-0.7835266) q[0];
rz(2.9339058) q[1];
sx q[1];
rz(-2.2238104) q[1];
sx q[1];
rz(1.9327365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9262919) q[0];
sx q[0];
rz(-2.4982105) q[0];
sx q[0];
rz(0.70633715) q[0];
x q[1];
rz(0.98548205) q[2];
sx q[2];
rz(-0.83067465) q[2];
sx q[2];
rz(0.58500803) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6912942) q[1];
sx q[1];
rz(-2.2992059) q[1];
sx q[1];
rz(-3.0888521) q[1];
rz(2.9840488) q[3];
sx q[3];
rz(-1.9374401) q[3];
sx q[3];
rz(0.025500209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8565159) q[2];
sx q[2];
rz(-0.74399844) q[2];
sx q[2];
rz(-0.12459717) q[2];
rz(-1.3397269) q[3];
sx q[3];
rz(-2.7281269) q[3];
sx q[3];
rz(1.193803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1272005) q[0];
sx q[0];
rz(-1.7008282) q[0];
sx q[0];
rz(0.097323962) q[0];
rz(-2.4864181) q[1];
sx q[1];
rz(-0.56071463) q[1];
sx q[1];
rz(-0.68897828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98039799) q[0];
sx q[0];
rz(-1.4965048) q[0];
sx q[0];
rz(2.0583389) q[0];
rz(-0.022975401) q[2];
sx q[2];
rz(-1.2189294) q[2];
sx q[2];
rz(-0.018154649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6072907) q[1];
sx q[1];
rz(-1.698415) q[1];
sx q[1];
rz(0.67903305) q[1];
x q[2];
rz(-2.9448255) q[3];
sx q[3];
rz(-2.4428058) q[3];
sx q[3];
rz(-0.014589498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9002306) q[2];
sx q[2];
rz(-2.1817544) q[2];
sx q[2];
rz(-2.0755365) q[2];
rz(-0.42029941) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(-1.458781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578388) q[0];
sx q[0];
rz(-2.3729615) q[0];
sx q[0];
rz(-0.88430697) q[0];
rz(-0.090944313) q[1];
sx q[1];
rz(-0.93946409) q[1];
sx q[1];
rz(1.1788064) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6323157) q[0];
sx q[0];
rz(-2.9686283) q[0];
sx q[0];
rz(-2.7559682) q[0];
x q[1];
rz(-2.2366899) q[2];
sx q[2];
rz(-1.5690815) q[2];
sx q[2];
rz(2.8752665) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43915054) q[1];
sx q[1];
rz(-1.0156173) q[1];
sx q[1];
rz(-1.3693643) q[1];
rz(-pi) q[2];
rz(-1.6903773) q[3];
sx q[3];
rz(-1.371583) q[3];
sx q[3];
rz(-1.3026893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3921597) q[2];
sx q[2];
rz(-1.9955111) q[2];
sx q[2];
rz(0.47449365) q[2];
rz(0.042304603) q[3];
sx q[3];
rz(-1.5174805) q[3];
sx q[3];
rz(-1.3206941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6831191) q[0];
sx q[0];
rz(-2.0404158) q[0];
sx q[0];
rz(-1.3018357) q[0];
rz(0.63109541) q[1];
sx q[1];
rz(-1.0429691) q[1];
sx q[1];
rz(-1.016681) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0250769) q[0];
sx q[0];
rz(-2.1943551) q[0];
sx q[0];
rz(-2.7843892) q[0];
rz(-pi) q[1];
rz(-1.0728371) q[2];
sx q[2];
rz(-0.11243056) q[2];
sx q[2];
rz(2.2635478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4415623) q[1];
sx q[1];
rz(-2.4228281) q[1];
sx q[1];
rz(1.5573386) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90908124) q[3];
sx q[3];
rz(-2.1611745) q[3];
sx q[3];
rz(2.8239241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3116698) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(2.1666918) q[2];
rz(-1.4026862) q[3];
sx q[3];
rz(-0.97167531) q[3];
sx q[3];
rz(-1.5506802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0065895157) q[0];
sx q[0];
rz(-2.4000744) q[0];
sx q[0];
rz(0.97493521) q[0];
rz(-1.4047752) q[1];
sx q[1];
rz(-1.732323) q[1];
sx q[1];
rz(2.7224532) q[1];
rz(-2.1640833) q[2];
sx q[2];
rz(-2.0089663) q[2];
sx q[2];
rz(0.51539863) q[2];
rz(-0.24079612) q[3];
sx q[3];
rz(-0.44178648) q[3];
sx q[3];
rz(1.0159258) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
