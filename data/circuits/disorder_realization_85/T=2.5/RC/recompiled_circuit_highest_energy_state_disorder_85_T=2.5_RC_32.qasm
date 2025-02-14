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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9250037) q[0];
sx q[0];
rz(-1.6356409) q[0];
sx q[0];
rz(-1.606604) q[0];
rz(-1.8433601) q[2];
sx q[2];
rz(-0.12530357) q[2];
sx q[2];
rz(1.7150777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3688226) q[1];
sx q[1];
rz(-1.5632679) q[1];
sx q[1];
rz(-2.9382985) q[1];
x q[2];
rz(1.3283861) q[3];
sx q[3];
rz(-0.42375733) q[3];
sx q[3];
rz(-1.7746944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8908995) q[2];
sx q[2];
rz(-0.74960274) q[2];
sx q[2];
rz(-3.0787943) q[2];
rz(-1.7167669) q[3];
sx q[3];
rz(-1.3751605) q[3];
sx q[3];
rz(0.73960251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0733136) q[0];
sx q[0];
rz(-0.24389076) q[0];
sx q[0];
rz(2.0252315) q[0];
rz(1.5721389) q[1];
sx q[1];
rz(-0.29466033) q[1];
sx q[1];
rz(-1.1192628) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.946944) q[0];
sx q[0];
rz(-2.4341321) q[0];
sx q[0];
rz(-0.0012673541) q[0];
rz(0.92932911) q[2];
sx q[2];
rz(-2.5615938) q[2];
sx q[2];
rz(-1.3652238) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.62323233) q[1];
sx q[1];
rz(-0.52930743) q[1];
sx q[1];
rz(-0.12959403) q[1];
x q[2];
rz(0.69208118) q[3];
sx q[3];
rz(-1.0861703) q[3];
sx q[3];
rz(-1.6177819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.08673) q[2];
sx q[2];
rz(-2.3708998) q[2];
sx q[2];
rz(1.1718303) q[2];
rz(-2.150676) q[3];
sx q[3];
rz(-2.1358229) q[3];
sx q[3];
rz(1.8876373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9988001) q[0];
sx q[0];
rz(-2.2642089) q[0];
sx q[0];
rz(0.073609322) q[0];
rz(-2.7403846) q[1];
sx q[1];
rz(-0.3708655) q[1];
sx q[1];
rz(0.71208316) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7382848) q[0];
sx q[0];
rz(-1.6845595) q[0];
sx q[0];
rz(3.0902181) q[0];
x q[1];
rz(-2.4855494) q[2];
sx q[2];
rz(-1.5696161) q[2];
sx q[2];
rz(-2.9805083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.449875) q[1];
sx q[1];
rz(-2.816201) q[1];
sx q[1];
rz(-2.8001333) q[1];
x q[2];
rz(2.1041342) q[3];
sx q[3];
rz(-1.6210307) q[3];
sx q[3];
rz(-0.929757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74505836) q[2];
sx q[2];
rz(-1.859954) q[2];
sx q[2];
rz(-0.38953951) q[2];
rz(0.39858308) q[3];
sx q[3];
rz(-1.7321209) q[3];
sx q[3];
rz(1.7194974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1280549) q[0];
sx q[0];
rz(-3.0792397) q[0];
sx q[0];
rz(-2.7259977) q[0];
rz(1.2526814) q[1];
sx q[1];
rz(-1.1377734) q[1];
sx q[1];
rz(-0.91928732) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.141826) q[0];
sx q[0];
rz(-1.5298889) q[0];
sx q[0];
rz(-1.4072988) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8785673) q[2];
sx q[2];
rz(-1.5490471) q[2];
sx q[2];
rz(1.1794752) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7722476) q[1];
sx q[1];
rz(-1.0076081) q[1];
sx q[1];
rz(2.289829) q[1];
rz(2.3727338) q[3];
sx q[3];
rz(-1.9058085) q[3];
sx q[3];
rz(1.5227536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2933942) q[2];
sx q[2];
rz(-1.5092809) q[2];
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
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070234805) q[0];
sx q[0];
rz(-0.6986343) q[0];
sx q[0];
rz(2.624766) q[0];
rz(-2.0961854) q[1];
sx q[1];
rz(-2.1357336) q[1];
sx q[1];
rz(2.7396835) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088083) q[0];
sx q[0];
rz(-2.0044745) q[0];
sx q[0];
rz(-2.3545803) q[0];
rz(-pi) q[1];
x q[1];
rz(1.044173) q[2];
sx q[2];
rz(-0.32794558) q[2];
sx q[2];
rz(-1.2931658) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7078122) q[1];
sx q[1];
rz(-1.3715327) q[1];
sx q[1];
rz(-0.0043745478) q[1];
x q[2];
rz(-1.4001717) q[3];
sx q[3];
rz(-2.5192891) q[3];
sx q[3];
rz(0.69392747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4154633) q[2];
sx q[2];
rz(-1.1795421) q[2];
sx q[2];
rz(-0.22280517) q[2];
rz(0.86083096) q[3];
sx q[3];
rz(-1.5102883) q[3];
sx q[3];
rz(1.2739325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.05805) q[0];
sx q[0];
rz(-1.8996176) q[0];
sx q[0];
rz(-0.85839957) q[0];
rz(-2.3421613) q[1];
sx q[1];
rz(-1.0170499) q[1];
sx q[1];
rz(1.3851091) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261968) q[0];
sx q[0];
rz(-1.3851278) q[0];
sx q[0];
rz(-1.8490922) q[0];
x q[1];
rz(1.8968514) q[2];
sx q[2];
rz(-1.6226124) q[2];
sx q[2];
rz(2.5246132) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0644324) q[1];
sx q[1];
rz(-0.61187498) q[1];
sx q[1];
rz(0.21122698) q[1];
x q[2];
rz(-2.0603745) q[3];
sx q[3];
rz(-1.7822232) q[3];
sx q[3];
rz(-1.2871773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73056209) q[2];
sx q[2];
rz(-2.4081814) q[2];
sx q[2];
rz(-2.9900271) q[2];
rz(-1.1847121) q[3];
sx q[3];
rz(-1.1368753) q[3];
sx q[3];
rz(2.241551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2043692) q[0];
sx q[0];
rz(-1.7974412) q[0];
sx q[0];
rz(-2.3580661) q[0];
rz(-2.9339058) q[1];
sx q[1];
rz(-2.2238104) q[1];
sx q[1];
rz(-1.9327365) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75654001) q[0];
sx q[0];
rz(-1.970743) q[0];
sx q[0];
rz(0.51838309) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83103387) q[2];
sx q[2];
rz(-1.990982) q[2];
sx q[2];
rz(-0.56545602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6912942) q[1];
sx q[1];
rz(-0.84238673) q[1];
sx q[1];
rz(0.052740514) q[1];
rz(-pi) q[2];
rz(-1.9416299) q[3];
sx q[3];
rz(-1.4238024) q[3];
sx q[3];
rz(1.539409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8565159) q[2];
sx q[2];
rz(-2.3975942) q[2];
sx q[2];
rz(0.12459717) q[2];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0143921) q[0];
sx q[0];
rz(-1.4407644) q[0];
sx q[0];
rz(-3.0442687) q[0];
rz(0.65517455) q[1];
sx q[1];
rz(-0.56071463) q[1];
sx q[1];
rz(2.4526144) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98039799) q[0];
sx q[0];
rz(-1.4965048) q[0];
sx q[0];
rz(-1.0832537) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1186173) q[2];
sx q[2];
rz(-1.9226632) q[2];
sx q[2];
rz(0.018154649) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.6072907) q[1];
sx q[1];
rz(-1.698415) q[1];
sx q[1];
rz(0.67903305) q[1];
rz(-0.19676713) q[3];
sx q[3];
rz(-0.69878687) q[3];
sx q[3];
rz(-0.014589498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9002306) q[2];
sx q[2];
rz(-2.1817544) q[2];
sx q[2];
rz(1.0660561) q[2];
rz(2.7212932) q[3];
sx q[3];
rz(-1.2371141) q[3];
sx q[3];
rz(1.6828116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083753839) q[0];
sx q[0];
rz(-2.3729615) q[0];
sx q[0];
rz(0.88430697) q[0];
rz(3.0506483) q[1];
sx q[1];
rz(-2.2021286) q[1];
sx q[1];
rz(1.9627862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0232078) q[0];
sx q[0];
rz(-1.7309446) q[0];
sx q[0];
rz(-1.636417) q[0];
rz(-pi) q[1];
rz(-1.5680204) q[2];
sx q[2];
rz(-0.66589543) q[2];
sx q[2];
rz(-1.3022873) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2388666) q[1];
sx q[1];
rz(-1.3999434) q[1];
sx q[1];
rz(-2.577223) q[1];
rz(-pi) q[2];
rz(2.6078958) q[3];
sx q[3];
rz(-0.23193905) q[3];
sx q[3];
rz(0.75702778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3921597) q[2];
sx q[2];
rz(-1.9955111) q[2];
sx q[2];
rz(0.47449365) q[2];
rz(-0.042304603) q[3];
sx q[3];
rz(-1.6241122) q[3];
sx q[3];
rz(-1.3206941) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45847356) q[0];
sx q[0];
rz(-1.1011769) q[0];
sx q[0];
rz(-1.3018357) q[0];
rz(-0.63109541) q[1];
sx q[1];
rz(-1.0429691) q[1];
sx q[1];
rz(1.016681) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0250769) q[0];
sx q[0];
rz(-2.1943551) q[0];
sx q[0];
rz(-0.35720346) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0687556) q[2];
sx q[2];
rz(-0.11243056) q[2];
sx q[2];
rz(0.87804483) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8606372) q[1];
sx q[1];
rz(-1.5619352) q[1];
sx q[1];
rz(0.85207664) q[1];
x q[2];
rz(-2.437463) q[3];
sx q[3];
rz(-1.0350772) q[3];
sx q[3];
rz(2.2975722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8299228) q[2];
sx q[2];
rz(-2.416478) q[2];
sx q[2];
rz(-0.97490087) q[2];
rz(1.4026862) q[3];
sx q[3];
rz(-2.1699173) q[3];
sx q[3];
rz(1.5909125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.97750935) q[2];
sx q[2];
rz(-2.0089663) q[2];
sx q[2];
rz(0.51539863) q[2];
rz(0.24079612) q[3];
sx q[3];
rz(-2.6998062) q[3];
sx q[3];
rz(-2.1256668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
