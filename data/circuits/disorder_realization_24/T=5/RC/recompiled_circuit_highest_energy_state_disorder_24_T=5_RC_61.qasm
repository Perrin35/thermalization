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
rz(-2.6744106) q[0];
sx q[0];
rz(-1.8888357) q[0];
sx q[0];
rz(2.3350265) q[0];
rz(-0.84713495) q[1];
sx q[1];
rz(5.7972941) q[1];
sx q[1];
rz(13.426933) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8019077) q[0];
sx q[0];
rz(-2.1640477) q[0];
sx q[0];
rz(0.16925933) q[0];
rz(0.2496232) q[2];
sx q[2];
rz(-2.2614517) q[2];
sx q[2];
rz(0.63930852) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1137915) q[1];
sx q[1];
rz(-1.5209201) q[1];
sx q[1];
rz(-2.8667843) q[1];
x q[2];
rz(0.87423726) q[3];
sx q[3];
rz(-0.69948602) q[3];
sx q[3];
rz(-1.4892088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90031558) q[2];
sx q[2];
rz(-0.52854717) q[2];
sx q[2];
rz(0.47145525) q[2];
rz(-2.6491162) q[3];
sx q[3];
rz(-1.9086647) q[3];
sx q[3];
rz(-1.2537778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8660698) q[0];
sx q[0];
rz(-2.7560784) q[0];
sx q[0];
rz(-0.27819124) q[0];
rz(-1.9506075) q[1];
sx q[1];
rz(-1.3326125) q[1];
sx q[1];
rz(1.8461548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0012521) q[0];
sx q[0];
rz(-2.3352726) q[0];
sx q[0];
rz(-0.6722404) q[0];
x q[1];
rz(-0.079600178) q[2];
sx q[2];
rz(-1.6222553) q[2];
sx q[2];
rz(-2.1934794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56486121) q[1];
sx q[1];
rz(-2.7612491) q[1];
sx q[1];
rz(1.5647792) q[1];
rz(-pi) q[2];
rz(2.2601028) q[3];
sx q[3];
rz(-2.3866598) q[3];
sx q[3];
rz(0.024294446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94464716) q[2];
sx q[2];
rz(-1.4713919) q[2];
sx q[2];
rz(0.50986457) q[2];
rz(1.4465796) q[3];
sx q[3];
rz(-2.6810985) q[3];
sx q[3];
rz(0.50486008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6079717) q[0];
sx q[0];
rz(-1.5190268) q[0];
sx q[0];
rz(2.1571889) q[0];
rz(3.1254752) q[1];
sx q[1];
rz(-2.0033483) q[1];
sx q[1];
rz(1.791753) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4181217) q[0];
sx q[0];
rz(-1.4051361) q[0];
sx q[0];
rz(0.70933527) q[0];
x q[1];
rz(1.0764112) q[2];
sx q[2];
rz(-2.5177023) q[2];
sx q[2];
rz(0.15984331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3193839) q[1];
sx q[1];
rz(-2.2693386) q[1];
sx q[1];
rz(-1.6234267) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0257019) q[3];
sx q[3];
rz(-2.3724764) q[3];
sx q[3];
rz(-0.18530857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4237889) q[2];
sx q[2];
rz(-0.70222792) q[2];
sx q[2];
rz(-0.88199893) q[2];
rz(1.1841904) q[3];
sx q[3];
rz(-2.2004674) q[3];
sx q[3];
rz(-0.73630303) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940359) q[0];
sx q[0];
rz(-1.2268257) q[0];
sx q[0];
rz(-0.077022821) q[0];
rz(0.19829622) q[1];
sx q[1];
rz(-1.6402596) q[1];
sx q[1];
rz(2.95453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27320751) q[0];
sx q[0];
rz(-1.7284231) q[0];
sx q[0];
rz(0.34021838) q[0];
rz(-pi) q[1];
rz(-0.18128975) q[2];
sx q[2];
rz(-1.9167056) q[2];
sx q[2];
rz(-0.51900253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.143972) q[1];
sx q[1];
rz(-1.629843) q[1];
sx q[1];
rz(2.4322492) q[1];
rz(-pi) q[2];
rz(1.2310394) q[3];
sx q[3];
rz(-1.8700784) q[3];
sx q[3];
rz(1.9055942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5483115) q[2];
sx q[2];
rz(-2.686794) q[2];
sx q[2];
rz(1.6853257) q[2];
rz(-0.71825394) q[3];
sx q[3];
rz(-1.3879489) q[3];
sx q[3];
rz(-2.3072402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2497571) q[0];
sx q[0];
rz(-1.3351853) q[0];
sx q[0];
rz(-2.7030113) q[0];
rz(2.7843685) q[1];
sx q[1];
rz(-1.8635187) q[1];
sx q[1];
rz(-1.2507778) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8880576) q[0];
sx q[0];
rz(-1.6707289) q[0];
sx q[0];
rz(-1.1424095) q[0];
x q[1];
rz(-3.0635034) q[2];
sx q[2];
rz(-1.8812582) q[2];
sx q[2];
rz(-3.1333528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4080271) q[1];
sx q[1];
rz(-1.1981315) q[1];
sx q[1];
rz(-1.9994451) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6762371) q[3];
sx q[3];
rz(-0.28959238) q[3];
sx q[3];
rz(0.48894879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.232051) q[2];
sx q[2];
rz(-2.6614058) q[2];
sx q[2];
rz(1.5117744) q[2];
rz(0.016228598) q[3];
sx q[3];
rz(-2.0461693) q[3];
sx q[3];
rz(2.3487263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4474354) q[0];
sx q[0];
rz(-1.8614391) q[0];
sx q[0];
rz(1.1496899) q[0];
rz(-0.42090526) q[1];
sx q[1];
rz(-1.4525388) q[1];
sx q[1];
rz(3.0973869) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71647787) q[0];
sx q[0];
rz(-1.2477861) q[0];
sx q[0];
rz(-2.4130505) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8755336) q[2];
sx q[2];
rz(-2.478577) q[2];
sx q[2];
rz(-0.11817486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.068983745) q[1];
sx q[1];
rz(-1.1812967) q[1];
sx q[1];
rz(-0.087444604) q[1];
x q[2];
rz(-0.44852528) q[3];
sx q[3];
rz(-0.57975804) q[3];
sx q[3];
rz(0.019252456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1786903) q[2];
sx q[2];
rz(-2.2447605) q[2];
sx q[2];
rz(2.1964591) q[2];
rz(1.6804228) q[3];
sx q[3];
rz(-1.9943359) q[3];
sx q[3];
rz(2.6233961) q[3];
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
rz(-pi/2) q[0];
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
rz(-1.9464251) q[0];
sx q[0];
rz(-0.31038809) q[0];
sx q[0];
rz(-0.84939605) q[0];
rz(1.7646344) q[1];
sx q[1];
rz(-0.57146776) q[1];
sx q[1];
rz(-1.8570541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81818578) q[0];
sx q[0];
rz(-1.0844829) q[0];
sx q[0];
rz(-1.125) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32676593) q[2];
sx q[2];
rz(-2.8180455) q[2];
sx q[2];
rz(1.7015333) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17749912) q[1];
sx q[1];
rz(-1.0455789) q[1];
sx q[1];
rz(-1.5296658) q[1];
rz(-2.2734145) q[3];
sx q[3];
rz(-2.4560611) q[3];
sx q[3];
rz(3.1044416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66594243) q[2];
sx q[2];
rz(-1.7721756) q[2];
sx q[2];
rz(-1.5744038) q[2];
rz(0.33268467) q[3];
sx q[3];
rz(-2.0761469) q[3];
sx q[3];
rz(-2.1596597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57629267) q[0];
sx q[0];
rz(-1.1197634) q[0];
sx q[0];
rz(2.0530307) q[0];
rz(2.2125878) q[1];
sx q[1];
rz(-1.6477511) q[1];
sx q[1];
rz(-0.30379024) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32919183) q[0];
sx q[0];
rz(-2.441933) q[0];
sx q[0];
rz(-0.60075642) q[0];
rz(-1.7342042) q[2];
sx q[2];
rz(-1.4532928) q[2];
sx q[2];
rz(2.8325956) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15081295) q[1];
sx q[1];
rz(-0.44027281) q[1];
sx q[1];
rz(-1.7899465) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63173826) q[3];
sx q[3];
rz(-2.4186385) q[3];
sx q[3];
rz(0.50448862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2527689) q[2];
sx q[2];
rz(-1.6656275) q[2];
sx q[2];
rz(-0.26407537) q[2];
rz(-1.1325599) q[3];
sx q[3];
rz(-2.8045636) q[3];
sx q[3];
rz(-1.0549217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3573414) q[0];
sx q[0];
rz(-0.34264523) q[0];
sx q[0];
rz(2.1270879) q[0];
rz(2.2134589) q[1];
sx q[1];
rz(-1.4200297) q[1];
sx q[1];
rz(2.6511505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2078932) q[0];
sx q[0];
rz(-1.5053147) q[0];
sx q[0];
rz(-1.4235953) q[0];
x q[1];
rz(0.3447475) q[2];
sx q[2];
rz(-1.9354068) q[2];
sx q[2];
rz(-2.5188336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7568276) q[1];
sx q[1];
rz(-0.30337983) q[1];
sx q[1];
rz(-1.7089173) q[1];
x q[2];
rz(-1.1808628) q[3];
sx q[3];
rz(-2.56757) q[3];
sx q[3];
rz(0.61709484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8162615) q[2];
sx q[2];
rz(-1.0518495) q[2];
sx q[2];
rz(3.0391147) q[2];
rz(-1.2517733) q[3];
sx q[3];
rz(-1.3636369) q[3];
sx q[3];
rz(1.6583091) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74806279) q[0];
sx q[0];
rz(-0.04638014) q[0];
sx q[0];
rz(1.8512132) q[0];
rz(-2.6875467) q[1];
sx q[1];
rz(-1.8171277) q[1];
sx q[1];
rz(-2.9581199) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.030122) q[0];
sx q[0];
rz(-0.18387499) q[0];
sx q[0];
rz(-1.4669815) q[0];
x q[1];
rz(-1.5558335) q[2];
sx q[2];
rz(-0.97285473) q[2];
sx q[2];
rz(2.9460965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.390229) q[1];
sx q[1];
rz(-1.5469264) q[1];
sx q[1];
rz(0.49225251) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0049263) q[3];
sx q[3];
rz(-1.4938032) q[3];
sx q[3];
rz(-2.1044097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4219249) q[2];
sx q[2];
rz(-2.0823961) q[2];
sx q[2];
rz(0.025040778) q[2];
rz(0.54343623) q[3];
sx q[3];
rz(-0.70356026) q[3];
sx q[3];
rz(2.8652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.339879) q[0];
sx q[0];
rz(-2.1262953) q[0];
sx q[0];
rz(-2.3332818) q[0];
rz(-2.8080151) q[1];
sx q[1];
rz(-1.7671276) q[1];
sx q[1];
rz(-0.50552013) q[1];
rz(-0.12814604) q[2];
sx q[2];
rz(-1.8983049) q[2];
sx q[2];
rz(-2.9111918) q[2];
rz(-0.89743817) q[3];
sx q[3];
rz(-1.031395) q[3];
sx q[3];
rz(-0.54678834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
