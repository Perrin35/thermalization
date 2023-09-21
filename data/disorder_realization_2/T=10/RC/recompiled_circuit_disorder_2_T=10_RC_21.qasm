OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(5.7603523) q[0];
sx q[0];
rz(8.8011959) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(2.6410988) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5252285) q[0];
sx q[0];
rz(-0.59983569) q[0];
sx q[0];
rz(-0.0068378011) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1262769) q[2];
sx q[2];
rz(-2.9657288) q[2];
sx q[2];
rz(-1.6989087) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.11195586) q[1];
sx q[1];
rz(-1.2458548) q[1];
sx q[1];
rz(-0.8128266) q[1];
x q[2];
rz(-2.9849103) q[3];
sx q[3];
rz(-1.1374258) q[3];
sx q[3];
rz(-0.79790786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.2228489) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-0.99213123) q[3];
sx q[3];
rz(-0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.37110776) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(1.3868015) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(2.4761377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9294445) q[0];
sx q[0];
rz(-1.11709) q[0];
sx q[0];
rz(-1.2544022) q[0];
rz(-pi) q[1];
rz(2.6639054) q[2];
sx q[2];
rz(-2.7825232) q[2];
sx q[2];
rz(1.2815086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3033777) q[1];
sx q[1];
rz(-1.8068552) q[1];
sx q[1];
rz(-2.6462376) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2802248) q[3];
sx q[3];
rz(-1.5966166) q[3];
sx q[3];
rz(0.65755075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5370496) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(-2.4501734) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-0.75876045) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.741515) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83126691) q[0];
sx q[0];
rz(-1.184549) q[0];
sx q[0];
rz(1.8589742) q[0];
rz(-2.5163469) q[2];
sx q[2];
rz(-2.9125104) q[2];
sx q[2];
rz(-0.89876995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6058265) q[1];
sx q[1];
rz(-1.2488135) q[1];
sx q[1];
rz(-0.95226007) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8758043) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(0.65762323) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(-1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(0.34805527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845683) q[0];
sx q[0];
rz(-0.94385249) q[0];
sx q[0];
rz(1.9355965) q[0];
rz(-pi) q[1];
rz(-1.1966755) q[2];
sx q[2];
rz(-1.4305563) q[2];
sx q[2];
rz(1.8082878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0749803) q[1];
sx q[1];
rz(-0.96251026) q[1];
sx q[1];
rz(0.10886701) q[1];
rz(-pi) q[2];
rz(0.16163687) q[3];
sx q[3];
rz(-1.703754) q[3];
sx q[3];
rz(2.0032351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-0.42281881) q[2];
rz(2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.2984498) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86622483) q[0];
sx q[0];
rz(-2.9242762) q[0];
sx q[0];
rz(2.2989681) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8680044) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(2.1538018) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4106771) q[1];
sx q[1];
rz(-1.2123322) q[1];
sx q[1];
rz(1.4645542) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88921806) q[3];
sx q[3];
rz(-0.71435706) q[3];
sx q[3];
rz(1.0970955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9412781) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(2.6679664) q[2];
rz(3.04223) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(0.9978869) q[0];
rz(2.2672794) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-0.46674892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022910718) q[0];
sx q[0];
rz(-0.70436275) q[0];
sx q[0];
rz(2.8055311) q[0];
rz(2.7141063) q[2];
sx q[2];
rz(-1.3689405) q[2];
sx q[2];
rz(1.0496548) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14086831) q[1];
sx q[1];
rz(-1.6182401) q[1];
sx q[1];
rz(0.15111698) q[1];
rz(-pi) q[2];
rz(1.407133) q[3];
sx q[3];
rz(-1.4758849) q[3];
sx q[3];
rz(-1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85990396) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5647282) q[2];
rz(2.5148897) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(0.41982857) q[0];
rz(-0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5931169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67328582) q[0];
sx q[0];
rz(-1.562403) q[0];
sx q[0];
rz(0.087455672) q[0];
rz(2.1599342) q[2];
sx q[2];
rz(-0.33827153) q[2];
sx q[2];
rz(2.8209518) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.067804) q[1];
sx q[1];
rz(-1.5556591) q[1];
sx q[1];
rz(3.036036) q[1];
x q[2];
rz(-1.1057304) q[3];
sx q[3];
rz(-1.8814058) q[3];
sx q[3];
rz(0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(2.5166683) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(-1.8274868) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02567357) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(-1.7220108) q[0];
rz(-0.74771379) q[2];
sx q[2];
rz(-1.3249825) q[2];
sx q[2];
rz(0.73033787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18409477) q[1];
sx q[1];
rz(-3.0093319) q[1];
sx q[1];
rz(1.1485419) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2460327) q[3];
sx q[3];
rz(-0.78403463) q[3];
sx q[3];
rz(1.3336381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4079995) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(-0.96380487) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(-1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.6864809) q[0];
sx q[0];
rz(-1.2458941) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(0.54661173) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0194191) q[0];
sx q[0];
rz(-2.008736) q[0];
sx q[0];
rz(2.7287448) q[0];
rz(-pi) q[1];
rz(-2.7466752) q[2];
sx q[2];
rz(-1.2262218) q[2];
sx q[2];
rz(-0.049023703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1281801) q[1];
sx q[1];
rz(-1.1798522) q[1];
sx q[1];
rz(0.18545111) q[1];
x q[2];
rz(-2.2253753) q[3];
sx q[3];
rz(-1.1715874) q[3];
sx q[3];
rz(0.44772128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0299915) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.4477504) q[2];
rz(2.0041806) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(-0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-2.4243673) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(2.4338914) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9428064) q[0];
sx q[0];
rz(-1.9137303) q[0];
sx q[0];
rz(-0.85533157) q[0];
rz(-pi) q[1];
rz(2.131358) q[2];
sx q[2];
rz(-0.26294225) q[2];
sx q[2];
rz(1.5514785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3860491) q[1];
sx q[1];
rz(-1.2554597) q[1];
sx q[1];
rz(0.83868933) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(-2.4952863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(2.2383402) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.306504) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(-2.3564561) q[2];
sx q[2];
rz(-2.5265836) q[2];
sx q[2];
rz(2.1267736) q[2];
rz(1.5296616) q[3];
sx q[3];
rz(-1.6251246) q[3];
sx q[3];
rz(2.3390935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
