OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(3.830885) q[0];
sx q[0];
rz(9.7552714) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7666727) q[0];
sx q[0];
rz(-1.8639038) q[0];
sx q[0];
rz(-0.93107443) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.032151392) q[2];
sx q[2];
rz(-1.8841779) q[2];
sx q[2];
rz(-1.9762447) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82422772) q[1];
sx q[1];
rz(-1.0365651) q[1];
sx q[1];
rz(2.6298916) q[1];
rz(-pi) q[2];
rz(0.96592824) q[3];
sx q[3];
rz(-0.60797193) q[3];
sx q[3];
rz(2.0917497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(1.5343792) q[2];
rz(-2.2065227) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222897) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(2.2252749) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5641862) q[0];
sx q[0];
rz(-2.2803377) q[0];
sx q[0];
rz(-2.2400411) q[0];
x q[1];
rz(2.3071438) q[2];
sx q[2];
rz(-1.6685899) q[2];
sx q[2];
rz(-0.18690878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.898145) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(1.2136202) q[1];
rz(-pi) q[2];
x q[2];
rz(1.968859) q[3];
sx q[3];
rz(-1.549984) q[3];
sx q[3];
rz(1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-0.82143482) q[2];
rz(0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-0.78330529) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538552) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(0.035765212) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(2.6170513) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7099972) q[0];
sx q[0];
rz(-1.4072197) q[0];
sx q[0];
rz(-0.030199108) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32326297) q[2];
sx q[2];
rz(-0.70463902) q[2];
sx q[2];
rz(1.3779674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2136038) q[1];
sx q[1];
rz(-0.065918006) q[1];
sx q[1];
rz(-0.35857486) q[1];
rz(-pi) q[2];
rz(-0.14822652) q[3];
sx q[3];
rz(-1.1550511) q[3];
sx q[3];
rz(1.9849329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77164578) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.6463722) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(-2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-0.088949732) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-2.451992) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(-0.50868209) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3454516) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(1.4622886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5213485) q[1];
sx q[1];
rz(-0.84940956) q[1];
sx q[1];
rz(-1.1299302) q[1];
rz(0.48456405) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(-1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(-2.518667) q[2];
rz(-1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(-1.4978283) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4125815) q[0];
sx q[0];
rz(-2.3238365) q[0];
sx q[0];
rz(1.9304995) q[0];
rz(-pi) q[1];
rz(-0.096502467) q[2];
sx q[2];
rz(-0.8000024) q[2];
sx q[2];
rz(0.9347136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51753804) q[1];
sx q[1];
rz(-2.2741286) q[1];
sx q[1];
rz(-1.5348977) q[1];
x q[2];
rz(1.4086401) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(-1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(-1.5636469) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(-2.5031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(-1.7061589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84076607) q[0];
sx q[0];
rz(-2.3593785) q[0];
sx q[0];
rz(-1.7813111) q[0];
rz(-pi) q[1];
rz(1.4749182) q[2];
sx q[2];
rz(-0.55170689) q[2];
sx q[2];
rz(-0.72223896) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3458851) q[1];
sx q[1];
rz(-2.7137623) q[1];
sx q[1];
rz(-1.2659628) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1141838) q[3];
sx q[3];
rz(-1.8084744) q[3];
sx q[3];
rz(2.4159367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0662213) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(-0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.28618318) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(1.2795992) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(2.1320027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855232) q[0];
sx q[0];
rz(-1.0970322) q[0];
sx q[0];
rz(1.1958634) q[0];
rz(-0.38452734) q[2];
sx q[2];
rz(-1.782801) q[2];
sx q[2];
rz(-0.42696135) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0040782) q[1];
sx q[1];
rz(-0.61811781) q[1];
sx q[1];
rz(-0.60933463) q[1];
x q[2];
rz(0.61543492) q[3];
sx q[3];
rz(-0.86343599) q[3];
sx q[3];
rz(-2.5497041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0916831) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(-2.696864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037576588) q[0];
sx q[0];
rz(-2.4009279) q[0];
sx q[0];
rz(0.91233493) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93034805) q[2];
sx q[2];
rz(-2.1313416) q[2];
sx q[2];
rz(-0.7330187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0900314) q[1];
sx q[1];
rz(-1.6981484) q[1];
sx q[1];
rz(-2.2584372) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29856155) q[3];
sx q[3];
rz(-2.7136554) q[3];
sx q[3];
rz(-3.0446133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.3172654) q[3];
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
rz(-pi) q[3];
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
rz(1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(-2.8544193) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(2.8093991) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6824326) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(2.0072323) q[0];
rz(-pi) q[1];
x q[1];
rz(3.095093) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(-2.7547714) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0041973) q[1];
sx q[1];
rz(-1.9837556) q[1];
sx q[1];
rz(-1.2910299) q[1];
rz(-pi) q[2];
rz(-1.9751777) q[3];
sx q[3];
rz(-2.0669193) q[3];
sx q[3];
rz(0.29508428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.7920114) q[3];
sx q[3];
rz(-0.65657842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.0861417) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(-2.0822051) q[0];
rz(-2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1232676) q[0];
sx q[0];
rz(-1.8013445) q[0];
sx q[0];
rz(-0.075320764) q[0];
rz(0.48760957) q[2];
sx q[2];
rz(-1.7874103) q[2];
sx q[2];
rz(-0.85607869) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49452457) q[1];
sx q[1];
rz(-1.8956603) q[1];
sx q[1];
rz(-1.4855794) q[1];
rz(-pi) q[2];
rz(0.20930807) q[3];
sx q[3];
rz(-2.6946687) q[3];
sx q[3];
rz(-2.350654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(-1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027325252) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(-1.04503) q[2];
sx q[2];
rz(-1.3617931) q[2];
sx q[2];
rz(1.1522273) q[2];
rz(2.1847235) q[3];
sx q[3];
rz(-1.7603684) q[3];
sx q[3];
rz(-2.6852222) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];