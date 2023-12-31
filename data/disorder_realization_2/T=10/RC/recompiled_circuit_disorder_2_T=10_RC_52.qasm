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
rz(-0.52283302) q[0];
sx q[0];
rz(-0.62358207) q[0];
rz(3.4317598) q[1];
sx q[1];
rz(5.5640339) q[1];
sx q[1];
rz(13.066864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60808027) q[0];
sx q[0];
rz(-0.97097662) q[0];
sx q[0];
rz(-1.5754726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093437336) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(1.1364394) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0296368) q[1];
sx q[1];
rz(-1.8957378) q[1];
sx q[1];
rz(0.8128266) q[1];
rz(-2.0088828) q[3];
sx q[3];
rz(-1.7128908) q[3];
sx q[3];
rz(-0.70664584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(-1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-0.97035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(-1.3868015) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5011713) q[0];
sx q[0];
rz(-1.8542395) q[0];
sx q[0];
rz(-2.6675176) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6639054) q[2];
sx q[2];
rz(-2.7825232) q[2];
sx q[2];
rz(1.2815086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0004955) q[1];
sx q[1];
rz(-0.54447237) q[1];
sx q[1];
rz(0.46846868) q[1];
rz(-1.5311702) q[3];
sx q[3];
rz(-2.431776) q[3];
sx q[3];
rz(2.1982847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.60454303) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(-1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(-0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(-1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.4000777) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62832075) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(-2.7404286) q[0];
x q[1];
rz(-0.18685762) q[2];
sx q[2];
rz(-1.7041022) q[2];
sx q[2];
rz(-1.2847628) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5249467) q[1];
sx q[1];
rz(-0.687462) q[1];
sx q[1];
rz(2.0928659) q[1];
rz(-1.2657884) q[3];
sx q[3];
rz(-1.7294356) q[3];
sx q[3];
rz(-0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(-0.65762323) q[2];
rz(-1.1714606) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79384971) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(1.6943278) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(0.34805527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570243) q[0];
sx q[0];
rz(-2.1977402) q[0];
sx q[0];
rz(-1.9355965) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1966755) q[2];
sx q[2];
rz(-1.7110363) q[2];
sx q[2];
rz(1.8082878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.066612331) q[1];
sx q[1];
rz(-0.96251026) q[1];
sx q[1];
rz(-0.10886701) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69339852) q[3];
sx q[3];
rz(-0.20892538) q[3];
sx q[3];
rz(2.8914176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7248914) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(2.7187738) q[2];
rz(-0.73741284) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(0.53043956) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(-1.2984498) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2753678) q[0];
sx q[0];
rz(-2.9242762) q[0];
sx q[0];
rz(0.84262459) q[0];
x q[1];
rz(0.60646306) q[2];
sx q[2];
rz(-1.3242553) q[2];
sx q[2];
rz(-2.390887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0188705) q[1];
sx q[1];
rz(-1.4713305) q[1];
sx q[1];
rz(-2.781267) q[1];
rz(-pi) q[2];
rz(-0.97814822) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(-0.076171906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(-3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6376003) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(0.9978869) q[0];
rz(0.87431327) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(-0.46674892) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8080374) q[0];
sx q[0];
rz(-1.7859965) q[0];
sx q[0];
rz(-2.4654885) q[0];
rz(-pi) q[1];
rz(1.7919962) q[2];
sx q[2];
rz(-1.1525407) q[2];
sx q[2];
rz(2.7115371) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0135865) q[1];
sx q[1];
rz(-0.15833536) q[1];
sx q[1];
rz(-0.30551417) q[1];
rz(-pi) q[2];
rz(-0.096188992) q[3];
sx q[3];
rz(-1.7337165) q[3];
sx q[3];
rz(0.067226203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85990396) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(1.5768645) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8108114) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(0.41982857) q[0];
rz(-2.9201674) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5484757) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2433462) q[0];
sx q[0];
rz(-1.6582489) q[0];
sx q[0];
rz(-1.5792219) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19303796) q[2];
sx q[2];
rz(-1.2912573) q[2];
sx q[2];
rz(-2.2045731) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35510264) q[1];
sx q[1];
rz(-3.0349602) q[1];
sx q[1];
rz(-0.14270466) q[1];
x q[2];
rz(2.1920131) q[3];
sx q[3];
rz(-0.55286828) q[3];
sx q[3];
rz(2.051193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(-1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.7109011) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-2.8523493) q[0];
rz(0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.8274868) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02567357) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(1.4195819) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35347519) q[2];
sx q[2];
rz(-2.3620053) q[2];
sx q[2];
rz(-2.0445063) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1738759) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(1.4500344) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8955599) q[3];
sx q[3];
rz(-0.78403463) q[3];
sx q[3];
rz(-1.3336381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(-1.2458941) q[0];
rz(0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(2.5949809) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36771691) q[0];
sx q[0];
rz(-1.1989294) q[0];
sx q[0];
rz(2.0433776) q[0];
rz(-pi) q[1];
x q[1];
rz(1.941628) q[2];
sx q[2];
rz(-1.2002581) q[2];
sx q[2];
rz(-1.3818936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6555772) q[1];
sx q[1];
rz(-1.3994819) q[1];
sx q[1];
rz(1.9678712) q[1];
x q[2];
rz(-2.176748) q[3];
sx q[3];
rz(-0.75111872) q[3];
sx q[3];
rz(-1.5918819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.6938422) q[2];
rz(2.0041806) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(2.4338914) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087698547) q[0];
sx q[0];
rz(-2.2367034) q[0];
sx q[0];
rz(0.44184394) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3466481) q[2];
sx q[2];
rz(-1.4321616) q[2];
sx q[2];
rz(0.56425205) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7555435) q[1];
sx q[1];
rz(-1.886133) q[1];
sx q[1];
rz(2.3029033) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0503065) q[3];
sx q[3];
rz(-2.4647053) q[3];
sx q[3];
rz(-2.6022079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(2.2383402) q[2];
rz(1.603027) q[3];
sx q[3];
rz(-0.86849803) q[3];
sx q[3];
rz(0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.306504) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-0.8159591) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(2.03394) q[2];
sx q[2];
rz(-1.1504428) q[2];
sx q[2];
rz(3.0124315) q[2];
rz(-1.5296616) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
