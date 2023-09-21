OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2064535) q[0];
sx q[0];
rz(-0.78092617) q[0];
sx q[0];
rz(2.934802) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(0.18016711) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59158303) q[0];
sx q[0];
rz(-1.2342493) q[0];
sx q[0];
rz(2.6495289) q[0];
rz(1.4397058) q[2];
sx q[2];
rz(-2.0799473) q[2];
sx q[2];
rz(-2.5703562) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.78046103) q[1];
sx q[1];
rz(-1.4968922) q[1];
sx q[1];
rz(-0.88688811) q[1];
rz(-1.7546685) q[3];
sx q[3];
rz(-2.2443218) q[3];
sx q[3];
rz(-2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(1.2940548) q[2];
rz(2.7358352) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(-2.7348203) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92293537) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(1.3719826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.820728) q[0];
sx q[0];
rz(-0.67824368) q[0];
sx q[0];
rz(-0.098537785) q[0];
rz(-pi) q[1];
rz(2.5955822) q[2];
sx q[2];
rz(-1.464932) q[2];
sx q[2];
rz(-0.0042303483) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5139211) q[1];
sx q[1];
rz(-1.6162335) q[1];
sx q[1];
rz(2.8915878) q[1];
rz(-0.45687859) q[3];
sx q[3];
rz(-0.51364726) q[3];
sx q[3];
rz(-1.9809686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(-2.1453693) q[2];
rz(1.0960724) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.4662194) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-2.0825785) q[0];
rz(1.9937218) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(2.0770729) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6305884) q[0];
sx q[0];
rz(-2.2017041) q[0];
sx q[0];
rz(2.0677807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9840165) q[2];
sx q[2];
rz(-1.1925979) q[2];
sx q[2];
rz(-2.6842897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3435622) q[1];
sx q[1];
rz(-2.8312632) q[1];
sx q[1];
rz(-1.6088435) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2952609) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(-0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(-2.0227382) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082821) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(1.0669473) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2571714) q[0];
sx q[0];
rz(-2.1002249) q[0];
sx q[0];
rz(0.04495312) q[0];
rz(2.8550451) q[2];
sx q[2];
rz(-2.8612125) q[2];
sx q[2];
rz(-2.7718411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49003562) q[1];
sx q[1];
rz(-2.2228096) q[1];
sx q[1];
rz(-1.055056) q[1];
rz(2.0526485) q[3];
sx q[3];
rz(-1.7715766) q[3];
sx q[3];
rz(-2.5364385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84714326) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(-0.28953826) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510058) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(2.3068413) q[0];
rz(1.897215) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(1.7117737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5127038) q[0];
sx q[0];
rz(-1.5733203) q[0];
sx q[0];
rz(-3.0789037) q[0];
rz(-1.3313815) q[2];
sx q[2];
rz(-1.2934226) q[2];
sx q[2];
rz(-2.036236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5809708) q[1];
sx q[1];
rz(-0.79172687) q[1];
sx q[1];
rz(0.03246275) q[1];
x q[2];
rz(1.3423213) q[3];
sx q[3];
rz(-2.4704774) q[3];
sx q[3];
rz(-1.445678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-3.0094106) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64602393) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(-1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-0.57055155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8342499) q[0];
sx q[0];
rz(-1.9295921) q[0];
sx q[0];
rz(0.3105727) q[0];
rz(0.678755) q[2];
sx q[2];
rz(-0.55870134) q[2];
sx q[2];
rz(-2.9847381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4482566) q[1];
sx q[1];
rz(-2.2168471) q[1];
sx q[1];
rz(-0.64530428) q[1];
x q[2];
rz(-1.5930575) q[3];
sx q[3];
rz(-2.2189757) q[3];
sx q[3];
rz(2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4346314) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(2.8213275) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657848) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(0.35807034) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(1.8310865) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983127) q[0];
sx q[0];
rz(-2.5710921) q[0];
sx q[0];
rz(0.35240726) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69775478) q[2];
sx q[2];
rz(-2.0332094) q[2];
sx q[2];
rz(-2.4042839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83289355) q[1];
sx q[1];
rz(-0.70736865) q[1];
sx q[1];
rz(-2.4127712) q[1];
x q[2];
rz(-2.2745423) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(-1.5501319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(-0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(1.1901201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1163568) q[0];
sx q[0];
rz(-1.802117) q[0];
sx q[0];
rz(3.0781834) q[0];
x q[1];
rz(1.6516079) q[2];
sx q[2];
rz(-1.5025856) q[2];
sx q[2];
rz(0.28997544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3548673) q[1];
sx q[1];
rz(-0.39199542) q[1];
sx q[1];
rz(-0.74704945) q[1];
rz(2.1281151) q[3];
sx q[3];
rz(-1.1790457) q[3];
sx q[3];
rz(-1.1607007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3999346) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(-1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25093108) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-2.9898306) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9230726) q[0];
sx q[0];
rz(-1.1306445) q[0];
sx q[0];
rz(1.7258304) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.62990909) q[2];
sx q[2];
rz(-1.7498778) q[2];
sx q[2];
rz(1.6023028) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18898447) q[1];
sx q[1];
rz(-0.70776716) q[1];
sx q[1];
rz(2.7680552) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2373958) q[3];
sx q[3];
rz(-0.70611806) q[3];
sx q[3];
rz(1.0977942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7342547) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(2.7091743) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.0712873) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(2.7684257) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(2.4180791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7454119) q[0];
sx q[0];
rz(-1.6999082) q[0];
sx q[0];
rz(-0.73208916) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4621387) q[2];
sx q[2];
rz(-1.6410769) q[2];
sx q[2];
rz(0.37014222) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3620421) q[1];
sx q[1];
rz(-0.83162809) q[1];
sx q[1];
rz(-1.9401624) q[1];
rz(0.017541842) q[3];
sx q[3];
rz(-2.8711257) q[3];
sx q[3];
rz(1.5568352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(-2.5881361) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(-1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2198467) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(-2.8593821) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(1.1214592) q[2];
sx q[2];
rz(-2.024414) q[2];
sx q[2];
rz(-1.450962) q[2];
rz(0.40209963) q[3];
sx q[3];
rz(-1.8310908) q[3];
sx q[3];
rz(2.6487333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
