OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93513918) q[0];
sx q[0];
rz(-2.3606665) q[0];
sx q[0];
rz(0.20679064) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(0.18016711) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500096) q[0];
sx q[0];
rz(-1.9073434) q[0];
sx q[0];
rz(2.6495289) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51282672) q[2];
sx q[2];
rz(-1.4564118) q[2];
sx q[2];
rz(-1.0637384) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3611316) q[1];
sx q[1];
rz(-1.4968922) q[1];
sx q[1];
rz(-2.2547045) q[1];
x q[2];
rz(0.22523017) q[3];
sx q[3];
rz(-0.69437829) q[3];
sx q[3];
rz(0.089689342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(1.8475378) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(-2.7348203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(-2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(-1.3719826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3208647) q[0];
sx q[0];
rz(-0.67824368) q[0];
sx q[0];
rz(-3.0430549) q[0];
rz(-2.9397474) q[2];
sx q[2];
rz(-0.55515528) q[2];
sx q[2];
rz(1.7472048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0731197) q[1];
sx q[1];
rz(-1.321055) q[1];
sx q[1];
rz(-1.6176893) q[1];
rz(-2.6847141) q[3];
sx q[3];
rz(-2.6279454) q[3];
sx q[3];
rz(1.160624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(-0.99622336) q[2];
rz(1.0960724) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4662194) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(1.0590142) q[0];
rz(1.9937218) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6305884) q[0];
sx q[0];
rz(-2.2017041) q[0];
sx q[0];
rz(-1.0738119) q[0];
rz(-pi) q[1];
rz(-2.9840165) q[2];
sx q[2];
rz(-1.9489947) q[2];
sx q[2];
rz(-2.6842897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73653446) q[1];
sx q[1];
rz(-1.5824123) q[1];
sx q[1];
rz(-1.2606773) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16626658) q[3];
sx q[3];
rz(-0.85321745) q[3];
sx q[3];
rz(1.0176413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(2.1030203) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(-1.0669473) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70908961) q[0];
sx q[0];
rz(-1.6095918) q[0];
sx q[0];
rz(1.0409271) q[0];
rz(-1.6520086) q[2];
sx q[2];
rz(-1.8394543) q[2];
sx q[2];
rz(3.0693698) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.74944118) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(0.72026003) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9158343) q[3];
sx q[3];
rz(-1.0994214) q[3];
sx q[3];
rz(-2.2798722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2944494) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(0.52044049) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(-0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(0.83475137) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(-1.429819) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833417) q[0];
sx q[0];
rz(-1.5081076) q[0];
sx q[0];
rz(1.5682674) q[0];
rz(-pi) q[1];
rz(2.8565035) q[2];
sx q[2];
rz(-1.8008917) q[2];
sx q[2];
rz(0.5321815) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1289542) q[1];
sx q[1];
rz(-1.5938938) q[1];
sx q[1];
rz(-2.3501293) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3423213) q[3];
sx q[3];
rz(-2.4704774) q[3];
sx q[3];
rz(1.6959147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4804068) q[2];
sx q[2];
rz(-1.0057665) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(-1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(2.5710411) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7084224) q[0];
sx q[0];
rz(-0.47009429) q[0];
sx q[0];
rz(-2.2545459) q[0];
rz(-1.19679) q[2];
sx q[2];
rz(-1.1454957) q[2];
sx q[2];
rz(-2.5380295) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30299444) q[1];
sx q[1];
rz(-1.0698776) q[1];
sx q[1];
rz(-0.8143199) q[1];
rz(-pi) q[2];
rz(0.029383226) q[3];
sx q[3];
rz(-0.64850649) q[3];
sx q[3];
rz(0.23119584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-2.3449507) q[2];
rz(0.32026511) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657848) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(2.7835223) q[0];
rz(2.8529196) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.3105062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3863556) q[0];
sx q[0];
rz(-2.1023395) q[0];
sx q[0];
rz(1.7887572) q[0];
rz(-pi) q[1];
rz(-2.481776) q[2];
sx q[2];
rz(-2.3265127) q[2];
sx q[2];
rz(1.8191402) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.83289355) q[1];
sx q[1];
rz(-2.434224) q[1];
sx q[1];
rz(-2.4127712) q[1];
rz(-pi) q[2];
rz(-0.60704105) q[3];
sx q[3];
rz(-0.97902521) q[3];
sx q[3];
rz(0.65601635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(-1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(2.8009801) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(-1.9514726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0252359) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(0.063409253) q[0];
x q[1];
rz(3.0731593) q[2];
sx q[2];
rz(-1.4901731) q[2];
sx q[2];
rz(1.8552519) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.076188033) q[1];
sx q[1];
rz(-1.8333865) q[1];
sx q[1];
rz(0.29448387) q[1];
rz(-pi) q[2];
rz(1.0134775) q[3];
sx q[3];
rz(-1.9625469) q[3];
sx q[3];
rz(1.980892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3999346) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(-0.41440543) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-0.32489052) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-2.9898306) q[0];
rz(1.4258619) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(2.192416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28578368) q[0];
sx q[0];
rz(-1.4306418) q[0];
sx q[0];
rz(2.6967718) q[0];
rz(-pi) q[1];
rz(2.8434535) q[2];
sx q[2];
rz(-0.65152822) q[2];
sx q[2];
rz(2.8704314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28724972) q[1];
sx q[1];
rz(-2.2209475) q[1];
sx q[1];
rz(1.2681505) q[1];
x q[2];
rz(-2.1612694) q[3];
sx q[3];
rz(-1.1579517) q[3];
sx q[3];
rz(-2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-0.43241832) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(2.4798415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-0.86078763) q[1];
sx q[1];
rz(-2.4180791) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3961807) q[0];
sx q[0];
rz(-1.6999082) q[0];
sx q[0];
rz(0.73208916) q[0];
rz(-2.4621387) q[2];
sx q[2];
rz(-1.5005158) q[2];
sx q[2];
rz(-2.7714504) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3620421) q[1];
sx q[1];
rz(-2.3099646) q[1];
sx q[1];
rz(1.9401624) q[1];
rz(-1.5659329) q[3];
sx q[3];
rz(-1.8412207) q[3];
sx q[3];
rz(1.566554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(-0.71183318) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.9217459) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(0.28221054) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(-0.49610207) q[2];
sx q[2];
rz(-1.169687) q[2];
sx q[2];
rz(-0.088427831) q[2];
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
