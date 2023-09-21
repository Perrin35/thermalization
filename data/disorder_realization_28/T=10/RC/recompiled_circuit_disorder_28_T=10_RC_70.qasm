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
rz(0.38987723) q[1];
sx q[1];
rz(-1.0607399) q[1];
sx q[1];
rz(-0.18016711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500096) q[0];
sx q[0];
rz(-1.2342493) q[0];
sx q[0];
rz(0.49206375) q[0];
rz(0.23001036) q[2];
sx q[2];
rz(-0.52431528) q[2];
sx q[2];
rz(2.8345248) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3611316) q[1];
sx q[1];
rz(-1.6447004) q[1];
sx q[1];
rz(0.88688811) q[1];
rz(-0.22523017) q[3];
sx q[3];
rz(-2.4472144) q[3];
sx q[3];
rz(-3.0519033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41123286) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(1.2940548) q[2];
rz(0.40575746) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(0.40677235) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92293537) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(0.12582114) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(1.3719826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94700891) q[0];
sx q[0];
rz(-2.2451375) q[0];
sx q[0];
rz(-1.4916923) q[0];
rz(0.20184529) q[2];
sx q[2];
rz(-0.55515528) q[2];
sx q[2];
rz(1.7472048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6276715) q[1];
sx q[1];
rz(-1.6162335) q[1];
sx q[1];
rz(0.25000484) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46866977) q[3];
sx q[3];
rz(-1.3523003) q[3];
sx q[3];
rz(0.0056497638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8877318) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(-2.1453693) q[2];
rz(2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4662194) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(-1.0590142) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-1.0645197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7721467) q[0];
sx q[0];
rz(-1.9659974) q[0];
sx q[0];
rz(0.69338436) q[0];
x q[1];
rz(2.9840165) q[2];
sx q[2];
rz(-1.9489947) q[2];
sx q[2];
rz(-0.45730293) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.83798422) q[1];
sx q[1];
rz(-1.260699) q[1];
sx q[1];
rz(-3.1293948) q[1];
rz(-pi) q[2];
rz(1.3833984) q[3];
sx q[3];
rz(-0.73323941) q[3];
sx q[3];
rz(0.76776615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(0.81937218) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(2.1030203) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(-1.0669473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70908961) q[0];
sx q[0];
rz(-1.6095918) q[0];
sx q[0];
rz(2.1006656) q[0];
x q[1];
rz(-1.6520086) q[2];
sx q[2];
rz(-1.8394543) q[2];
sx q[2];
rz(3.0693698) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2414788) q[1];
sx q[1];
rz(-0.80726868) q[1];
sx q[1];
rz(2.5680053) q[1];
x q[2];
rz(-1.9846109) q[3];
sx q[3];
rz(-2.6226351) q[3];
sx q[3];
rz(1.3299691) q[3];
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
rz(-1.8013022) q[3];
sx q[3];
rz(0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510058) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(-2.3068413) q[0];
rz(1.897215) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.429819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5127038) q[0];
sx q[0];
rz(-1.5682724) q[0];
sx q[0];
rz(-0.062688962) q[0];
rz(1.8102112) q[2];
sx q[2];
rz(-1.2934226) q[2];
sx q[2];
rz(-2.036236) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5347848) q[1];
sx q[1];
rz(-0.77960289) q[1];
sx q[1];
rz(1.603655) q[1];
rz(-pi) q[2];
rz(0.17794869) q[3];
sx q[3];
rz(-0.92015172) q[3];
sx q[3];
rz(-1.4072756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(1.2832114) q[2];
rz(3.0094106) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(-1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(2.5710411) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7084224) q[0];
sx q[0];
rz(-2.6714984) q[0];
sx q[0];
rz(2.2545459) q[0];
x q[1];
rz(-1.19679) q[2];
sx q[2];
rz(-1.9960969) q[2];
sx q[2];
rz(-0.60356319) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3444896) q[1];
sx q[1];
rz(-0.87901607) q[1];
sx q[1];
rz(2.2441191) q[1];
x q[2];
rz(1.5485351) q[3];
sx q[3];
rz(-2.2189757) q[3];
sx q[3];
rz(2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(-2.3449507) q[2];
rz(-0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(1.8626574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(-2.7835223) q[0];
rz(-0.2886731) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(-1.8310865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0689439) q[0];
sx q[0];
rz(-1.7582969) q[0];
sx q[0];
rz(2.5995273) q[0];
x q[1];
rz(2.481776) q[2];
sx q[2];
rz(-0.81507999) q[2];
sx q[2];
rz(1.8191402) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83289355) q[1];
sx q[1];
rz(-2.434224) q[1];
sx q[1];
rz(0.72882148) q[1];
x q[2];
rz(2.2745423) q[3];
sx q[3];
rz(-2.3208445) q[3];
sx q[3];
rz(-1.5501319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(-0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(-1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26577935) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.9514726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2954138) q[0];
sx q[0];
rz(-0.23970397) q[0];
sx q[0];
rz(1.3079877) q[0];
rz(-1.6516079) q[2];
sx q[2];
rz(-1.6390071) q[2];
sx q[2];
rz(-2.8516172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0654046) q[1];
sx q[1];
rz(-1.8333865) q[1];
sx q[1];
rz(0.29448387) q[1];
rz(-0.45300071) q[3];
sx q[3];
rz(-1.060033) q[3];
sx q[3];
rz(0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74165806) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(2.7271872) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(-1.4258619) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(0.94917667) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.21852) q[0];
sx q[0];
rz(-2.0109482) q[0];
sx q[0];
rz(-1.7258304) q[0];
rz(-pi) q[1];
rz(0.29813913) q[2];
sx q[2];
rz(-2.4900644) q[2];
sx q[2];
rz(2.8704314) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4703428) q[1];
sx q[1];
rz(-1.3312695) q[1];
sx q[1];
rz(-2.4688979) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90419681) q[3];
sx q[3];
rz(-2.4354746) q[3];
sx q[3];
rz(2.0437984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-0.43241832) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-2.4798415) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(0.7235136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0821738) q[0];
sx q[0];
rz(-2.2954303) q[0];
sx q[0];
rz(1.7436149) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0300272) q[2];
sx q[2];
rz(-2.4590883) q[2];
sx q[2];
rz(-2.0276558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8404624) q[1];
sx q[1];
rz(-2.3311619) q[1];
sx q[1];
rz(2.7644972) q[1];
rz(-pi) q[2];
rz(1.5659329) q[3];
sx q[3];
rz(-1.300372) q[3];
sx q[3];
rz(-1.5750386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.92131203) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(-0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217459) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(-2.8593821) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(-1.1214592) q[2];
sx q[2];
rz(-1.1171787) q[2];
sx q[2];
rz(1.6906307) q[2];
rz(-1.8525193) q[3];
sx q[3];
rz(-1.9586133) q[3];
sx q[3];
rz(-1.9546399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
