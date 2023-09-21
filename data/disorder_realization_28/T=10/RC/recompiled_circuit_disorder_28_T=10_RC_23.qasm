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
rz(3.9225188) q[0];
sx q[0];
rz(9.6315686) q[0];
rz(0.38987723) q[1];
sx q[1];
rz(-1.0607399) q[1];
sx q[1];
rz(-0.18016711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5500096) q[0];
sx q[0];
rz(-1.9073434) q[0];
sx q[0];
rz(0.49206375) q[0];
rz(-pi) q[1];
rz(-2.6287659) q[2];
sx q[2];
rz(-1.6851808) q[2];
sx q[2];
rz(2.0778542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85045056) q[1];
sx q[1];
rz(-0.88911118) q[1];
sx q[1];
rz(0.095231685) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68182919) q[3];
sx q[3];
rz(-1.7141984) q[3];
sx q[3];
rz(-1.6553866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(-1.2940548) q[2];
rz(0.40575746) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92293537) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(0.12582114) q[0];
rz(0.80548349) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(1.7696101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94700891) q[0];
sx q[0];
rz(-2.2451375) q[0];
sx q[0];
rz(1.4916923) q[0];
x q[1];
rz(0.20184529) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(1.3943878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6276715) q[1];
sx q[1];
rz(-1.6162335) q[1];
sx q[1];
rz(2.8915878) q[1];
x q[2];
rz(-1.814718) q[3];
sx q[3];
rz(-2.0274649) q[3];
sx q[3];
rz(-1.4671385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(-2.1453693) q[2];
rz(1.0960724) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(1.0590142) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(-2.0770729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36944593) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-0.69338436) q[0];
rz(1.9532922) q[2];
sx q[2];
rz(-1.4244392) q[2];
sx q[2];
rz(-1.0548897) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7980305) q[1];
sx q[1];
rz(-2.8312632) q[1];
sx q[1];
rz(-1.6088435) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2952609) q[3];
sx q[3];
rz(-1.4457821) q[3];
sx q[3];
rz(-0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(2.8783669) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(2.3222205) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.3942962) q[0];
rz(2.1030203) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(2.0746453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88442125) q[0];
sx q[0];
rz(-1.0413678) q[0];
sx q[0];
rz(-3.0966395) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8720886) q[2];
sx q[2];
rz(-1.4925033) q[2];
sx q[2];
rz(1.6646202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74944118) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(-2.4213326) q[1];
x q[2];
rz(0.22575836) q[3];
sx q[3];
rz(-1.0994214) q[3];
sx q[3];
rz(2.2798722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84714326) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(-2.8520544) q[2];
rz(-0.52044049) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(-2.3068413) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(1.7117737) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833417) q[0];
sx q[0];
rz(-1.6334851) q[0];
sx q[0];
rz(-1.5682674) q[0];
rz(-1.3313815) q[2];
sx q[2];
rz(-1.2934226) q[2];
sx q[2];
rz(-2.036236) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6068078) q[1];
sx q[1];
rz(-0.77960289) q[1];
sx q[1];
rz(1.5379377) q[1];
x q[2];
rz(0.91246446) q[3];
sx q[3];
rz(-1.7121127) q[3];
sx q[3];
rz(-0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(0.13218203) q[3];
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
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(-1.5006789) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(0.57055155) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30734277) q[0];
sx q[0];
rz(-1.9295921) q[0];
sx q[0];
rz(0.3105727) q[0];
rz(-pi) q[1];
rz(-0.45285593) q[2];
sx q[2];
rz(-1.9100683) q[2];
sx q[2];
rz(1.1277744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3444896) q[1];
sx q[1];
rz(-0.87901607) q[1];
sx q[1];
rz(2.2441191) q[1];
rz(-pi) q[2];
rz(-1.5930575) q[3];
sx q[3];
rz(-0.92261693) q[3];
sx q[3];
rz(0.19433403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4346314) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(-2.3449507) q[2];
rz(0.32026511) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(1.8626574) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(-2.7835223) q[0];
rz(2.8529196) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(1.8310865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072648777) q[0];
sx q[0];
rz(-1.3832958) q[0];
sx q[0];
rz(-0.54206538) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99408044) q[2];
sx q[2];
rz(-2.1834282) q[2];
sx q[2];
rz(2.6660369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.14168921) q[1];
sx q[1];
rz(-2.018376) q[1];
sx q[1];
rz(-0.56772851) q[1];
rz(-2.2745423) q[3];
sx q[3];
rz(-2.3208445) q[3];
sx q[3];
rz(-1.5914608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(-0.99747783) q[2];
rz(-2.5618662) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(-1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(-2.8009801) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(1.1901201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2954138) q[0];
sx q[0];
rz(-2.9018887) q[0];
sx q[0];
rz(1.3079877) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0731593) q[2];
sx q[2];
rz(-1.4901731) q[2];
sx q[2];
rz(1.8552519) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.076188033) q[1];
sx q[1];
rz(-1.3082062) q[1];
sx q[1];
rz(-0.29448387) q[1];
rz(-0.90772273) q[3];
sx q[3];
rz(-0.66909664) q[3];
sx q[3];
rz(0.95975403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74165806) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(-0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(-2.192416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.855809) q[0];
sx q[0];
rz(-1.7109509) q[0];
sx q[0];
rz(-2.6967718) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8434535) q[2];
sx q[2];
rz(-0.65152822) q[2];
sx q[2];
rz(0.27116129) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9526082) q[1];
sx q[1];
rz(-2.4338255) q[1];
sx q[1];
rz(2.7680552) q[1];
rz(-pi) q[2];
rz(2.1612694) q[3];
sx q[3];
rz(-1.983641) q[3];
sx q[3];
rz(1.0124029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7342547) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(-1.6010823) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-2.4798415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(2.7684257) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(2.4180791) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3961807) q[0];
sx q[0];
rz(-1.6999082) q[0];
sx q[0];
rz(2.4095035) q[0];
x q[1];
rz(0.67945393) q[2];
sx q[2];
rz(-1.5005158) q[2];
sx q[2];
rz(-2.7714504) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3011303) q[1];
sx q[1];
rz(-0.81043078) q[1];
sx q[1];
rz(0.37709548) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1240508) q[3];
sx q[3];
rz(-0.27046698) q[3];
sx q[3];
rz(1.5568352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92131203) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(0.55345654) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.2198467) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(-0.28221054) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(2.0201335) q[2];
sx q[2];
rz(-1.1171787) q[2];
sx q[2];
rz(1.6906307) q[2];
rz(-1.2890733) q[3];
sx q[3];
rz(-1.1829794) q[3];
sx q[3];
rz(1.1869528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];