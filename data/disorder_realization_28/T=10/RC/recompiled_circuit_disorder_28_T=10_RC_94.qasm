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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59158303) q[0];
sx q[0];
rz(-1.2342493) q[0];
sx q[0];
rz(-0.49206375) q[0];
rz(-1.4397058) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(0.57123643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.78046103) q[1];
sx q[1];
rz(-1.4968922) q[1];
sx q[1];
rz(-0.88688811) q[1];
x q[2];
rz(2.4597635) q[3];
sx q[3];
rz(-1.4273943) q[3];
sx q[3];
rz(1.6553866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41123286) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(-1.2940548) q[2];
rz(2.7358352) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(-0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2186573) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(-0.12582114) q[0];
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
rz(-2.3208647) q[0];
sx q[0];
rz(-0.67824368) q[0];
sx q[0];
rz(-0.098537785) q[0];
rz(-1.4470909) q[2];
sx q[2];
rz(-2.1134085) q[2];
sx q[2];
rz(-1.6306842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.068472915) q[1];
sx q[1];
rz(-1.321055) q[1];
sx q[1];
rz(1.6176893) q[1];
x q[2];
rz(-1.814718) q[3];
sx q[3];
rz(-2.0274649) q[3];
sx q[3];
rz(1.6744542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(-0.99622336) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(-0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4662194) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(2.0825785) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(1.0645197) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7721467) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-2.4482083) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15757615) q[2];
sx q[2];
rz(-1.1925979) q[2];
sx q[2];
rz(-0.45730293) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4050582) q[1];
sx q[1];
rz(-1.5824123) q[1];
sx q[1];
rz(-1.8809153) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2952609) q[3];
sx q[3];
rz(-1.4457821) q[3];
sx q[3];
rz(0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(2.8783669) q[2];
rz(-1.1188544) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(1.0669473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79558668) q[0];
sx q[0];
rz(-0.53115244) q[0];
sx q[0];
rz(-1.4941494) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8550451) q[2];
sx q[2];
rz(-2.8612125) q[2];
sx q[2];
rz(2.7718411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3921515) q[1];
sx q[1];
rz(-1.9736119) q[1];
sx q[1];
rz(2.4213326) q[1];
rz(1.0889441) q[3];
sx q[3];
rz(-1.370016) q[3];
sx q[3];
rz(0.60515412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2944494) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(-0.28953826) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(-0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(0.83475137) q[0];
rz(-1.897215) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.7117737) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5127038) q[0];
sx q[0];
rz(-1.5733203) q[0];
sx q[0];
rz(3.0789037) q[0];
rz(-pi) q[1];
rz(-1.3313815) q[2];
sx q[2];
rz(-1.8481701) q[2];
sx q[2];
rz(2.036236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1289542) q[1];
sx q[1];
rz(-1.5938938) q[1];
sx q[1];
rz(0.79146339) q[1];
rz(-pi) q[2];
x q[2];
rz(0.91246446) q[3];
sx q[3];
rz(-1.4294799) q[3];
sx q[3];
rz(0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(-2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64602393) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-2.6640889) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-2.5710411) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7659089) q[0];
sx q[0];
rz(-1.8609957) q[0];
sx q[0];
rz(1.9460088) q[0];
x q[1];
rz(2.4628377) q[2];
sx q[2];
rz(-2.5828913) q[2];
sx q[2];
rz(0.15685454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3444896) q[1];
sx q[1];
rz(-2.2625766) q[1];
sx q[1];
rz(-0.89747353) q[1];
x q[2];
rz(1.5930575) q[3];
sx q[3];
rz(-0.92261693) q[3];
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
rz(-1.0501477) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(-0.32026511) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(1.8626574) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(-2.7835223) q[0];
rz(-2.8529196) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.8310865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.755237) q[0];
sx q[0];
rz(-1.0392531) q[0];
sx q[0];
rz(1.3528354) q[0];
rz(-2.1475122) q[2];
sx q[2];
rz(-0.95816441) q[2];
sx q[2];
rz(0.47555579) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3086991) q[1];
sx q[1];
rz(-2.434224) q[1];
sx q[1];
rz(-2.4127712) q[1];
rz(-pi) q[2];
rz(2.2566124) q[3];
sx q[3];
rz(-2.0639869) q[3];
sx q[3];
rz(-2.5964338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8075809) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(2.1441148) q[2];
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
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1163568) q[0];
sx q[0];
rz(-1.802117) q[0];
sx q[0];
rz(0.063409253) q[0];
rz(1.6516079) q[2];
sx q[2];
rz(-1.6390071) q[2];
sx q[2];
rz(2.8516172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.076188033) q[1];
sx q[1];
rz(-1.8333865) q[1];
sx q[1];
rz(0.29448387) q[1];
rz(-pi) q[2];
rz(-2.6885919) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3999346) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-0.41440543) q[2];
rz(-1.7587781) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(-1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(0.94917667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5715036) q[0];
sx q[0];
rz(-2.6766258) q[0];
sx q[0];
rz(0.31682195) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7911712) q[2];
sx q[2];
rz(-0.95249635) q[2];
sx q[2];
rz(0.097620336) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4703428) q[1];
sx q[1];
rz(-1.3312695) q[1];
sx q[1];
rz(2.4688979) q[1];
x q[2];
rz(-0.90419681) q[3];
sx q[3];
rz(-0.70611806) q[3];
sx q[3];
rz(1.0977942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-2.7091743) q[2];
rz(-1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(-2.4798415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0703053) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(-2.7684257) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-0.7235136) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3961807) q[0];
sx q[0];
rz(-1.4416845) q[0];
sx q[0];
rz(-0.73208916) q[0];
rz(-pi) q[1];
rz(-1.6610442) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(1.8842763) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.605466) q[1];
sx q[1];
rz(-1.3007174) q[1];
sx q[1];
rz(2.3675766) q[1];
rz(0.27042737) q[3];
sx q[3];
rz(-1.575483) q[3];
sx q[3];
rz(0.0029431012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(-2.5881361) q[2];
rz(-2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(1.0420943) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(2.8593821) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(-0.49610207) q[2];
sx q[2];
rz(-1.169687) q[2];
sx q[2];
rz(-0.088427831) q[2];
rz(-2.5440352) q[3];
sx q[3];
rz(-2.6664824) q[3];
sx q[3];
rz(0.53371724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];