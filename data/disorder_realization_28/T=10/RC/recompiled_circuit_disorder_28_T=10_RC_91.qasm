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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500096) q[0];
sx q[0];
rz(-1.2342493) q[0];
sx q[0];
rz(-2.6495289) q[0];
x q[1];
rz(1.7018868) q[2];
sx q[2];
rz(-2.0799473) q[2];
sx q[2];
rz(-0.57123643) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.441592) q[1];
sx q[1];
rz(-2.4543426) q[1];
sx q[1];
rz(-1.4541461) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22523017) q[3];
sx q[3];
rz(-2.4472144) q[3];
sx q[3];
rz(-3.0519033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(-1.8475378) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(-0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-3.0157715) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(1.7696101) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1945837) q[0];
sx q[0];
rz(-0.8964552) q[0];
sx q[0];
rz(1.6499004) q[0];
x q[1];
rz(2.9397474) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(-1.3943878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0224277) q[1];
sx q[1];
rz(-0.25401527) q[1];
sx q[1];
rz(-0.18175061) q[1];
rz(-pi) q[2];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.3523003) q[3];
sx q[3];
rz(0.0056497638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4662194) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(2.0825785) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7721467) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-0.69338436) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1946482) q[2];
sx q[2];
rz(-2.7333439) q[2];
sx q[2];
rz(2.2779235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.7980305) q[1];
sx q[1];
rz(-2.8312632) q[1];
sx q[1];
rz(1.6088435) q[1];
x q[2];
rz(-0.84633175) q[3];
sx q[3];
rz(-1.4457821) q[3];
sx q[3];
rz(0.66305977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(-0.26322571) q[2];
rz(2.0227382) q[3];
sx q[3];
rz(-1.2759195) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(-1.0385723) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(-2.0746453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.346006) q[0];
sx q[0];
rz(-2.6104402) q[0];
sx q[0];
rz(1.6474433) q[0];
x q[1];
rz(2.8720886) q[2];
sx q[2];
rz(-1.6490893) q[2];
sx q[2];
rz(-1.4769725) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.651557) q[1];
sx q[1];
rz(-2.2228096) q[1];
sx q[1];
rz(-2.0865366) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1569818) q[3];
sx q[3];
rz(-2.6226351) q[3];
sx q[3];
rz(-1.3299691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2944494) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(-2.8520544) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.897215) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(1.429819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0833417) q[0];
sx q[0];
rz(-1.6334851) q[0];
sx q[0];
rz(-1.5682674) q[0];
x q[1];
rz(-0.28508913) q[2];
sx q[2];
rz(-1.3407009) q[2];
sx q[2];
rz(2.6094112) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5809708) q[1];
sx q[1];
rz(-0.79172687) q[1];
sx q[1];
rz(0.03246275) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17794869) q[3];
sx q[3];
rz(-2.2214409) q[3];
sx q[3];
rz(1.4072756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66118583) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(1.2832114) q[2];
rz(0.13218203) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64602393) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(2.6640889) q[0];
rz(-1.6409138) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(0.57055155) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8342499) q[0];
sx q[0];
rz(-1.2120005) q[0];
sx q[0];
rz(2.8310199) q[0];
rz(-pi) q[1];
rz(-2.6887367) q[2];
sx q[2];
rz(-1.9100683) q[2];
sx q[2];
rz(-1.1277744) q[2];
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
rz(1.5930575) q[3];
sx q[3];
rz(-0.92261693) q[3];
sx q[3];
rz(2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4346314) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(2.8213275) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(1.8626574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657848) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(-0.35807034) q[0];
rz(-0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.3105062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7983127) q[0];
sx q[0];
rz(-2.5710921) q[0];
sx q[0];
rz(2.7891854) q[0];
rz(-pi) q[1];
rz(2.4438379) q[2];
sx q[2];
rz(-1.1083833) q[2];
sx q[2];
rz(-2.4042839) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6984182) q[1];
sx q[1];
rz(-2.0768754) q[1];
sx q[1];
rz(-2.0884104) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5345516) q[3];
sx q[3];
rz(-0.97902521) q[3];
sx q[3];
rz(2.4855763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(2.5618662) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(-1.4355481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(-1.2365201) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(-1.9514726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6105881) q[0];
sx q[0];
rz(-1.6325145) q[0];
sx q[0];
rz(-1.3390263) q[0];
x q[1];
rz(-2.2731407) q[2];
sx q[2];
rz(-0.10570279) q[2];
sx q[2];
rz(-2.5603574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5684143) q[1];
sx q[1];
rz(-1.8548994) q[1];
sx q[1];
rz(1.8446288) q[1];
rz(-pi) q[2];
rz(2.2338699) q[3];
sx q[3];
rz(-0.66909664) q[3];
sx q[3];
rz(-2.1818386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.74165806) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(0.41440543) q[2];
rz(-1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(2.8167021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(-1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-2.192416) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28578368) q[0];
sx q[0];
rz(-1.4306418) q[0];
sx q[0];
rz(-2.6967718) q[0];
rz(-2.8434535) q[2];
sx q[2];
rz(-2.4900644) q[2];
sx q[2];
rz(-0.27116129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9526082) q[1];
sx q[1];
rz(-0.70776716) q[1];
sx q[1];
rz(0.37353746) q[1];
x q[2];
rz(-0.90419681) q[3];
sx q[3];
rz(-2.4354746) q[3];
sx q[3];
rz(-1.0977942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-0.43241832) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(-0.37316698) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(-2.4180791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31691566) q[0];
sx q[0];
rz(-2.4002889) q[0];
sx q[0];
rz(-2.9497428) q[0];
x q[1];
rz(-3.0300272) q[2];
sx q[2];
rz(-2.4590883) q[2];
sx q[2];
rz(1.1139368) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8404624) q[1];
sx q[1];
rz(-0.81043078) q[1];
sx q[1];
rz(0.37709548) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5756597) q[3];
sx q[3];
rz(-1.300372) q[3];
sx q[3];
rz(-1.5750386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92131203) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(-0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198467) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
rz(2.8593821) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(-2.0201335) q[2];
sx q[2];
rz(-2.024414) q[2];
sx q[2];
rz(-1.450962) q[2];
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
