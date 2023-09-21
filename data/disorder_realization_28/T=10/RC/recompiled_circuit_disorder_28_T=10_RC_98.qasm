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
rz(-2.9614255) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3375815) q[0];
sx q[0];
rz(-1.1085701) q[0];
sx q[0];
rz(-1.1929212) q[0];
x q[1];
rz(1.4397058) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(2.5703562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85045056) q[1];
sx q[1];
rz(-2.2524815) q[1];
sx q[1];
rz(0.095231685) q[1];
rz(-pi) q[2];
rz(0.68182919) q[3];
sx q[3];
rz(-1.4273943) q[3];
sx q[3];
rz(-1.6553866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92293537) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(3.0157715) q[0];
rz(2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(-1.7696101) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208647) q[0];
sx q[0];
rz(-2.463349) q[0];
sx q[0];
rz(0.098537785) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20184529) q[2];
sx q[2];
rz(-0.55515528) q[2];
sx q[2];
rz(-1.7472048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5139211) q[1];
sx q[1];
rz(-1.6162335) q[1];
sx q[1];
rz(0.25000484) q[1];
rz(-pi) q[2];
rz(0.45687859) q[3];
sx q[3];
rz(-0.51364726) q[3];
sx q[3];
rz(-1.160624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.25386086) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(-0.99622336) q[2];
rz(-1.0960724) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(1.0590142) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(1.0645197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6305884) q[0];
sx q[0];
rz(-2.2017041) q[0];
sx q[0];
rz(-2.0677807) q[0];
rz(1.9532922) q[2];
sx q[2];
rz(-1.4244392) q[2];
sx q[2];
rz(2.086703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7980305) q[1];
sx q[1];
rz(-0.31032944) q[1];
sx q[1];
rz(-1.5327492) q[1];
rz(1.3833984) q[3];
sx q[3];
rz(-2.4083532) q[3];
sx q[3];
rz(-0.76776615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5082821) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(-1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.4515406) q[1];
sx q[1];
rz(2.0746453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.346006) q[0];
sx q[0];
rz(-0.53115244) q[0];
sx q[0];
rz(-1.6474433) q[0];
rz(1.6520086) q[2];
sx q[2];
rz(-1.8394543) q[2];
sx q[2];
rz(-3.0693698) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3921515) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(-0.72026003) q[1];
rz(-pi) q[2];
rz(-1.9846109) q[3];
sx q[3];
rz(-0.51895751) q[3];
sx q[3];
rz(1.8116236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.84714326) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4905869) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(0.83475137) q[0];
rz(-1.897215) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(-1.429819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017905047) q[0];
sx q[0];
rz(-3.078853) q[0];
sx q[0];
rz(-0.040266589) q[0];
rz(-0.69447563) q[2];
sx q[2];
rz(-2.777213) q[2];
sx q[2];
rz(-0.37730544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5606219) q[1];
sx q[1];
rz(-0.79172687) q[1];
sx q[1];
rz(-0.03246275) q[1];
rz(-pi) q[2];
rz(0.17794869) q[3];
sx q[3];
rz(-0.92015172) q[3];
sx q[3];
rz(1.734317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4804068) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(1.8583813) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(-0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(1.6409138) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-2.5710411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7659089) q[0];
sx q[0];
rz(-1.280597) q[0];
sx q[0];
rz(1.9460088) q[0];
rz(-pi) q[1];
rz(-0.45285593) q[2];
sx q[2];
rz(-1.9100683) q[2];
sx q[2];
rz(1.1277744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4482566) q[1];
sx q[1];
rz(-0.92474557) q[1];
sx q[1];
rz(2.4962884) q[1];
rz(-pi) q[2];
rz(1.5485351) q[3];
sx q[3];
rz(-2.2189757) q[3];
sx q[3];
rz(-0.19433403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4346314) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(-0.32026511) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(-1.8626574) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0657848) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(-0.35807034) q[0];
rz(-2.8529196) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(1.3105062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072648777) q[0];
sx q[0];
rz(-1.3832958) q[0];
sx q[0];
rz(-0.54206538) q[0];
rz(-2.481776) q[2];
sx q[2];
rz(-2.3265127) q[2];
sx q[2];
rz(-1.3224524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6984182) q[1];
sx q[1];
rz(-1.0647173) q[1];
sx q[1];
rz(-2.0884104) q[1];
rz(-pi) q[2];
rz(-2.2566124) q[3];
sx q[3];
rz(-2.0639869) q[3];
sx q[3];
rz(-0.54515884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(-2.1441148) q[2];
rz(-2.5618662) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(1.7060446) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26577935) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(0.34061256) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(-1.1901201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0252359) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(0.063409253) q[0];
rz(-1.6516079) q[2];
sx q[2];
rz(-1.5025856) q[2];
sx q[2];
rz(-0.28997544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.076188033) q[1];
sx q[1];
rz(-1.8333865) q[1];
sx q[1];
rz(2.8471088) q[1];
x q[2];
rz(-2.1281151) q[3];
sx q[3];
rz(-1.1790457) q[3];
sx q[3];
rz(-1.980892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74165806) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(1.7587781) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5700891) q[0];
sx q[0];
rz(-0.46496689) q[0];
sx q[0];
rz(-0.31682195) q[0];
rz(2.5116836) q[2];
sx q[2];
rz(-1.7498778) q[2];
sx q[2];
rz(-1.5392898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6712499) q[1];
sx q[1];
rz(-1.3312695) q[1];
sx q[1];
rz(-2.4688979) q[1];
rz(-pi) q[2];
rz(-0.90419681) q[3];
sx q[3];
rz(-2.4354746) q[3];
sx q[3];
rz(-1.0977942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(-2.7091743) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-2.4798415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0712873) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(-2.7684257) q[0];
rz(-2.7846653) q[1];
sx q[1];
rz(-0.86078763) q[1];
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
rz(-1.6999082) q[0];
sx q[0];
rz(0.73208916) q[0];
rz(-pi) q[1];
rz(-3.0300272) q[2];
sx q[2];
rz(-0.68250436) q[2];
sx q[2];
rz(2.0276558) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.605466) q[1];
sx q[1];
rz(-1.3007174) q[1];
sx q[1];
rz(-0.77401604) q[1];
x q[2];
rz(-2.8711653) q[3];
sx q[3];
rz(-1.575483) q[3];
sx q[3];
rz(0.0029431012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(0.7278022) q[2];
sx q[2];
rz(-2.5143378) q[2];
sx q[2];
rz(0.85744748) q[2];
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