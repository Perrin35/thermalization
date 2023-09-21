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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102403) q[0];
sx q[0];
rz(-2.5533479) q[0];
sx q[0];
rz(-2.504185) q[0];
rz(-pi) q[1];
rz(1.7018868) q[2];
sx q[2];
rz(-2.0799473) q[2];
sx q[2];
rz(2.5703562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78046103) q[1];
sx q[1];
rz(-1.6447004) q[1];
sx q[1];
rz(-0.88688811) q[1];
rz(-pi) q[2];
rz(2.4597635) q[3];
sx q[3];
rz(-1.7141984) q[3];
sx q[3];
rz(1.486206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(-1.8475378) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.5016705) q[3];
sx q[3];
rz(0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.92293537) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(1.7696101) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4683511) q[0];
sx q[0];
rz(-1.5090319) q[0];
sx q[0];
rz(0.67586918) q[0];
rz(-pi) q[1];
rz(1.6945018) q[2];
sx q[2];
rz(-1.0281841) q[2];
sx q[2];
rz(1.6306842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0731197) q[1];
sx q[1];
rz(-1.8205376) q[1];
sx q[1];
rz(-1.6176893) q[1];
x q[2];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.7892924) q[3];
sx q[3];
rz(3.1359429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8877318) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(-2.1453693) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(2.0825785) q[0];
rz(1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(1.0645197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76751906) q[0];
sx q[0];
rz(-2.3600183) q[0];
sx q[0];
rz(2.563345) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1946482) q[2];
sx q[2];
rz(-0.40824879) q[2];
sx q[2];
rz(-0.86366913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83798422) q[1];
sx q[1];
rz(-1.260699) q[1];
sx q[1];
rz(-0.012197818) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7581942) q[3];
sx q[3];
rz(-0.73323941) q[3];
sx q[3];
rz(-0.76776615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.069783) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(1.3942962) q[0];
rz(-1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(2.0746453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.346006) q[0];
sx q[0];
rz(-2.6104402) q[0];
sx q[0];
rz(-1.4941494) q[0];
x q[1];
rz(-2.8550451) q[2];
sx q[2];
rz(-0.28038014) q[2];
sx q[2];
rz(-2.7718411) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9001138) q[1];
sx q[1];
rz(-2.334324) q[1];
sx q[1];
rz(2.5680053) q[1];
rz(-2.9158343) q[3];
sx q[3];
rz(-1.0994214) q[3];
sx q[3];
rz(2.2798722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84714326) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(-2.8520544) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(-0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6510058) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(2.3068413) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(-1.429819) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5127038) q[0];
sx q[0];
rz(-1.5682724) q[0];
sx q[0];
rz(-0.062688962) q[0];
rz(-2.8565035) q[2];
sx q[2];
rz(-1.3407009) q[2];
sx q[2];
rz(0.5321815) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5606219) q[1];
sx q[1];
rz(-0.79172687) q[1];
sx q[1];
rz(3.1091299) q[1];
rz(-2.963644) q[3];
sx q[3];
rz(-0.92015172) q[3];
sx q[3];
rz(-1.4072756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66118583) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-2.81288) q[3];
sx q[3];
rz(-0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64602393) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-2.6640889) q[0];
rz(1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(0.57055155) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8342499) q[0];
sx q[0];
rz(-1.9295921) q[0];
sx q[0];
rz(2.8310199) q[0];
rz(-1.19679) q[2];
sx q[2];
rz(-1.1454957) q[2];
sx q[2];
rz(-2.5380295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3444896) q[1];
sx q[1];
rz(-2.2625766) q[1];
sx q[1];
rz(2.2441191) q[1];
rz(-pi) q[2];
rz(0.64829867) q[3];
sx q[3];
rz(-1.5885421) q[3];
sx q[3];
rz(1.363021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70696124) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(0.79664191) q[2];
rz(2.8213275) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(-1.8626574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657848) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(-0.35807034) q[0];
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
rz(1.7983127) q[0];
sx q[0];
rz(-0.57050059) q[0];
sx q[0];
rz(0.35240726) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1475122) q[2];
sx q[2];
rz(-2.1834282) q[2];
sx q[2];
rz(-2.6660369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9999034) q[1];
sx q[1];
rz(-2.018376) q[1];
sx q[1];
rz(0.56772851) q[1];
rz(0.86705039) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(1.5914608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(-2.1441148) q[2];
rz(2.5618662) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(2.8009801) q[0];
rz(-1.2365201) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.1901201) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8461788) q[0];
sx q[0];
rz(-0.23970397) q[0];
sx q[0];
rz(-1.833605) q[0];
x q[1];
rz(-1.6516079) q[2];
sx q[2];
rz(-1.5025856) q[2];
sx q[2];
rz(2.8516172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.076188033) q[1];
sx q[1];
rz(-1.3082062) q[1];
sx q[1];
rz(-2.8471088) q[1];
rz(-pi) q[2];
rz(-0.45300071) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(-0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3999346) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(-1.3828145) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(-0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(-1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(0.94917667) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5700891) q[0];
sx q[0];
rz(-2.6766258) q[0];
sx q[0];
rz(-0.31682195) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5116836) q[2];
sx q[2];
rz(-1.3917149) q[2];
sx q[2];
rz(1.6023028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.28724972) q[1];
sx q[1];
rz(-2.2209475) q[1];
sx q[1];
rz(-1.2681505) q[1];
rz(0.98032326) q[3];
sx q[3];
rz(-1.983641) q[3];
sx q[3];
rz(-1.0124029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40733797) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(2.7091743) q[2];
rz(-1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(-2.4798415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.0712873) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(2.7684257) q[0];
rz(0.35692731) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(0.7235136) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31691566) q[0];
sx q[0];
rz(-0.74130374) q[0];
sx q[0];
rz(2.9497428) q[0];
x q[1];
rz(-1.6610442) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(-1.2573164) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.77955058) q[1];
sx q[1];
rz(-0.83162809) q[1];
sx q[1];
rz(1.9401624) q[1];
rz(-pi) q[2];
rz(3.1240508) q[3];
sx q[3];
rz(-0.27046698) q[3];
sx q[3];
rz(-1.5847575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(0.55345654) q[2];
rz(0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-2.6454906) q[2];
sx q[2];
rz(-1.9719057) q[2];
sx q[2];
rz(3.0531648) q[2];
rz(0.59755748) q[3];
sx q[3];
rz(-2.6664824) q[3];
sx q[3];
rz(0.53371724) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
