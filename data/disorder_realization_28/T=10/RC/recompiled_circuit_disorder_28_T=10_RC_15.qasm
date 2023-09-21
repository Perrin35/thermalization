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
rz(-0.20679064) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(0.18016711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313523) q[0];
sx q[0];
rz(-2.5533479) q[0];
sx q[0];
rz(-2.504185) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23001036) q[2];
sx q[2];
rz(-2.6172774) q[2];
sx q[2];
rz(-0.3070679) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78046103) q[1];
sx q[1];
rz(-1.4968922) q[1];
sx q[1];
rz(-0.88688811) q[1];
x q[2];
rz(1.3869242) q[3];
sx q[3];
rz(-0.89727083) q[3];
sx q[3];
rz(-0.20007381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(1.2940548) q[2];
rz(0.40575746) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(-0.40677235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(2.3361092) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(-1.7696101) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208647) q[0];
sx q[0];
rz(-2.463349) q[0];
sx q[0];
rz(-3.0430549) q[0];
rz(-pi) q[1];
rz(0.20184529) q[2];
sx q[2];
rz(-2.5864374) q[2];
sx q[2];
rz(-1.7472048) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0224277) q[1];
sx q[1];
rz(-0.25401527) q[1];
sx q[1];
rz(-2.959842) q[1];
x q[2];
rz(2.6729229) q[3];
sx q[3];
rz(-1.3523003) q[3];
sx q[3];
rz(3.1359429) q[3];
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
rz(-2.0455202) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(2.2912099) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-2.0825785) q[0];
rz(-1.1478708) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-2.0770729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5110042) q[0];
sx q[0];
rz(-2.2017041) q[0];
sx q[0];
rz(1.0738119) q[0];
x q[1];
rz(1.9532922) q[2];
sx q[2];
rz(-1.7171535) q[2];
sx q[2];
rz(1.0548897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7980305) q[1];
sx q[1];
rz(-2.8312632) q[1];
sx q[1];
rz(1.5327492) q[1];
rz(0.84633175) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(-2.4785329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0718096) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(-2.8783669) q[2];
rz(2.0227382) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(-2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(1.0669473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70908961) q[0];
sx q[0];
rz(-1.5320008) q[0];
sx q[0];
rz(1.0409271) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26950403) q[2];
sx q[2];
rz(-1.4925033) q[2];
sx q[2];
rz(1.6646202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74944118) q[1];
sx q[1];
rz(-1.9736119) q[1];
sx q[1];
rz(-2.4213326) q[1];
rz(2.9158343) q[3];
sx q[3];
rz(-1.0994214) q[3];
sx q[3];
rz(0.86172047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905869) q[0];
sx q[0];
rz(-0.0087954272) q[0];
sx q[0];
rz(-2.3068413) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(-1.429819) q[1];
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
x q[1];
rz(0.69447563) q[2];
sx q[2];
rz(-2.777213) q[2];
sx q[2];
rz(-2.7642872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5606219) q[1];
sx q[1];
rz(-2.3498658) q[1];
sx q[1];
rz(0.03246275) q[1];
rz(-0.17794869) q[3];
sx q[3];
rz(-2.2214409) q[3];
sx q[3];
rz(1.734317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4804068) q[2];
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
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(-2.6640889) q[0];
rz(-1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(2.5710411) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43317023) q[0];
sx q[0];
rz(-0.47009429) q[0];
sx q[0];
rz(-0.8870468) q[0];
rz(1.19679) q[2];
sx q[2];
rz(-1.1454957) q[2];
sx q[2];
rz(-0.60356319) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30299444) q[1];
sx q[1];
rz(-2.0717151) q[1];
sx q[1];
rz(-0.8143199) q[1];
x q[2];
rz(-3.1122094) q[3];
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
rz(-2.091445) q[2];
sx q[2];
rz(-0.79664191) q[2];
rz(-2.8213275) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(1.8626574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(0.35807034) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(-1.3105062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34328) q[0];
sx q[0];
rz(-0.57050059) q[0];
sx q[0];
rz(-0.35240726) q[0];
rz(0.69775478) q[2];
sx q[2];
rz(-2.0332094) q[2];
sx q[2];
rz(0.73730872) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6984182) q[1];
sx q[1];
rz(-1.0647173) q[1];
sx q[1];
rz(2.0884104) q[1];
rz(2.2566124) q[3];
sx q[3];
rz(-1.0776057) q[3];
sx q[3];
rz(2.5964338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(-0.99747783) q[2];
rz(0.57972646) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(-1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(-1.9050725) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(-1.9514726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0252359) q[0];
sx q[0];
rz(-1.802117) q[0];
sx q[0];
rz(3.0781834) q[0];
rz(-0.868452) q[2];
sx q[2];
rz(-0.10570279) q[2];
sx q[2];
rz(-0.58123523) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.076188033) q[1];
sx q[1];
rz(-1.8333865) q[1];
sx q[1];
rz(-0.29448387) q[1];
rz(0.45300071) q[3];
sx q[3];
rz(-1.060033) q[3];
sx q[3];
rz(2.9651027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.74165806) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(0.41440543) q[2];
rz(-1.3828145) q[3];
sx q[3];
rz(-1.4004935) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25093108) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(-1.4258619) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(0.94917667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230726) q[0];
sx q[0];
rz(-2.0109482) q[0];
sx q[0];
rz(-1.4157622) q[0];
rz(-1.3504215) q[2];
sx q[2];
rz(-2.1890963) q[2];
sx q[2];
rz(3.0439723) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6712499) q[1];
sx q[1];
rz(-1.8103231) q[1];
sx q[1];
rz(0.67269477) q[1];
x q[2];
rz(-0.48524951) q[3];
sx q[3];
rz(-2.1059548) q[3];
sx q[3];
rz(2.8458965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7342547) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(1.5405103) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(2.4798415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(-0.35692731) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-0.7235136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7454119) q[0];
sx q[0];
rz(-1.6999082) q[0];
sx q[0];
rz(0.73208916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6610442) q[2];
sx q[2];
rz(-2.2482578) q[2];
sx q[2];
rz(1.2573164) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3620421) q[1];
sx q[1];
rz(-2.3099646) q[1];
sx q[1];
rz(-1.9401624) q[1];
rz(-pi) q[2];
rz(-0.27042737) q[3];
sx q[3];
rz(-1.575483) q[3];
sx q[3];
rz(3.1386496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.7582515) q[2];
sx q[2];
rz(0.55345654) q[2];
rz(2.4297595) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(-2.0994983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
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
rz(0.7278022) q[2];
sx q[2];
rz(-2.5143378) q[2];
sx q[2];
rz(0.85744748) q[2];
rz(-0.40209963) q[3];
sx q[3];
rz(-1.3105018) q[3];
sx q[3];
rz(-0.49285938) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
