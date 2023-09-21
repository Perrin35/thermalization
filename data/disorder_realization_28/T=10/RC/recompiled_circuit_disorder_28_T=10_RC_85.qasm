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
rz(-2.5500096) q[0];
sx q[0];
rz(-1.9073434) q[0];
sx q[0];
rz(0.49206375) q[0];
x q[1];
rz(1.4397058) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(2.5703562) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78046103) q[1];
sx q[1];
rz(-1.6447004) q[1];
sx q[1];
rz(0.88688811) q[1];
x q[2];
rz(-1.3869242) q[3];
sx q[3];
rz(-2.2443218) q[3];
sx q[3];
rz(2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7303598) q[2];
sx q[2];
rz(-2.2683472) q[2];
sx q[2];
rz(-1.2940548) q[2];
rz(0.40575746) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(2.7348203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3719826) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94700891) q[0];
sx q[0];
rz(-2.2451375) q[0];
sx q[0];
rz(1.4916923) q[0];
rz(-pi) q[1];
rz(-0.54601045) q[2];
sx q[2];
rz(-1.6766607) q[2];
sx q[2];
rz(-3.1373623) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6276715) q[1];
sx q[1];
rz(-1.6162335) q[1];
sx q[1];
rz(-0.25000484) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6847141) q[3];
sx q[3];
rz(-0.51364726) q[3];
sx q[3];
rz(1.9809686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(-1.0960724) q[3];
sx q[3];
rz(-1.4527861) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.6753733) q[0];
sx q[0];
rz(-0.96005625) q[0];
sx q[0];
rz(-2.0825785) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-2.4025326) q[1];
sx q[1];
rz(1.0645197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36944593) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-0.69338436) q[0];
rz(-pi) q[1];
rz(-1.9532922) q[2];
sx q[2];
rz(-1.4244392) q[2];
sx q[2];
rz(1.0548897) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.83798422) q[1];
sx q[1];
rz(-1.8808937) q[1];
sx q[1];
rz(-0.012197818) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9753261) q[3];
sx q[3];
rz(-2.2883752) q[3];
sx q[3];
rz(2.1239514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.069783) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(2.8783669) q[2];
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
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333106) q[0];
sx q[0];
rz(-3.0968956) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(-2.1030203) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(1.0669473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.432503) q[0];
sx q[0];
rz(-1.5320008) q[0];
sx q[0];
rz(-2.1006656) q[0];
rz(2.8720886) q[2];
sx q[2];
rz(-1.6490893) q[2];
sx q[2];
rz(-1.4769725) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9001138) q[1];
sx q[1];
rz(-2.334324) q[1];
sx q[1];
rz(-0.57358731) q[1];
x q[2];
rz(-2.0526485) q[3];
sx q[3];
rz(-1.7715766) q[3];
sx q[3];
rz(2.5364385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.3501945) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
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
rz(-1.8908187) q[1];
sx q[1];
rz(-1.7117737) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058250931) q[0];
sx q[0];
rz(-1.5081076) q[0];
sx q[0];
rz(-1.5682674) q[0];
x q[1];
rz(-2.8565035) q[2];
sx q[2];
rz(-1.3407009) q[2];
sx q[2];
rz(0.5321815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5606219) q[1];
sx q[1];
rz(-2.3498658) q[1];
sx q[1];
rz(-3.1091299) q[1];
x q[2];
rz(0.17794869) q[3];
sx q[3];
rz(-2.2214409) q[3];
sx q[3];
rz(-1.734317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(1.2832114) q[2];
rz(0.13218203) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(2.8018518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4955687) q[0];
sx q[0];
rz(-3.0806354) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(0.57055155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3756838) q[0];
sx q[0];
rz(-1.8609957) q[0];
sx q[0];
rz(-1.1955839) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6887367) q[2];
sx q[2];
rz(-1.9100683) q[2];
sx q[2];
rz(2.0138182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3444896) q[1];
sx q[1];
rz(-0.87901607) q[1];
sx q[1];
rz(-0.89747353) q[1];
rz(-pi) q[2];
rz(3.1122094) q[3];
sx q[3];
rz(-0.64850649) q[3];
sx q[3];
rz(2.9103968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4346314) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(0.79664191) q[2];
rz(2.8213275) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.2789352) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-0.2886731) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.3105062) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863556) q[0];
sx q[0];
rz(-2.1023395) q[0];
sx q[0];
rz(-1.7887572) q[0];
rz(-0.69775478) q[2];
sx q[2];
rz(-1.1083833) q[2];
sx q[2];
rz(-2.4042839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14168921) q[1];
sx q[1];
rz(-1.1232166) q[1];
sx q[1];
rz(2.5738641) q[1];
x q[2];
rz(2.2745423) q[3];
sx q[3];
rz(-2.3208445) q[3];
sx q[3];
rz(-1.5501319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(2.5618662) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(-1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.3842936) q[1];
sx q[1];
rz(-1.1901201) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1163568) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(-0.063409253) q[0];
rz(-pi) q[1];
rz(-2.2731407) q[2];
sx q[2];
rz(-3.0358899) q[2];
sx q[2];
rz(-0.58123523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0654046) q[1];
sx q[1];
rz(-1.3082062) q[1];
sx q[1];
rz(2.8471088) q[1];
rz(-pi) q[2];
rz(2.6885919) q[3];
sx q[3];
rz(-1.060033) q[3];
sx q[3];
rz(0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906616) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(-0.94917667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.21852) q[0];
sx q[0];
rz(-2.0109482) q[0];
sx q[0];
rz(-1.7258304) q[0];
x q[1];
rz(0.62990909) q[2];
sx q[2];
rz(-1.3917149) q[2];
sx q[2];
rz(1.6023028) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8543429) q[1];
sx q[1];
rz(-2.2209475) q[1];
sx q[1];
rz(-1.8734422) q[1];
rz(0.98032326) q[3];
sx q[3];
rz(-1.1579517) q[3];
sx q[3];
rz(-2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7342547) q[2];
sx q[2];
rz(-0.39525017) q[2];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0703053) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(2.7684257) q[0];
rz(2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-0.7235136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3961807) q[0];
sx q[0];
rz(-1.4416845) q[0];
sx q[0];
rz(0.73208916) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0300272) q[2];
sx q[2];
rz(-0.68250436) q[2];
sx q[2];
rz(1.1139368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3620421) q[1];
sx q[1];
rz(-0.83162809) q[1];
sx q[1];
rz(1.9401624) q[1];
rz(-pi) q[2];
rz(3.1240508) q[3];
sx q[3];
rz(-2.8711257) q[3];
sx q[3];
rz(-1.5568352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2202806) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(2.5881361) q[2];
rz(2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2198467) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
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
rz(-2.739493) q[3];
sx q[3];
rz(-1.8310908) q[3];
sx q[3];
rz(2.6487333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
