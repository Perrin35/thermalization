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
rz(-2.7517154) q[1];
sx q[1];
rz(1.0607399) q[1];
sx q[1];
rz(12.386204) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102403) q[0];
sx q[0];
rz(-2.5533479) q[0];
sx q[0];
rz(-2.504185) q[0];
rz(-0.51282672) q[2];
sx q[2];
rz(-1.6851808) q[2];
sx q[2];
rz(-2.0778542) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3611316) q[1];
sx q[1];
rz(-1.4968922) q[1];
sx q[1];
rz(2.2547045) q[1];
x q[2];
rz(-1.7546685) q[3];
sx q[3];
rz(-0.89727083) q[3];
sx q[3];
rz(2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41123286) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
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
rz(-0.80635726) q[1];
sx q[1];
rz(-1.7696101) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.820728) q[0];
sx q[0];
rz(-2.463349) q[0];
sx q[0];
rz(-0.098537785) q[0];
rz(2.5955822) q[2];
sx q[2];
rz(-1.6766607) q[2];
sx q[2];
rz(-3.1373623) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0731197) q[1];
sx q[1];
rz(-1.8205376) q[1];
sx q[1];
rz(-1.6176893) q[1];
rz(-pi) q[2];
rz(2.6847141) q[3];
sx q[3];
rz(-2.6279454) q[3];
sx q[3];
rz(1.9809686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(2.1453693) q[2];
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
rz(-pi/2) q[2];
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
rz(-2.0770729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3740736) q[0];
sx q[0];
rz(-2.3600183) q[0];
sx q[0];
rz(2.563345) q[0];
rz(-pi) q[1];
rz(1.1883005) q[2];
sx q[2];
rz(-1.4244392) q[2];
sx q[2];
rz(-2.086703) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3435622) q[1];
sx q[1];
rz(-0.31032944) q[1];
sx q[1];
rz(1.5327492) q[1];
rz(-pi) q[2];
rz(2.9753261) q[3];
sx q[3];
rz(-0.85321745) q[3];
sx q[3];
rz(-1.0176413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.069783) q[2];
sx q[2];
rz(-1.1867384) q[2];
sx q[2];
rz(-0.26322571) q[2];
rz(-1.1188544) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(-0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082821) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.3942962) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(1.0669473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70908961) q[0];
sx q[0];
rz(-1.6095918) q[0];
sx q[0];
rz(1.0409271) q[0];
rz(-2.8720886) q[2];
sx q[2];
rz(-1.6490893) q[2];
sx q[2];
rz(-1.6646202) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3921515) q[1];
sx q[1];
rz(-1.9736119) q[1];
sx q[1];
rz(0.72026003) q[1];
rz(-pi) q[2];
rz(-1.9846109) q[3];
sx q[3];
rz(-2.6226351) q[3];
sx q[3];
rz(-1.8116236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(2.8520544) q[2];
rz(0.52044049) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(2.382544) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.2507739) q[1];
sx q[1];
rz(1.7117737) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058250931) q[0];
sx q[0];
rz(-1.6334851) q[0];
sx q[0];
rz(-1.5682674) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.447117) q[2];
sx q[2];
rz(-2.777213) q[2];
sx q[2];
rz(-2.7642872) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6068078) q[1];
sx q[1];
rz(-2.3619898) q[1];
sx q[1];
rz(1.5379377) q[1];
x q[2];
rz(2.963644) q[3];
sx q[3];
rz(-2.2214409) q[3];
sx q[3];
rz(1.734317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66118583) q[2];
sx q[2];
rz(-1.0057665) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(3.0094106) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64602393) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(-1.5006789) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-2.5710411) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8342499) q[0];
sx q[0];
rz(-1.9295921) q[0];
sx q[0];
rz(-0.3105727) q[0];
rz(-2.4628377) q[2];
sx q[2];
rz(-0.55870134) q[2];
sx q[2];
rz(0.15685454) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4482566) q[1];
sx q[1];
rz(-2.2168471) q[1];
sx q[1];
rz(-2.4962884) q[1];
rz(-0.029383226) q[3];
sx q[3];
rz(-0.64850649) q[3];
sx q[3];
rz(-0.23119584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4346314) q[2];
sx q[2];
rz(-2.091445) q[2];
sx q[2];
rz(2.3449507) q[2];
rz(2.8213275) q[3];
sx q[3];
rz(-1.0864778) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(0.35807034) q[0];
rz(2.8529196) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(1.8310865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34328) q[0];
sx q[0];
rz(-2.5710921) q[0];
sx q[0];
rz(0.35240726) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69775478) q[2];
sx q[2];
rz(-1.1083833) q[2];
sx q[2];
rz(-0.73730872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14168921) q[1];
sx q[1];
rz(-1.1232166) q[1];
sx q[1];
rz(-0.56772851) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5345516) q[3];
sx q[3];
rz(-2.1625674) q[3];
sx q[3];
rz(2.4855763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8075809) q[2];
sx q[2];
rz(-0.58471218) q[2];
sx q[2];
rz(2.1441148) q[2];
rz(2.5618662) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.4900603) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(1.2365201) q[1];
sx q[1];
rz(-0.75729901) q[1];
sx q[1];
rz(-1.9514726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1163568) q[0];
sx q[0];
rz(-1.3394757) q[0];
sx q[0];
rz(0.063409253) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.868452) q[2];
sx q[2];
rz(-3.0358899) q[2];
sx q[2];
rz(0.58123523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78672532) q[1];
sx q[1];
rz(-2.7495972) q[1];
sx q[1];
rz(2.3945432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6885919) q[3];
sx q[3];
rz(-2.0815597) q[3];
sx q[3];
rz(-0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74165806) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(1.7587781) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(-0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25093108) q[0];
sx q[0];
rz(-2.119976) q[0];
sx q[0];
rz(-0.1517621) q[0];
rz(1.7157308) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(0.94917667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28578368) q[0];
sx q[0];
rz(-1.7109509) q[0];
sx q[0];
rz(0.44482081) q[0];
x q[1];
rz(0.62990909) q[2];
sx q[2];
rz(-1.3917149) q[2];
sx q[2];
rz(-1.5392898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28724972) q[1];
sx q[1];
rz(-0.92064511) q[1];
sx q[1];
rz(1.2681505) q[1];
rz(-pi) q[2];
rz(-0.48524951) q[3];
sx q[3];
rz(-1.0356379) q[3];
sx q[3];
rz(-2.8458965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(-2.7091743) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0712873) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(-2.7684257) q[0];
rz(0.35692731) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(-0.7235136) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3961807) q[0];
sx q[0];
rz(-1.4416845) q[0];
sx q[0];
rz(2.4095035) q[0];
rz(-pi) q[1];
rz(-1.4805484) q[2];
sx q[2];
rz(-0.89333488) q[2];
sx q[2];
rz(-1.8842763) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.605466) q[1];
sx q[1];
rz(-1.3007174) q[1];
sx q[1];
rz(-2.3675766) q[1];
rz(-pi) q[2];
x q[2];
rz(0.017541842) q[3];
sx q[3];
rz(-0.27046698) q[3];
sx q[3];
rz(1.5847575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92131203) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-2.5881361) q[2];
rz(2.4297595) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9217459) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(-0.28221054) q[1];
sx q[1];
rz(-2.3574775) q[1];
sx q[1];
rz(-2.2025253) q[1];
rz(0.49610207) q[2];
sx q[2];
rz(-1.9719057) q[2];
sx q[2];
rz(3.0531648) q[2];
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
