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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313523) q[0];
sx q[0];
rz(-0.58824476) q[0];
sx q[0];
rz(-0.63740762) q[0];
rz(-pi) q[1];
rz(1.7018868) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(-2.5703562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85045056) q[1];
sx q[1];
rz(-2.2524815) q[1];
sx q[1];
rz(-0.095231685) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4597635) q[3];
sx q[3];
rz(-1.7141984) q[3];
sx q[3];
rz(1.486206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(1.8475378) q[2];
rz(-2.7358352) q[3];
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
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92293537) q[0];
sx q[0];
rz(-1.0359534) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(1.3719826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67324154) q[0];
sx q[0];
rz(-1.5090319) q[0];
sx q[0];
rz(-0.67586918) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20184529) q[2];
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
rz(-3.0224277) q[1];
sx q[1];
rz(-2.8875774) q[1];
sx q[1];
rz(2.959842) q[1];
rz(-2.6729229) q[3];
sx q[3];
rz(-1.7892924) q[3];
sx q[3];
rz(3.1359429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25386086) q[2];
sx q[2];
rz(-2.4527206) q[2];
sx q[2];
rz(0.99622336) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(-0.85038275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-2.0825785) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-1.0645197) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6305884) q[0];
sx q[0];
rz(-2.2017041) q[0];
sx q[0];
rz(-1.0738119) q[0];
rz(-pi) q[1];
rz(0.15757615) q[2];
sx q[2];
rz(-1.1925979) q[2];
sx q[2];
rz(2.6842897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.73653446) q[1];
sx q[1];
rz(-1.5824123) q[1];
sx q[1];
rz(1.2606773) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9753261) q[3];
sx q[3];
rz(-2.2883752) q[3];
sx q[3];
rz(-2.1239514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.069783) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(0.26322571) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(-0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.5082821) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(-1.3942962) q[0];
rz(-1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(2.0746453) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2571714) q[0];
sx q[0];
rz(-2.1002249) q[0];
sx q[0];
rz(-3.0966395) q[0];
rz(-1.6520086) q[2];
sx q[2];
rz(-1.8394543) q[2];
sx q[2];
rz(3.0693698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2414788) q[1];
sx q[1];
rz(-0.80726868) q[1];
sx q[1];
rz(-2.5680053) q[1];
rz(1.1569818) q[3];
sx q[3];
rz(-2.6226351) q[3];
sx q[3];
rz(1.3299691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84714326) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(-0.28953826) q[2];
rz(-2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(-2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905869) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(-0.83475137) q[0];
rz(1.897215) q[1];
sx q[1];
rz(-1.8908187) q[1];
sx q[1];
rz(1.429819) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017905047) q[0];
sx q[0];
rz(-0.062739685) q[0];
sx q[0];
rz(0.040266589) q[0];
x q[1];
rz(0.28508913) q[2];
sx q[2];
rz(-1.8008917) q[2];
sx q[2];
rz(-0.5321815) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1289542) q[1];
sx q[1];
rz(-1.5938938) q[1];
sx q[1];
rz(2.3501293) q[1];
x q[2];
rz(-2.2291282) q[3];
sx q[3];
rz(-1.7121127) q[3];
sx q[3];
rz(-0.055012881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66118583) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(-1.2832114) q[2];
rz(-3.0094106) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(-0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64602393) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(-0.47750372) q[0];
rz(-1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(2.5710411) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3756838) q[0];
sx q[0];
rz(-1.8609957) q[0];
sx q[0];
rz(1.9460088) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6887367) q[2];
sx q[2];
rz(-1.2315244) q[2];
sx q[2];
rz(1.1277744) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3444896) q[1];
sx q[1];
rz(-0.87901607) q[1];
sx q[1];
rz(-0.89747353) q[1];
rz(-pi) q[2];
rz(2.493294) q[3];
sx q[3];
rz(-1.5530506) q[3];
sx q[3];
rz(1.363021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(-2.7835223) q[0];
rz(2.8529196) q[1];
sx q[1];
rz(-2.6117548) q[1];
sx q[1];
rz(-1.3105062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3863556) q[0];
sx q[0];
rz(-1.0392531) q[0];
sx q[0];
rz(-1.3528354) q[0];
x q[1];
rz(2.1475122) q[2];
sx q[2];
rz(-2.1834282) q[2];
sx q[2];
rz(0.47555579) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3086991) q[1];
sx q[1];
rz(-2.434224) q[1];
sx q[1];
rz(0.72882148) q[1];
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
rz(-pi) q[1];
rz(0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(-2.1441148) q[2];
rz(-2.5618662) q[3];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26577935) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(-0.34061256) q[0];
rz(-1.2365201) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.1901201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2954138) q[0];
sx q[0];
rz(-2.9018887) q[0];
sx q[0];
rz(-1.3079877) q[0];
rz(1.6516079) q[2];
sx q[2];
rz(-1.5025856) q[2];
sx q[2];
rz(0.28997544) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3548673) q[1];
sx q[1];
rz(-0.39199542) q[1];
sx q[1];
rz(0.74704945) q[1];
rz(2.6885919) q[3];
sx q[3];
rz(-1.060033) q[3];
sx q[3];
rz(0.1764899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3999346) q[2];
sx q[2];
rz(-2.2592893) q[2];
sx q[2];
rz(-0.41440543) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25093108) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(0.1517621) q[0];
rz(-1.4258619) q[1];
sx q[1];
rz(-1.7646004) q[1];
sx q[1];
rz(2.192416) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9230726) q[0];
sx q[0];
rz(-2.0109482) q[0];
sx q[0];
rz(1.4157622) q[0];
rz(-pi) q[1];
rz(0.29813913) q[2];
sx q[2];
rz(-2.4900644) q[2];
sx q[2];
rz(2.8704314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8543429) q[1];
sx q[1];
rz(-0.92064511) q[1];
sx q[1];
rz(-1.2681505) q[1];
rz(2.1612694) q[3];
sx q[3];
rz(-1.1579517) q[3];
sx q[3];
rz(2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.40733797) q[2];
sx q[2];
rz(-0.39525017) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.5099022) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0703053) q[0];
sx q[0];
rz(-2.8685331) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(0.35692731) q[1];
sx q[1];
rz(-0.86078763) q[1];
sx q[1];
rz(2.4180791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7454119) q[0];
sx q[0];
rz(-1.6999082) q[0];
sx q[0];
rz(0.73208916) q[0];
rz(-pi) q[1];
rz(-3.0300272) q[2];
sx q[2];
rz(-2.4590883) q[2];
sx q[2];
rz(-2.0276558) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3620421) q[1];
sx q[1];
rz(-0.83162809) q[1];
sx q[1];
rz(-1.9401624) q[1];
rz(-pi) q[2];
rz(-0.017541842) q[3];
sx q[3];
rz(-2.8711257) q[3];
sx q[3];
rz(-1.5568352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.92131203) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-0.55345654) q[2];
rz(2.4297595) q[3];
sx q[3];
rz(-0.40829855) q[3];
sx q[3];
rz(1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.2198467) q[0];
sx q[0];
rz(-2.5384359) q[0];
sx q[0];
rz(0.13721101) q[0];
rz(0.28221054) q[1];
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
