OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.763971) q[0];
sx q[0];
rz(-0.42855898) q[0];
sx q[0];
rz(-2.9052486) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(-0.16170391) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16888389) q[0];
sx q[0];
rz(-2.5492269) q[0];
sx q[0];
rz(-1.374872) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5125844) q[2];
sx q[2];
rz(-1.7663284) q[2];
sx q[2];
rz(-2.9511676) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3907416) q[1];
sx q[1];
rz(-0.027719434) q[1];
sx q[1];
rz(1.8230468) q[1];
x q[2];
rz(1.7908279) q[3];
sx q[3];
rz(-1.9373978) q[3];
sx q[3];
rz(-2.0442968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24254313) q[2];
sx q[2];
rz(-0.87729064) q[2];
sx q[2];
rz(0.38929942) q[2];
rz(-0.23046514) q[3];
sx q[3];
rz(-3.1233628) q[3];
sx q[3];
rz(2.5267595) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56735754) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(1.6472598) q[0];
rz(1.5859454) q[1];
sx q[1];
rz(-0.21265282) q[1];
sx q[1];
rz(2.0071323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581438) q[0];
sx q[0];
rz(-0.76144816) q[0];
sx q[0];
rz(-2.7694031) q[0];
x q[1];
rz(-2.0087162) q[2];
sx q[2];
rz(-0.88405245) q[2];
sx q[2];
rz(1.1816292) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4006869) q[1];
sx q[1];
rz(-2.4016651) q[1];
sx q[1];
rz(0.94018971) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5364356) q[3];
sx q[3];
rz(-0.89563771) q[3];
sx q[3];
rz(1.9534115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.050934164) q[2];
sx q[2];
rz(-2.2218573) q[2];
sx q[2];
rz(1.301379) q[2];
rz(-2.0672412) q[3];
sx q[3];
rz(-2.8315872) q[3];
sx q[3];
rz(1.5955135) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12304561) q[0];
sx q[0];
rz(-0.30236852) q[0];
sx q[0];
rz(-2.5240335) q[0];
rz(1.1005719) q[1];
sx q[1];
rz(-0.01958422) q[1];
sx q[1];
rz(2.7380131) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36369187) q[0];
sx q[0];
rz(-3.0263623) q[0];
sx q[0];
rz(-1.4794502) q[0];
rz(-pi) q[1];
rz(1.4694655) q[2];
sx q[2];
rz(-1.0472385) q[2];
sx q[2];
rz(-3.0803185) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61768764) q[1];
sx q[1];
rz(-0.20719658) q[1];
sx q[1];
rz(2.2788384) q[1];
rz(-0.37637122) q[3];
sx q[3];
rz(-2.5507002) q[3];
sx q[3];
rz(-1.8666238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5201716) q[2];
sx q[2];
rz(-1.4474063) q[2];
sx q[2];
rz(2.3604895) q[2];
rz(0.33629867) q[3];
sx q[3];
rz(-1.4399485) q[3];
sx q[3];
rz(-2.4380016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6853365) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(-1.9235032) q[0];
rz(1.7958027) q[1];
sx q[1];
rz(-0.0082052611) q[1];
sx q[1];
rz(1.9075314) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0722375) q[0];
sx q[0];
rz(-1.5610106) q[0];
sx q[0];
rz(2.9517677) q[0];
rz(-2.7379964) q[2];
sx q[2];
rz(-2.4880313) q[2];
sx q[2];
rz(0.60774481) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.36653301) q[1];
sx q[1];
rz(-1.2121736) q[1];
sx q[1];
rz(0.26945646) q[1];
rz(-0.0099700516) q[3];
sx q[3];
rz(-1.5640576) q[3];
sx q[3];
rz(-2.2611195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9464843) q[2];
sx q[2];
rz(-2.7184964) q[2];
sx q[2];
rz(1.1484324) q[2];
rz(2.3871683) q[3];
sx q[3];
rz(-1.9781338) q[3];
sx q[3];
rz(0.26328009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.36963439) q[0];
sx q[0];
rz(-3.0535871) q[0];
sx q[0];
rz(2.5391286) q[0];
rz(-0.98214904) q[1];
sx q[1];
rz(-3.1352477) q[1];
sx q[1];
rz(-0.51012653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54464357) q[0];
sx q[0];
rz(-1.7644797) q[0];
sx q[0];
rz(2.9525443) q[0];
x q[1];
rz(2.192239) q[2];
sx q[2];
rz(-2.992625) q[2];
sx q[2];
rz(-1.1412652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1044961) q[1];
sx q[1];
rz(-1.949661) q[1];
sx q[1];
rz(-2.563743) q[1];
rz(2.3272133) q[3];
sx q[3];
rz(-1.571221) q[3];
sx q[3];
rz(0.18691508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0977352) q[2];
sx q[2];
rz(-2.0243702) q[2];
sx q[2];
rz(2.1065693) q[2];
rz(1.467661) q[3];
sx q[3];
rz(-2.7799455) q[3];
sx q[3];
rz(-0.58290946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3007091) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(-1.0915407) q[0];
rz(-2.7409399) q[1];
sx q[1];
rz(-3.1392097) q[1];
sx q[1];
rz(-2.5750652) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3728722) q[0];
sx q[0];
rz(-2.3864288) q[0];
sx q[0];
rz(-1.7674957) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69366535) q[2];
sx q[2];
rz(-2.1595507) q[2];
sx q[2];
rz(0.47364571) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5623746) q[1];
sx q[1];
rz(-1.7953582) q[1];
sx q[1];
rz(-2.5689305) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4584389) q[3];
sx q[3];
rz(-1.0010202) q[3];
sx q[3];
rz(0.36142413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94622508) q[2];
sx q[2];
rz(-0.75104284) q[2];
sx q[2];
rz(-1.6132149) q[2];
rz(-2.1402806) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(-0.24054578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840851) q[0];
sx q[0];
rz(-0.79428285) q[0];
sx q[0];
rz(-1.1450144) q[0];
rz(1.6192294) q[1];
sx q[1];
rz(-3.122819) q[1];
sx q[1];
rz(-1.9782664) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32864654) q[0];
sx q[0];
rz(-0.29314679) q[0];
sx q[0];
rz(1.4100083) q[0];
rz(-pi) q[1];
rz(-1.7567891) q[2];
sx q[2];
rz(-2.459068) q[2];
sx q[2];
rz(0.96772099) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.12194828) q[1];
sx q[1];
rz(-2.5202529) q[1];
sx q[1];
rz(1.3430194) q[1];
x q[2];
rz(2.9807253) q[3];
sx q[3];
rz(-2.5013574) q[3];
sx q[3];
rz(1.7750027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5160949) q[2];
sx q[2];
rz(-2.1558546) q[2];
sx q[2];
rz(2.7682448) q[2];
rz(0.18800023) q[3];
sx q[3];
rz(-1.5865654) q[3];
sx q[3];
rz(1.9787623) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9173376) q[0];
sx q[0];
rz(-0.52977109) q[0];
sx q[0];
rz(2.8944471) q[0];
rz(-2.9180134) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(1.6252958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2151536) q[0];
sx q[0];
rz(-1.5064539) q[0];
sx q[0];
rz(-1.6872726) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4319576) q[2];
sx q[2];
rz(-0.81474761) q[2];
sx q[2];
rz(0.88569631) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5471955) q[1];
sx q[1];
rz(-1.5548348) q[1];
sx q[1];
rz(2.7475824) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0706625) q[3];
sx q[3];
rz(-2.5462357) q[3];
sx q[3];
rz(-1.8728674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.197) q[2];
sx q[2];
rz(-2.4304515) q[2];
sx q[2];
rz(-1.5468583) q[2];
rz(-2.366015) q[3];
sx q[3];
rz(-1.2838793) q[3];
sx q[3];
rz(-1.6790793) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32770661) q[0];
sx q[0];
rz(-0.96868181) q[0];
sx q[0];
rz(2.0342597) q[0];
rz(-2.2164717) q[1];
sx q[1];
rz(-3.1396301) q[1];
sx q[1];
rz(2.3855239) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7118581) q[0];
sx q[0];
rz(-1.2639771) q[0];
sx q[0];
rz(-1.0064678) q[0];
rz(-1.5816566) q[2];
sx q[2];
rz(-2.6107222) q[2];
sx q[2];
rz(-2.6185991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7489657) q[1];
sx q[1];
rz(-2.4973329) q[1];
sx q[1];
rz(0.49796748) q[1];
rz(-pi) q[2];
rz(-0.1875137) q[3];
sx q[3];
rz(-0.57996677) q[3];
sx q[3];
rz(0.43646508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6206616) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(-1.0010285) q[2];
rz(-1.3641317) q[3];
sx q[3];
rz(-2.2083211) q[3];
sx q[3];
rz(-0.53396839) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022024632) q[0];
sx q[0];
rz(-1.7689995) q[0];
sx q[0];
rz(-0.46510988) q[0];
rz(-1.3348835) q[1];
sx q[1];
rz(-0.3723793) q[1];
sx q[1];
rz(1.575527) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50349277) q[0];
sx q[0];
rz(-1.800312) q[0];
sx q[0];
rz(-3.0969308) q[0];
x q[1];
rz(-0.31120531) q[2];
sx q[2];
rz(-2.0796144) q[2];
sx q[2];
rz(-2.6151163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82619691) q[1];
sx q[1];
rz(-3.1387355) q[1];
sx q[1];
rz(-0.14551659) q[1];
x q[2];
rz(1.6080721) q[3];
sx q[3];
rz(-1.1916984) q[3];
sx q[3];
rz(2.9994734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89432013) q[2];
sx q[2];
rz(-3.0970116) q[2];
sx q[2];
rz(1.0856005) q[2];
rz(-1.9191437) q[3];
sx q[3];
rz(-0.51210755) q[3];
sx q[3];
rz(1.2781757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6923675) q[0];
sx q[0];
rz(-1.4341555) q[0];
sx q[0];
rz(1.7644802) q[0];
rz(-1.5832681) q[1];
sx q[1];
rz(-0.91455864) q[1];
sx q[1];
rz(0.22462489) q[1];
rz(-1.575532) q[2];
sx q[2];
rz(-1.4700459) q[2];
sx q[2];
rz(-2.8503304) q[2];
rz(-2.0381773) q[3];
sx q[3];
rz(-0.49020185) q[3];
sx q[3];
rz(1.310631) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
