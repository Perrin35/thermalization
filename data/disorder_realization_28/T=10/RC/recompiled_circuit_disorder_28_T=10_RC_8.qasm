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
rz(1.0607399) q[1];
sx q[1];
rz(12.386204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59158303) q[0];
sx q[0];
rz(-1.2342493) q[0];
sx q[0];
rz(0.49206375) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9115823) q[2];
sx q[2];
rz(-0.52431528) q[2];
sx q[2];
rz(2.8345248) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3611316) q[1];
sx q[1];
rz(-1.6447004) q[1];
sx q[1];
rz(2.2547045) q[1];
rz(-pi) q[2];
rz(2.4597635) q[3];
sx q[3];
rz(-1.7141984) q[3];
sx q[3];
rz(1.486206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7303598) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(-1.2940548) q[2];
rz(-2.7358352) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(2.7348203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92293537) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(3.0157715) q[0];
rz(0.80548349) q[1];
sx q[1];
rz(-2.3352354) q[1];
sx q[1];
rz(-1.7696101) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94700891) q[0];
sx q[0];
rz(-0.8964552) q[0];
sx q[0];
rz(-1.4916923) q[0];
rz(0.54601045) q[2];
sx q[2];
rz(-1.464932) q[2];
sx q[2];
rz(0.0042303483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.068472915) q[1];
sx q[1];
rz(-1.321055) q[1];
sx q[1];
rz(-1.5239034) q[1];
rz(1.3268746) q[3];
sx q[3];
rz(-2.0274649) q[3];
sx q[3];
rz(-1.4671385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(-0.99622336) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.4662194) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(-2.0825785) q[0];
rz(-1.9937218) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(-1.0645197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3740736) q[0];
sx q[0];
rz(-0.78157434) q[0];
sx q[0];
rz(0.57824761) q[0];
x q[1];
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
rz(0.7980305) q[1];
sx q[1];
rz(-0.31032944) q[1];
sx q[1];
rz(1.5327492) q[1];
rz(-pi) q[2];
rz(0.84633175) q[3];
sx q[3];
rz(-1.6958106) q[3];
sx q[3];
rz(-2.4785329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0718096) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(0.26322571) q[2];
rz(-1.1188544) q[3];
sx q[3];
rz(-1.8656732) q[3];
sx q[3];
rz(-0.81937218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082821) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.3942962) q[0];
rz(-1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(-1.0669473) q[1];
x q[2];
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
rz(1.6520086) q[2];
sx q[2];
rz(-1.8394543) q[2];
sx q[2];
rz(0.072222885) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3921515) q[1];
sx q[1];
rz(-1.1679808) q[1];
sx q[1];
rz(-2.4213326) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9846109) q[3];
sx q[3];
rz(-0.51895751) q[3];
sx q[3];
rz(-1.8116236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2944494) q[2];
sx q[2];
rz(-1.7913982) q[2];
sx q[2];
rz(0.28953826) q[2];
rz(-2.6211522) q[3];
sx q[3];
rz(-1.8013022) q[3];
sx q[3];
rz(0.7590487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-2.3068413) q[0];
rz(-1.2443776) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(-1.429819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058250931) q[0];
sx q[0];
rz(-1.5081076) q[0];
sx q[0];
rz(1.5682674) q[0];
rz(-pi) q[1];
rz(-1.8102112) q[2];
sx q[2];
rz(-1.8481701) q[2];
sx q[2];
rz(1.1053567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5809708) q[1];
sx q[1];
rz(-0.79172687) q[1];
sx q[1];
rz(0.03246275) q[1];
rz(-pi) q[2];
rz(-0.17794869) q[3];
sx q[3];
rz(-0.92015172) q[3];
sx q[3];
rz(1.4072756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4955687) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(1.5006789) q[1];
sx q[1];
rz(-1.6093107) q[1];
sx q[1];
rz(0.57055155) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30734277) q[0];
sx q[0];
rz(-1.2120005) q[0];
sx q[0];
rz(-2.8310199) q[0];
x q[1];
rz(0.678755) q[2];
sx q[2];
rz(-0.55870134) q[2];
sx q[2];
rz(0.15685454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30299444) q[1];
sx q[1];
rz(-2.0717151) q[1];
sx q[1];
rz(2.3272728) q[1];
rz(-pi) q[2];
rz(0.64829867) q[3];
sx q[3];
rz(-1.5885421) q[3];
sx q[3];
rz(1.363021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4346314) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(-2.3449507) q[2];
rz(-0.32026511) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-2.1053574) q[0];
sx q[0];
rz(0.35807034) q[0];
rz(2.8529196) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(-1.8310865) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0689439) q[0];
sx q[0];
rz(-1.7582969) q[0];
sx q[0];
rz(-0.54206538) q[0];
rz(-pi) q[1];
rz(0.65981664) q[2];
sx q[2];
rz(-2.3265127) q[2];
sx q[2];
rz(1.8191402) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3086991) q[1];
sx q[1];
rz(-0.70736865) q[1];
sx q[1];
rz(-2.4127712) q[1];
x q[2];
rz(-0.60704105) q[3];
sx q[3];
rz(-0.97902521) q[3];
sx q[3];
rz(0.65601635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8075809) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(-2.5618662) q[3];
sx q[3];
rz(-0.72967356) q[3];
sx q[3];
rz(-1.4355481) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(-2.8009801) q[0];
rz(1.9050725) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.1901201) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0252359) q[0];
sx q[0];
rz(-1.802117) q[0];
sx q[0];
rz(0.063409253) q[0];
rz(1.4899848) q[2];
sx q[2];
rz(-1.6390071) q[2];
sx q[2];
rz(0.28997544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5684143) q[1];
sx q[1];
rz(-1.8548994) q[1];
sx q[1];
rz(1.8446288) q[1];
x q[2];
rz(2.6885919) q[3];
sx q[3];
rz(-1.060033) q[3];
sx q[3];
rz(-2.9651027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74165806) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(2.7271872) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(0.32489052) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.855809) q[0];
sx q[0];
rz(-1.7109509) q[0];
sx q[0];
rz(-2.6967718) q[0];
rz(-1.7911712) q[2];
sx q[2];
rz(-0.95249635) q[2];
sx q[2];
rz(3.0439723) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8543429) q[1];
sx q[1];
rz(-0.92064511) q[1];
sx q[1];
rz(1.2681505) q[1];
rz(0.48524951) q[3];
sx q[3];
rz(-2.1059548) q[3];
sx q[3];
rz(0.29569611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(1.6010823) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(-0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
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
rz(2.4180791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0821738) q[0];
sx q[0];
rz(-0.84616236) q[0];
sx q[0];
rz(1.7436149) q[0];
rz(2.4621387) q[2];
sx q[2];
rz(-1.6410769) q[2];
sx q[2];
rz(-2.7714504) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53612667) q[1];
sx q[1];
rz(-1.8408753) q[1];
sx q[1];
rz(0.77401604) q[1];
rz(1.5756597) q[3];
sx q[3];
rz(-1.8412207) q[3];
sx q[3];
rz(-1.5750386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.92131203) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(-2.5881361) q[2];
rz(-0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(-1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
