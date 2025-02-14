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
rz(0.23634401) q[0];
rz(1.4689245) q[1];
sx q[1];
rz(-1.5552712) q[1];
sx q[1];
rz(0.16170391) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40370052) q[0];
sx q[0];
rz(-0.99125112) q[0];
sx q[0];
rz(3.011322) q[0];
rz(-pi) q[1];
rz(-1.6290083) q[2];
sx q[2];
rz(-1.3752642) q[2];
sx q[2];
rz(-0.19042507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0031944) q[1];
sx q[1];
rz(-1.5439543) q[1];
sx q[1];
rz(0.0069199847) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51689619) q[3];
sx q[3];
rz(-2.7166043) q[3];
sx q[3];
rz(1.6551415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8990495) q[2];
sx q[2];
rz(-2.264302) q[2];
sx q[2];
rz(-2.7522932) q[2];
rz(-0.23046514) q[3];
sx q[3];
rz(-0.018229818) q[3];
sx q[3];
rz(-2.5267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56735754) q[0];
sx q[0];
rz(-0.95139545) q[0];
sx q[0];
rz(1.6472598) q[0];
rz(1.5556473) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(2.0071323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0527549) q[0];
sx q[0];
rz(-0.87273926) q[0];
sx q[0];
rz(-1.9044756) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4058543) q[2];
sx q[2];
rz(-1.904907) q[2];
sx q[2];
rz(-0.10057893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74090577) q[1];
sx q[1];
rz(-2.4016651) q[1];
sx q[1];
rz(-2.2014029) q[1];
x q[2];
rz(-1.0026917) q[3];
sx q[3];
rz(-2.3062996) q[3];
sx q[3];
rz(-2.7138674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.050934164) q[2];
sx q[2];
rz(-0.91973534) q[2];
sx q[2];
rz(1.301379) q[2];
rz(2.0672412) q[3];
sx q[3];
rz(-2.8315872) q[3];
sx q[3];
rz(-1.5955135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018547) q[0];
sx q[0];
rz(-0.30236852) q[0];
sx q[0];
rz(2.5240335) q[0];
rz(2.0410208) q[1];
sx q[1];
rz(-3.1220084) q[1];
sx q[1];
rz(-0.40357959) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27173938) q[0];
sx q[0];
rz(-1.4560486) q[0];
sx q[0];
rz(-0.010557584) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4694655) q[2];
sx q[2];
rz(-2.0943542) q[2];
sx q[2];
rz(0.061274139) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8858635) q[1];
sx q[1];
rz(-1.4366062) q[1];
sx q[1];
rz(-1.7291452) q[1];
rz(1.3290492) q[3];
sx q[3];
rz(-1.0261593) q[3];
sx q[3];
rz(-2.310809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6214211) q[2];
sx q[2];
rz(-1.4474063) q[2];
sx q[2];
rz(-0.78110313) q[2];
rz(0.33629867) q[3];
sx q[3];
rz(-1.7016442) q[3];
sx q[3];
rz(-0.70359105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562562) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(1.2180895) q[0];
rz(1.7958027) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(-1.9075314) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.06935519) q[0];
sx q[0];
rz(-1.580582) q[0];
sx q[0];
rz(2.9517677) q[0];
rz(0.40359621) q[2];
sx q[2];
rz(-2.4880313) q[2];
sx q[2];
rz(0.60774481) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36653301) q[1];
sx q[1];
rz(-1.2121736) q[1];
sx q[1];
rz(-0.26945646) q[1];
rz(-pi) q[2];
rz(2.5472121) q[3];
sx q[3];
rz(-0.012033741) q[3];
sx q[3];
rz(-1.2846701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1951083) q[2];
sx q[2];
rz(-0.4230963) q[2];
sx q[2];
rz(1.9931603) q[2];
rz(-0.75442433) q[3];
sx q[3];
rz(-1.9781338) q[3];
sx q[3];
rz(-2.8783126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7719583) q[0];
sx q[0];
rz(-0.088005528) q[0];
sx q[0];
rz(-2.5391286) q[0];
rz(2.1594436) q[1];
sx q[1];
rz(-0.0063449675) q[1];
sx q[1];
rz(-2.6314661) q[1];
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
rz(-pi) q[1];
rz(-2.192239) q[2];
sx q[2];
rz(-2.992625) q[2];
sx q[2];
rz(-2.0003275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.229631) q[1];
sx q[1];
rz(-2.103064) q[1];
sx q[1];
rz(-1.1271354) q[1];
rz(-pi) q[2];
rz(-1.5714151) q[3];
sx q[3];
rz(-0.75641707) q[3];
sx q[3];
rz(-1.3843313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0977352) q[2];
sx q[2];
rz(-2.0243702) q[2];
sx q[2];
rz(-2.1065693) q[2];
rz(-1.6739316) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(-2.5586832) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408836) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(-1.0915407) q[0];
rz(2.7409399) q[1];
sx q[1];
rz(-0.0023829208) q[1];
sx q[1];
rz(0.56652743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5015903) q[0];
sx q[0];
rz(-2.3079607) q[0];
sx q[0];
rz(0.18192525) q[0];
rz(-pi) q[1];
rz(0.80711694) q[2];
sx q[2];
rz(-0.87701488) q[2];
sx q[2];
rz(1.6859695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57921806) q[1];
sx q[1];
rz(-1.3462344) q[1];
sx q[1];
rz(-0.57266219) q[1];
x q[2];
rz(-2.9683365) q[3];
sx q[3];
rz(-0.57954407) q[3];
sx q[3];
rz(-0.15523191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.94622508) q[2];
sx q[2];
rz(-2.3905498) q[2];
sx q[2];
rz(1.6132149) q[2];
rz(2.1402806) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(0.24054578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2840851) q[0];
sx q[0];
rz(-2.3473098) q[0];
sx q[0];
rz(1.9965782) q[0];
rz(-1.5223632) q[1];
sx q[1];
rz(-3.122819) q[1];
sx q[1];
rz(-1.9782664) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0534819) q[0];
sx q[0];
rz(-1.5245175) q[0];
sx q[0];
rz(1.2812216) q[0];
x q[1];
rz(-2.9923963) q[2];
sx q[2];
rz(-0.90221221) q[2];
sx q[2];
rz(-0.72982349) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2625573) q[1];
sx q[1];
rz(-1.702629) q[1];
sx q[1];
rz(0.96178949) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6895503) q[3];
sx q[3];
rz(-0.94014478) q[3];
sx q[3];
rz(-1.5753559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6254977) q[2];
sx q[2];
rz(-2.1558546) q[2];
sx q[2];
rz(2.7682448) q[2];
rz(2.9535924) q[3];
sx q[3];
rz(-1.5865654) q[3];
sx q[3];
rz(1.1628304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.2242551) q[0];
sx q[0];
rz(-2.6118216) q[0];
sx q[0];
rz(2.8944471) q[0];
rz(0.22357926) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(-1.5162969) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36316582) q[0];
sx q[0];
rz(-1.6870305) q[0];
sx q[0];
rz(3.0768125) q[0];
rz(-pi) q[1];
rz(-2.7096351) q[2];
sx q[2];
rz(-0.81474761) q[2];
sx q[2];
rz(-0.88569631) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1113562) q[1];
sx q[1];
rz(-1.9647536) q[1];
sx q[1];
rz(-1.5880821) q[1];
rz(-0.59418847) q[3];
sx q[3];
rz(-1.5310413) q[3];
sx q[3];
rz(0.36082855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94459263) q[2];
sx q[2];
rz(-2.4304515) q[2];
sx q[2];
rz(1.5468583) q[2];
rz(2.366015) q[3];
sx q[3];
rz(-1.2838793) q[3];
sx q[3];
rz(1.6790793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32770661) q[0];
sx q[0];
rz(-0.96868181) q[0];
sx q[0];
rz(-2.0342597) q[0];
rz(-2.2164717) q[1];
sx q[1];
rz(-3.1396301) q[1];
sx q[1];
rz(-0.75606871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8116279) q[0];
sx q[0];
rz(-1.0357619) q[0];
sx q[0];
rz(-2.7828548) q[0];
rz(-pi) q[1];
rz(-2.101641) q[2];
sx q[2];
rz(-1.565298) q[2];
sx q[2];
rz(1.0571684) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9097164) q[1];
sx q[1];
rz(-1.2798339) q[1];
sx q[1];
rz(-2.5582486) q[1];
rz(-pi) q[2];
rz(-1.4492726) q[3];
sx q[3];
rz(-1.0022707) q[3];
sx q[3];
rz(-2.4820676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5209311) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(-2.1405641) q[2];
rz(-1.3641317) q[3];
sx q[3];
rz(-2.2083211) q[3];
sx q[3];
rz(2.6076243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.022024632) q[0];
sx q[0];
rz(-1.3725932) q[0];
sx q[0];
rz(0.46510988) q[0];
rz(-1.3348835) q[1];
sx q[1];
rz(-2.7692134) q[1];
sx q[1];
rz(-1.575527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69746419) q[0];
sx q[0];
rz(-0.23374548) q[0];
sx q[0];
rz(1.759619) q[0];
rz(1.0687629) q[2];
sx q[2];
rz(-0.58922592) q[2];
sx q[2];
rz(3.0844944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82619691) q[1];
sx q[1];
rz(-3.1387355) q[1];
sx q[1];
rz(2.9960761) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0483143) q[3];
sx q[3];
rz(-2.7607548) q[3];
sx q[3];
rz(-0.24254984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89432013) q[2];
sx q[2];
rz(-0.044581052) q[2];
sx q[2];
rz(-1.0856005) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6923675) q[0];
sx q[0];
rz(-1.7074371) q[0];
sx q[0];
rz(-1.3771124) q[0];
rz(1.5583246) q[1];
sx q[1];
rz(-0.91455864) q[1];
sx q[1];
rz(0.22462489) q[1];
rz(-0.046810016) q[2];
sx q[2];
rz(-3.0407314) q[2];
sx q[2];
rz(-2.8032816) q[2];
rz(1.1261945) q[3];
sx q[3];
rz(-1.3570519) q[3];
sx q[3];
rz(-2.9828664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
