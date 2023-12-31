OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8310175) q[0];
sx q[0];
rz(-2.0895045) q[0];
sx q[0];
rz(-1.6488099) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(-1.0746497) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2079017) q[0];
sx q[0];
rz(-1.5655193) q[0];
sx q[0];
rz(-0.018215608) q[0];
x q[1];
rz(-2.3518402) q[2];
sx q[2];
rz(-1.4573922) q[2];
sx q[2];
rz(1.5838768) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5684811) q[1];
sx q[1];
rz(-0.1460349) q[1];
sx q[1];
rz(-2.5558429) q[1];
rz(-pi) q[2];
rz(0.53332897) q[3];
sx q[3];
rz(-1.1733857) q[3];
sx q[3];
rz(-1.0226137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.5011903) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.879724) q[0];
sx q[0];
rz(-0.12517087) q[0];
sx q[0];
rz(-2.9414888) q[0];
rz(-pi) q[1];
rz(1.00374) q[2];
sx q[2];
rz(-2.6008285) q[2];
sx q[2];
rz(2.8218249) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7565631) q[1];
sx q[1];
rz(-0.90365138) q[1];
sx q[1];
rz(-1.3515616) q[1];
x q[2];
rz(1.9618481) q[3];
sx q[3];
rz(-1.5763596) q[3];
sx q[3];
rz(-0.55913505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8555277) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(-1.9747915) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.5037781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8381739) q[0];
sx q[0];
rz(-1.5056599) q[0];
sx q[0];
rz(1.6584048) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9861228) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(-0.090678064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5895264) q[1];
sx q[1];
rz(-2.8376736) q[1];
sx q[1];
rz(0.1429096) q[1];
rz(-3.010473) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(-2.2281447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.1733615) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(2.0898576) q[0];
rz(-0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.9925041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968397) q[0];
sx q[0];
rz(-1.4595932) q[0];
sx q[0];
rz(2.8993594) q[0];
x q[1];
rz(-0.43196584) q[2];
sx q[2];
rz(-2.1199806) q[2];
sx q[2];
rz(0.84749046) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9149575) q[1];
sx q[1];
rz(-1.2775704) q[1];
sx q[1];
rz(-2.3073767) q[1];
rz(-pi) q[2];
rz(2.3991688) q[3];
sx q[3];
rz(-1.8291049) q[3];
sx q[3];
rz(-1.0148259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(-1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.4887061) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-1.0523798) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-0.88422424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0852172) q[0];
sx q[0];
rz(-1.3803122) q[0];
sx q[0];
rz(-2.9604244) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3166979) q[2];
sx q[2];
rz(-0.56606228) q[2];
sx q[2];
rz(2.538946) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65518889) q[1];
sx q[1];
rz(-0.15895325) q[1];
sx q[1];
rz(0.43227355) q[1];
rz(0.98975011) q[3];
sx q[3];
rz(-0.61468609) q[3];
sx q[3];
rz(-0.97782545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(-0.70619839) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(0.96907369) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1554402) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(-0.040779671) q[0];
rz(-pi) q[1];
rz(1.5029807) q[2];
sx q[2];
rz(-0.8578476) q[2];
sx q[2];
rz(0.80578795) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59086696) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(-0.73901891) q[1];
rz(-pi) q[2];
rz(1.045741) q[3];
sx q[3];
rz(-1.1610371) q[3];
sx q[3];
rz(-2.4214782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51450729) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(-1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.4356027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.766151) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(-0.96813162) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1232386) q[2];
sx q[2];
rz(-1.5510501) q[2];
sx q[2];
rz(0.28723082) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.095852921) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(2.2687001) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0649879) q[3];
sx q[3];
rz(-1.1899663) q[3];
sx q[3];
rz(-1.3169469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4009565) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62194659) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(-1.4343065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.442996) q[2];
sx q[2];
rz(-1.5489849) q[2];
sx q[2];
rz(-0.70805659) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8088088) q[1];
sx q[1];
rz(-2.60751) q[1];
sx q[1];
rz(-0.10314718) q[1];
x q[2];
rz(0.41994628) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(-1.9672293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0662971) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-0.31059206) q[0];
rz(-0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(-0.92393595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5101178) q[0];
sx q[0];
rz(-0.81561136) q[0];
sx q[0];
rz(2.5328013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.02248259) q[2];
sx q[2];
rz(-2.6528931) q[2];
sx q[2];
rz(1.894001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9388705) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(2.2321781) q[1];
x q[2];
rz(0.95122533) q[3];
sx q[3];
rz(-0.95622548) q[3];
sx q[3];
rz(0.28704498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(-1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-0.038169233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036457688) q[0];
sx q[0];
rz(-1.5974701) q[0];
sx q[0];
rz(2.1795576) q[0];
rz(-pi) q[1];
rz(1.9534692) q[2];
sx q[2];
rz(-1.4854991) q[2];
sx q[2];
rz(1.485699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.888615) q[1];
sx q[1];
rz(-1.9441609) q[1];
sx q[1];
rz(-3.0967767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2323654) q[3];
sx q[3];
rz(-2.0709166) q[3];
sx q[3];
rz(-1.3487032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.9434628) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(0.23588022) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772298) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(1.4981131) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(2.6565068) q[2];
sx q[2];
rz(-0.1242287) q[2];
sx q[2];
rz(2.6835174) q[2];
rz(-2.2723972) q[3];
sx q[3];
rz(-1.8164608) q[3];
sx q[3];
rz(1.7432004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
