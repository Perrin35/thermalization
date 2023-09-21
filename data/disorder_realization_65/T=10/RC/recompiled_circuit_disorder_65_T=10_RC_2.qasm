OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(-1.0746497) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2225392) q[0];
sx q[0];
rz(-0.018964501) q[0];
sx q[0];
rz(-2.8595964) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9825748) q[2];
sx q[2];
rz(-2.3454925) q[2];
sx q[2];
rz(-0.12479347) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5831456) q[1];
sx q[1];
rz(-1.6513283) q[1];
sx q[1];
rz(-3.0196378) q[1];
x q[2];
rz(-2.4514276) q[3];
sx q[3];
rz(-2.4881722) q[3];
sx q[3];
rz(-1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.82912123) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.6281698) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(-0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.5011903) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.507507) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(-3.0188942) q[0];
rz(-pi) q[1];
rz(-1.00374) q[2];
sx q[2];
rz(-0.54076414) q[2];
sx q[2];
rz(2.8218249) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.048763927) q[1];
sx q[1];
rz(-1.7424912) q[1];
sx q[1];
rz(-2.462639) q[1];
x q[2];
rz(3.1355751) q[3];
sx q[3];
rz(-1.179751) q[3];
sx q[3];
rz(-1.0139549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(-1.1668011) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(1.0789543) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.6378145) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5119748) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(-2.2114121) q[0];
rz(-pi) q[1];
rz(0.90942803) q[2];
sx q[2];
rz(-2.4436908) q[2];
sx q[2];
rz(2.1991889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7392002) q[1];
sx q[1];
rz(-1.270073) q[1];
sx q[1];
rz(1.6154358) q[1];
x q[2];
rz(1.4012666) q[3];
sx q[3];
rz(-2.4781514) q[3];
sx q[3];
rz(-2.0142114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7109795) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-2.553489) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96823111) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(1.9925041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34654348) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(-1.4562777) q[0];
rz(-2.7096268) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(0.84749046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.541084) q[1];
sx q[1];
rz(-0.87227548) q[1];
sx q[1];
rz(0.38703106) q[1];
x q[2];
rz(-2.7690414) q[3];
sx q[3];
rz(-2.3636892) q[3];
sx q[3];
rz(2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(-1.8956005) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4887061) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(-1.9309689) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592504) q[0];
sx q[0];
rz(-1.3929402) q[0];
sx q[0];
rz(1.3772208) q[0];
rz(-pi) q[1];
rz(2.3166979) q[2];
sx q[2];
rz(-0.56606228) q[2];
sx q[2];
rz(-0.60264665) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9235232) q[1];
sx q[1];
rz(-1.7150208) q[1];
sx q[1];
rz(-1.5037392) q[1];
x q[2];
rz(-2.1518425) q[3];
sx q[3];
rz(-2.5269066) q[3];
sx q[3];
rz(0.97782545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9050682) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(-1.0926251) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(0.96907369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7965308) q[0];
sx q[0];
rz(-3.025265) q[0];
sx q[0];
rz(-1.2141114) q[0];
rz(3.063383) q[2];
sx q[2];
rz(-2.4259896) q[2];
sx q[2];
rz(-2.2323334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6094738) q[1];
sx q[1];
rz(-1.9874959) q[1];
sx q[1];
rz(-2.6344224) q[1];
x q[2];
rz(-3*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(-0.43760854) q[2];
rz(0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(1.3597885) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.4356027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3754417) q[0];
sx q[0];
rz(-1.0784082) q[0];
sx q[0];
rz(-0.96813162) q[0];
x q[1];
rz(1.5905459) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(-1.8576647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.095852921) q[1];
sx q[1];
rz(-2.6624655) q[1];
sx q[1];
rz(0.87289255) q[1];
rz(0.42922677) q[3];
sx q[3];
rz(-2.0373404) q[3];
sx q[3];
rz(-0.45688094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.7763604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6909957) q[0];
sx q[0];
rz(-2.9857675) q[0];
sx q[0];
rz(2.0777006) q[0];
rz(1.442996) q[2];
sx q[2];
rz(-1.5926077) q[2];
sx q[2];
rz(-2.4335361) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33278388) q[1];
sx q[1];
rz(-2.60751) q[1];
sx q[1];
rz(3.0384455) q[1];
x q[2];
rz(2.7435999) q[3];
sx q[3];
rz(-0.45106217) q[3];
sx q[3];
rz(-3.1068902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0662971) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-0.42262849) q[2];
rz(-0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(0.2376093) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(-0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(2.2176567) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63147488) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(-0.60879137) q[0];
x q[1];
rz(3.1191101) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(1.894001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.82211557) q[1];
sx q[1];
rz(-1.176782) q[1];
sx q[1];
rz(1.0072717) q[1];
rz(-pi) q[2];
rz(-0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(-0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.49736398) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(-0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-3.1034234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725806) q[0];
sx q[0];
rz(-0.60927143) q[0];
sx q[0];
rz(1.5241745) q[0];
rz(-3.0496809) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(-3.0222169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8074179) q[1];
sx q[1];
rz(-1.6125229) q[1];
sx q[1];
rz(-1.9445022) q[1];
x q[2];
rz(0.60572259) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(0.13525087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.9434628) q[2];
rz(1.0839869) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643628) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.4981131) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-1.6289564) q[2];
sx q[2];
rz(-1.6806316) q[2];
sx q[2];
rz(-3.1113839) q[2];
rz(1.2002767) q[3];
sx q[3];
rz(-2.4051718) q[3];
sx q[3];
rz(-2.6889599) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
