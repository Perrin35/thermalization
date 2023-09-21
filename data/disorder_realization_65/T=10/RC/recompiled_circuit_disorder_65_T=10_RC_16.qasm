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
rz(2.066943) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.933691) q[0];
sx q[0];
rz(-1.5760734) q[0];
sx q[0];
rz(0.018215608) q[0];
rz(-0.78975241) q[2];
sx q[2];
rz(-1.6842004) q[2];
sx q[2];
rz(-1.5577158) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97779146) q[1];
sx q[1];
rz(-1.692354) q[1];
sx q[1];
rz(-1.6519283) q[1];
rz(-pi) q[2];
rz(0.69016506) q[3];
sx q[3];
rz(-2.4881722) q[3];
sx q[3];
rz(-1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.5134229) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(-0.2127969) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(1.6404023) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0813635) q[0];
sx q[0];
rz(-1.4481359) q[0];
sx q[0];
rz(-1.5457904) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.00374) q[2];
sx q[2];
rz(-2.6008285) q[2];
sx q[2];
rz(-2.8218249) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3850296) q[1];
sx q[1];
rz(-2.2379413) q[1];
sx q[1];
rz(1.7900311) q[1];
rz(-pi) q[2];
rz(-1.1797446) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(-2.5824576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-2.2272026) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(-2.0626383) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(-1.5037781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62961783) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(2.2114121) q[0];
x q[1];
rz(0.9861228) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(-0.090678064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1551731) q[1];
sx q[1];
rz(-1.6134312) q[1];
sx q[1];
rz(-0.30100545) q[1];
rz(-pi) q[2];
x q[2];
rz(1.740326) q[3];
sx q[3];
rz(-0.66344122) q[3];
sx q[3];
rz(1.1273813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(2.553489) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(-2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.1490885) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950492) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(1.6853149) q[0];
rz(-pi) q[1];
rz(-0.97781424) q[2];
sx q[2];
rz(-1.2056418) q[2];
sx q[2];
rz(0.95945537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2266352) q[1];
sx q[1];
rz(-1.8640222) q[1];
sx q[1];
rz(-2.3073767) q[1];
rz(-0.37255128) q[3];
sx q[3];
rz(-0.77790341) q[3];
sx q[3];
rz(-0.82749623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2465308) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(-0.66876283) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4887061) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(0.88422424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28750557) q[0];
sx q[0];
rz(-2.8794718) q[0];
sx q[0];
rz(-2.322312) q[0];
x q[1];
rz(2.0074559) q[2];
sx q[2];
rz(-1.9433125) q[2];
sx q[2];
rz(2.8357752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7985178) q[1];
sx q[1];
rz(-1.5044364) q[1];
sx q[1];
rz(0.14454486) q[1];
x q[2];
rz(0.98975011) q[3];
sx q[3];
rz(-2.5269066) q[3];
sx q[3];
rz(-2.1637672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23652442) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(-2.0489676) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(-2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-0.96907369) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1554402) q[0];
sx q[0];
rz(-1.4618205) q[0];
sx q[0];
rz(-0.040779671) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71408748) q[2];
sx q[2];
rz(-1.519515) q[2];
sx q[2];
rz(2.4209792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8816996) q[1];
sx q[1];
rz(-2.0310146) q[1];
sx q[1];
rz(-2.0395181) q[1];
rz(-pi) q[2];
x q[2];
rz(8*pi/11) q[3];
sx q[3];
rz(-2.487605) q[3];
sx q[3];
rz(1.6884782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51450729) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(-0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3095734) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(-1.4490022) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.4356027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5098269) q[0];
sx q[0];
rz(-1.0477715) q[0];
sx q[0];
rz(2.564389) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82201634) q[2];
sx q[2];
rz(-3.1146345) q[2];
sx q[2];
rz(-2.6798623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0457397) q[1];
sx q[1];
rz(-2.6624655) q[1];
sx q[1];
rz(2.2687001) q[1];
x q[2];
rz(0.88019754) q[3];
sx q[3];
rz(-2.5185563) q[3];
sx q[3];
rz(0.33720371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4009565) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(0.17399542) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85132861) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(-2.6538387) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4505969) q[0];
sx q[0];
rz(-0.15582514) q[0];
sx q[0];
rz(1.063892) q[0];
rz(-pi) q[1];
rz(-1.4012785) q[2];
sx q[2];
rz(-0.12963824) q[2];
sx q[2];
rz(2.4469751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8088088) q[1];
sx q[1];
rz(-2.60751) q[1];
sx q[1];
rz(-0.10314718) q[1];
rz(-pi) q[2];
rz(-0.39799277) q[3];
sx q[3];
rz(-2.6905305) q[3];
sx q[3];
rz(3.1068902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(-0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(2.2176567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16238427) q[0];
sx q[0];
rz(-0.93063336) q[0];
sx q[0];
rz(-2.1167273) q[0];
rz(2.6529979) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(-0.34305629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98665806) q[1];
sx q[1];
rz(-2.0866052) q[1];
sx q[1];
rz(2.6845279) q[1];
rz(-pi) q[2];
rz(0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.3941992) q[2];
rz(-1.9723069) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.753153) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-3.1034234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.569012) q[0];
sx q[0];
rz(-0.60927143) q[0];
sx q[0];
rz(1.6174181) q[0];
rz(-pi) q[1];
rz(-1.3456887) q[2];
sx q[2];
rz(-0.39160608) q[2];
sx q[2];
rz(-3.0181146) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0109509) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(-1.4569267) q[1];
rz(-pi) q[2];
rz(-2.2977703) q[3];
sx q[3];
rz(-2.3355964) q[3];
sx q[3];
rz(-2.3674915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(-1.0839869) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.1643628) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(1.6434796) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(-3.0315728) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
rz(2.2723972) q[3];
sx q[3];
rz(-1.3251318) q[3];
sx q[3];
rz(-1.3983923) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
