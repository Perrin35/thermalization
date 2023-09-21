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
rz(1.6488099) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(5.5881349) q[1];
sx q[1];
rz(10.499428) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91905347) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(0.28199621) q[0];
rz(-pi) q[1];
rz(-2.9825748) q[2];
sx q[2];
rz(-0.79610014) q[2];
sx q[2];
rz(0.12479347) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5831456) q[1];
sx q[1];
rz(-1.4902643) q[1];
sx q[1];
rz(-3.0196378) q[1];
rz(-pi) q[2];
rz(-0.69016506) q[3];
sx q[3];
rz(-0.65342045) q[3];
sx q[3];
rz(2.0131468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82912123) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.2233541) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(1.6404023) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26186865) q[0];
sx q[0];
rz(-3.0164218) q[0];
sx q[0];
rz(-0.20010389) q[0];
rz(-2.1378527) q[2];
sx q[2];
rz(-2.6008285) q[2];
sx q[2];
rz(2.8218249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0928287) q[1];
sx q[1];
rz(-1.3991014) q[1];
sx q[1];
rz(-2.462639) q[1];
rz(1.5853911) q[3];
sx q[3];
rz(-2.7505034) q[3];
sx q[3];
rz(0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(1.0789543) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.5037781) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30341879) q[0];
sx q[0];
rz(-1.6359328) q[0];
sx q[0];
rz(1.6584048) q[0];
x q[1];
rz(-0.47567993) q[2];
sx q[2];
rz(-2.1026346) q[2];
sx q[2];
rz(-1.4059517) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5895264) q[1];
sx q[1];
rz(-2.8376736) q[1];
sx q[1];
rz(0.1429096) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13111968) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-2.6617699) q[2];
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
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96823111) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(-2.5412718) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(-1.9925041) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968397) q[0];
sx q[0];
rz(-1.6819994) q[0];
sx q[0];
rz(-2.8993594) q[0];
rz(-2.1707702) q[2];
sx q[2];
rz(-2.4568825) q[2];
sx q[2];
rz(3.0175356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60050868) q[1];
sx q[1];
rz(-2.2693172) q[1];
sx q[1];
rz(2.7545616) q[1];
x q[2];
rz(2.7690414) q[3];
sx q[3];
rz(-0.77790341) q[3];
sx q[3];
rz(2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-1.0523798) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(0.88422424) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8540871) q[0];
sx q[0];
rz(-0.26212087) q[0];
sx q[0];
rz(-0.8192807) q[0];
rz(-pi) q[1];
rz(2.7344633) q[2];
sx q[2];
rz(-1.1659157) q[2];
sx q[2];
rz(-1.4332353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.65518889) q[1];
sx q[1];
rz(-2.9826394) q[1];
sx q[1];
rz(-0.43227355) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.103881) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(-2.0562293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1554402) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(-0.040779671) q[0];
x q[1];
rz(-1.6386119) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(2.3358047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5507257) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(0.73901891) q[1];
rz(8*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(-1.6884782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51450729) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.3597885) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(-1.4356027) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5098269) q[0];
sx q[0];
rz(-2.0938211) q[0];
sx q[0];
rz(-2.564389) q[0];
x q[1];
rz(0.82201634) q[2];
sx q[2];
rz(-3.1146345) q[2];
sx q[2];
rz(-2.6798623) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3064733) q[1];
sx q[1];
rz(-1.2700348) q[1];
sx q[1];
rz(1.1919828) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0766047) q[3];
sx q[3];
rz(-1.1899663) q[3];
sx q[3];
rz(1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4009565) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(-2.1334355) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(3.0086349) q[0];
rz(0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.7763604) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4505969) q[0];
sx q[0];
rz(-0.15582514) q[0];
sx q[0];
rz(-2.0777006) q[0];
rz(-pi) q[1];
x q[1];
rz(3.119602) q[2];
sx q[2];
rz(-1.4430265) q[2];
sx q[2];
rz(2.2760504) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9284968) q[1];
sx q[1];
rz(-2.1017385) q[1];
sx q[1];
rz(1.6316158) q[1];
rz(1.385231) q[3];
sx q[3];
rz(-1.1572596) q[3];
sx q[3];
rz(2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-0.42262849) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3850708) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(0.71682741) q[0];
x q[1];
rz(-1.5588435) q[2];
sx q[2];
rz(-2.0593615) q[2];
sx q[2];
rz(1.2221297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.98665806) q[1];
sx q[1];
rz(-1.0549874) q[1];
sx q[1];
rz(2.6845279) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-2.101427) q[3];
sx q[3];
rz(0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(-0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(0.038169233) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036457688) q[0];
sx q[0];
rz(-1.5441226) q[0];
sx q[0];
rz(-0.96203502) q[0];
rz(-0.091911749) q[2];
sx q[2];
rz(-1.9520063) q[2];
sx q[2];
rz(0.11937571) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3341748) q[1];
sx q[1];
rz(-1.6125229) q[1];
sx q[1];
rz(1.1970904) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2323654) q[3];
sx q[3];
rz(-2.0709166) q[3];
sx q[3];
rz(-1.7928894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(-1.4981131) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(0.11001982) q[2];
sx q[2];
rz(-1.6286055) q[2];
sx q[2];
rz(1.5946228) q[2];
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
