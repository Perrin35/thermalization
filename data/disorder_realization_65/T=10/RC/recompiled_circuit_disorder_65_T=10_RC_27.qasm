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
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(-1.0746497) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2225392) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(-0.28199621) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9825748) q[2];
sx q[2];
rz(-0.79610014) q[2];
sx q[2];
rz(3.0167992) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5731116) q[1];
sx q[1];
rz(-2.9955578) q[1];
sx q[1];
rz(-2.5558429) q[1];
rz(-pi) q[2];
rz(-0.69016506) q[3];
sx q[3];
rz(-0.65342045) q[3];
sx q[3];
rz(-1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.6281698) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(-1.6404023) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0813635) q[0];
sx q[0];
rz(-1.4481359) q[0];
sx q[0];
rz(-1.5958022) q[0];
rz(-pi) q[1];
rz(-2.03962) q[2];
sx q[2];
rz(-1.8509682) q[2];
sx q[2];
rz(-2.3902992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0928287) q[1];
sx q[1];
rz(-1.3991014) q[1];
sx q[1];
rz(-0.6789536) q[1];
x q[2];
rz(-1.9618481) q[3];
sx q[3];
rz(-1.5763596) q[3];
sx q[3];
rz(-2.5824576) q[3];
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
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-2.0626383) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.6378145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2730946) q[0];
sx q[0];
rz(-1.4833741) q[0];
sx q[0];
rz(3.0762061) q[0];
rz(-pi) q[1];
rz(-2.1554699) q[2];
sx q[2];
rz(-1.976527) q[2];
sx q[2];
rz(-3.0509146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1551731) q[1];
sx q[1];
rz(-1.5281614) q[1];
sx q[1];
rz(0.30100545) q[1];
x q[2];
rz(-0.13111968) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(2.2281447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(2.3977996) q[2];
rz(-0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(1.051735) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.9925041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34654348) q[0];
sx q[0];
rz(-1.8115037) q[0];
sx q[0];
rz(-1.4562777) q[0];
rz(-2.1637784) q[2];
sx q[2];
rz(-1.2056418) q[2];
sx q[2];
rz(-0.95945537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1059873) q[1];
sx q[1];
rz(-2.3590901) q[1];
sx q[1];
rz(1.1483907) q[1];
rz(-pi) q[2];
rz(-1.9150919) q[3];
sx q[3];
rz(-2.2831884) q[3];
sx q[3];
rz(2.815849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8540871) q[0];
sx q[0];
rz(-2.8794718) q[0];
sx q[0];
rz(-0.8192807) q[0];
x q[1];
rz(0.82489478) q[2];
sx q[2];
rz(-2.5755304) q[2];
sx q[2];
rz(-0.60264665) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4864038) q[1];
sx q[1];
rz(-0.15895325) q[1];
sx q[1];
rz(-0.43227355) q[1];
rz(-1.0377117) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(-1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23652442) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(-1.3177692) q[0];
rz(-0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(2.172519) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7965308) q[0];
sx q[0];
rz(-0.11632761) q[0];
sx q[0];
rz(1.9274812) q[0];
x q[1];
rz(1.6386119) q[2];
sx q[2];
rz(-0.8578476) q[2];
sx q[2];
rz(-0.80578795) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8816996) q[1];
sx q[1];
rz(-1.1105781) q[1];
sx q[1];
rz(-2.0395181) q[1];
rz(-pi) q[2];
rz(-0.46521516) q[3];
sx q[3];
rz(-1.0930659) q[3];
sx q[3];
rz(0.623869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51450729) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(-0.43760854) q[2];
rz(-3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095734) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(1.3597885) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.4356027) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5933701) q[0];
sx q[0];
rz(-2.3832294) q[0];
sx q[0];
rz(-2.3286657) q[0];
x q[1];
rz(-0.82201634) q[2];
sx q[2];
rz(-0.026958131) q[2];
sx q[2];
rz(0.4617304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0457397) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(-0.87289255) q[1];
x q[2];
rz(-2.0766047) q[3];
sx q[3];
rz(-1.9516264) q[3];
sx q[3];
rz(-1.8246458) q[3];
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
rz(0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85132861) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.7763604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62194659) q[0];
sx q[0];
rz(-1.4953817) q[0];
sx q[0];
rz(1.4343065) q[0];
rz(1.6985967) q[2];
sx q[2];
rz(-1.5489849) q[2];
sx q[2];
rz(-2.4335361) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3268765) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(-2.609842) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7216464) q[3];
sx q[3];
rz(-1.4010324) q[3];
sx q[3];
rz(-1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-0.42262849) q[2];
rz(-0.95831174) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(-0.2376093) q[3];
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
rz(3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(-0.31059206) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16238427) q[0];
sx q[0];
rz(-2.2109593) q[0];
sx q[0];
rz(1.0248653) q[0];
x q[1];
rz(2.6529979) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(-0.34305629) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20272217) q[1];
sx q[1];
rz(-0.67516967) q[1];
sx q[1];
rz(-0.90941456) q[1];
rz(-pi) q[2];
rz(-2.4531035) q[3];
sx q[3];
rz(-0.8430891) q[3];
sx q[3];
rz(-1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6258433) q[0];
sx q[0];
rz(-0.96228296) q[0];
sx q[0];
rz(3.1090816) q[0];
rz(-0.091911749) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(3.0222169) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8074179) q[1];
sx q[1];
rz(-1.6125229) q[1];
sx q[1];
rz(-1.1970904) q[1];
rz(-0.60572259) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(3.0063418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(-1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(-1.6434796) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(2.6565068) q[2];
sx q[2];
rz(-0.1242287) q[2];
sx q[2];
rz(2.6835174) q[2];
rz(0.31717832) q[3];
sx q[3];
rz(-2.2472897) q[3];
sx q[3];
rz(-0.030285611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
