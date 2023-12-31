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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63700919) q[0];
sx q[0];
rz(-1.5890117) q[0];
sx q[0];
rz(1.5760742) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4104112) q[2];
sx q[2];
rz(-2.3540902) q[2];
sx q[2];
rz(-0.10057848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5731116) q[1];
sx q[1];
rz(-0.1460349) q[1];
sx q[1];
rz(-2.5558429) q[1];
rz(-1.1172469) q[3];
sx q[3];
rz(-1.082886) q[3];
sx q[3];
rz(2.818056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.6281698) q[2];
rz(1.1734022) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(1.2233541) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.879724) q[0];
sx q[0];
rz(-3.0164218) q[0];
sx q[0];
rz(-2.9414888) q[0];
rz(2.03962) q[2];
sx q[2];
rz(-1.2906244) q[2];
sx q[2];
rz(0.75129347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7565631) q[1];
sx q[1];
rz(-0.90365138) q[1];
sx q[1];
rz(1.7900311) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9618481) q[3];
sx q[3];
rz(-1.5763596) q[3];
sx q[3];
rz(0.55913505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(-1.5141053) q[2];
rz(-1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062155) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(1.6378145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2730946) q[0];
sx q[0];
rz(-1.6582186) q[0];
sx q[0];
rz(-0.065386535) q[0];
x q[1];
rz(-0.90942803) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(-0.94240377) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1551731) q[1];
sx q[1];
rz(-1.5281614) q[1];
sx q[1];
rz(2.8405872) q[1];
rz(1.740326) q[3];
sx q[3];
rz(-0.66344122) q[3];
sx q[3];
rz(1.1273813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96823111) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(-0.60032088) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.9925041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950492) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(-1.4562777) q[0];
rz(-pi) q[1];
rz(0.97082246) q[2];
sx q[2];
rz(-2.4568825) q[2];
sx q[2];
rz(3.0175356) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.035605343) q[1];
sx q[1];
rz(-2.3590901) q[1];
sx q[1];
rz(1.993202) q[1];
rz(-pi) q[2];
rz(-1.2265008) q[3];
sx q[3];
rz(-2.2831884) q[3];
sx q[3];
rz(0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(1.8956005) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(-3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4887061) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(1.0523798) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-0.88422424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0852172) q[0];
sx q[0];
rz(-1.7612805) q[0];
sx q[0];
rz(0.18116829) q[0];
rz(-pi) q[1];
rz(-0.82489478) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21806949) q[1];
sx q[1];
rz(-1.7150208) q[1];
sx q[1];
rz(1.5037392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36966427) q[3];
sx q[3];
rz(-2.0737994) q[3];
sx q[3];
rz(-2.8408185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23652442) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(-1.0926251) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5613865) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(-1.6798621) q[0];
x q[1];
rz(2.4275052) q[2];
sx q[2];
rz(-1.6220777) q[2];
sx q[2];
rz(2.4209792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59086696) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(-0.73901891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6763775) q[3];
sx q[3];
rz(-1.0930659) q[3];
sx q[3];
rz(-2.5177237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(-3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095734) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(1.4490022) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(1.70599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6317658) q[0];
sx q[0];
rz(-2.0938211) q[0];
sx q[0];
rz(2.564389) q[0];
rz(-1.5905459) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(1.8576647) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85305271) q[1];
sx q[1];
rz(-1.9318046) q[1];
sx q[1];
rz(2.8193874) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0649879) q[3];
sx q[3];
rz(-1.1899663) q[3];
sx q[3];
rz(1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(3.0086349) q[0];
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
rz(2.5196461) q[0];
sx q[0];
rz(-1.4953817) q[0];
sx q[0];
rz(1.7072862) q[0];
rz(-pi) q[1];
rz(-3.119602) q[2];
sx q[2];
rz(-1.6985661) q[2];
sx q[2];
rz(-0.86554229) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3268765) q[1];
sx q[1];
rz(-1.5183581) q[1];
sx q[1];
rz(-2.609842) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7563617) q[3];
sx q[3];
rz(-1.1572596) q[3];
sx q[3];
rz(2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(-2.7189642) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(0.2376093) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(0.31059206) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(2.2176567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9792084) q[0];
sx q[0];
rz(-2.2109593) q[0];
sx q[0];
rz(2.1167273) q[0];
rz(3.1191101) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(1.894001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98665806) q[1];
sx q[1];
rz(-1.0549874) q[1];
sx q[1];
rz(-0.45706473) q[1];
rz(-pi) q[2];
rz(0.95122533) q[3];
sx q[3];
rz(-2.1853672) q[3];
sx q[3];
rz(-0.28704498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-0.038169233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036457688) q[0];
sx q[0];
rz(-1.5974701) q[0];
sx q[0];
rz(-0.96203502) q[0];
rz(-pi) q[1];
rz(-1.795904) q[2];
sx q[2];
rz(-0.39160608) q[2];
sx q[2];
rz(3.0181146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1306418) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(-1.4569267) q[1];
rz(-2.5358701) q[3];
sx q[3];
rz(-2.1401792) q[3];
sx q[3];
rz(-0.13525087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.9434628) q[2];
rz(-1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-0.23588022) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643628) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(-1.6434796) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-0.11001982) q[2];
sx q[2];
rz(-1.5129871) q[2];
sx q[2];
rz(-1.5469698) q[2];
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
