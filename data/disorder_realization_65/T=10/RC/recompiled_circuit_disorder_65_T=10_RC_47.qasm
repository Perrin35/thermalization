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
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.933691) q[0];
sx q[0];
rz(-1.5760734) q[0];
sx q[0];
rz(0.018215608) q[0];
x q[1];
rz(1.4104112) q[2];
sx q[2];
rz(-2.3540902) q[2];
sx q[2];
rz(-3.0410142) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5731116) q[1];
sx q[1];
rz(-2.9955578) q[1];
sx q[1];
rz(-0.58574974) q[1];
rz(-pi) q[2];
rz(2.4514276) q[3];
sx q[3];
rz(-0.65342045) q[3];
sx q[3];
rz(-1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82912123) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.6281698) q[2];
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
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728977) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5630961) q[1];
sx q[1];
rz(-1.5011903) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.507507) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(-0.12269845) q[0];
rz(-pi) q[1];
rz(-1.00374) q[2];
sx q[2];
rz(-0.54076414) q[2];
sx q[2];
rz(2.8218249) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7565631) q[1];
sx q[1];
rz(-2.2379413) q[1];
sx q[1];
rz(-1.3515616) q[1];
x q[2];
rz(-1.9618481) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(2.5824576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(-1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(2.2272026) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(1.0789543) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(1.6378145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5119748) q[0];
sx q[0];
rz(-0.10911988) q[0];
sx q[0];
rz(-0.93018053) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2321646) q[2];
sx q[2];
rz(-0.69790188) q[2];
sx q[2];
rz(-0.94240377) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.40239247) q[1];
sx q[1];
rz(-1.8715197) q[1];
sx q[1];
rz(1.5261569) q[1];
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
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(0.7437931) q[2];
rz(-2.8841694) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-0.58810365) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(-0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.9925041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34654348) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(1.6853149) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97082246) q[2];
sx q[2];
rz(-0.68471013) q[2];
sx q[2];
rz(3.0175356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.541084) q[1];
sx q[1];
rz(-0.87227548) q[1];
sx q[1];
rz(0.38703106) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7690414) q[3];
sx q[3];
rz(-2.3636892) q[3];
sx q[3];
rz(-2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2465308) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(2.4728298) q[2];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-1.724842) q[1];
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
rz(0.82489478) q[2];
sx q[2];
rz(-2.5755304) q[2];
sx q[2];
rz(2.538946) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7985178) q[1];
sx q[1];
rz(-1.5044364) q[1];
sx q[1];
rz(-2.9970478) q[1];
rz(2.103881) q[3];
sx q[3];
rz(-1.8928877) q[3];
sx q[3];
rz(-2.0562293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9050682) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(-1.0926251) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4051751) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(1.3177692) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(0.96907369) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3450619) q[0];
sx q[0];
rz(-0.11632761) q[0];
sx q[0];
rz(1.9274812) q[0];
rz(-pi) q[1];
rz(0.078209608) q[2];
sx q[2];
rz(-0.71560301) q[2];
sx q[2];
rz(-2.2323334) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.25989306) q[1];
sx q[1];
rz(-1.1105781) q[1];
sx q[1];
rz(-1.1020745) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0958517) q[3];
sx q[3];
rz(-1.1610371) q[3];
sx q[3];
rz(0.72011442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(-0.43760854) q[2];
rz(-0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(-1.4356027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.5905459) q[2];
sx q[2];
rz(-1.5891468) q[2];
sx q[2];
rz(1.2839279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85305271) q[1];
sx q[1];
rz(-1.9318046) q[1];
sx q[1];
rz(2.8193874) q[1];
rz(-pi) q[2];
rz(2.0766047) q[3];
sx q[3];
rz(-1.1899663) q[3];
sx q[3];
rz(-1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4009565) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-0.1329578) q[0];
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
rz(0.62194659) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(1.7072862) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7403142) q[2];
sx q[2];
rz(-0.12963824) q[2];
sx q[2];
rz(0.6946176) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8088088) q[1];
sx q[1];
rz(-0.53408264) q[1];
sx q[1];
rz(-0.10314718) q[1];
rz(-pi) q[2];
rz(-1.7563617) q[3];
sx q[3];
rz(-1.1572596) q[3];
sx q[3];
rz(-0.47172771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(0.95831174) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(2.9039834) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-0.31059206) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-0.92393595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9792084) q[0];
sx q[0];
rz(-2.2109593) q[0];
sx q[0];
rz(1.0248653) q[0];
x q[1];
rz(-2.6529979) q[2];
sx q[2];
rz(-1.5813507) q[2];
sx q[2];
rz(0.34305629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1549346) q[1];
sx q[1];
rz(-1.0549874) q[1];
sx q[1];
rz(0.45706473) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(-2.2480272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.3941992) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7681463) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(-0.038169233) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.105135) q[0];
sx q[0];
rz(-1.5974701) q[0];
sx q[0];
rz(0.96203502) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3456887) q[2];
sx q[2];
rz(-2.7499866) q[2];
sx q[2];
rz(0.12347808) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.888615) q[1];
sx q[1];
rz(-1.1974317) q[1];
sx q[1];
rz(-3.0967767) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60572259) q[3];
sx q[3];
rz(-2.1401792) q[3];
sx q[3];
rz(-3.0063418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
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
rz(-pi/2) q[0];
x q[2];
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
rz(2.6565068) q[2];
sx q[2];
rz(-0.1242287) q[2];
sx q[2];
rz(2.6835174) q[2];
rz(-0.86919541) q[3];
sx q[3];
rz(-1.3251318) q[3];
sx q[3];
rz(-1.3983923) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
