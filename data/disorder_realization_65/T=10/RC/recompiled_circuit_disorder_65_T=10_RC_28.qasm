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
rz(-2.2812023) q[1];
sx q[1];
rz(5.5881349) q[1];
sx q[1];
rz(10.499428) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91905347) q[0];
sx q[0];
rz(-0.018964501) q[0];
sx q[0];
rz(2.8595964) q[0];
rz(2.3518402) q[2];
sx q[2];
rz(-1.4573922) q[2];
sx q[2];
rz(-1.5838768) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5731116) q[1];
sx q[1];
rz(-0.1460349) q[1];
sx q[1];
rz(-0.58574974) q[1];
x q[2];
rz(-0.69016506) q[3];
sx q[3];
rz(-2.4881722) q[3];
sx q[3];
rz(1.1284459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.6281698) q[2];
rz(1.9681905) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(-2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5728977) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(1.9182385) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(-1.6404023) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6340856) q[0];
sx q[0];
rz(-1.5956143) q[0];
sx q[0];
rz(3.0188942) q[0];
rz(-2.03962) q[2];
sx q[2];
rz(-1.2906244) q[2];
sx q[2];
rz(2.3902992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4109107) q[1];
sx q[1];
rz(-2.4446179) q[1];
sx q[1];
rz(-0.26941401) q[1];
rz(1.5853911) q[3];
sx q[3];
rz(-0.39108927) q[3];
sx q[3];
rz(-0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(-1.5141053) q[2];
rz(-1.1668011) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-2.0626383) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(-1.6378145) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5119748) q[0];
sx q[0];
rz(-3.0324728) q[0];
sx q[0];
rz(-0.93018053) q[0];
rz(-2.2321646) q[2];
sx q[2];
rz(-2.4436908) q[2];
sx q[2];
rz(2.1991889) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7392002) q[1];
sx q[1];
rz(-1.8715197) q[1];
sx q[1];
rz(-1.6154358) q[1];
rz(-pi) q[2];
rz(0.13111968) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(-2.2281447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-0.58810365) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96823111) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3455032) q[0];
sx q[0];
rz(-0.26608276) q[0];
sx q[0];
rz(0.43568133) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1707702) q[2];
sx q[2];
rz(-0.68471013) q[2];
sx q[2];
rz(0.1240571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1059873) q[1];
sx q[1];
rz(-0.78250256) q[1];
sx q[1];
rz(1.1483907) q[1];
x q[2];
rz(1.2265008) q[3];
sx q[3];
rz(-2.2831884) q[3];
sx q[3];
rz(2.815849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(0.66876283) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(-1.0523798) q[0];
rz(-1.9309689) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8540871) q[0];
sx q[0];
rz(-2.8794718) q[0];
sx q[0];
rz(0.8192807) q[0];
rz(-pi) q[1];
rz(0.82489478) q[2];
sx q[2];
rz(-0.56606228) q[2];
sx q[2];
rz(-2.538946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.21806949) q[1];
sx q[1];
rz(-1.4265718) q[1];
sx q[1];
rz(-1.5037392) q[1];
rz(-pi) q[2];
rz(2.103881) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(-1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.3177692) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9861525) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(3.100813) q[0];
rz(1.6386119) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(-2.3358047) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5507257) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(0.73901891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0958517) q[3];
sx q[3];
rz(-1.1610371) q[3];
sx q[3];
rz(2.4214782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(0.058241025) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83201927) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(-1.3597885) q[0];
rz(1.4490022) q[1];
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
rz(0.3754417) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(-2.173461) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5510467) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(-1.2839279) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0457397) q[1];
sx q[1];
rz(-0.47912712) q[1];
sx q[1];
rz(2.2687001) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7123659) q[3];
sx q[3];
rz(-2.0373404) q[3];
sx q[3];
rz(-2.6847117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(0.17399542) q[2];
rz(-2.1334355) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85132861) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(0.48775396) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62194659) q[0];
sx q[0];
rz(-1.4953817) q[0];
sx q[0];
rz(1.4343065) q[0];
x q[1];
rz(1.7403142) q[2];
sx q[2];
rz(-0.12963824) q[2];
sx q[2];
rz(2.4469751) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9284968) q[1];
sx q[1];
rz(-1.0398541) q[1];
sx q[1];
rz(1.6316158) q[1];
x q[2];
rz(-0.39799277) q[3];
sx q[3];
rz(-0.45106217) q[3];
sx q[3];
rz(0.034702452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(2.1832809) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(0.92393595) q[1];
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
rz(3.1191101) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(-1.2475916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9388705) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(0.90941456) q[1];
x q[2];
rz(2.4531035) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(-1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(0.53846255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(-1.753153) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(0.038169233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5725806) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(-1.5241745) q[0];
rz(-pi) q[1];
rz(1.3456887) q[2];
sx q[2];
rz(-0.39160608) q[2];
sx q[2];
rz(3.0181146) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0109509) q[1];
sx q[1];
rz(-0.37591939) q[1];
sx q[1];
rz(1.6846659) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.1756211) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(-1.1981298) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
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
rz(1.5126363) q[2];
sx q[2];
rz(-1.6806316) q[2];
sx q[2];
rz(-3.1113839) q[2];
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
