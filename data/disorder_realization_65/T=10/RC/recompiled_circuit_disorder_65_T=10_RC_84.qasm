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
rz(2.0895045) q[0];
sx q[0];
rz(10.917561) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(2.066943) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91905347) q[0];
sx q[0];
rz(-0.018964501) q[0];
sx q[0];
rz(2.8595964) q[0];
rz(0.78975241) q[2];
sx q[2];
rz(-1.4573922) q[2];
sx q[2];
rz(1.5838768) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1638012) q[1];
sx q[1];
rz(-1.692354) q[1];
sx q[1];
rz(1.4896643) q[1];
rz(-pi) q[2];
rz(2.0243458) q[3];
sx q[3];
rz(-1.082886) q[3];
sx q[3];
rz(2.818056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82912123) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.5134229) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
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
rz(-1.5784966) q[1];
sx q[1];
rz(-1.6404023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6340856) q[0];
sx q[0];
rz(-1.5459783) q[0];
sx q[0];
rz(-3.0188942) q[0];
rz(-2.8295849) q[2];
sx q[2];
rz(-1.1216251) q[2];
sx q[2];
rz(0.95865196) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4109107) q[1];
sx q[1];
rz(-0.69697471) q[1];
sx q[1];
rz(-2.8721786) q[1];
rz(-pi) q[2];
rz(-1.5853911) q[3];
sx q[3];
rz(-0.39108927) q[3];
sx q[3];
rz(0.99816834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(1.5141053) q[2];
rz(1.1668011) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(0.91439009) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(-2.0626383) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.6378145) q[1];
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
rz(-pi) q[1];
rz(-0.90942803) q[2];
sx q[2];
rz(-2.4436908) q[2];
sx q[2];
rz(0.94240377) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5895264) q[1];
sx q[1];
rz(-2.8376736) q[1];
sx q[1];
rz(2.9986831) q[1];
rz(-pi) q[2];
rz(1.4012666) q[3];
sx q[3];
rz(-0.66344122) q[3];
sx q[3];
rz(-1.1273813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7109795) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(2.3977996) q[2];
rz(2.8841694) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96823111) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(-1.9925041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7950492) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(1.4562777) q[0];
x q[1];
rz(2.1637784) q[2];
sx q[2];
rz(-1.9359509) q[2];
sx q[2];
rz(2.1821373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.541084) q[1];
sx q[1];
rz(-2.2693172) q[1];
sx q[1];
rz(2.7545616) q[1];
x q[2];
rz(-0.74242384) q[3];
sx q[3];
rz(-1.8291049) q[3];
sx q[3];
rz(-1.0148259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-0.66876283) q[2];
rz(1.8956005) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(2.2573684) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592504) q[0];
sx q[0];
rz(-1.7486524) q[0];
sx q[0];
rz(1.7643719) q[0];
x q[1];
rz(-0.40712936) q[2];
sx q[2];
rz(-1.9756769) q[2];
sx q[2];
rz(-1.7083573) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7985178) q[1];
sx q[1];
rz(-1.6371562) q[1];
sx q[1];
rz(0.14454486) q[1];
rz(-pi) q[2];
rz(0.36966427) q[3];
sx q[3];
rz(-2.0737994) q[3];
sx q[3];
rz(0.30077416) q[3];
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
rz(2.0489676) q[2];
rz(-0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.3177692) q[0];
rz(2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-2.172519) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1554402) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(3.100813) q[0];
rz(0.078209608) q[2];
sx q[2];
rz(-0.71560301) q[2];
sx q[2];
rz(-2.2323334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.25989306) q[1];
sx q[1];
rz(-2.0310146) q[1];
sx q[1];
rz(-2.0395181) q[1];
rz(1.045741) q[3];
sx q[3];
rz(-1.9805555) q[3];
sx q[3];
rz(-0.72011442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51450729) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(-2.7039841) q[2];
rz(-3.0833516) q[3];
sx q[3];
rz(-0.85844675) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(-1.7818041) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(1.4356027) q[1];
rz(-pi/2) q[2];
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
rz(-pi) q[1];
rz(-1.5905459) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(-1.2839279) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3064733) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.9496099) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(-2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-2.984916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(-1.7763604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93850219) q[0];
sx q[0];
rz(-1.4346968) q[0];
sx q[0];
rz(0.076119856) q[0];
rz(-3.119602) q[2];
sx q[2];
rz(-1.6985661) q[2];
sx q[2];
rz(2.2760504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8088088) q[1];
sx q[1];
rz(-2.60751) q[1];
sx q[1];
rz(-0.10314718) q[1];
rz(-2.7216464) q[3];
sx q[3];
rz(-1.7405602) q[3];
sx q[3];
rz(1.9672293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-1.5035166) q[1];
sx q[1];
rz(0.92393595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9792084) q[0];
sx q[0];
rz(-0.93063336) q[0];
sx q[0];
rz(2.1167273) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1191101) q[2];
sx q[2];
rz(-2.6528931) q[2];
sx q[2];
rz(-1.2475916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3194771) q[1];
sx q[1];
rz(-1.9648106) q[1];
sx q[1];
rz(-2.134321) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95122533) q[3];
sx q[3];
rz(-0.95622548) q[3];
sx q[3];
rz(2.8545477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49736398) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(-1.3941992) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37344638) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-3.1034234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6258433) q[0];
sx q[0];
rz(-2.1793097) q[0];
sx q[0];
rz(-0.032511062) q[0];
x q[1];
rz(-0.091911749) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(-0.11937571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0109509) q[1];
sx q[1];
rz(-0.37591939) q[1];
sx q[1];
rz(1.6846659) q[1];
rz(-pi) q[2];
rz(-0.60572259) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(-0.13525087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(-1.9434628) q[2];
rz(2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9772298) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(-1.4981131) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(-1.5126363) q[2];
sx q[2];
rz(-1.460961) q[2];
sx q[2];
rz(0.030208781) q[2];
rz(-0.31717832) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];