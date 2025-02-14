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
rz(1.3801112) q[0];
sx q[0];
rz(-1.5910281) q[0];
sx q[0];
rz(-0.38183364) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.763321) q[0];
sx q[0];
rz(-2.1594783) q[0];
sx q[0];
rz(2.4886542) q[0];
x q[1];
rz(-1.5906232) q[2];
sx q[2];
rz(-1.2006491) q[2];
sx q[2];
rz(1.4408979) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0142483) q[1];
sx q[1];
rz(-0.44522983) q[1];
sx q[1];
rz(2.8365652) q[1];
x q[2];
rz(1.0409058) q[3];
sx q[3];
rz(-1.232551) q[3];
sx q[3];
rz(1.6774981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(1.5748242) q[2];
rz(-1.1664248) q[3];
sx q[3];
rz(-1.0650485) q[3];
sx q[3];
rz(-0.69275698) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-0.80802149) q[0];
sx q[0];
rz(-1.1760733) q[0];
rz(0.067497079) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(2.8228021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1677396) q[0];
sx q[0];
rz(-1.8813846) q[0];
sx q[0];
rz(0.74790252) q[0];
rz(2.6574357) q[2];
sx q[2];
rz(-2.6645497) q[2];
sx q[2];
rz(-2.8271528) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31736703) q[1];
sx q[1];
rz(-1.1928802) q[1];
sx q[1];
rz(2.1721858) q[1];
rz(-pi) q[2];
x q[2];
rz(0.12572581) q[3];
sx q[3];
rz(-1.2866402) q[3];
sx q[3];
rz(1.8234258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0789644) q[2];
sx q[2];
rz(-2.2603409) q[2];
sx q[2];
rz(0.47743615) q[2];
rz(-1.3580953) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(3.131598) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7497712) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(2.0023477) q[0];
rz(-0.40621743) q[1];
sx q[1];
rz(-1.258305) q[1];
sx q[1];
rz(-0.67799813) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11230532) q[0];
sx q[0];
rz(-0.59354085) q[0];
sx q[0];
rz(-1.4112006) q[0];
rz(-0.35786104) q[2];
sx q[2];
rz(-2.1659019) q[2];
sx q[2];
rz(0.33085631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78732508) q[1];
sx q[1];
rz(-2.4128822) q[1];
sx q[1];
rz(2.4142152) q[1];
x q[2];
rz(1.0150681) q[3];
sx q[3];
rz(-0.58150269) q[3];
sx q[3];
rz(2.8648928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4262126) q[2];
sx q[2];
rz(-1.5072301) q[2];
sx q[2];
rz(1.0769843) q[2];
rz(-1.2601323) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(-2.8857359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28427163) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(0.76048365) q[0];
rz(2.7898232) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(2.9913091) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12211299) q[0];
sx q[0];
rz(-1.5727336) q[0];
sx q[0];
rz(-1.5674855) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9530832) q[2];
sx q[2];
rz(-1.4103544) q[2];
sx q[2];
rz(-0.40942243) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1898415) q[1];
sx q[1];
rz(-0.69165666) q[1];
sx q[1];
rz(2.1053949) q[1];
rz(0.49154776) q[3];
sx q[3];
rz(-0.50954362) q[3];
sx q[3];
rz(0.59922892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1111697) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(0.58491582) q[2];
rz(3.0758744) q[3];
sx q[3];
rz(-1.1394371) q[3];
sx q[3];
rz(-1.0544624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86972648) q[0];
sx q[0];
rz(-1.2249185) q[0];
sx q[0];
rz(-2.4140893) q[0];
rz(-1.8456521) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(-2.4268699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6520157) q[0];
sx q[0];
rz(-1.6418818) q[0];
sx q[0];
rz(0.22146225) q[0];
x q[1];
rz(-2.7834846) q[2];
sx q[2];
rz(-0.3198959) q[2];
sx q[2];
rz(1.7826155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89940577) q[1];
sx q[1];
rz(-1.9493628) q[1];
sx q[1];
rz(2.942571) q[1];
rz(-pi) q[2];
rz(-2.1398224) q[3];
sx q[3];
rz(-2.4845893) q[3];
sx q[3];
rz(-3.0083619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99594816) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(1.145251) q[2];
rz(-0.20259419) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(-0.32874671) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8458493) q[0];
sx q[0];
rz(-2.8491617) q[0];
sx q[0];
rz(-0.16786815) q[0];
rz(2.5766418) q[1];
sx q[1];
rz(-2.0940557) q[1];
sx q[1];
rz(1.8126743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8362373) q[0];
sx q[0];
rz(-2.1026045) q[0];
sx q[0];
rz(2.1270069) q[0];
x q[1];
rz(-1.6199048) q[2];
sx q[2];
rz(-1.430776) q[2];
sx q[2];
rz(0.55632178) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4327287) q[1];
sx q[1];
rz(-1.8058746) q[1];
sx q[1];
rz(2.8793596) q[1];
rz(-pi) q[2];
rz(-1.331948) q[3];
sx q[3];
rz(-1.788056) q[3];
sx q[3];
rz(-0.9700194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7557175) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(2.7937549) q[2];
rz(-2.8268585) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(-2.0033964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9921853) q[0];
sx q[0];
rz(-1.6353761) q[0];
sx q[0];
rz(0.95112479) q[0];
rz(-2.8490207) q[1];
sx q[1];
rz(-0.89076275) q[1];
sx q[1];
rz(2.1975885) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68444618) q[0];
sx q[0];
rz(-1.5862042) q[0];
sx q[0];
rz(1.2627361) q[0];
x q[1];
rz(-0.79156178) q[2];
sx q[2];
rz(-2.5063609) q[2];
sx q[2];
rz(2.5286365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85691014) q[1];
sx q[1];
rz(-1.7439242) q[1];
sx q[1];
rz(-0.022927479) q[1];
rz(-pi) q[2];
rz(0.24686019) q[3];
sx q[3];
rz(-1.180397) q[3];
sx q[3];
rz(2.2228783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(1.4415119) q[2];
rz(-3.1169543) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(2.1520481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330419) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(-3.0035875) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.6777918) q[1];
sx q[1];
rz(-2.2053351) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9247583) q[0];
sx q[0];
rz(-2.6568597) q[0];
sx q[0];
rz(-2.4527571) q[0];
x q[1];
rz(0.9844043) q[2];
sx q[2];
rz(-0.83293646) q[2];
sx q[2];
rz(3.0393103) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4605751) q[1];
sx q[1];
rz(-1.5854757) q[1];
sx q[1];
rz(-1.3880066) q[1];
rz(-pi) q[2];
rz(-2.4560087) q[3];
sx q[3];
rz(-1.3635677) q[3];
sx q[3];
rz(-2.384546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.998698) q[2];
sx q[2];
rz(-1.2724718) q[2];
sx q[2];
rz(-1.3188837) q[2];
rz(0.12254347) q[3];
sx q[3];
rz(-0.69869852) q[3];
sx q[3];
rz(-2.7114649) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1634624) q[0];
sx q[0];
rz(-2.6977111) q[0];
sx q[0];
rz(2.7769856) q[0];
rz(1.7568024) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(2.2273831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.171998) q[0];
sx q[0];
rz(-1.1437609) q[0];
sx q[0];
rz(-2.9946505) q[0];
x q[1];
rz(-1.6169102) q[2];
sx q[2];
rz(-0.83674016) q[2];
sx q[2];
rz(3.0457599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5678775) q[1];
sx q[1];
rz(-0.93522954) q[1];
sx q[1];
rz(-1.741841) q[1];
x q[2];
rz(1.7201603) q[3];
sx q[3];
rz(-1.0348667) q[3];
sx q[3];
rz(-1.5045297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.30283516) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(1.5128822) q[2];
rz(2.5527939) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47422472) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(-2.8571416) q[0];
rz(-2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(0.42025748) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66061879) q[0];
sx q[0];
rz(-1.7848564) q[0];
sx q[0];
rz(-1.1461036) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3510397) q[2];
sx q[2];
rz(-1.9091064) q[2];
sx q[2];
rz(-1.4684791) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6984562) q[1];
sx q[1];
rz(-1.4520169) q[1];
sx q[1];
rz(-0.08929792) q[1];
rz(2.5758366) q[3];
sx q[3];
rz(-1.8744191) q[3];
sx q[3];
rz(-2.2591801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9922716) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(0.9672271) q[2];
rz(1.017717) q[3];
sx q[3];
rz(-1.8571721) q[3];
sx q[3];
rz(0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.16348542) q[0];
sx q[0];
rz(-1.0657943) q[0];
sx q[0];
rz(-0.97070538) q[0];
rz(3.0178487) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(-0.32044784) q[2];
sx q[2];
rz(-1.6978227) q[2];
sx q[2];
rz(2.6200079) q[2];
rz(-1.5354034) q[3];
sx q[3];
rz(-1.2589068) q[3];
sx q[3];
rz(1.3990228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
