OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1239399) q[0];
sx q[0];
rz(-0.91683638) q[0];
sx q[0];
rz(2.7066948) q[0];
rz(-1.544156) q[1];
sx q[1];
rz(-0.51900744) q[1];
sx q[1];
rz(-2.5595698) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25044826) q[0];
sx q[0];
rz(-1.8316557) q[0];
sx q[0];
rz(2.9820739) q[0];
rz(-pi) q[1];
rz(-0.93272722) q[2];
sx q[2];
rz(-1.9550394) q[2];
sx q[2];
rz(-2.8926433) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7705045) q[1];
sx q[1];
rz(-2.6367715) q[1];
sx q[1];
rz(1.0434898) q[1];
rz(-2.1277027) q[3];
sx q[3];
rz(-0.32809533) q[3];
sx q[3];
rz(1.2893639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5883098) q[2];
sx q[2];
rz(-1.2574235) q[2];
sx q[2];
rz(2.8498939) q[2];
rz(2.6907673) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(-2.5341471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.70158231) q[0];
sx q[0];
rz(-2.1511183) q[0];
sx q[0];
rz(1.095358) q[0];
rz(1.9502684) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(2.2390168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13690925) q[0];
sx q[0];
rz(-0.97494128) q[0];
sx q[0];
rz(-2.3803902) q[0];
x q[1];
rz(2.8475855) q[2];
sx q[2];
rz(-2.2180811) q[2];
sx q[2];
rz(0.33366007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2486156) q[1];
sx q[1];
rz(-1.4505523) q[1];
sx q[1];
rz(-0.29386947) q[1];
x q[2];
rz(-0.91530494) q[3];
sx q[3];
rz(-0.73507959) q[3];
sx q[3];
rz(0.32441586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1878745) q[2];
sx q[2];
rz(-2.1672921) q[2];
sx q[2];
rz(-3.0226959) q[2];
rz(-0.64905727) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(-1.4547179) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.651799) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(-0.5157665) q[0];
rz(2.0491397) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(1.2752424) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3608801) q[0];
sx q[0];
rz(-0.045340625) q[0];
sx q[0];
rz(-0.41098292) q[0];
x q[1];
rz(2.5162773) q[2];
sx q[2];
rz(-0.38039243) q[2];
sx q[2];
rz(-1.8521295) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0787492) q[1];
sx q[1];
rz(-0.8991407) q[1];
sx q[1];
rz(0.25154227) q[1];
rz(-pi) q[2];
rz(0.79408349) q[3];
sx q[3];
rz(-1.6424254) q[3];
sx q[3];
rz(-0.41251999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.065993) q[2];
sx q[2];
rz(-1.8241355) q[2];
sx q[2];
rz(-1.1598178) q[2];
rz(2.6521111) q[3];
sx q[3];
rz(-1.5533841) q[3];
sx q[3];
rz(2.421853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8664261) q[0];
sx q[0];
rz(-0.32298276) q[0];
sx q[0];
rz(-2.6599356) q[0];
rz(0.9306759) q[1];
sx q[1];
rz(-1.296867) q[1];
sx q[1];
rz(-0.01893386) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1279432) q[0];
sx q[0];
rz(-0.054220323) q[0];
sx q[0];
rz(-1.9677866) q[0];
x q[1];
rz(1.5140947) q[2];
sx q[2];
rz(-0.66326521) q[2];
sx q[2];
rz(0.5723638) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1065797) q[1];
sx q[1];
rz(-2.0845319) q[1];
sx q[1];
rz(0.68961838) q[1];
x q[2];
rz(2.3778524) q[3];
sx q[3];
rz(-0.9943878) q[3];
sx q[3];
rz(1.1820205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21861741) q[2];
sx q[2];
rz(-0.3564035) q[2];
sx q[2];
rz(0.74529988) q[2];
rz(-0.37799147) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(-0.88855827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22076631) q[0];
sx q[0];
rz(-0.31196088) q[0];
sx q[0];
rz(-1.1908603) q[0];
rz(2.6530755) q[1];
sx q[1];
rz(-2.2256475) q[1];
sx q[1];
rz(2.3337505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69881781) q[0];
sx q[0];
rz(-1.4801868) q[0];
sx q[0];
rz(1.242799) q[0];
x q[1];
rz(-2.0338397) q[2];
sx q[2];
rz(-1.3658189) q[2];
sx q[2];
rz(0.092627138) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.399209) q[1];
sx q[1];
rz(-1.1157554) q[1];
sx q[1];
rz(1.4728404) q[1];
x q[2];
rz(-2.9414953) q[3];
sx q[3];
rz(-0.95733023) q[3];
sx q[3];
rz(2.4870949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0025582) q[2];
sx q[2];
rz(-2.148874) q[2];
sx q[2];
rz(-2.9774418) q[2];
rz(0.74603355) q[3];
sx q[3];
rz(-1.7697216) q[3];
sx q[3];
rz(3.0677838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1604851) q[0];
sx q[0];
rz(-1.2657413) q[0];
sx q[0];
rz(0.72738457) q[0];
rz(2.6630867) q[1];
sx q[1];
rz(-1.6886657) q[1];
sx q[1];
rz(-0.7368288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19664431) q[0];
sx q[0];
rz(-1.6262615) q[0];
sx q[0];
rz(-1.7086298) q[0];
x q[1];
rz(1.9987891) q[2];
sx q[2];
rz(-0.30747947) q[2];
sx q[2];
rz(-2.3124419) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9068272) q[1];
sx q[1];
rz(-1.2052457) q[1];
sx q[1];
rz(-0.91616456) q[1];
rz(-0.013929587) q[3];
sx q[3];
rz(-1.1831814) q[3];
sx q[3];
rz(-2.882021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9937146) q[2];
sx q[2];
rz(-1.0241877) q[2];
sx q[2];
rz(2.4564339) q[2];
rz(-1.0088751) q[3];
sx q[3];
rz(-2.4515371) q[3];
sx q[3];
rz(2.8907997) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20550263) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(-1.093338) q[0];
rz(3.1178442) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(2.645983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4706866) q[0];
sx q[0];
rz(-1.5829594) q[0];
sx q[0];
rz(-0.06605102) q[0];
x q[1];
rz(1.3643814) q[2];
sx q[2];
rz(-1.9999747) q[2];
sx q[2];
rz(2.6292554) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62140761) q[1];
sx q[1];
rz(-1.2222792) q[1];
sx q[1];
rz(-1.7640616) q[1];
x q[2];
rz(-1.0823332) q[3];
sx q[3];
rz(-1.9562072) q[3];
sx q[3];
rz(-3.1311991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89639837) q[2];
sx q[2];
rz(-0.79309016) q[2];
sx q[2];
rz(2.8450361) q[2];
rz(-2.631393) q[3];
sx q[3];
rz(-1.5254131) q[3];
sx q[3];
rz(0.27142522) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940014) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(1.9872794) q[0];
rz(-1.1555903) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(-1.4591699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4302368) q[0];
sx q[0];
rz(-0.85604307) q[0];
sx q[0];
rz(-1.6011333) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5695249) q[2];
sx q[2];
rz(-2.4172815) q[2];
sx q[2];
rz(1.6189761) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.3698813) q[1];
sx q[1];
rz(-2.071131) q[1];
sx q[1];
rz(2.9369563) q[1];
rz(-1.0351712) q[3];
sx q[3];
rz(-1.0058306) q[3];
sx q[3];
rz(-1.21711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.83361721) q[2];
sx q[2];
rz(-0.54215446) q[2];
sx q[2];
rz(-1.3758434) q[2];
rz(0.086325072) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(-1.8858006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1983222) q[0];
sx q[0];
rz(-2.4477796) q[0];
sx q[0];
rz(2.2721403) q[0];
rz(-2.4122639) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(-0.23390153) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1764404) q[0];
sx q[0];
rz(-1.6185221) q[0];
sx q[0];
rz(-0.41712572) q[0];
rz(2.9314165) q[2];
sx q[2];
rz(-1.4769722) q[2];
sx q[2];
rz(1.2503586) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6944551) q[1];
sx q[1];
rz(-1.4233973) q[1];
sx q[1];
rz(2.6958392) q[1];
rz(-1.2083268) q[3];
sx q[3];
rz(-2.3911282) q[3];
sx q[3];
rz(-0.01736162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0539315) q[2];
sx q[2];
rz(-1.6475185) q[2];
sx q[2];
rz(1.9653448) q[2];
rz(-0.67115274) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(1.5161071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43435708) q[0];
sx q[0];
rz(-0.86152995) q[0];
sx q[0];
rz(-1.0435411) q[0];
rz(-0.097298233) q[1];
sx q[1];
rz(-1.0568591) q[1];
sx q[1];
rz(-1.1204488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53172511) q[0];
sx q[0];
rz(-1.6311444) q[0];
sx q[0];
rz(0.97452314) q[0];
rz(-pi) q[1];
rz(-0.78586833) q[2];
sx q[2];
rz(-0.95986569) q[2];
sx q[2];
rz(-2.9402006) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35353684) q[1];
sx q[1];
rz(-2.9573943) q[1];
sx q[1];
rz(-0.26923979) q[1];
rz(-1.2213732) q[3];
sx q[3];
rz(-0.50822483) q[3];
sx q[3];
rz(2.1247381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6843159) q[2];
sx q[2];
rz(-2.5280759) q[2];
sx q[2];
rz(1.3555917) q[2];
rz(1.5654303) q[3];
sx q[3];
rz(-1.0836982) q[3];
sx q[3];
rz(-2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044357) q[0];
sx q[0];
rz(-2.5813527) q[0];
sx q[0];
rz(-0.78975633) q[0];
rz(-2.4850028) q[1];
sx q[1];
rz(-1.1142535) q[1];
sx q[1];
rz(2.5744892) q[1];
rz(0.91083679) q[2];
sx q[2];
rz(-1.177236) q[2];
sx q[2];
rz(0.81753035) q[2];
rz(-0.58224596) q[3];
sx q[3];
rz(-2.7979294) q[3];
sx q[3];
rz(2.9002849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
