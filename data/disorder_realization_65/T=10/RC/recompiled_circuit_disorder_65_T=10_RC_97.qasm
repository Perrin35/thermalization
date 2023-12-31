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
rz(2.066943) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5045835) q[0];
sx q[0];
rz(-1.5890117) q[0];
sx q[0];
rz(1.5760742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3518402) q[2];
sx q[2];
rz(-1.4573922) q[2];
sx q[2];
rz(1.5577158) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.97779146) q[1];
sx q[1];
rz(-1.692354) q[1];
sx q[1];
rz(-1.6519283) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1172469) q[3];
sx q[3];
rz(-1.082886) q[3];
sx q[3];
rz(0.32353668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82912123) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(1.6281698) q[2];
rz(1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-0.16076316) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.5011903) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26186865) q[0];
sx q[0];
rz(-0.12517087) q[0];
sx q[0];
rz(-2.9414888) q[0];
x q[1];
rz(-1.00374) q[2];
sx q[2];
rz(-2.6008285) q[2];
sx q[2];
rz(-2.8218249) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0928287) q[1];
sx q[1];
rz(-1.7424912) q[1];
sx q[1];
rz(-2.462639) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9618481) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(-2.5824576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(-1.9619933) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.5037781) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8684981) q[0];
sx q[0];
rz(-1.4833741) q[0];
sx q[0];
rz(0.065386535) q[0];
rz(-2.6659127) q[2];
sx q[2];
rz(-2.1026346) q[2];
sx q[2];
rz(1.4059517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5895264) q[1];
sx q[1];
rz(-0.30391903) q[1];
sx q[1];
rz(-2.9986831) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.740326) q[3];
sx q[3];
rz(-2.4781514) q[3];
sx q[3];
rz(-2.0142114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(0.7437931) q[2];
rz(0.25742325) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(0.58810365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(0.60032088) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1968397) q[0];
sx q[0];
rz(-1.6819994) q[0];
sx q[0];
rz(0.24223321) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7096268) q[2];
sx q[2];
rz(-2.1199806) q[2];
sx q[2];
rz(-2.2941022) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1059873) q[1];
sx q[1];
rz(-2.3590901) q[1];
sx q[1];
rz(-1.1483907) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2265008) q[3];
sx q[3];
rz(-2.2831884) q[3];
sx q[3];
rz(-0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2465308) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(1.2459922) q[3];
sx q[3];
rz(-2.1462838) q[3];
sx q[3];
rz(-0.016383735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.65288654) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(1.0523798) q[0];
rz(-1.2106238) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563755) q[0];
sx q[0];
rz(-1.7612805) q[0];
sx q[0];
rz(2.9604244) q[0];
x q[1];
rz(-1.1341368) q[2];
sx q[2];
rz(-1.1982802) q[2];
sx q[2];
rz(-2.8357752) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9235232) q[1];
sx q[1];
rz(-1.4265718) q[1];
sx q[1];
rz(1.6378535) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98975011) q[3];
sx q[3];
rz(-0.61468609) q[3];
sx q[3];
rz(-0.97782545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-2.4464567) q[2];
sx q[2];
rz(1.0926251) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(-2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(-1.8238235) q[0];
rz(-2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(-0.96907369) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1554402) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(3.100813) q[0];
x q[1];
rz(1.6386119) q[2];
sx q[2];
rz(-0.8578476) q[2];
sx q[2];
rz(2.3358047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5507257) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(-2.4025737) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.045741) q[3];
sx q[3];
rz(-1.1610371) q[3];
sx q[3];
rz(2.4214782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(0.43760854) q[2];
rz(-3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.5313914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(1.6925905) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(-1.70599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5098269) q[0];
sx q[0];
rz(-1.0477715) q[0];
sx q[0];
rz(-2.564389) q[0];
rz(-1.5905459) q[2];
sx q[2];
rz(-1.5891468) q[2];
sx q[2];
rz(1.2839279) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3064733) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.1919828) q[1];
x q[2];
rz(1.0649879) q[3];
sx q[3];
rz(-1.9516264) q[3];
sx q[3];
rz(-1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4009565) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(2.1334355) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.290264) q[0];
sx q[0];
rz(-2.2785318) q[0];
sx q[0];
rz(0.1329578) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-0.98633352) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62194659) q[0];
sx q[0];
rz(-1.6462109) q[0];
sx q[0];
rz(1.7072862) q[0];
rz(-pi) q[1];
rz(0.021990701) q[2];
sx q[2];
rz(-1.6985661) q[2];
sx q[2];
rz(2.2760504) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3268765) q[1];
sx q[1];
rz(-1.5183581) q[1];
sx q[1];
rz(-2.609842) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7563617) q[3];
sx q[3];
rz(-1.1572596) q[3];
sx q[3];
rz(0.47172771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.075295538) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(2.1832809) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(2.8310006) q[0];
rz(-2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(2.2176567) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5101178) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(0.60879137) q[0];
rz(0.02248259) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(1.2475916) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.20272217) q[1];
sx q[1];
rz(-0.67516967) q[1];
sx q[1];
rz(-2.2321781) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4531035) q[3];
sx q[3];
rz(-2.2985035) q[3];
sx q[3];
rz(1.9635995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49736398) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(1.1692858) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(-0.038169233) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725806) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(1.5241745) q[0];
rz(-pi) q[1];
rz(3.0496809) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(-0.11937571) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3341748) q[1];
sx q[1];
rz(-1.5290698) q[1];
sx q[1];
rz(1.9445022) q[1];
rz(0.90922728) q[3];
sx q[3];
rz(-1.070676) q[3];
sx q[3];
rz(1.3487032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.9434628) q[2];
rz(1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(-1.5126363) q[2];
sx q[2];
rz(-1.460961) q[2];
sx q[2];
rz(0.030208781) q[2];
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
