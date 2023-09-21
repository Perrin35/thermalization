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
rz(-2.3518402) q[2];
sx q[2];
rz(-1.6842004) q[2];
sx q[2];
rz(1.5577158) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5731116) q[1];
sx q[1];
rz(-2.9955578) q[1];
sx q[1];
rz(-0.58574974) q[1];
rz(2.4514276) q[3];
sx q[3];
rz(-0.65342045) q[3];
sx q[3];
rz(2.0131468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-0.99779904) q[2];
sx q[2];
rz(-1.5134229) q[2];
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
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56869498) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.6404023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0813635) q[0];
sx q[0];
rz(-1.6934568) q[0];
sx q[0];
rz(-1.5958022) q[0];
rz(1.00374) q[2];
sx q[2];
rz(-2.6008285) q[2];
sx q[2];
rz(-0.31976779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.048763927) q[1];
sx q[1];
rz(-1.7424912) q[1];
sx q[1];
rz(-0.6789536) q[1];
rz(-0.0060175671) q[3];
sx q[3];
rz(-1.179751) q[3];
sx q[3];
rz(2.1276377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.286065) q[2];
sx q[2];
rz(-0.76670402) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1353772) q[0];
sx q[0];
rz(-2.0897946) q[0];
sx q[0];
rz(2.0626383) q[0];
rz(1.1795993) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(1.6378145) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2730946) q[0];
sx q[0];
rz(-1.4833741) q[0];
sx q[0];
rz(-3.0762061) q[0];
rz(-2.1554699) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(-0.090678064) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9864195) q[1];
sx q[1];
rz(-1.6134312) q[1];
sx q[1];
rz(-2.8405872) q[1];
x q[2];
rz(0.13111968) q[3];
sx q[3];
rz(-2.2230806) q[3];
sx q[3];
rz(0.91344792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-2.3977996) q[2];
rz(-2.8841694) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96823111) q[0];
sx q[0];
rz(-0.1156725) q[0];
sx q[0];
rz(2.0898576) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-2.5225263) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7950492) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(1.6853149) q[0];
rz(-pi) q[1];
rz(2.7096268) q[2];
sx q[2];
rz(-2.1199806) q[2];
sx q[2];
rz(-2.2941022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1059873) q[1];
sx q[1];
rz(-0.78250256) q[1];
sx q[1];
rz(-1.1483907) q[1];
rz(0.37255128) q[3];
sx q[3];
rz(-2.3636892) q[3];
sx q[3];
rz(2.3140964) q[3];
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
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4887061) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(-1.2106238) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(2.2573684) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.592504) q[0];
sx q[0];
rz(-1.3929402) q[0];
sx q[0];
rz(-1.7643719) q[0];
rz(-pi) q[1];
rz(-2.7344633) q[2];
sx q[2];
rz(-1.9756769) q[2];
sx q[2];
rz(-1.4332353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.65518889) q[1];
sx q[1];
rz(-2.9826394) q[1];
sx q[1];
rz(0.43227355) q[1];
x q[2];
rz(-2.7719284) q[3];
sx q[3];
rz(-1.0677932) q[3];
sx q[3];
rz(-0.30077416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(-2.8975899) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7364175) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(1.3177692) q[0];
rz(0.70619839) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(2.172519) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5613865) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(-1.6798621) q[0];
rz(-pi) q[1];
rz(1.5029807) q[2];
sx q[2];
rz(-2.2837451) q[2];
sx q[2];
rz(-0.80578795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59086696) q[1];
sx q[1];
rz(-0.64462763) q[1];
sx q[1];
rz(2.4025737) q[1];
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
x q[1];
rz(-0.51450729) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(2.7039841) q[2];
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
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83201927) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(-1.3597885) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.4356027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3754417) q[0];
sx q[0];
rz(-2.0631844) q[0];
sx q[0];
rz(2.173461) q[0];
rz(-pi) q[1];
rz(-1.5905459) q[2];
sx q[2];
rz(-1.5524459) q[2];
sx q[2];
rz(-1.2839279) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83511931) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.9496099) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88019754) q[3];
sx q[3];
rz(-0.62303632) q[3];
sx q[3];
rz(-2.8043889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.2601968) q[2];
sx q[2];
rz(-0.17399542) q[2];
rz(-1.0081572) q[3];
sx q[3];
rz(-2.5139136) q[3];
sx q[3];
rz(-0.15667668) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.290264) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-3.0086349) q[0];
rz(2.6538387) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.7763604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2030905) q[0];
sx q[0];
rz(-1.4346968) q[0];
sx q[0];
rz(-0.076119856) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6985967) q[2];
sx q[2];
rz(-1.5489849) q[2];
sx q[2];
rz(-2.4335361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9284968) q[1];
sx q[1];
rz(-2.1017385) q[1];
sx q[1];
rz(-1.5099768) q[1];
rz(-pi) q[2];
rz(-1.385231) q[3];
sx q[3];
rz(-1.1572596) q[3];
sx q[3];
rz(-2.6698649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0662971) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(2.7189642) q[2];
rz(2.1832809) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(-2.8310006) q[0];
rz(2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(0.92393595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7565219) q[0];
sx q[0];
rz(-1.1413045) q[0];
sx q[0];
rz(-0.71682741) q[0];
rz(0.02248259) q[2];
sx q[2];
rz(-2.6528931) q[2];
sx q[2];
rz(1.894001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20272217) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(-0.90941456) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1903673) q[3];
sx q[3];
rz(-0.95622548) q[3];
sx q[3];
rz(0.28704498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6442287) q[2];
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
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-0.099530749) q[0];
sx q[0];
rz(-1.3884397) q[0];
rz(3.1346079) q[1];
sx q[1];
rz(-2.3313637) q[1];
sx q[1];
rz(3.1034234) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6258433) q[0];
sx q[0];
rz(-0.96228296) q[0];
sx q[0];
rz(0.032511062) q[0];
rz(0.091911749) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(0.11937571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0109509) q[1];
sx q[1];
rz(-0.37591939) q[1];
sx q[1];
rz(-1.6846659) q[1];
rz(-pi) q[2];
rz(2.2977703) q[3];
sx q[3];
rz(-2.3355964) q[3];
sx q[3];
rz(-0.7741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1756211) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(-2.0576058) q[3];
sx q[3];
rz(-2.4626648) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(1.4981131) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-0.48508587) q[2];
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